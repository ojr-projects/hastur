# %% prep_selex.py: HT-SELEX inputs (Jolma et al. 2013, ENA PRJEB3289) for the dev/test splits
"""
Per TF library (TF, barcode, batch), as in Comsequent's benchmarks/phase4_selex/40_selex_run.py:
  fg.fa = the final selection cycle (first N_FG reads without N), bg.fa = the cycle-0 (ZeroCycle) library with the same
  barcode (first N_BG reads; BG_MULT x larger than Comsequent's matched setting, as the theory recommends).
Reads are STREAMED (curl | gunzip, stopped after N reads) rather than downloading whole fastq files.
truth.json: every JASPAR 2024 CORE vertebrate (redundant) matrix whose name matches the TF: consensus strings
(Comsequent's rule: base >= 0.5, two-base letter >= 0.75) and the matrices themselves (for PWM similarity).
Usage (box): python data/prep_selex.py dev [max_libraries]     env: SELEX_WORK (~/data/selex), SELEX_N_FG, SELEX_N_BG
"""
import concurrent.futures as cf
import gzip
import json
import os
import subprocess
import sys
import threading

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
SPLIT = os.path.join(HERE, '..', 'splits', 'selex_split.tsv')
WORK = os.path.expanduser(os.environ.get('SELEX_WORK', '~/data/selex'))
DL = os.path.join(WORK, '_downloads')
N_FG = int(os.environ.get('SELEX_N_FG', '20000'))
N_BG = int(os.environ.get('SELEX_N_BG', '100000'))
ENA = ('https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB3289&result=read_run'
       '&fields=run_accession,sample_alias,fastq_bytes,fastq_ftp,read_count&format=tsv&limit=0')
JASPAR = 'https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_redundant_pfms_meme.txt'
PAIR = {frozenset('AG'): 'R', frozenset('CT'): 'Y', frozenset('CG'): 'S', frozenset('AT'): 'W', frozenset('GT'): 'K',
        frozenset('AC'): 'M'}
_SEM = threading.BoundedSemaphore(3)          # EBI refuses many parallel connections


def libraries():
    path = os.path.join(DL, 'ena_runs.tsv')
    if not os.path.exists(path):
        subprocess.run(['curl', '-sSL', '--retry', '5', '-o', path, ENA], check=True)
    d = pd.read_csv(path, sep='\t')
    a = d.sample_alias.str.split('_')
    d['tf'], d['barcode'], d['cycle'] = a.str[0], a.str[1], a.str[-1].astype(int)
    d['batch'] = a.apply(lambda x: '_'.join(x[2:-1]))
    zero = d[d.tf == 'ZeroCycle'].drop_duplicates('barcode').set_index('barcode')
    rows = []
    for (tf, bc, batch), g in d[d.tf != 'ZeroCycle'].groupby(['tf', 'barcode', 'batch']):
        if bc not in zero.index:
            continue
        top = g.sort_values('cycle').iloc[-1]
        rows.append(dict(lib=f'{tf}_{bc}_{batch}', tf=tf, barcode=bc, cycle=int(top.cycle), fg_run=top.run_accession,
                         fg_url='https://' + top.fastq_ftp.split(';')[0], bg_run=zero.loc[bc, 'run_accession'],
                         bg_url='https://' + zero.loc[bc, 'fastq_ftp'].split(';')[0]))
    return pd.DataFrame(rows)


def consensus(rows):
    s = []
    for r in rows:
        o = np.argsort(r)[::-1]
        if r[o[0]] >= 0.5:
            s.append('ACGT'[o[0]])
        elif r[o[0]] + r[o[1]] >= 0.75:
            s.append(PAIR[frozenset('ACGT'[o[0]] + 'ACGT'[o[1]])])
        else:
            s.append('N')
    return ''.join(s).strip('N')


def jaspar():
    path = os.path.join(DL, 'jaspar2024_vertebrates_redundant.meme')
    if not os.path.exists(path):
        subprocess.run(['curl', '-sSL', '--retry', '5', '-o', path, JASPAR], check=True)
    out, name, rows, rd = {}, None, [], False
    for line in open(path):
        if line.startswith('MOTIF'):
            if name and rows:
                out.setdefault(name, []).append(rows)
            f = line.split()
            name, rows, rd = (f[2] if len(f) > 2 else f[1]).upper(), [], False
        elif line.startswith('letter-probability'):
            rd = True
        elif rd:
            p = line.split()
            if len(p) == 4:
                rows.append([float(x) for x in p])
            elif rows:
                rd = False
    if name and rows:
        out.setdefault(name, []).append(rows)
    return out


def fetch(url, run, n, out):
    """First n reads without N from a gzipped fastq URL, streamed."""
    if os.path.exists(out):
        return
    with _SEM:
        p = subprocess.Popen(['curl', '-sSL', '--retry', '10', '--retry-connrefused', url], stdout=subprocess.PIPE)
        k = 0
        with gzip.open(p.stdout, 'rt') as fh, open(out + '.tmp', 'w') as fo:
            try:
                for i, line in enumerate(fh):
                    if i % 4 == 1:
                        s = line.strip().upper()
                        if s and 'N' not in s:
                            fo.write(f'>r{k}\n{s}\n')
                            k += 1
                            if k >= n:
                                break
            except (EOFError, OSError):
                pass
        p.kill()
        p.wait()
    if k < min(n, 1000):
        raise RuntimeError(f'{run}: only {k} reads')
    os.replace(out + '.tmp', out)


def build(r, J):
    d = os.path.join(WORK, r.lib)
    os.makedirs(d, exist_ok=True)
    if os.path.exists(os.path.join(d, 'truth.json')):
        return f'{r.lib} done'
    mats = J.get(r.tf.upper(), [])
    if not mats:
        return f'{r.lib} skipped: no JASPAR matrix for {r.tf}'
    fetch(r.fg_url, r.fg_run, N_FG, os.path.join(d, 'fg.fa'))
    bgc = os.path.join(DL, f'{r.bg_run}_{N_BG}.fa')      # cycle-0 libraries are shared across TFs
    fetch(r.bg_url, r.bg_run, N_BG, bgc)
    if not os.path.exists(os.path.join(d, 'bg.fa')):
        os.symlink(bgc, os.path.join(d, 'bg.fa'))
    json.dump(dict(dataset=r.lib, tf=r.tf, cycle=r.cycle, truth=[consensus(m) or 'N' for m in mats],
                   pwms=[np.asarray(m).tolist() for m in mats]), open(os.path.join(d, 'truth.json'), 'w'))
    return f'{r.lib} ok'


if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv) > 1 else 'dev'
    cap = int(sys.argv[2]) if len(sys.argv) > 2 else 10 ** 9
    os.makedirs(DL, exist_ok=True)
    sp = pd.read_csv(SPLIT, sep='\t')
    libs = libraries()
    libs = libs[libs.tf.isin(set(sp[sp.split == which].tf))].sort_values('lib')
    print(f'{len(libs)} {which} libraries for {libs.tf.nunique()} TFs; building up to {cap}', flush=True)
    libs = libs.head(cap)
    J = jaspar()
    with cf.ThreadPoolExecutor(6) as ex:
        for f in cf.as_completed([ex.submit(build, r, J) for r in libs.itertuples()]):
            try:
                print(f.result(), flush=True)
            except Exception as e:  # noqa: BLE001
                print('FAILED', repr(e), flush=True)
