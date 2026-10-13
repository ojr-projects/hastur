# %% prep_chip.py: ChIP inputs for the dev/test splits (splits/chip_split.tsv)
"""
Builds, per factor: fg.fa = top TOP conservative-IDR summits ±HALF bp (ENCODE GRCh38), bg.fa = GC-matched genomic
windows not overlapping any peak, BG_MULT times as many as the foreground.

The sequence logic is copied from Comsequent's benchmarks/phase5_negctrl/50_chip_negctrl.py (build), with the same
RNG seed. The background sampler is sequential, so the first n_fg background sequences are exactly Comsequent's P5
background (bg1x.fa is written separately for like-for-like comparison with P5).
Usage (box): python data/prep_chip.py [dev|test|all]     env: CHIP_WORK (default ~/data/chip), CHIP_BG_MULT (10)
Each FASTA has a BED alongside (same windows, same order; 0-based half-open). For a unit built before the BEDs
existed, the windows are recomputed, checked identical to the existing FASTAs, and only the BEDs are written.
"""
import bisect
import concurrent.futures as cf
import gzip
import multiprocessing
import json
import os
import re
import subprocess
import sys
import urllib.request

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
SPLIT = os.environ.get('CHIP_SPLIT', os.path.join(HERE, '..', 'splits', 'chip_split.tsv'))   # test2: chip_test2_split.tsv
WORK = os.path.expanduser(os.environ.get('CHIP_WORK', '~/data/chip'))
DATA = os.path.join(WORK, '_downloads')
BG_MULT = int(os.environ.get('CHIP_BG_MULT', '10'))
MASKED = os.environ.get('CHIP_MASKED') == '1'   # soft-masked repeats -> N in the written sequences (same windows)
TOP, HALF = 5000, 100
GENOME_URL = 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz'


def get(path):
    req = urllib.request.Request('https://www.encodeproject.org' + path, headers={'Accept': 'application/json'})
    with urllib.request.urlopen(req, timeout=120) as r:
        return json.load(r)


def peaks_file(factor, cell):
    q = (f'/search/?type=File&output_type=conservative+IDR+thresholded+peaks&assembly=GRCh38&file_format=bed'
         f'&status=released&target.label={factor}&biosample_ontology.term_name={cell}&limit=all&frame=object')
    try:
        g = get(q)['@graph']
    except urllib.error.HTTPError as e:     # ENCODE returns 404 for an empty search
        if e.code == 404:
            g = []
        else:
            raise
    if not g:
        raise LookupError(f'no conservative IDR peaks for {factor} in {cell}')
    f = sorted(g, key=lambda x: x.get('date_created', ''))[-1]
    out = os.path.join(DATA, f'{factor}_{f["accession"]}.bed.gz')
    if not os.path.exists(out):
        subprocess.run(['curl', '-sSL', '--retry', '5', '-o', out + '.part', 'https://www.encodeproject.org' + f['href']],
                       check=True)
        os.replace(out + '.part', out)
    return out, f['accession']


def genome():
    fa = os.path.join(DATA, 'hg38.fa.gz')
    if not os.path.exists(fa):
        subprocess.run(['curl', '-sSL', '--retry', '5', '-o', fa + '.part', GENOME_URL], check=True)
        os.replace(fa + '.part', fa)
    seqs, name, buf = {}, None, []
    with gzip.open(fa, 'rt') as fh:
        for line in fh:
            if line.startswith('>'):
                if name:
                    seqs[name] = ''.join(buf) if MASKED else ''.join(buf).upper()
                name, buf = line[1:].split()[0], []
            elif name and re.fullmatch(r'chr[0-9XY]+', name):
                buf.append(line.strip())
    if name and buf:
        seqs[name] = ''.join(buf) if MASKED else ''.join(buf).upper()
    return {k: v for k, v in seqs.items() if re.fullmatch(r'chr[0-9XY]+', k)}


def gc(s):
    return (s.count('G') + s.count('C')) / max(1, len(s))


def mask(s):
    """Soft-masked (lowercase) bases -> N; the caller decides on s.upper(), so the windows are those of the unmasked build."""
    return re.sub('[acgtn]', 'N', s)


def write_bed(path, coords):
    with open(path + '.tmp', 'w') as fh:
        for i, (c, a, b) in enumerate(coords):
            fh.write(f'{c}\t{a}\t{b}\ts{i}\n')
    os.replace(path + '.tmp', path)


def read_fa(path):
    return [l.rstrip('\n') for l in open(path) if not l.startswith('>')]


def write_fa(path, seqs):
    with open(path + '.tmp', 'w') as fh:
        for i, s in enumerate(seqs):
            fh.write(f'>s{i}\n{s}\n')
    os.replace(path + '.tmp', path)


def build(factor, cell, G):
    d = os.path.join(WORK, factor)
    os.makedirs(d, exist_ok=True)
    fg, bg, bg1 = (os.path.join(d, f) for f in ('fg.fa', 'bg.fa', 'bg1x.fa'))
    beds = [os.path.join(d, f) for f in ('fg.bed', 'bg.bed', 'bg1x.bed')]
    if os.path.exists(bg1) and all(os.path.exists(b) for b in beds):
        return 'done'
    backfill = os.path.exists(bg1)
    ref = os.environ.get('CHIP_REF')     # rebuild from the recorded peak file (not ENCODE's latest)
    if ref and os.path.exists(os.path.join(ref, factor, 'info.json')):
        acc = json.load(open(os.path.join(ref, factor, 'info.json')))['peaks']
        bed = os.path.join(DATA, f'{factor}_{acc}.bed.gz')
    else:
        bed, acc = peaks_file(factor, cell)
    p = pd.read_csv(bed, sep='\t', header=None)
    p = p[p[0].isin(G)].sort_values(6, ascending=False)
    occupied = {}
    for r in p.itertuples():
        occupied.setdefault(r[1], []).append((r[2], r[3]))
    # O(log n) overlap test: peaks sorted by start, prefix maximum of their ends (same decisions as a linear scan)
    occ_idx = {}
    for c, iv in occupied.items():
        iv = sorted(iv)
        occ_idx[c] = ([a for a, _ in iv], np.maximum.accumulate([b for _, b in iv]))

    def overlaps(c, st, en):
        if c not in occ_idx:
            return False
        starts, pmax = occ_idx[c]
        i = bisect.bisect_left(starts, en)
        return i > 0 and pmax[i - 1] > st
    top, top_u, top_c = [], [], []
    for r in p.itertuples():
        c = r[2] + r[10]
        s0 = G[r[1]][c - HALF:c + HALF + 1]
        s = s0.upper()
        if len(s) == 2 * HALF + 1 and 'N' not in s:
            top.append(mask(s0) if MASKED else s)
            top_u.append(s)
            top_c.append((r[1], c - HALF, c + HALF + 1))
        if len(top) >= TOP:
            break
    rng = np.random.default_rng(20261009)                  # same seed as Comsequent P5
    bins = np.clip((np.array([gc(s) for s in top_u]) * 50).astype(int), 0, 49)   # decisions on unmasked windows
    need1 = np.bincount(bins, minlength=50)
    out, firstn, out_c = [], None, []
    chroms = list(G)
    lens = np.array([len(G[c]) for c in chroms], float)
    for rep in range(BG_MULT):                             # rep 0 reproduces P5's background exactly
        have, tries = np.zeros(50, int), 0
        while have.sum() < need1.sum() and tries < 200 * TOP:
            tries += 1
            c = chroms[rng.choice(len(chroms), p=lens / lens.sum())]
            st = int(rng.integers(0, len(G[c]) - 2 * HALF - 1))
            s0 = G[c][st:st + 2 * HALF + 1]
            s = s0.upper()
            if 'N' in s or overlaps(c, st, st + 2 * HALF + 1):
                continue
            b = min(49, int(gc(s) * 50))
            if have[b] < need1[b]:
                out.append(mask(s0) if MASKED else s)
                out_c.append((c, st, st + 2 * HALF + 1))
                have[b] += 1
        if rep == 0:
            firstn = list(out)
    if backfill:
        for path, ss in ((fg, top), (bg, out), (bg1, firstn)):
            if read_fa(path) != ss:
                raise RuntimeError(f'rebuilt windows differ from {path}; BED not written')
    else:
        write_fa(fg, top)
        write_fa(bg, out)
        write_fa(bg1, firstn)
    write_bed(beds[0], top_c)
    write_bed(beds[1], out_c)
    write_bed(beds[2], out_c[:len(firstn)])
    if backfill:
        return 'BED backfilled (windows identical)'
    json.dump(dict(factor=factor, cell=cell, peaks=acc, n_fg=len(top), n_bg=len(out), n_bg1=len(firstn)),
              open(os.path.join(d, 'info.json'), 'w'))
    return f'{len(top)} fg, {len(out)} bg'


GENOME = None


def build_fork(factor, cell):
    return build(factor, cell, GENOME)


if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv) > 1 else 'dev'
    sp = pd.read_csv(SPLIT, sep='\t')
    if which != 'all':
        sp = sp[sp.split == which]
    os.makedirs(DATA, exist_ok=True)
    print('loading hg38', flush=True)
    G = genome()
    print(f'{len(G)} chromosomes', flush=True)
    globals()['GENOME'] = G   # inherited by forked workers (copy-on-write)
    with cf.ProcessPoolExecutor(int(os.environ.get('CHIP_WORKERS', '24')),
                                mp_context=multiprocessing.get_context('fork')) as ex:
        futs = {ex.submit(build_fork, r.factor, r.cell): r.factor for r in sp.itertuples()}
        for f in cf.as_completed(futs):
            try:
                print(futs[f], f.result(), flush=True)
            except Exception as e:  # noqa: BLE001
                print(futs[f], 'FAILED', repr(e), flush=True)
