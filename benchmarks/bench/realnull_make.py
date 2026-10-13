# %% realnull_make.py: real-genome motif-free pairs for FR1 (splits/protocol_FR1.md, nulls 2 and 3)
"""
<root>/<pair>/{fg.fa, bg.fa, truth.json(null=True)}; 20 replicates per kind; background 10x the foreground.
  human_exch_n{500,2000}   two disjoint random sets of hg38 windows (201 bp, no N)
  human_gc_n{500,2000}     fg: random windows with GC in [0.55, 0.65]; bg: random genome-wide windows
  <genome>_{intergenic,coding}_n500   exchangeable splits of random 100-bp windows lying wholly inside intergenic
                           (between annotated genes) or coding (CDS) intervals, for C. jejuni (AT-rich), E. coli and
                           B. cenocepacia (GC-rich)
Usage: python bench/realnull_make.py <root> [human|bacteria|all]     env SEED_OFFSET (confirmation: 5000)
fg.bed / bg.bed sit alongside (same windows and order). For a pair built before the BEDs existed, the windows are
recomputed, checked identical to the existing FASTAs, and only the BEDs are written.
"""
import gzip
import json
import os
import re
import sys

import numpy as np

HG38 = os.path.expanduser(os.environ.get('HG38', '~/data/chip/_downloads/hg38.fa.gz'))
MASKED = os.environ.get('MASKED') == '1'    # soft-masked repeats -> N in the written sequences (same windows)
GEN = '/home/oliver/code/comsequent/data/genomes'
BACT = {'cjejuni': 'Campylobacter_jejuni_RM1221/ncbi_dataset/data/GCF_000011865.1',
        'ecoli': 'Escherichia_coli_str_K-12_substr._MG1655/ncbi_dataset/data/GCF_000005845.2',
        'bceno': 'Burkholderia_cenocepacia_J2315/ncbi_dataset/data/GCF_000009485.1'}
REPS, BGX = 20, 10
OFF = int(os.environ.get('SEED_OFFSET', '0'))   # FR1 confirmation nulls: SEED_OFFSET=5000 (fresh replicates)


def read_fasta(path, keep=None):
    op = gzip.open if path.endswith('.gz') else open
    seqs, name, buf = {}, None, []
    with op(path, 'rt') as fh:
        for line in fh:
            if line.startswith('>'):
                if name and (keep is None or keep(name)):
                    seqs[name] = ''.join(buf) if MASKED else ''.join(buf).upper()
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name and (keep is None or keep(name)):
        seqs[name] = ''.join(buf) if MASKED else ''.join(buf).upper()
    return seqs


def done(d):
    return all(os.path.exists(os.path.join(d, f)) for f in ('truth.json', 'fg.bed', 'bg.bed'))


def write(d, fg, bg, kind, rep):
    """fg, bg: (seqs, keys, coords) from sample_windows."""
    backfill = os.path.exists(os.path.join(d, 'truth.json'))
    os.makedirs(d, exist_ok=True)
    for nm, (ss, _, cc) in (('fg', fg), ('bg', bg)):
        if backfill:
            old = [l.rstrip('\n') for l in open(os.path.join(d, nm + '.fa')) if not l.startswith('>')]
            if old != ss:
                raise RuntimeError(f'rebuilt windows differ from {d}/{nm}.fa; BED not written')
        else:
            with open(os.path.join(d, nm + '.fa'), 'w') as fh:
                for i, s in enumerate(ss):
                    fh.write(f'>s{i}\n{s}\n')
        with open(os.path.join(d, nm + '.bed'), 'w') as fh:
            for i, (c, a, b) in enumerate(cc):
                fh.write(f'{c}\t{a}\t{b}\ts{i}\n')
    if not backfill:
        json.dump(dict(dataset=os.path.basename(d), kind=kind, rep=rep, null=True, truth=[]),
                  open(os.path.join(d, 'truth.json'), 'w'))


def sample_windows(rng, intervals, seqs, n, w, cond=None, exclude=None):
    """n windows of width w drawn uniformly from (chrom, start, end) intervals; no N; optional GC condition."""
    lens = np.array([e - s - w for c, s, e in intervals], float)
    ok = lens > 0
    iv = [intervals[i] for i in np.flatnonzero(ok)]
    p = lens[ok] / lens[ok].sum()
    out, keys, coords = [], set(), []
    while len(out) < n:
        c, s, e = iv[rng.choice(len(iv), p=p)]
        st = int(rng.integers(s, e - w))
        if (c, st) in keys or (exclude is not None and (c, st) in exclude):
            continue
        x0 = seqs[c][st:st + w]
        x = x0.upper()                   # decisions on the unmasked window, so masked builds keep the same windows
        if 'N' in x or len(x) < w:
            continue
        if cond is not None and not cond(x):
            continue
        keys.add((c, st))
        coords.append((c, st, st + w))
        out.append(re.sub('[acgtn]', 'N', x0) if MASKED else x)
    return out, keys, coords


if __name__ == '__main__':
    root, which = sys.argv[1], (sys.argv[2] if len(sys.argv) > 2 else 'all')
    if which in ('human', 'all'):
        G = read_fasta(HG38, keep=lambda n: re.fullmatch(r'chr[0-9XY]+', n) is not None)
        iv = [(c, 0, len(s)) for c, s in G.items()]
        gc = lambda x: 0.55 <= (x.count('G') + x.count('C')) / len(x) <= 0.65
        for n in (500, 2000):
            for rep in range(REPS):
                rng = np.random.default_rng(770000 + n + rep + OFF)
                d = os.path.join(root, f'human_exch_n{n}__r{rep}')
                if not done(d):
                    fg = sample_windows(rng, iv, G, n, 201)
                    bg = sample_windows(rng, iv, G, BGX * n, 201, exclude=fg[1])
                    write(d, fg, bg, f'human_exch_n{n}', rep)
                d = os.path.join(root, f'human_gc_n{n}__r{rep}')
                if not done(d):
                    fg = sample_windows(rng, iv, G, n, 201, cond=gc)
                    bg = sample_windows(rng, iv, G, BGX * n, 201, exclude=fg[1])
                    write(d, fg, bg, f'human_gc_n{n}', rep)
            print('human', n, flush=True)
        del G
    if which in ('bacteria', 'all'):
        for g, rel in BACT.items():
            base = os.path.join(GEN, rel)
            fna = [f for f in os.listdir(base) if f.endswith('.fna')][0]
            S = read_fasta(os.path.join(base, fna))
            cds = {c: [] for c in S}
            genes = {c: [] for c in S}
            for line in open(os.path.join(base, 'genomic.gff')):
                if line.startswith('#'):
                    continue
                f = line.rstrip('\n').split('\t')
                if len(f) < 9 or f[0] not in S:
                    continue
                if f[2] == 'CDS':
                    cds[f[0]].append((int(f[3]) - 1, int(f[4])))
                if f[2] in ('gene', 'pseudogene'):
                    genes[f[0]].append((int(f[3]) - 1, int(f[4])))
            coding = [(c, a, b) for c, ab in cds.items() for a, b in ab]
            inter = []
            for c, ab in genes.items():
                ab = sorted(ab)
                end = 0
                for a, b in ab:
                    if a > end:
                        inter.append((c, end, a))
                    end = max(end, b)
                if end < len(S[c]):
                    inter.append((c, end, len(S[c])))
            for region, ivs in (('intergenic', inter), ('coding', coding)):
                for rep in range(REPS):
                    d = os.path.join(root, f'{g}_{region}_n500__r{rep}')
                    if done(d):
                        continue
                    rng = np.random.default_rng(880000 + rep + 31 * len(g) + len(region) + OFF)
                    fg = sample_windows(rng, ivs, S, 500, 100)
                    bg = sample_windows(rng, ivs, S, BGX * 500, 100, exclude=fg[1])
                    write(d, fg, bg, f'{g}_{region}', rep)
            print(g, flush=True)
