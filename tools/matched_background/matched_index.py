# %% matched_index.py: hg38 covariate index for PROMETHEUS Mode A (docs/PROMETHEUS_PLAN.md, C1 and step 1)
"""
Tiles the genome into W-bp windows (step W) and records, per window, the covariates the matched sets are drawn on:
  gc        G+C fraction over non-N bases
  cpg_oe    CpG observed/expected = n_CG * n_valid / (n_C * n_G)   (0 when n_C or n_G is 0)
  rep       soft-masked fraction (RepeatMasker + TRF lowercase in the UCSC assembly)
  rep_class RepeatMasker repClass at the window centre (rmsk.txt.gz; 0 = none)
  tss       distance from the window centre to the nearest annotated TSS (GENCODE basic transcripts)
  n_valid   non-N bases
Raw covariates only: the buckets used for matching are registered with the Mode A round, never tuned here.
The same function annotates arbitrary BED regions (fg sets), so fg and index covariates are computed identically.

Usage:
  python matched_index.py index    <hg38.fa.gz> <rmsk.txt.gz> <gencode.gtf.gz> <out_dir> [W=201]
  python matched_index.py annotate <index_dir> <regions.bed> <out.tsv>
"""
import gzip
import json
import os
import re
import sys

import numpy as np

CHROM = re.compile(r'chr[0-9XY]+')
LUT_GC = np.zeros(256, np.uint8)
LUT_N = np.zeros(256, np.uint8)
LUT_LOW = np.zeros(256, np.uint8)
LUT_C = np.zeros(256, np.uint8)
LUT_G = np.zeros(256, np.uint8)
for ch in b'GCgc':
    LUT_GC[ch] = 1
for ch in b'Nn':
    LUT_N[ch] = 1
for ch in b'acgtn':
    LUT_LOW[ch] = 1
for ch in b'Cc':
    LUT_C[ch] = 1
for ch in b'Gg':
    LUT_G[ch] = 1


def read_genome(path):
    seqs, name, buf = {}, None, []
    with gzip.open(path, 'rb') as fh:
        for line in fh:
            if line.startswith(b'>'):
                if name and CHROM.fullmatch(name):
                    seqs[name] = np.frombuffer(b''.join(buf), np.uint8)
                name, buf = line[1:].split()[0].decode(), []
            elif name and CHROM.fullmatch(name):
                buf.append(line.rstrip())
    if name and CHROM.fullmatch(name):
        seqs[name] = np.frombuffer(b''.join(buf), np.uint8)
    return seqs


def read_tss(gtf):
    tss = {}
    with gzip.open(gtf, 'rt') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.split('\t', 8)
            if f[2] != 'transcript' or not CHROM.fullmatch(f[0]):
                continue
            tss.setdefault(f[0], []).append(int(f[3]) - 1 if f[6] == '+' else int(f[4]) - 1)
    return {c: np.unique(np.array(v, np.int64)) for c, v in tss.items()}


def read_rmsk(path):
    """{chrom: (starts, ends, class_ids)} sorted by start; class names indexed from 1."""
    names, rows = {}, {}
    with gzip.open(path, 'rt') as fh:
        for line in fh:
            f = line.split('\t')
            c = f[5]
            if not CHROM.fullmatch(c):
                continue
            k = names.setdefault(f[11], len(names) + 1)
            rows.setdefault(c, []).append((int(f[6]), int(f[7]), k))
    out = {}
    for c, r in rows.items():
        a = np.array(sorted(r), np.int64)
        out[c] = (a[:, 0], a[:, 1], a[:, 2].astype(np.uint8))
    return out, [n for n, _ in sorted(names.items(), key=lambda x: x[1])]


def class_at(rm, c, pos):
    """repClass covering each position (the latest-starting repeat that covers it; 0 = none)."""
    if c not in rm:
        return np.zeros(len(pos), np.uint8)
    st, en, k = rm[c]
    i = np.searchsorted(st, pos, side='right') - 1
    ok = (i >= 0) & (en[np.maximum(i, 0)] > pos)
    return np.where(ok, k[np.maximum(i, 0)], 0).astype(np.uint8)


def tss_dist(tss, c, pos):
    if c not in tss:
        return np.full(len(pos), 2 ** 31 - 1, np.int64)
    t = tss[c]
    i = np.searchsorted(t, pos)
    lo = np.abs(pos - t[np.clip(i - 1, 0, len(t) - 1)])
    hi = np.abs(t[np.clip(i, 0, len(t) - 1)] - pos)
    return np.minimum(lo, hi)


def covariates(seq, starts, ends):
    """Per-region covariates from cumulative sums over a chromosome's byte array (any region widths)."""
    def cs(x):
        return np.concatenate([[0], np.cumsum(x, dtype=np.int64)])
    n = cs(LUT_N[seq])
    g = cs(LUT_GC[seq])
    low = cs(LUT_LOW[seq])
    cc = cs(LUT_C[seq])
    gg = cs(LUT_G[seq])
    cg = cs(np.concatenate([LUT_C[seq[:-1]] & LUT_G[seq[1:]], [0]]).astype(np.uint8))
    s, e = np.asarray(starts, np.int64), np.asarray(ends, np.int64)
    nv = (e - s) - (n[e] - n[s])
    nc, ng = cc[e] - cc[s], gg[e] - gg[s]
    ncg = cg[e - 1] - cg[s]                       # CG dinucleotides starting inside [s, e-1)
    with np.errstate(divide='ignore', invalid='ignore'):
        gc = np.where(nv > 0, (g[e] - g[s]) / nv, 0.0)
        oe = np.where((nc > 0) & (ng > 0), ncg * nv / (nc * ng), 0.0)
        rep = np.where(e > s, (low[e] - low[s]) / (e - s), 0.0)
    return nv, gc, oe, rep


def build_index(fa, rmsk, gtf, out, W=201):
    os.makedirs(out, exist_ok=True)
    G = read_genome(fa)
    tss = read_tss(gtf)
    rm, classes = read_rmsk(rmsk)
    chroms = sorted(G, key=lambda c: (len(c) > 5, c))
    cols = {k: [] for k in ('chrom', 'start', 'n_valid', 'gc', 'cpg_oe', 'rep', 'rep_class', 'tss')}
    for ci, c in enumerate(chroms):
        L = len(G[c])
        st = np.arange(0, L - W + 1, W, dtype=np.int64)
        nv, gc, oe, rep = covariates(G[c], st, st + W)
        mid = st + W // 2
        cols['chrom'].append(np.full(len(st), ci, np.uint8))
        cols['start'].append(st.astype(np.uint32))
        cols['n_valid'].append(nv.astype(np.uint16))
        cols['gc'].append(gc.astype(np.float32))
        cols['cpg_oe'].append(oe.astype(np.float32))
        cols['rep'].append(rep.astype(np.float32))
        cols['rep_class'].append(class_at(rm, c, mid))
        cols['tss'].append(np.minimum(tss_dist(tss, c, mid), 2 ** 31 - 1).astype(np.int32))
        print(c, len(st), flush=True)
    np.savez(os.path.join(out, 'windows.npz'), **{k: np.concatenate(v) for k, v in cols.items()})
    np.savez(os.path.join(out, 'tss.npz'), **tss)
    rs = {}
    for c, (a, b, k) in rm.items():
        rs[c + '_s'], rs[c + '_e'], rs[c + '_k'] = a, b, k
    np.savez(os.path.join(out, 'rmsk.npz'), **rs)
    json.dump(dict(W=W, chroms=chroms, rep_classes=['none'] + classes, genome=os.path.abspath(fa),
                   rmsk=os.path.abspath(rmsk), gtf=os.path.abspath(gtf)),
              open(os.path.join(out, 'index.json'), 'w'), indent=1)


def annotate(index_dir, bed, out):
    meta = json.load(open(os.path.join(index_dir, 'index.json')))
    G = read_genome(meta['genome'])
    tss = dict(np.load(os.path.join(index_dir, 'tss.npz')))
    r = np.load(os.path.join(index_dir, 'rmsk.npz'))
    rm = {c: (r[c + '_s'], r[c + '_e'], r[c + '_k']) for c in meta['chroms'] if c + '_s' in r}
    rows = [l.rstrip('\n').split('\t') for l in open(bed) if l.strip()]
    res = [None] * len(rows)
    by = {}
    for i, f in enumerate(rows):
        by.setdefault(f[0], []).append(i)
    for c, idx in by.items():
        st = np.array([int(rows[i][1]) for i in idx])
        en = np.array([int(rows[i][2]) for i in idx])
        nv, gc, oe, rep = covariates(G[c], st, en)
        mid = (st + en) // 2
        k = class_at(rm, c, mid)
        d = tss_dist(tss, c, mid)
        for j, i in enumerate(idx):
            res[i] = (nv[j], gc[j], oe[j], rep[j], meta['rep_classes'][k[j]], d[j])
    with open(out, 'w') as fh:
        fh.write('name\tchrom\tstart\tend\tn_valid\tgc\tcpg_oe\trep\trep_class\ttss\n')
        for f, x in zip(rows, res):
            fh.write(f'{f[3] if len(f) > 3 else "."}\t{f[0]}\t{f[1]}\t{f[2]}\t{x[0]}\t{x[1]:.4f}\t{x[2]:.4f}\t'
                     f'{x[3]:.4f}\t{x[4]}\t{x[5]}\n')


if __name__ == '__main__':
    if sys.argv[1] == 'index':
        build_index(*sys.argv[2:6], W=int(sys.argv[6]) if len(sys.argv) > 6 else 201)
    elif sys.argv[1] == 'annotate':
        annotate(*sys.argv[2:5])
    else:
        sys.exit(__doc__)
