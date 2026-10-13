# %% match_sets.py: Mode A matched backgrounds and strata (PROMETHEUS C1; docs/PROMETHEUS_PLAN.md §3)
"""
For every fg region (BED, 201 bp), draws K genome tiles from the hg38 covariate index (matched_index.py) in the same
covariate bucket. Writes, into <out_dir>:
  fg.fa        the fg sequences (re-extracted from the genome by coordinate, same order)
  bg.fa        K matched tiles per fg region, in fg order
  bg.bed       their coordinates
  strata.txt   one stratum id per sequence, fg then bg (input to `prometheus --strata`)
  match.json   buckets used, coarsening events, counts
Buckets (fixed here, registered with the Mode A round; informed by step 2b, where unmasked false calls on
GC-selected nulls sat in SINE copies and masked ones in GC-rich unique sequence):
  GC        0.05-wide bins
  CpG o/e   < 0.2, 0.2-0.4, 0.4-0.6, >= 0.6
  repeat    soft-masked share < 0.1, 0.1-0.5, >= 0.5; plus SINE vs other at the centre when >= 0.5
  TSS       (with --tss) < 1 kb, 1-5 kb, 5-50 kb, >= 50 kb
A bucket with fewer than 4K free tiles is coarsened by dropping TSS, then CpG, then the repeat class, then halving
the GC resolution (logged per region). Tiles overlapping any fg region or any --exclude interval are never drawn,
and no tile is drawn twice.
Usage: python match_sets.py <index_dir> <fg.bed> <out_dir> [--k 10] [--tss] [--exclude peaks.bed[.gz]] [--masked]
       [--seed a,b,...]
"""
import argparse
import gzip
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import matched_index as mi  # noqa: E402

W = 201


def bins(gc, oe, rep, rcls_sine, tss, level, use_tss):
    """Integer bucket keys at coarsening level 0..4 (arrays in, array out)."""
    gres = 0.05 if level < 4 else 0.10
    g = np.clip((gc / gres).astype(np.int64), 0, 20)
    o = np.digitize(oe, [0.2, 0.4, 0.6])
    r = np.digitize(rep, [0.1, 0.5])
    rs = np.where(r == 2, rcls_sine.astype(np.int64), 0)
    t = np.digitize(tss, [1000, 5000, 50000])
    if not use_tss or level >= 1:
        t = np.zeros_like(t)
    if level >= 2:
        o = np.zeros_like(o)
    if level >= 3:
        rs = np.zeros_like(rs)
    return (((g * 4 + o) * 3 + r) * 2 + rs) * 4 + t


def read_bed(p):
    op = gzip.open if p.endswith('.gz') else open
    with op(p, 'rt') as fh:
        return [(f[0], int(f[1]), int(f[2])) for f in (l.split('\t') for l in fh) if len(f) >= 3 and not f[0].startswith('#')]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('index_dir')
    ap.add_argument('fg_bed')
    ap.add_argument('out_dir')
    ap.add_argument('--k', type=int, default=10)
    ap.add_argument('--tss', action='store_true')
    ap.add_argument('--exclude', default='')
    ap.add_argument('--masked', action='store_true')
    ap.add_argument('--seed', default='9100')
    ap.add_argument('--genome', default=None, help='preloaded genome (internal)')
    a = ap.parse_args()
    run(a.index_dir, a.fg_bed, a.out_dir, a.k, a.tss, a.exclude, a.masked, [int(x) for x in a.seed.split(',')])


_CACHE = {}


def load_index(index_dir):
    if index_dir not in _CACHE:
        meta = json.load(open(os.path.join(index_dir, 'index.json')))
        z = dict(np.load(os.path.join(index_dir, 'windows.npz')))
        sine = np.array([n.startswith('SINE') for n in meta['rep_classes']])
        z['sine'] = sine[z['rep_class']]
        chrom_off = {}
        for ci, c in enumerate(meta['chroms']):
            idx = np.flatnonzero(z['chrom'] == ci)
            chrom_off[c] = (int(idx[0]), len(idx)) if len(idx) else (0, 0)
        tss = dict(np.load(os.path.join(index_dir, 'tss.npz')))
        rz = np.load(os.path.join(index_dir, 'rmsk.npz'))
        rm = {c: (rz[c + '_s'], rz[c + '_e'], rz[c + '_k']) for c in meta['chroms'] if c + '_s' in rz}
        G = mi.read_genome(meta['genome'])
        _CACHE[index_dir] = (meta, z, chrom_off, tss, rm, sine, G)
    return _CACHE[index_dir]


_KEYS = {}


def sorted_keys(index_dir, use_tss):
    """Per coarsening level: bucket keys of every full-length tile, sorted, with the tile ids (cached)."""
    if (index_dir, use_tss) not in _KEYS:
        z = load_index(index_dir)[1]
        tiles = np.flatnonzero(z['n_valid'] == W)
        out = []
        for lev in range(5):
            tk = bins(z['gc'][tiles], z['cpg_oe'][tiles], z['rep'][tiles], z['sine'][tiles], z['tss'][tiles], lev, use_tss)
            order = np.argsort(tk, kind='stable')
            out.append((tk[order], tiles[order]))
        _KEYS[(index_dir, use_tss)] = out
    return _KEYS[(index_dir, use_tss)]


def run(index_dir, fg_bed, out_dir, k, use_tss, exclude, masked, seed):
    meta, z, chrom_off, tss, rm, sine, G = load_index(index_dir)
    os.makedirs(out_dir, exist_ok=True)
    fg = read_bed(fg_bed)
    # fg covariates, computed exactly as the index's
    n = len(fg)
    gc, oe, rep, rs, td = (np.zeros(n) for _ in range(5))
    by = {}
    for i, (c, s, e) in enumerate(fg):
        by.setdefault(c, []).append(i)
    for c, idx in by.items():
        st = np.array([fg[i][1] for i in idx])
        en = np.array([fg[i][2] for i in idx])
        # per-region slices (identical values; whole-chromosome cumulative sums cost ~12 GB per worker on chr1)
        cov = [mi.covariates(G[c][a:b], np.array([0]), np.array([b - a])) for a, b in zip(st, en)]
        g_ = np.array([x[1][0] for x in cov])
        o_ = np.array([x[2][0] for x in cov])
        r_ = np.array([x[3][0] for x in cov])
        mid = (st + en) // 2
        cls = mi.class_at(rm, c, mid)
        gc[idx], oe[idx], rep[idx] = g_, o_, r_
        rs[idx] = sine[cls]
        td[idx] = mi.tss_dist(tss, c, mid)
    # usable tiles: full length, no N, not overlapping fg or excluded intervals
    ok = z['n_valid'] == W
    for c, s, e in fg + (read_bed(exclude) if exclude else []):
        if c not in chrom_off:
            continue
        off, cnt = chrom_off[c]
        a, b = max(0, s // W), min(cnt, (e - 1) // W + 1)
        ok[off + a:off + b] = False
    rng = np.random.default_rng(seed)
    used = ~ok                                   # excluded and used tiles are never drawn
    keys_by_level = sorted_keys(index_dir, use_tss)
    fk = [bins(gc, oe, rep, rs, td, lev, use_tss) for lev in range(5)]
    bg_idx, strata_fg, strata_bg, coarsen, short = [], [], [], {}, 0
    for i in range(n):
        for lev in range(5):
            sk, st = keys_by_level[lev]
            lo, hi = np.searchsorted(sk, fk[lev][i]), np.searchsorted(sk, fk[lev][i], side='right')
            cand = st[lo:hi]
            if (hi - lo) - used[cand].sum() >= 4 * k or lev == 4:
                break
        if lev:
            coarsen[lev] = coarsen.get(lev, 0) + 1
        pick = []
        tries = 0
        while len(pick) < k and tries < 50 * k and len(cand):
            j = cand[rng.integers(len(cand))]
            tries += 1
            if not used[j]:
                used[j] = True
                pick.append(j)
        short += len(pick) < k
        sid = f'L{lev}_{int(fk[lev][i])}'
        strata_fg.append(sid)
        strata_bg += [sid] * len(pick)
        bg_idx += pick
    chroms = meta['chroms']

    def seq(c, s, e):
        x = G[c][s:e].tobytes().decode()
        return ''.join('N' if ch in 'acgtn' else ch for ch in x) if masked else x.upper()
    with open(os.path.join(out_dir, 'fg.fa'), 'w') as fh:
        for i, (c, s, e) in enumerate(fg):
            fh.write(f'>f{i}\n{seq(c, s, e)}\n')
    with open(os.path.join(out_dir, 'bg.fa'), 'w') as fh, open(os.path.join(out_dir, 'bg.bed'), 'w') as fb:
        for j, t in enumerate(bg_idx):
            c = chroms[z['chrom'][t]]
            s = int(z['start'][t])
            fh.write(f'>b{j}\n{seq(c, s, s + W)}\n')
            fb.write(f'{c}\t{s}\t{s + W}\tb{j}\n')
    with open(os.path.join(out_dir, 'strata.txt'), 'w') as fh:
        fh.write('\n'.join(strata_fg + strata_bg) + '\n')
    json.dump(dict(n_fg=n, n_bg=len(bg_idx), k=k, tss=use_tss, masked=masked, seed=seed,
                   coarsened={str(kk): v for kk, v in sorted(coarsen.items())}, n_strata=len(set(strata_fg)),
                   short=int(short)),
              open(os.path.join(out_dir, 'match.json'), 'w'))
    return n, len(bg_idx), coarsen


if __name__ == '__main__':
    main()
