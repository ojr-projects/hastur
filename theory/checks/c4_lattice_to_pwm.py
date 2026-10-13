# %% c4_lattice_to_pwm.py: can lattice detection + all-data EM represent soft motifs at small N? (THEORY.md sections 3, 6)
"""
Mechanism check for the reconciled pipeline on PWM-class simulations (iid background, 100 bp, both strands, control
set of equal size, one site per carrying sequence as in Comsequent's shared_sim):

  ranking   rank of the first motif-consistent call among all canonical patterns for
              EXACT_k5   Comsequent-style contrast z over every gapped 5-of-13 pattern (495 shapes x 1024)
              EXACT_L8   contrast z over all contiguous 8-mers
              SOFT_L8    the same 8-mer table smoothed by a per-position letter kernel (match 1, mismatch b) by mode
                         products: a soft count for every consensus at once (THEORY.md 3.5)
  represent PWM similarity to the truth (mean over truth columns of 1 - total variation, best offset/strand) for
              LAT+EM     ZOOPS EM on ALL fg sequences, started from each of the top 10 EXACT_k5 seeds (anchored
                         frame counts, U1-style), best final model chosen by fg-vs-control log-likelihood gain
              LAT8+EM    the same EM started from the top 10 seeds of SOFT_L8 (the full-length shape table)
              *_seed1    EM from the single top-ranked seed: the motif is chosen by the calibrated detection statistic
                         and EM supplies only the representation (THEORY.md 7: decouple detection from estimation)
              ORACLE+EM  the same EM started from the true PWM (ceiling for this EM)
              STREME     STREME 5.5.9 top motif, fg vs the same control
            (each also as best of its top 3: STREME's rank 1 is sometimes a non-significant junk motif)
Run: ~/micromamba/envs/nanomotif/bin/python theory/checks/c4_lattice_to_pwm.py [reps]
"""
import itertools
import os
import re
import subprocess
import sys
import tempfile
import time

import numpy as np
import pandas as pd
from numpy.lib.stride_tricks import sliding_window_view

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
STREME = os.path.expanduser('~/micromamba/envs/meme/bin/streme')
MOTIFS = os.path.expanduser('~/micromamba/envs/homer/share/homer/data/knownTFs/vertebrates/known.motifs')
LSEQ = 100
Q = np.full(4, 0.25)
REPS = int(sys.argv[1]) if len(sys.argv) > 1 else 5
IC_LEVELS = {'low': 0.7, 'mid': 1.2, 'high': 1.7}


# %% Motifs and simulation
def pc_from_ic(ic):
    lo, hi = 0.25, 0.99999
    for _ in range(80):
        q = (lo + hi) / 2
        v = 2 + q * np.log2(q) + (1 - q) * np.log2((1 - q) / 3)
        lo, hi = (q, hi) if v < ic else (lo, q)
    return q


def sim_pwm(L, level, seed=0):
    r = np.random.default_rng(seed)
    pc = pc_from_ic(IC_LEVELS[level])
    m = np.full((L, 4), (1 - pc) / 3)
    m[np.arange(L), r.integers(0, 4, L)] = pc
    return m


def homer_pwm(prefix):
    name, rows = None, []
    for line in open(MOTIFS):
        if line.startswith('>'):
            if name and name.startswith(prefix):
                break
            name, rows = line.split('\t')[1], []
        elif line.strip():
            rows.append([float(v) for v in line.split()])
    m = np.clip(np.array(rows), 1e-3, None)
    return m / m.sum(axis=1, keepdims=True)


def simulate(n, pi, pwm, seed):
    r = np.random.default_rng(seed)
    fg = r.integers(0, 4, (n, LSEQ), dtype=np.int64)
    bg = r.integers(0, 4, (n, LSEQ), dtype=np.int64)
    L = len(pwm)
    cdf = np.cumsum(pwm, axis=1)
    for i in np.flatnonzero(r.random(n) < pi):
        site = (r.random((L, 1)) > cdf).sum(axis=1)
        if r.random() < 0.5:
            site = 3 - site[::-1]
        t = r.integers(0, LSEQ - L + 1)
        fg[i, t:t + L] = site
    return fg, bg


def rc(X):
    return 3 - X[:, ::-1]


# %% Lattice counting (both strands)
def codes(X, offs):
    span = offs[-1] + 1
    T = X.shape[1] - span + 1
    c = np.zeros((X.shape[0], T), dtype=np.int64)
    for o in offs:
        c = c * 4 + X[:, o:o + T]
    return c


def count(X, offs):
    k = len(offs)
    return (np.bincount(codes(X, offs).ravel(), minlength=4 ** k)
            + np.bincount(codes(rc(X), offs).ravel(), minlength=4 ** k)).astype(float)


def lattice_k5(fg, bg, k=5, w=13):
    """Contrast z for every gapped k-of-w pattern; returns list of (z, offsets, letters) for canonical patterns."""
    out = []
    for rest in itertools.combinations(range(1, w), k - 1):
        offs = (0,) + rest
        F, Bc = count(fg, offs), count(bg, offs)
        z = (F - Bc) / np.sqrt(F + Bc + 1)
        span = offs[-1]
        roffs = tuple(sorted(span - o for o in offs))
        digits = (np.arange(4 ** k)[:, None] // 4 ** np.arange(k - 1, -1, -1)) % 4      # letters per offset
        # reverse complement pattern: letter at offset o becomes 3-letter at span-o
        rdig = (3 - digits)[:, ::-1]
        rcode = (rdig * 4 ** np.arange(k - 1, -1, -1)).sum(axis=1)
        keep = (roffs > offs) | ((roffs == offs) & (rcode >= np.arange(4 ** k)))
        idx = np.flatnonzero(keep & (z > 2))
        out += [(z[i], offs, digits[i]) for i in idx]
    out.sort(key=lambda t: -t[0])
    return out


def consistent(offs, letters, pwm, thr=0.4):
    """Every specified (offset, letter) lies inside the motif at a column giving the letter p >= thr (either strand)."""
    L = len(pwm)
    for m in (pwm, pwm[::-1, ::-1]):
        for shift in range(-offs[-1], L):
            if all(0 <= o + shift < L and m[o + shift, b] >= thr for o, b in zip(offs, letters)):
                return True
    return False


def mode_smooth(Tn, K, L):
    T = Tn.reshape((4,) * L)
    for j in range(L):
        T = np.moveaxis(np.tensordot(K, T, axes=([1], [j])), 0, j)
    return T.ravel()


def lmer_ranks(fg, bg, pwm, L=8, b_soft=(0.05, 0.15)):
    """Rank of the best motif-consistent consensus L-mer under exact and soft-kernel contrast z."""
    offs = tuple(range(L))
    F, Bc = count(fg, offs), count(bg, offs)
    digits = (np.arange(4 ** L)[:, None] // 4 ** np.arange(L - 1, -1, -1)) % 4
    # consistent L-mers: every letter has p >= 0.4 in the motif at some alignment (motif may be longer than L)
    ok = np.zeros(4 ** L, bool)
    Lm = len(pwm)
    for m in (pwm, pwm[::-1, ::-1]):
        for shift in range(0, max(Lm - L, 0) + 1):
            sub = m[shift:shift + L]
            if len(sub) < L:
                continue
            ok |= (sub[np.arange(L), digits] >= 0.4).all(axis=1)
    res = {}
    stats = {'EXACT_L8': (F - Bc) / np.sqrt(F + Bc + 1)}
    for b in b_soft:
        K = (1 - b) * np.eye(4) + b * np.ones((4, 4))
        K2 = K ** 2
        TF, TB = mode_smooth(F, K, L), mode_smooth(Bc, K, L)
        V = mode_smooth(F, K2, L) + mode_smooth(Bc, K2, L)
        stats[f'SOFT_L8_b{b}'] = (TF - TB) / np.sqrt(V + 1)
    rcode = ((3 - digits)[:, ::-1] * 4 ** np.arange(L - 1, -1, -1)).sum(axis=1)
    canon = np.arange(4 ** L) <= rcode
    seeds = None
    for nm, z in stats.items():
        order = np.argsort(-z)
        order = order[canon[order]]           # both strands are counted: keep one of each reverse-complement pair
        first = np.flatnonzero(ok[order])
        res[nm] = int(first[0] + 1) if len(first) else 10 ** 6
        if nm == f'SOFT_L8_b{b_soft[0]}':
            seeds = [digits[i] for i in order[:10]]
    return res, seeds


# %% ZOOPS EM on all sequences (both strands)
def windows(X, W):
    v = sliding_window_view(X, W, axis=1)                       # N x T x W
    return np.concatenate([v, 3 - v[:, :, ::-1]], axis=1)     # N x 2T x W


def zoops_ll(Xw, pwm, gamma):
    lp = np.log(pwm) - np.log(Q)[None, :]
    S = lp[np.arange(pwm.shape[0])[None, None, :], Xw].sum(axis=2)    # N x 2T log LR
    n2t = S.shape[1]
    lse = np.logaddexp.reduce(S, axis=1)
    ll = np.logaddexp(np.log(1 - gamma), np.log(gamma / n2t) + lse)
    return S, ll


def em(X, init, iters=60, gamma=0.3, pseudo=0.25):
    W = len(init)
    Xw = windows(X, W)
    pwm = init.copy()
    for _ in range(iters):
        S, ll = zoops_ll(Xw, pwm, gamma)
        n2t = S.shape[1]
        post = np.exp(np.log(gamma / n2t) + S - ll[:, None])        # N x 2T
        cnt = np.zeros((W, 4))
        for j in range(W):
            cnt[j] = np.bincount(Xw[:, :, j].ravel(), weights=post.ravel(), minlength=4)
        new = (cnt + pseudo * Q) / (cnt.sum(axis=1, keepdims=True) + pseudo)
        gamma = float(np.clip(post.sum() / len(X), 0.02, 0.98))
        if np.abs(new - pwm).max() < 1e-4:
            pwm = new
            break
        pwm = new
    return pwm, gamma


def trim(pwm, min_ic=0.3):
    ic = 2 + (pwm * np.log2(pwm)).sum(axis=1)
    keep = np.flatnonzero(ic >= min_ic)
    return pwm[keep[0]:keep[-1] + 1] if len(keep) else pwm


def anchored_init(fg, bg, offs, letters, flank=4):
    """U1-style frame PWM: column counts in a frame around every occurrence of the seed, fg minus bg."""
    span = offs[-1] + 1
    W = span + 2 * flank
    cols = []
    for X in (fg, bg):
        Xw = windows(X, W)                                    # N x 2T x W (both strands)
        m = np.ones(Xw.shape[:2], bool)
        for o, b in zip(offs, letters):
            m &= Xw[:, :, flank + o] == b
        cnt = np.stack([np.bincount(Xw[:, :, j][m], minlength=4) for j in range(W)])
        cols.append(cnt.astype(float))
    ex = np.clip(cols[0] - cols[1], 0, None) + 0.5
    return ex / ex.sum(axis=1, keepdims=True)


def contrast_gain(fg, bg, pwm, gamma):
    _, l1 = zoops_ll(windows(fg, len(pwm)), pwm, gamma)
    _, l0 = zoops_ll(windows(bg, len(pwm)), pwm, gamma)
    return l1.sum() - l0.sum()


def similarity(learned, truth):
    """Mean over truth columns of 1 - TV distance to the aligned learned column (uniform where not covered)."""
    L, W = len(truth), len(learned)
    best = 0.0
    for m in (learned, learned[::-1, ::-1]):
        for off in range(-W + 1, L):
            s = 0.0
            for j in range(L):
                i = j - off
                col = m[i] if 0 <= i < W else Q
                s += 1 - 0.5 * np.abs(col - truth[j]).sum()
            best = max(best, s / L)
    return best


def null_similarity(truth):
    return np.mean([1 - 0.5 * np.abs(Q - c).sum() for c in truth])


# %% STREME
def write_fa(path, X):
    with open(path, 'w') as f:
        for i, s in enumerate(X):
            f.write(f'>s{i}\n' + ''.join('ACGT'[c] for c in s) + '\n')


def run_streme(fg, bg):
    with tempfile.TemporaryDirectory() as d:
        write_fa(f'{d}/fg.fa', fg)
        write_fa(f'{d}/bg.fa', bg)
        t0 = time.time()
        out = subprocess.run([STREME, '--p', f'{d}/fg.fa', '--n', f'{d}/bg.fa', '--dna', '--nmotifs', '3', '--text'],
                             capture_output=True, text=True, errors='replace', timeout=1800).stdout
        dt = time.time() - t0
    ms = re.findall(r'letter-probability matrix:[^\n]*\n((?:\s*[0-9.]+\s+[0-9.]+\s+[0-9.]+\s+[0-9.]+\s*\n)+)', out)
    pwms = [np.array([[float(v) for v in line.split()] for line in m.strip().split('\n')]) for m in ms]
    return pwms, dt


# %% Run
CONDS = [('pwm_L8_mid', 500, 0.1), ('pwm_L8_mid', 500, 0.3), ('pwm_L8_mid', 2000, 0.05), ('pwm_L8_mid', 2000, 0.1),
         ('pwm_L12_low', 2000, 0.1), ('pwm_L12_low', 2000, 0.3),
         ('NRF1', 2000, 0.1), ('REST', 2000, 0.1), ('CTCF(', 2000, 0.1), ('Sox2', 2000, 0.1), ('p53(', 2000, 0.1)]
if os.environ.get('C4_ONLY'):          # e.g. C4_ONLY=3 runs one condition (timing / debugging)
    CONDS = [CONDS[int(i)] for i in os.environ['C4_ONLY'].split(',')]
TAG = os.environ.get('C4_TAG', '')
rows = []
for name, n, pi in CONDS:
    pwm = sim_pwm(int(name[5:].split('_')[0]), name.split('_')[-1]) if name.startswith('pwm_') else homer_pwm(name)
    for rep in range(REPS):
        fg, bg = simulate(n, pi, pwm, seed=1000 * rep + n + int(pi * 1000))
        row = dict(motif=name, n=n, pi=pi, rep=rep, L=len(pwm), null_sim=null_similarity(pwm))
        t0 = time.time()
        lat = lattice_k5(fg, bg)
        row['t_lattice'] = time.time() - t0
        cons_ranks = [i for i, (z, o, l) in enumerate(lat[:5000]) if consistent(o, l, pwm)]
        row['rank_EXACT_k5'] = cons_ranks[0] + 1 if cons_ranks else 10 ** 6
        ranks8, seeds8 = lmer_ranks(fg, bg, pwm)
        row.update({f'rank_{k}': v for k, v in ranks8.items()})
        t0 = time.time()
        fits = []
        for z, o, l in lat[:10]:
            init = anchored_init(fg, bg, o, l)
            p, g = em(fg, init)
            p = trim(p)
            fits.append((contrast_gain(fg, bg, p, g), p))
        row['sim_LAT_EM_seed1'] = similarity(fits[0][1], pwm)      # EM from the top-ranked seed (selection by detection z)
        fits.sort(key=lambda t: -t[0])
        row['t_lat_em'] = time.time() - t0
        row['sim_LAT_EM'] = similarity(fits[0][1], pwm)
        row['sim_LAT_EM_best3'] = max(similarity(p, pwm) for _, p in fits[:3])
        row['W_LAT_EM'] = len(fits[0][1])
        t0 = time.time()
        fits8 = []
        for d in seeds8:
            p, g = em(fg, anchored_init(fg, bg, tuple(range(8)), d))
            p = trim(p)
            fits8.append((contrast_gain(fg, bg, p, g), p))
        row['sim_LAT8_EM_seed1'] = similarity(fits8[0][1], pwm)
        fits8.sort(key=lambda t: -t[0])
        row['t_lat8_em'] = time.time() - t0
        row['sim_LAT8_EM'] = similarity(fits8[0][1], pwm)
        row['sim_LAT8_EM_best3'] = max(similarity(p, pwm) for _, p in fits8[:3])
        po, go = em(fg, np.vstack([np.full((2, 4), .25), pwm, np.full((2, 4), .25)]))
        row['sim_ORACLE_EM'] = similarity(trim(po), pwm)
        sp, dt = run_streme(fg, bg)
        row['sim_STREME'] = similarity(sp[0], pwm) if sp else np.nan
        row['sim_STREME_best3'] = max(similarity(p, pwm) for p in sp[:3]) if sp else np.nan
        row['t_streme'] = dt
        rows.append(row)
        print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in row.items()}, flush=True)
df = pd.DataFrame(rows)
df.to_csv(os.path.join(OUT, f'c4_lattice_to_pwm{TAG}.tsv'), sep='\t', index=False)

# %% Summary
g = df.groupby(['motif', 'n', 'pi'], sort=False)
summ = g.agg(L=('L', 'first'), null=('null_sim', 'first'),
             k5_top10=('rank_EXACT_k5', lambda r: np.mean(r <= 10)),
             L8_top10=('rank_EXACT_L8', lambda r: np.mean(r <= 10)),
             soft05_top10=('rank_SOFT_L8_b0.05', lambda r: np.mean(r <= 10)),
             soft15_top10=('rank_SOFT_L8_b0.15', lambda r: np.mean(r <= 10)),
             LAT_EM=('sim_LAT_EM', 'mean'), LAT_EM3=('sim_LAT_EM_best3', 'mean'),
             LAT8_EM=('sim_LAT8_EM', 'mean'), LAT8_EM3=('sim_LAT8_EM_best3', 'mean'),
             LAT_EM_s1=('sim_LAT_EM_seed1', 'mean'), LAT8_EM_s1=('sim_LAT8_EM_seed1', 'mean'), ORACLE_EM=('sim_ORACLE_EM', 'mean'),
             STREME=('sim_STREME', 'mean'), STREME3=('sim_STREME_best3', 'mean'),
             t_lat=('t_lattice', 'median'), t_em=('t_lat_em', 'median'), t_streme=('t_streme', 'median'))
print(summ.round(2).to_string())
