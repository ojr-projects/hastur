# %% c5_composition_ranking.py: composition drift vs motif ranking (THEORY.md section 5)
"""
A foreground of CpG-island-like composition (GC ~0.6, CpG not depleted) against a genome-like background (GC ~0.41,
CpG depleted), with a planted PWM motif. All contiguous 8-mers are ranked by:

  ENR      contrast z of raw counts (fg vs bg)                                   (Comsequent's most sensitive detector)
  MKV_m    double-centred contrast: log(O/E_m) in fg minus the same in bg, where E_m is the order-m Markov expectation
           fitted to each set ITSELF (from its own (m+1)-mer counts); z with Poisson SEs.  m = 1, 2, 3
           Composition captured by an order-m chain cancels in each set separately, so a fg/bg difference in
           composition (here a dinucleotide shift) cannot rank; a motif is structure of order > m+1 and survives.
Background sizes: equal to the foreground, and 10x (counting is cheap; the contrast variance falls by 1 + N_fg/N_bg).
Readouts: rank of the first motif-consistent 8-mer; share of the top 10 that are composition calls (GC >= 0.75 or
>= 2 CpG); calls with z > 5 on a motif-free pair with the same composition shift (false discoveries).
Run: ~/micromamba/envs/nanomotif/bin/python theory/checks/c5_composition_ranking.py
"""
import os

import numpy as np
import pandas as pd

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
LSEQ, L = 100, 8


def chain(freq, cg_ratio):
    """Order-1 transition matrix with stationary-ish base frequencies and P(G|C) scaled by cg_ratio."""
    P = np.tile(freq, (4, 1))
    P[1, 2] *= cg_ratio
    return P / P.sum(axis=1, keepdims=True)


BG = chain(np.array([0.295, 0.205, 0.205, 0.295]), 0.25)      # genome-like: CpG depleted
FG = chain(np.array([0.19, 0.31, 0.31, 0.19]), 0.95)          # CpG-island-like


def sample_chain(P, n, rng):
    X = np.empty((n, LSEQ), dtype=np.int64)
    X[:, 0] = rng.integers(0, 4, n)
    cdf = np.cumsum(P, axis=1)
    for t in range(1, LSEQ):
        u = rng.random(n)
        X[:, t] = (u[:, None] > cdf[X[:, t - 1]]).sum(axis=1)
    return X


def plant(X, pwm, pi, rng):
    X = X.copy()
    cdf = np.cumsum(pwm, axis=1)
    for i in np.flatnonzero(rng.random(len(X)) < pi):
        site = (rng.random((len(pwm), 1)) > cdf).sum(axis=1)
        if rng.random() < 0.5:
            site = 3 - site[::-1]
        t = rng.integers(0, LSEQ - len(pwm) + 1)
        X[i, t:t + len(pwm)] = site
    return X


def kcount(X, k):
    T = X.shape[1] - k + 1
    out = np.zeros(4 ** k)
    for Y in (X, 3 - X[:, ::-1]):
        c = np.zeros((len(Y), T), dtype=np.int64)
        for o in range(k):
            c = c * 4 + Y[:, o:o + T]
        out += np.bincount(c.ravel(), minlength=4 ** k)
    return out


def markov_expect(X, m):
    """Order-m Markov expectation of every 8-mer from the set's own (m+1)- and m-mer counts (both strands)."""
    cm1, cm = kcount(X, m + 1) + 0.5, kcount(X, m) + 0.5
    n8 = kcount(X, L).sum()
    digits = (np.arange(4 ** L)[:, None] // 4 ** np.arange(L - 1, -1, -1)) % 4
    def code(cols):
        return (digits[:, cols] * 4 ** np.arange(len(cols) - 1, -1, -1)).sum(axis=1)
    logE = np.log(cm1[code(list(range(m + 1)))] / cm1.sum())
    for i in range(m + 1, L):
        logE += np.log(cm1[code(list(range(i - m, i + 1)))]) - np.log(cm[code(list(range(i - m, i)))])
    return n8 * np.exp(logE)


def stats(fg, bg):
    F, B = kcount(fg, L), kcount(bg, L)
    r = F.sum() / B.sum()                          # fg/bg window ratio (a larger background lowers the contrast variance)
    out = {'ENR': (F - r * B) / np.sqrt(F + r ** 2 * B + 1)}
    for m in (1, 2, 3):
        Ef, Eb = markov_expect(fg, m), markov_expect(bg, m)
        d = np.log((F + 0.5) / Ef) - np.log((B + 0.5) / Eb)
        out[f'MKV_{m}'] = d / np.sqrt(1 / (F + 0.5) + 1 / (B + 0.5))
    return out


digits = (np.arange(4 ** L)[:, None] // 4 ** np.arange(L - 1, -1, -1)) % 4
gc = np.isin(digits, [1, 2]).mean(axis=1)
cpg = ((digits[:, :-1] == 1) & (digits[:, 1:] == 2)).sum(axis=1)
is_comp = (gc >= 0.75) | (cpg >= 2)
rcode = ((3 - digits)[:, ::-1] * 4 ** np.arange(L - 1, -1, -1)).sum(axis=1)
canon = np.arange(4 ** L) <= rcode


def pc_from_ic(ic):
    lo, hi = 0.25, 0.99999
    for _ in range(80):
        q = (lo + hi) / 2
        v = 2 + q * np.log2(q) + (1 - q) * np.log2((1 - q) / 3)
        lo, hi = (q, hi) if v < ic else (lo, q)
    return q


def motif(seed):
    r = np.random.default_rng(seed)
    pc = pc_from_ic(1.2)
    m = np.full((L, 4), (1 - pc) / 3)
    m[np.arange(L), r.integers(0, 4, L)] = pc
    return m


rows = []
for rep in range(5):
    rng = np.random.default_rng(rep)
    pwm = motif(100 + rep)
    ok = np.zeros(4 ** L, bool)
    for m in (pwm, pwm[::-1, ::-1]):
        ok |= (m[np.arange(L), digits] >= 0.4).all(axis=1)
    for n, pi, bgx in [(2000, 0.1, 1), (2000, 0.3, 1), (10000, 0.05, 1), (2000, 0.1, 10), (2000, 0.3, 10)]:
        fg = plant(sample_chain(FG, n, rng), pwm, pi, rng)
        bg = sample_chain(BG, n * bgx, rng)
        null_fg = sample_chain(FG, n, rng)
        st, st0 = stats(fg, bg), stats(null_fg, bg)
        for nm in st:
            z = st[nm]
            order = np.argsort(-z)
            order = order[canon[order]]
            first = np.flatnonzero(ok[order])
            rows.append(dict(rep=rep, n=n, pi=pi, bg_mult=bgx, stat=nm, rank_truth=int(first[0] + 1) if len(first) else 10 ** 6,
                             top10_comp=float(is_comp[order[:10]].mean()),
                             null_calls_z5=int((st0[nm][canon] > 5).sum())))
df = pd.DataFrame(rows)
df.to_csv(os.path.join(OUT, 'c5_composition_ranking.tsv'), sep='\t', index=False)
summ = df.groupby(['n', 'pi', 'bg_mult', 'stat']).agg(truth_top1=('rank_truth', lambda r: np.mean(r == 1)),
                                           truth_top10=('rank_truth', lambda r: np.mean(r <= 10)),
                                           median_rank=('rank_truth', 'median'),
                                           top10_composition=('top10_comp', 'mean'),
                                           motif_free_calls_z5=('null_calls_z5', 'mean'))
print(summ.round(2).to_string())
