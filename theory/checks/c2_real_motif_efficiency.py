# %% c2_real_motif_efficiency.py: where the dilution gap comes from, on 472 real vertebrate TF motifs (THEORY.md section 3)
"""
Oracle detection rates for nested detector classes, computed exactly from each motif's columns.

A detector is a window set A; with N windows, a per-window site rate pi, site mass P(A) and background mass Q(A),
the count z-score is z^2 = N * pi^2 P(A)^2 / (pi P(A) + Q(A)) (Poisson, known background). The information rate
r = max over the class of pi^2 P^2 / (pi P + Q) gives the windows needed for z = 5: N = 25 / r.

Classes (each contains the one before it, so the rate ratios decompose the gap):
  EXACT_k  conjunction of the consensus letter at <= k positions            (Comsequent lattice cell)
  IUPAC_k  conjunction of letter sets (top 1/2/3 letters by p/q) at <= k     (degenerate lattice, U3)
  SOFT_k   likelihood-ratio threshold on the k most informative positions   (rank-1 weighting of one shape table)
  SOFT_L   likelihood-ratio threshold on all columns                        (oracle PWM scan; Neyman-Pearson bound)

Gap factors: letter = EXACT/IUPAC, weighting = IUPAC/SOFT_k, truncation = SOFT_k/SOFT_L (all as data multiples).
Run: ~/micromamba/envs/nanomotif/bin/python theory/checks/c2_real_motif_efficiency.py   (C2_LIB=jaspar for JASPAR)
"""
import os
import re

import numpy as np
import pandas as pd

OUT = os.path.join(os.path.dirname(__file__), '..', 'results')
LIB = os.environ.get('C2_LIB', 'homer')     # 'homer': 472 HOMER ChIP-derived matrices; 'jaspar': 1,293 JASPAR matrices
MOTIFS = os.path.expanduser('~/micromamba/envs/homer/share/homer/data/knownTFs/vertebrates/' +
                            ('known.motifs' if LIB == 'homer' else 'all.motifs'))
Q = np.array([0.295, 0.205, 0.205, 0.295])        # human-like background (GC 0.41)
PIS = [1e-2, 1e-3, 1e-4]                           # per-window site rates: SELEX-like, ChIP-like, sparse
BIN = 0.02                                         # log-scale bin width (nats) for the DP and the score grid


def read_homer(path):
    out, name, rows = [], None, []
    for line in open(path):
        if line.startswith('>'):
            if name:
                out.append((name, np.array(rows)))
            f = line[1:].rstrip('\n').split('\t')
            name, rows = f[1], []
        elif line.strip():
            rows.append([float(v) for v in line.split()])
    out.append((name, np.array(rows)))
    res = []
    for nm, m in out:
        m = np.clip(m, 1e-3, None)
        res.append((nm, m / m.sum(axis=1, keepdims=True)))
    return res


def rate(P, Qm, pi):
    """max over frontier points of pi^2 P^2 / (pi P + Q)."""
    return float(np.max(pi ** 2 * P ** 2 / (pi * P + Qm)))


def conj_frontier(options, k):
    """DP over positions: each position offers (skip) or one of its options (logp, logq). Returns arrays of (P, Q) for
    the best logP per (#specified <= k, logQ bin)."""
    nb = int(np.ceil(k * 8 / BIN)) + 2
    best = np.full((k + 1, nb), -np.inf)
    best[0, 0] = 0.0
    for opts in options:
        new = best.copy()
        for lp, lq in opts:
            sh = int(round(-lq / BIN))
            if sh >= nb:
                continue
            cand = np.full_like(best, -np.inf)
            cand[1:, sh:] = best[:-1, :nb - sh] + lp
            new = np.maximum(new, cand)
        best = new
    lqs = -np.arange(nb) * BIN
    m = np.isfinite(best[1:])
    P = np.exp(best[1:][m])
    Qv = np.exp(np.broadcast_to(lqs, best[1:].shape)[m])
    return P, Qv


def soft_frontier(cols):
    """Exact score distributions of the log-likelihood-ratio sum over the given columns, under the site model and
    the background; returns tail masses P(s >= t), Q(s >= t) over all thresholds t (the Neyman-Pearson frontier)."""
    lr = np.log(cols / Q)                                # L x 4
    idx = np.round(lr / BIN).astype(int)
    lo = idx.min(axis=1).sum()
    hi = idx.max(axis=1).sum()
    n = hi - lo + 1
    dp_p = np.zeros(n)
    dp_q = np.zeros(n)
    off = 0
    dp_p[0] = dp_q[0] = 1.0
    width = 1
    for j in range(len(cols)):
        mn = idx[j].min()
        new_w = width + idx[j].max() - mn
        np_, nq = np.zeros(new_w), np.zeros(new_w)
        for b in range(4):
            s = idx[j, b] - mn
            np_[s:s + width] += dp_p[:width] * cols[j, b]
            nq[s:s + width] += dp_q[:width] * Q[b]
        dp_p[:new_w], dp_q[:new_w] = np_, nq
        width = new_w
    Pt = np.cumsum(dp_p[:width][::-1])[::-1]
    Qt = np.cumsum(dp_q[:width][::-1])[::-1]
    return Pt, Qt


def column_stats(c):
    ratio = c / Q
    chi2 = float((c ** 2 / Q).sum() - 1)
    kl = float((c * np.log(ratio)).sum())
    top = np.argsort(-ratio)
    kappa1 = float(c[top[0]] ** 2 / Q[top[0]] / (chi2 + 1))
    ic = float(2 + (c * np.log2(c)).sum())
    return chi2, kl, kappa1, ic, float(c.max())


def analyse(name, m, ks=(5, 7)):
    L = len(m)
    stats = np.array([column_stats(c) for c in m])
    ratio = m / Q
    order = np.argsort(-ratio, axis=1)
    exact_opts = [[(np.log(m[j, order[j, 0]]), np.log(Q[order[j, 0]]))] for j in range(L)]
    iupac_opts = []
    for j in range(L):
        o = []
        for s in (1, 2, 3):
            B = order[j, :s]
            o.append((np.log(m[j, B].sum()), np.log(Q[B].sum())))
        iupac_opts.append(o)
    row = dict(motif=name, L=L, IC=stats[:, 3].sum(), mean_chi2=stats[:, 0].mean(), mean_kappa1=stats[:, 2].mean(),
               frac_soft_cols=float((stats[:, 4] < 0.6).mean()))
    Ps, Qs = soft_frontier(m)
    by_info = np.argsort(-stats[:, 1])
    fr = {}
    for k in ks:
        fr[f'EXACT_{k}'] = conj_frontier(exact_opts, k)
        fr[f'IUPAC_{k}'] = conj_frontier(iupac_opts, k)
        fr[f'SOFT_{k}'] = soft_frontier(m[np.sort(by_info[:min(k, L)])])
    fr['SOFT_L'] = (Ps, Qs)
    for pi in PIS:
        for key, (P, Qv) in fr.items():
            row[f'{key}@{pi:g}'] = 25 / rate(P, Qv, pi)          # windows needed for z = 5
    return row


# %% Real motifs + crisp references
motifs = read_homer(MOTIFS)
if LIB == 'jaspar':
    motifs = [(nm, m) for nm, m in motifs if nm.endswith('/Jaspar')]
crisp = {'GATC (Dam)': 'GATC', 'GANTC (CcrM)': 'GANTC', 'CCWGG (Dcm)': 'CCWGG', 'EcoKI': 'AACNNNNNNGTGC',
         'GCNNNANNNNGC': 'GCNNNANNNNGC'}
iupac = {'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'W': 'AT', 'N': 'ACGT'}
for nm, pat in crisp.items():
    m = np.array([[1.0 if b in iupac[ch] else 1e-4 for b in 'ACGT'] for ch in pat])
    motifs.append(('CRISP:' + nm, m / m.sum(axis=1, keepdims=True)))

rows = [analyse(nm, m) for nm, m in motifs]
df = pd.DataFrame(rows)
df.to_csv(os.path.join(OUT, 'c2_real_motif_efficiency' + ('' if LIB == 'homer' else '_' + LIB) + '.tsv'), sep='\t', index=False)

# %% Summaries
real = df[~df.motif.str.startswith('CRISP')]
print(f'{len(real)} real motifs; median L {real.L.median():.0f}, median total IC {real.IC.median():.1f} bits, '
      f'median per-column kappa1 (single-letter share of chi2+1) {real.mean_kappa1.median():.2f}')
for pi in PIS:
    s = lambda a, b: (real[f'{a}@{pi:g}'] / real[f'{b}@{pi:g}'])
    print(f'\n== per-window site rate pi = {pi:g}: data multiple vs oracle PWM (median [IQR]; share > 2x; share > 10x)')
    for k in (5, 7):
        for a in (f'EXACT_{k}', f'IUPAC_{k}', f'SOFT_{k}'):
            r = s(a, 'SOFT_L')
            print(f'  {a:9s}: {r.median():6.2f} [{r.quantile(.25):.2f}-{r.quantile(.75):.2f}]  '
                  f'>2x {np.mean(r > 2):.2f}  >10x {np.mean(r > 10):.2f}')
        lt, wt, tr = s(f'EXACT_{k}', f'IUPAC_{k}'), s(f'IUPAC_{k}', f'SOFT_{k}'), s(f'SOFT_{k}', 'SOFT_L')
        print(f'  k={k} decomposition (median): letter {lt.median():.2f}x, weighting {wt.median():.2f}x, '
              f'truncation {tr.median():.2f}x; geometric-mean shares of log gap: '
              f'{np.log(lt).mean() / np.log(s(f"EXACT_{k}", "SOFT_L")).mean():.2f} / '
              f'{np.log(wt).mean() / np.log(s(f"EXACT_{k}", "SOFT_L")).mean():.2f} / '
              f'{np.log(tr).mean() / np.log(s(f"EXACT_{k}", "SOFT_L")).mean():.2f}')
    c = df[df.motif.str.startswith('CRISP')]
    print('  crisp references EXACT_7/SOFT_L:', ', '.join(f'{m.split(":")[1]} {v:.2f}' for m, v in
                                                     zip(c.motif, c[f'EXACT_7@{pi:g}'] / c[f'SOFT_L@{pi:g}'])))

# Named motifs Comsequent missed on ENCODE ChIP (P5) vs ones it found
print('\n== named motifs, pi = 1e-3: EXACT_7 / IUPAC_7 / SOFT_7 data multiples vs oracle PWM')
for pat in ['REST', 'Sp1', 'ZNF143', 'NRF1', 'CTCF\\(', 'GABPA', 'YY1', 'Max', 'USF1', 'Elk1', 'Atf1', 'FOXA1',
            'GATA3', 'CEBPB', 'NFY', 'HNF4a', 'AP-1', 'TEAD', 'Oct4', 'Sox2', 'PU.1\\(', 'p53\\(', 'NFkB-p65\\(']:
    d = real[real.motif.str.contains('^' + pat, regex=True)]
    if len(d):
        r = d.iloc[0]
        print(f"  {r.motif[:40]:40s} L {r.L:2d} IC {r.IC:5.1f}: "
              f"{r['EXACT_7@0.001'] / r['SOFT_L@0.001']:6.1f} / {r['IUPAC_7@0.001'] / r['SOFT_L@0.001']:5.1f} / "
              f"{r['SOFT_7@0.001'] / r['SOFT_L@0.001']:5.1f}")

# Feature correlations of the EXACT_7 gap
g = np.log(real['EXACT_7@0.001'] / real['SOFT_L@0.001'])
for f in ['L', 'IC', 'mean_chi2', 'mean_kappa1', 'frac_soft_cols']:
    print(f'  Spearman(log EXACT_7 gap, {f}) = {pd.Series(g).corr(real[f], method="spearman"):+.2f}')
