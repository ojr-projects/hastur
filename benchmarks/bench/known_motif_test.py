# %% known_motif_test.py: known-motif test mode, quick look on a promoter screen (ROADMAP_v2 §4 item 2; C2 alternative)
"""
Each bank motif is tested as a whole PWM: per region, the best log-odds window score (both strands, uniform
background, +0.01 pseudocount); statistic = Welch t (fg vs bg) of those best scores, as a z. Family-wise control:
label-permutation max-T over the bank (N_PERM permutations of the pooled regions, max over motifs), Gumbel fit by
moments as F2, alpha 0.05 (the known family's whole budget in this quick look). The region scores are computed once,
because every unit draws fg and bg from the same region population; each unit is then index arithmetic.
Leave-own-out: optional per-unit exclusion of the module's annotated TFs (truth.json 'truths').
Usage: python bench/known_motif_test.py <units_root> <bank.meme> <out.tsv> [--loo]
"""
import json
import math
import os
import re
import sys

import numpy as np

N_PERM, ALPHA = 200, 0.05
IDX = {c: i for i, c in enumerate('ACGT')}


def read_meme(p):
    txt = open(p).read()
    return {n: np.array([[float(x) for x in l.split()] for l in body.strip().split('\n')])
            for n, body in re.findall(r'MOTIF (\S+)\nletter-probability[^\n]*\n((?:[0-9. ]+\n)+)', txt)}


def read_fa(p):
    out, name = {}, None
    for line in open(p):
        if line.startswith('>'):
            name = line[1:].strip()
        else:
            out[name] = line.strip()
    return out


def best_scores(seqs, P):
    S = np.log2((P / P.sum(axis=1, keepdims=True) + 0.01) / 0.26)
    W = len(S)
    out = np.full(len(seqs), -np.inf)
    for i, s in enumerate(seqs):
        x = np.array([IDX.get(c, -1) for c in s])
        n = len(x) - W + 1
        if n <= 0:
            continue
        best = -np.inf
        for Q in (S, S[::-1, ::-1]):
            sc = np.zeros(n)
            for j in range(W):
                col = x[j:j + n]
                sc += np.where(col >= 0, Q[j][np.clip(col, 0, 3)], 0.0)
            best = max(best, sc.max())
        out[i] = best
    return out


def welch_z(X, m):
    """z of Welch t for fg mask m vs the rest, for every motif row of X (motifs x regions)."""
    a, b = X[:, m], X[:, ~m]
    va, vb = a.var(axis=1, ddof=1) / a.shape[1], b.var(axis=1, ddof=1) / b.shape[1]
    return (a.mean(axis=1) - b.mean(axis=1)) / np.sqrt(va + vb + 1e-12)


def gumbel_thr(mx):
    mean, sd = np.mean(mx), np.std(mx, ddof=1)
    beta = sd * math.sqrt(6) / math.pi
    mu = mean - 0.5772 * beta
    return mu, beta


def main():
    root, bankp, out = sys.argv[1:4]
    loo = '--loo' in sys.argv
    bank = read_meme(bankp)
    names = sorted(bank)
    units = sorted(u for u in os.listdir(root) if os.path.exists(os.path.join(root, u, 'fg.fa')))
    u0 = os.path.join(root, units[0])
    regions = {**read_fa(os.path.join(u0, 'fg.fa')), **read_fa(os.path.join(u0, 'bg.fa'))}
    ids = sorted(regions, key=lambda r: int(r[1:]))
    pos = {r: i for i, r in enumerate(ids)}
    X = np.array([best_scores([regions[r] for r in ids], bank[n]) for n in names])
    print(f'{len(names)} motifs x {len(ids)} regions scored', flush=True)
    rows = []
    for u in units:
        d = os.path.join(root, u)
        j = json.load(open(os.path.join(d, 'truth.json')))
        fg = [pos[l[1:].strip()] for l in open(os.path.join(d, 'fg.fa')) if l.startswith('>')]
        m = np.zeros(len(ids), bool)
        m[fg] = True
        keep = np.ones(len(names), bool)
        if loo and not j.get('null'):
            own = {t['tf'] for t in j.get('truths', [])}
            keep = np.array([n not in own for n in names])
        elif loo and j.get('null'):
            own = set()
        z = welch_z(X[keep], m)
        rng = np.random.default_rng([4242, len(fg), sum(fg) % 100003])
        mx = []
        for _ in range(N_PERM):
            pm = np.zeros(len(ids), bool)
            pm[rng.choice(len(ids), len(fg), replace=False)] = True
            mx.append(welch_z(X[keep], pm).max())
        mu, beta = gumbel_thr(mx)
        kn = [n for n, k in zip(names, keep) if k]
        for n, zz in sorted(zip(kn, z), key=lambda x: -x[1])[:5]:
            p_adj = min(1.0, -math.expm1(-math.exp(-(zz - mu) / beta)))
            rows.append(dict(unit=u, null=bool(j.get('null')), module=j.get('module'), motif=n, z=round(float(zz), 3),
                             p_adj=p_adj, sig=p_adj <= ALPHA,
                             truth=','.join(t['tf'] for t in j.get('truths', [])) if not j.get('null') else ''))
    import pandas as pd
    df = pd.DataFrame(rows)
    df.to_csv(out, sep='\t', index=False)
    n, p = df[df.null], df[~df.null]
    fc = n.groupby('unit').sig.sum()
    print(f'null runs {fc.size}: false calls/run {fc.mean():.3f}, runs with any {(fc > 0).mean():.3f}')
    p = p.assign(hit=[mo in (t.split(',') if t else []) for mo, t in zip(p.motif, p.truth)])
    ev = p[p.truth != ''].unit.unique()
    h = p[p.sig & p.hit].unit.unique()
    print(f'modules with an annotated bank TF: {len(ev)}; own regulator called (sig): {len(h)}')
    print(f'native calls on all modules {int(p.sig.sum())}; sig calls not the annotated TF {int((p.sig & ~p.hit).sum())}')


if __name__ == '__main__':
    main()
