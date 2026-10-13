# %% sim_make.py: the simulation dev panel (ROADMAP D0)
"""
Writes <root>/<dataset>/{fg.fa, bg.fa, truth.json} for every cell x replicate of the panel. 100-bp sequences, one
site per carrying sequence, random strand and position. Backgrounds: iid uniform, or a CpG-drift pair (fg
CpG-island-like order-1 chain, bg genome-like CpG-depleted chain). Background size BG_MULT x foreground for every tool.
Usage: python bench/sim_make.py <root> [reps]
"""
import json
import os
import sys

import numpy as np

HOMER_LIB = os.path.expanduser('~/micromamba/envs/homer/share/homer/data/knownTFs/vertebrates/known.motifs')
LSEQ, BG_MULT = 100, 5
IC = {'low': 0.7, 'mid': 1.2, 'high': 1.7}
CRISP = {'GCNNNANNNNGC': 'GCNNNANNNNGC', 'EcoKI': 'AACNNNNNNGTGC', 'GATC': 'GATC'}
# cell: (motif, N, pi, background)
PANEL = [('pwm_L8_mid', 500, 0.1, 'iid'), ('pwm_L8_mid', 2000, 0.05, 'iid'), ('pwm_L8_mid', 2000, 0.1, 'iid'),
         ('pwm_L12_low', 2000, 0.1, 'iid'), ('REST', 2000, 0.1, 'iid'), ('NRF1', 2000, 0.1, 'iid'),
         ('p53(', 2000, 0.1, 'iid'), ('GCNNNANNNNGC', 2000, 0.05, 'iid'), ('EcoKI', 2000, 0.05, 'iid'),
         ('GATC', 2000, 0.05, 'iid'), ('pwm_L8_mid', 2000, 0.1, 'cpg')]


def pc_from_ic(ic):
    lo, hi = 0.25, 0.99999
    for _ in range(80):
        q = (lo + hi) / 2
        v = 2 + q * np.log2(q) + (1 - q) * np.log2((1 - q) / 3)
        lo, hi = (q, hi) if v < ic else (lo, q)
    return q


def motif_pwm(name, rep):
    if name.startswith('pwm_L'):
        L, lvl = name[5:].split('_')
        r = np.random.default_rng(1000 + rep)
        pc = pc_from_ic(IC[lvl])
        m = np.full((int(L), 4), (1 - pc) / 3)
        m[np.arange(int(L)), r.integers(0, 4, int(L))] = pc
        return m
    if name in CRISP:
        iu = {'A': [1, 0, 0, 0], 'C': [0, 1, 0, 0], 'G': [0, 0, 1, 0], 'T': [0, 0, 0, 1], 'N': [1, 1, 1, 1]}
        m = np.array([iu[c] for c in CRISP[name]], float)
        return m / m.sum(axis=1, keepdims=True)
    name_, rows = None, []
    for line in open(HOMER_LIB):
        if line.startswith('>'):
            if name_ and name_.startswith(name):
                break
            name_, rows = line.split('\t')[1], []
        elif line.strip():
            rows.append([float(v) for v in line.split()])
    m = np.clip(np.array(rows), 1e-3, None)
    return m / m.sum(axis=1, keepdims=True)


def consensus(m):
    if (m.max(axis=1) < 0.3).all():
        return 'N'
    rev = {'AC': 'M', 'AG': 'R', 'AT': 'W', 'CG': 'S', 'CT': 'Y', 'GT': 'K'}
    out = []
    for col in m:
        o = np.argsort(-col)
        if col[o[0]] >= 0.6:
            out.append('ACGT'[o[0]])
        elif col[o[0]] + col[o[1]] >= 0.8:
            out.append(rev[''.join(sorted('ACGT'[o[0]] + 'ACGT'[o[1]]))])
        else:
            out.append('N')
    return ''.join(out).strip('N')


def chain(freq, cg):
    P = np.tile(freq, (4, 1))
    P[1, 2] *= cg
    return P / P.sum(axis=1, keepdims=True)


def sample(kind, n, rng, fg_side):
    if kind == 'iid':
        return rng.integers(0, 4, (n, LSEQ))
    P = chain(np.array([0.19, 0.31, 0.31, 0.19]), 0.95) if fg_side else chain(np.array([0.295, 0.205, 0.205, 0.295]), 0.25)
    X = np.empty((n, LSEQ), np.int64)
    X[:, 0] = rng.integers(0, 4, n)
    cdf = np.cumsum(P, axis=1)
    for t in range(1, LSEQ):
        X[:, t] = (rng.random(n)[:, None] > cdf[X[:, t - 1]]).sum(axis=1)
    return X


def plant(X, pwm, pi, rng):
    X = X.copy()
    L = len(pwm)
    cdf = np.cumsum(pwm, axis=1)
    for i in np.flatnonzero(rng.random(len(X)) < pi):
        site = (rng.random((L, 1)) > cdf).sum(axis=1)
        if rng.random() < 0.5:
            site = 3 - site[::-1]
        t = rng.integers(0, LSEQ - L + 1)
        X[i, t:t + L] = site
    return X


def write_fa(path, X):
    with open(path, 'w') as fh:
        for i, s in enumerate(X):
            fh.write(f'>s{i}\n' + ''.join('ACGT'[c] for c in s) + '\n')


if __name__ == '__main__':
    root = sys.argv[1]
    reps = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    for name, n, pi, kind in PANEL:
        for rep in range(reps):
            did = f"{name.rstrip('(')}__n{n}__pi{pi}__{kind}__r{rep}"
            d = os.path.join(root, did)
            if os.path.exists(os.path.join(d, 'truth.json')):
                continue
            os.makedirs(d, exist_ok=True)
            rng = np.random.default_rng(7919 * rep + n + int(pi * 1e4) + len(name) * 31 + len(kind) * 101)
            pwm = motif_pwm(name, rep)
            write_fa(os.path.join(d, 'fg.fa'), plant(sample(kind, n, rng, True), pwm, pi, rng))
            write_fa(os.path.join(d, 'bg.fa'), sample(kind, n * BG_MULT, rng, False))
            json.dump(dict(dataset=did, motif=name.rstrip('('), n=n, pi=pi, background=kind, rep=rep,
                           truth=[consensus(pwm)], pwm=pwm.tolist()), open(os.path.join(d, 'truth.json'), 'w'))
            print(did, consensus(pwm), flush=True)
