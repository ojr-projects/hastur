# %% c8_p1_chance_floor.py: chance floor of the P1 methylation lenient rule (Comsequent results, read-only)
"""
Scores Comsequent's anchored-walk top 10 and Nanomotif's calls against letter-shuffled truths with P1's own
`contains_truth` rule (imported verbatim from 14_meth_evaluate.py), to check that the methylation win is not a chance
artefact of the matching rule. Also reports Nanomotif's list length (its lenient check has no top-10 cap).
Run: ~/micromamba/envs/nanomotif/bin/python theory/checks/c8_p1_chance_floor.py
"""
import glob
import os

import numpy as np
import pandas as pd

P1 = '/home/oliver/code/comsequent/benchmarks/phase1_methylation'
src = open(os.path.join(P1, '14_meth_evaluate.py')).read()
exec(src[src.index('LENIENT_EXTRA = 2'):src.index("if 'anchored_walk_calls' in d:")])   # contains_truth, lenient_found

d = pd.concat(pd.read_csv(f, sep='\t') for f in glob.glob(os.path.join(P1, 'results/parts/meth_downsample_r10m4/*.tsv')))
rng = np.random.default_rng(0)


def shuf(t):
    letters = [c for c in t if c != 'N']
    rng.shuffle(letters)
    it = iter(letters)
    return ''.join(c if c == 'N' else next(it) for c in t)


real_c = np.mean([lenient_found(c, m, 10) for c, m in zip(d.anchored_walk_calls, d.motif)])
real_n = np.mean([lenient_found(c, m) for c, m in zip(d.nanomotif_calls, d.motif)])
ch_c, ch_n = [], []
for c, nc, m in zip(d.anchored_walk_calls, d.nanomotif_calls, d.motif):
    sh = [s for s in (shuf(m) for _ in range(5)) if s != m]
    if sh:
        ch_c.append(np.mean([lenient_found(c, s, 10) for s in sh]))
        ch_n.append(np.mean([lenient_found(nc, s) for s in sh]))
n_calls = d.nanomotif_calls.astype(str).str.split().str.len()
print(f'{len(d)} motif x run rows (meth_downsample_r10m4)')
print(f'Comsequent anchored walk, top 10: lenient {real_c:.2f}, shuffled-truth chance {np.mean(ch_c):.3f}')
print(f'Nanomotif, all calls:             lenient {real_n:.2f}, shuffled-truth chance {np.mean(ch_n):.3f}')
print(f'Nanomotif calls per run: median {n_calls.median():.0f}, max {n_calls.max():.0f}')
