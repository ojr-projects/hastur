# %% null_make.py: motif-free calibration pairs (DEV_INSIGHTS 3.7: calibrate before reading any power figure)
"""
<root>/<pair>/{fg.fa, bg.fa, truth.json(null=True)}; N = 2000 x 100 bp, background 5x, 10 replicates each of
  iid      fg and bg iid uniform
  markov   fg and bg from the same genome-like order-1 chain (CpG depleted)
  cpgshift fg CpG-island-like chain, bg genome-like chain (a composition difference, no motif)
Readout (bench/null_eval.py): calls each tool flags as significant; HASTUR's should be ~0.05 per run at most.
Usage: python bench/null_make.py <root> [seed_offset]
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sim_make import chain, write_fa  # noqa: E402

LSEQ, N, BGX = 100, 2000, 5
REPS = int(os.environ.get("NULL_REPS", "10"))
GENOME = chain(np.array([0.295, 0.205, 0.205, 0.295]), 0.25)
ISLAND = chain(np.array([0.19, 0.31, 0.31, 0.19]), 0.95)


BURN = 100   # stationary start: without it the uniform first base is a non-stationary artefact that alone produced
             # every residual false call under the CpG shift (fixed 2026-10-10; panels made before carry the artefact)


def sample(P, n, rng):
    if P is None:
        return rng.integers(0, 4, (n, LSEQ))
    X = np.empty((n, LSEQ + BURN), np.int64)
    X[:, 0] = rng.integers(0, 4, n)
    cdf = np.cumsum(P, axis=1)
    for t in range(1, LSEQ + BURN):
        X[:, t] = (rng.random(n)[:, None] > cdf[X[:, t - 1]]).sum(axis=1)
    return X[:, BURN:]


if __name__ == '__main__':
    root = sys.argv[1]
    seed0 = int(sys.argv[2]) if len(sys.argv) > 2 else 0      # LN4 confirmation panel: seed offset 900000
    for kind, (Pf, Pb) in {'iid': (None, None), 'markov': (GENOME, GENOME), 'cpgshift': (ISLAND, GENOME)}.items():
        for rep in range(REPS):
            d = os.path.join(root, f'null_{kind}__r{rep}')
            if os.path.exists(os.path.join(d, 'truth.json')):
                continue
            os.makedirs(d, exist_ok=True)
            rng = np.random.default_rng([424242, seed0, rep, len(kind)])   # tuple seed (2026-10-11 fix; the earlier additive seed made offsets 10,000 apart share 10 of 20 replicates)
            write_fa(os.path.join(d, 'fg.fa'), sample(Pf, N, rng))
            write_fa(os.path.join(d, 'bg.fa'), sample(Pb, N * BGX, rng))
            json.dump(dict(dataset=os.path.basename(d), kind=kind, rep=rep, null=True, truth=[]),
                      open(os.path.join(d, 'truth.json'), 'w'))
            print(d, flush=True)
