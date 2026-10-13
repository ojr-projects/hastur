# %% repro.py: replicate reproducibility of significant calls (FR1 also-reported; ROADMAP D8)
"""
Low-N panel: each factor x N has three independent subsamples (r0, r1, r2). For every call a tool reports as
significant at its native threshold (HASTUR: Yellow Sign; STREME default run: E <= 0.05 or S <= 0.05), ask whether a
matching significant call (symmetric matrix similarity >= SIM) appears in each other replicate.
Reported per tool and N: significant calls per run, the share of them reproduced in at least one / both other
replicates, and the share of factor x N cells whose truth-matching motif is significant in all three.
Usage: python bench/repro.py <lowN root> <truth_map.json> [tools=hlat_a3,streme]
"""
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'eval'))
from frontier import scored_calls  # noqa: E402
from score import pwm_sim  # noqa: E402

SIM = 0.8


def sym(a, b):
    return max(pwm_sim(a, b), pwm_sim(b, a))


def sig_calls(d, tool):
    c = scored_calls(d, tool)
    if c is None:
        return None
    if tool.startswith('hlat'):
        sig = pd.read_csv(os.path.join(d, tool, 'hlat.tsv'), sep='\t').sig.astype(bool).tolist() if c else []
        return [m for (s, m), g in zip(c, sig) if g]
    return [m for s, m in c if s >= -np.log10(0.05) - 1e-9]


def main():
    root, tmap = sys.argv[1], json.load(open(sys.argv[2]))
    tools = (sys.argv[3] if len(sys.argv) > 3 else 'hlat_a3,streme').split(',')
    cells = {}
    for ds in sorted(os.listdir(root)):
        f, n, r = ds.split('__')
        cells.setdefault((f, n[1:]), []).append(ds)
    rows, cell_rows = [], []
    for (f, n), reps in cells.items():
        for tool in tools:
            S = {r: sig_calls(os.path.join(root, r), tool) for r in reps}
            if any(v is None for v in S.values()) or len(S) < 3:
                continue
            for r, calls in S.items():
                others = [S[o] for o in S if o != r]
                for m in calls:
                    hit = [any(sym(m, x) >= SIM for x in oc) for oc in others]
                    rows.append(dict(factor=f, n=int(n), tool=tool, rep=r, any_other=any(hit), both_other=all(hit)))
            if f in tmap:
                T = np.array(tmap[f]['pwm'])
                found = [any(pwm_sim(m, T) >= 0.85 for m in S[r]) for r in reps]
                cell_rows.append(dict(factor=f, n=int(n), tool=tool, truth_sig_all3=all(found), truth_sig_any=any(found),
                                      n_sig_mean=np.mean([len(S[r]) for r in reps])))
    d, c = pd.DataFrame(rows), pd.DataFrame(cell_rows)
    a = d.groupby(['n', 'tool']).agg(sig_calls=('any_other', 'size'), repro_any=('any_other', 'mean'),
                                     repro_both=('both_other', 'mean')).round(3)
    a['sig_per_run'] = (a.sig_calls / d.groupby(['n', 'tool']).rep.nunique()).round(2)
    print(a.to_string())
    print('\ntruth-matching motif significant in all three replicates / in any (matrix-truth factors):')
    print(c.groupby(['n', 'tool'])[['truth_sig_all3', 'truth_sig_any', 'n_sig_mean']].mean().round(3).to_string())


if __name__ == '__main__':
    main()
