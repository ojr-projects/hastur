# %% regulon_report.py: D1 readouts 1-4 (splits/protocol_D1.md, amendment 1), applied mechanically
"""
Usage: python bench/regulon_report.py <d1_calls.tsv> <regulon_oracle.tsv>
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from regulon_eval import thresholds  # noqa: E402

HL = os.environ.get('HL_TOOL', 'hlat_a3')   # LN7: the HASTUR arm to report

RATES = (0.05, 0.1, 0.5)
rng = np.random.default_rng(0)
d = pd.read_csv(sys.argv[1], sep='\t')
tau, EV = thresholds(sys.argv[2])
d['runtime'] = pd.to_numeric(d.runtime, errors='coerce')
units = d.drop_duplicates(['unit', 'tool'])
pos = d[~d.null]
nul = d[d.null]
n_null = nul.groupby('tool').unit.nunique()
L = pos.dropna(subset=['tf']).drop_duplicates('unit').set_index('unit')
site_len = {u: None for u in L.index}


def per_regulon(g, mask):
    """unit -> 1 if any call passing `mask` hits; every evaluable regulon present (0 when no output)."""
    h = g[mask & g.hit.fillna(False).astype(bool)].groupby('unit').size().gt(0)
    return pd.Series({u: int(h.get(u, False)) for u in EV})


def boot(x, B=10000):
    x = np.asarray(x, float)
    m = rng.choice(x, size=(B, len(x))).mean(axis=1)
    return np.quantile(m, 0.025), np.quantile(m, 0.975)


print(f'units per tool: {units.groupby("tool").unit.nunique().to_dict()}; null runs: {n_null.to_dict()}')
print(f'STREME S-mode units (no hold-out; amendment 1): {units[units.s_mode == True].groupby("tool").unit.nunique().to_dict()}')
print(f'evaluable regulons: {len(EV)}')

print('\n== 1. native operating points')
for tool in (HL, 'streme'):
    p, n = pos[pos.tool == tool], nul[nul.tool == tool]
    if p.empty or n.empty:
        continue
    nat = p.native.fillna(False).astype(bool)
    hits = per_regulon(p, nat)
    ch = p[nat & p.chance_hit.fillna(False).astype(bool)].groupby('unit').size().gt(0)
    chance = np.mean([bool(ch.get(u, False)) for u in EV])
    fc = n.groupby('unit').native.apply(lambda s: s.fillna(False).astype(bool).sum())
    calls34 = int(nat.sum())
    ne = int(p[nat & ~p.unit.isin(EV)].shape[0])
    lo, hi = boot(hits.values)
    fdr = fc.mean() * L.shape[0] / max(calls34, 1)
    print(f'{tool:8s} sensitivity {hits.mean():.3f} [{lo:.3f},{hi:.3f}] ({hits.sum()}/{len(EV)}) | false calls/null run '
          f'{fc.mean():.3f} (runs with any: {(fc > 0).mean():.3f}) | native calls on 34 regulons {calls34} '
          f'(non-evaluable: {ne}) | est. screen FDR {fdr:.3f} | mismatched-truth chance {chance:.3f}')
    print(f'         missed: {" ".join(u for u in EV if not hits[u])}')

print('\n== 2. frontier (sensitivity on evaluable regulons vs false calls per null run)')
best = {}
for tool in (HL, 'streme_all', 'homer'):
    p, n = pos[pos.tool == tool], nul[nul.tool == tool]
    if p.empty or n.empty:
        continue
    grid = np.unique(np.concatenate([p.s.dropna().values, n.s.dropna().values, [np.inf]]))
    ns = n.groupby('unit').s.apply(lambda s: np.sort(s.dropna().values))
    curve = []
    for t in grid:
        fc = np.mean([np.sum(v >= t) for v in ns]) if len(ns) else 0.0
        curve.append((t, fc))
    line = [f'{tool:10s}']
    for r in RATES:
        ok = [t for t, fc in curve if fc <= r]
        t = min(ok)
        h = per_regulon(p, p.s >= t)
        lo, hi = boot(h.values)
        best[(tool, r)] = h
        line.append(f'@{r}: {h.mean():.3f} [{lo:.3f},{hi:.3f}] (t={t:.2f})')
    print(' | '.join(line))
for r in RATES:
    if (HL, r) in best and ('streme_all', r) in best:
        x = (best[(HL, r)] - best[('streme_all', r)]).values
        lo, hi = boot(x)
        print(f'HASTUR - STREME sensitivity @ {r} false calls/run: {x.mean():+.3f} [{lo:+.3f},{hi:+.3f}]')

print('\n== 3. cost (runtime s over all units)')
print(units.groupby('tool').runtime.agg(['sum', 'median', 'max']).round(1).to_string())

print('\n== 4. secondary (descriptive)')
# site length from the oracle-free truth file is not in the calls table; take it from the unit name -> regulon_make
tl = {}
root = os.path.dirname(os.path.abspath(sys.argv[2]))
for u in EV:
    p = os.path.join(root, 'regulon', u, 'truth.json')
    if os.path.exists(p):
        import json
        tl[u] = json.load(open(p))['site_len']
if tl:
    short = [u for u in EV if tl.get(u, 99) <= 16]
    long_ = [u for u in EV if tl.get(u, 0) >= 17]
    for tool in (HL, 'streme'):
        p = pos[pos.tool == tool]
        if p.empty:
            continue
        h = per_regulon(p, p.native.fillna(False).astype(bool))
        print(f'{tool:8s} native sensitivity: sites <= 16 bp {h[short].mean():.3f} (n={len(short)}) | '
              f'>= 17 bp {h[long_].mean():.3f} (n={len(long_)})')
    for tool in (HL, 'streme_all', 'homer'):
        if (tool, 0.1) in best:
            h = best[(tool, 0.1)]
            print(f'{tool:10s} @0.1 sensitivity: <= 16 bp {h[short].mean():.3f} | >= 17 bp {h[long_].mean():.3f}')
for tool in (HL, 'streme', 'streme_all', 'homer'):
    p = pos[(pos.tool == tool) & pos.unit.isin(EV) & pos.hit.fillna(False).astype(bool)]
    if not p.empty:
        b = p.sort_values('sites').groupby('unit').tail(1)
        print(f'{tool:10s} best hitting call per regulon: sites {b.sites.mean():.3f}, P5 {b.p5.astype(float).mean():.3f}, '
              f'PWM sim {b.sim.mean():.3f} (n={len(b)})')
