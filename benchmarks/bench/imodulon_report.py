# %% imodulon_report.py: E1 readouts 1-5 (splits/protocol_E1.md), applied mechanically
"""
Usage: python bench/imodulon_report.py <e1_calls.tsv> <e1_oracle.tsv> <units_root>
"""
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from imodulon_eval import thresholds  # noqa: E402

HL = 'hlat_f2r'
NATIVE = (HL, 'streme', 'meme', 'xstreme')
FRONT = (HL, 'streme_all', 'homer', 'meme')
RATES = (0.05, 0.1, 0.5)
NOT_D1 = {'OxyR', 'EvgA', 'GadW', 'GalR', 'GalS', 'GntR', 'TyrR'}
rng = np.random.default_rng(0)
d = pd.read_csv(sys.argv[1], sep='\t')
tau, ev = thresholds(sys.argv[2])
root = sys.argv[3]
EM = sorted({m for m, _ in ev})
mods = pd.read_csv(os.path.join(root, 'modules.tsv'), sep='\t')
inc = mods[mods.included].set_index('safe')
d['runtime'] = pd.to_numeric(d.runtime, errors='coerce')
units = d.drop_duplicates(['unit', 'tool'])
pos, nul = d[~d.null], d[d.null]
n_mod = inc.shape[0]


def per_module(g, mask, which=EM):
    h = g[mask & g.hit.fillna(False).astype(bool)].groupby('unit').size().gt(0)
    return pd.Series({u: int(h.get(u, False)) for u in which})


def boot(x, B=10000):
    x = np.asarray(x, float)
    m = rng.choice(x, size=(B, len(x))).mean(axis=1)
    return np.quantile(m, 0.025), np.quantile(m, 0.975)


print(f'units per tool: {units.groupby("tool").unit.nunique().to_dict()}')
print(f'null runs per tool: {nul.groupby("tool").unit.nunique().to_dict()}')
print(f'STREME S-mode units: {units[units.s_mode == True].groupby("tool").unit.nunique().to_dict()}')
print(f'modules {n_mod}; evaluable {len(EM)}')

print('\n== 1. native operating points')
for tool in NATIVE:
    p, n = pos[pos.tool == tool], nul[nul.tool == tool]
    if p.empty or n.empty:
        continue
    nat = p.native.fillna(False).astype(bool)
    hits = per_module(p, nat)
    ch = p[nat & p.chance_hit.fillna(False).astype(bool)].groupby('unit').size().gt(0)
    chance = np.mean([bool(ch.get(u, False)) for u in EM])
    fc = n.groupby('unit').native.apply(lambda s: s.fillna(False).astype(bool).sum())
    calls = int(nat.sum())
    ne = int(p[nat & ~p.unit.isin(EM)].shape[0])
    lo, hi = boot(hits.values)
    fdr = fc.mean() * n_mod / max(calls, 1)
    print(f'{tool:8s} sensitivity {hits.mean():.3f} [{lo:.3f},{hi:.3f}] ({hits.sum()}/{len(EM)}) | false calls/null run '
          f'{fc.mean():.3f} (runs with any {(fc > 0).mean():.3f}, n={len(fc)}) | native calls on {n_mod} modules {calls} '
          f'(non-evaluable {ne}) | est. screen FDR {fdr:.3f} | mismatched-truth chance {chance:.3f}')

print('\n== 2. frontier (sensitivity on evaluable modules vs false calls per null run)')
best = {}
for tool in FRONT:
    p, n = pos[pos.tool == tool], nul[nul.tool == tool]
    if p.empty or n.empty:
        continue
    grid = np.unique(np.concatenate([p.s.dropna().values, n.s.dropna().values, [np.inf]]))
    ns = [np.sort(v.dropna().values) for _, v in n.groupby('unit').s]
    nrun = n.unit.nunique()
    allv = np.sort(np.concatenate(ns)) if ns else np.array([])
    line = [f'{tool:10s}']
    for r in RATES:
        fcs = (len(allv) - np.searchsorted(allv, grid, side='left')) / nrun
        t = grid[np.argmax(fcs <= r)]
        h = per_module(p, p.s >= t)
        lo, hi = boot(h.values)
        best[(tool, r)] = h
        line.append(f'@{r}: {h.mean():.3f} [{lo:.3f},{hi:.3f}] (t={t:.2f})')
    print(' | '.join(line))
for other in FRONT[1:]:
    for r in RATES:
        if (HL, r) in best and (other, r) in best:
            x = (best[(HL, r)] - best[(other, r)]).values
            lo, hi = boot(x)
            print(f'HASTUR - {other:10s} @ {r}: {x.mean():+.3f} [{lo:+.3f},{hi:+.3f}]')

print('\n== 3. native calls on modules without a cis truth (descriptive)')
nc = inc[~inc.index.isin(EM)]
for tool in NATIVE:
    p = pos[(pos.tool == tool) & pos.unit.isin(nc.index)]
    if p.empty:
        continue
    k = p.groupby('unit').native.apply(lambda s: s.fillna(False).astype(bool).sum()).reindex(nc.index, fill_value=0)
    by = k.groupby(nc.category).mean().round(2).to_dict()
    print(f'{tool:8s} calls per module {k.mean():.2f}; modules with none {(k == 0).mean():.2f}; by category {by}')

print('\n== 4. cost (runtime s over all units)')
print(units.groupby('tool').runtime.agg(['sum', 'median', 'max']).round(1).to_string())

print('\n== 5. strata (descriptive)')
size = inc.loc[EM, 'n_regions']
sl = {}
for m in EM:
    j = json.load(open(os.path.join(root, m, 'truth.json')))
    sl[m] = max(t['site_len'] for t in j['truths'] if (m, t['tf']) in ev)
fresh = [m for m in EM if any(t in NOT_D1 for t in inc.loc[m, 'truth_tfs'].split(','))]
strata = {'5-9 regions': [m for m in EM if size[m] < 10], '10-19': [m for m in EM if 10 <= size[m] < 20],
          '>=20': [m for m in EM if size[m] >= 20], 'site<=16': [m for m in EM if sl[m] <= 16],
          'site>=17': [m for m in EM if sl[m] >= 17], 'not-D1 TFs': fresh}
for tool in NATIVE:
    p = pos[pos.tool == tool]
    if p.empty:
        continue
    h = per_module(p, p.native.fillna(False).astype(bool))
    print(f'{tool:8s} native: ' + ' | '.join(f'{k} {h[v].mean():.2f} (n={len(v)})' for k, v in strata.items()))
for tool in FRONT:
    if (tool, 0.1) in best:
        h = best[(tool, 0.1)]
        print(f'{tool:10s} @0.1: ' + ' | '.join(f'{k} {h[v].mean():.2f}' for k, v in strata.items()))
print('\nper module (native hit by tool; @0.1 hit for frontier tools)')
tab = pd.DataFrame({t: per_module(pos[pos.tool == t], pos[pos.tool == t].native.fillna(False).astype(bool))
                    for t in NATIVE if not pos[pos.tool == t].empty})
for t in FRONT:
    if (t, 0.1) in best:
        tab[t + '@0.1'] = best[(t, 0.1)]
tab.insert(0, 'regions', size)
tab.insert(1, 'truth', [inc.loc[m, 'truth_tfs'] for m in tab.index])
print(tab.to_string())
