# %% c4_summarize.py: merge c4 part files and print the THEORY.md section 6.3 table
"""Run after c4_lattice_to_pwm.py (one part file per condition): ~/micromamba/envs/nanomotif/bin/python theory/checks/c4_summarize.py"""
import glob
import os

import numpy as np
import pandas as pd

RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
parts = sorted(glob.glob(os.path.join(RES, 'c4_lattice_to_pwm_part*.tsv')))
df = pd.concat(pd.read_csv(p, sep='\t') for p in parts)
df.to_csv(os.path.join(RES, 'c4_lattice_to_pwm.tsv'), sep='\t', index=False)

order = [('pwm_L8_mid', 500, 0.1), ('pwm_L8_mid', 500, 0.3), ('pwm_L8_mid', 2000, 0.05), ('pwm_L8_mid', 2000, 0.1),
         ('pwm_L12_low', 2000, 0.1), ('pwm_L12_low', 2000, 0.3), ('NRF1', 2000, 0.1), ('REST', 2000, 0.1),
         ('CTCF(', 2000, 0.1), ('Sox2', 2000, 0.1), ('p53(', 2000, 0.1)]
top10 = lambda s: 'n/a' if (s >= 10 ** 6).all() else f'{np.mean(s <= 10):.1f}'   # n/a: no motif-consistent pattern exists
rows = []
for m, n, pi in order:
    d = df[(df.motif == m) & (df.n == n) & (df.pi == pi)]
    if d.empty:
        continue
    f = lambda c: f'{d[c].mean():.2f}'
    rows.append({'motif': m.rstrip('('), 'N': n, 'π': pi, 'L': int(d.L.iloc[0]), 'reps': len(d),
                 'top10 k5': top10(d.rank_EXACT_k5), 'top10 L8': top10(d.rank_EXACT_L8),
                 'top10 soft L8': top10(d['rank_SOFT_L8_b0.05']),
                 'LAT+EM (k5 seeds, by gain)': f('sim_LAT_EM'), 'LAT8+EM top seed': f('sim_LAT8_EM_seed1'),
                 'LAT8+EM by gain': f('sim_LAT8_EM'), 'LAT+EM top k5 seed': f('sim_LAT_EM_seed1'),
                 'ORACLE+EM': f('sim_ORACLE_EM'), 'STREME top 1': f('sim_STREME'), 'STREME best of 3': f('sim_STREME_best3'),
                 'flat': f'{d.null_sim.iloc[0]:.2f}'})
t = pd.DataFrame(rows)
print('| ' + ' | '.join(t.columns) + ' |')
print('|' + '---|' * len(t.columns))
for _, r in t.iterrows():
    print('| ' + ' | '.join(str(v) for v in r.values) + ' |')

# Paired comparisons over all datasets
print()
for a, b in [('sim_LAT8_EM_seed1', 'sim_STREME'), ('sim_LAT8_EM_seed1', 'sim_STREME_best3'), ('sim_LAT_EM', 'sim_STREME'),
             ('sim_LAT8_EM_seed1', 'sim_LAT8_EM'), ('sim_LAT8_EM_seed1', 'sim_ORACLE_EM')]:
    dd = df[a] - df[b]
    print(f'{a} - {b}: mean {dd.mean():+.3f}, wins {np.mean(dd > 0.02):.2f}, losses {np.mean(dd < -0.02):.2f}')
print('median runtime per dataset (s): lattice', df.t_lattice.median().round(1), '| EM from 10 k5 seeds',
      df.t_lat_em.median().round(1), '| EM from 10 L8 seeds', df.t_lat8_em.median().round(1), '| STREME', df.t_streme.median().round(1))
