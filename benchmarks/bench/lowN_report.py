# %% lowN_report.py: apply the registered LN1 decision rules (splits/protocol_lowN.md) mechanically
"""
Input: the evaluate.py table for the low-N panel (and optionally the N = 5,000 dev table).
Prints per N and tool: located (contains) top-10 / top-1, P5 top-10, best PWM similarity, each with chance; then the
paired differences HASTUR arm - STREME and centred - plain over units (factor x replicate), with 95% bootstrap CIs,
and the breakpoint verdicts. Everything is repeated without REST (beyond STREME/HOMER default widths).
Usage: python bench/lowN_report.py <lowN_eval.tsv>
"""
import sys

import numpy as np
import pandas as pd

B = 10000
rng = np.random.default_rng(0)


def boot_ci(x):
    x = np.asarray(x, float)
    m = rng.choice(x, size=(B, len(x)), replace=True).mean(axis=1)
    return x.mean(), np.quantile(m, 0.025), np.quantile(m, 0.975)


def paired(d, a, b, col='hit10_contains'):
    p = d.pivot_table(index='dataset', columns='tool', values=col, aggfunc='first')
    p = p.dropna(subset=[a, b])
    return boot_ci(p[a].astype(float) - p[b].astype(float)), len(p)


def report(d, label):
    print(f'\n######## {label}')
    s = d.groupby(['n', 'tool']).agg(units=('dataset', 'size'), located10=('hit10_contains', 'mean'),
                                     located1=('hit1_contains', 'mean'), ch_loc=('chance10_contains', 'mean'),
                                     p5_10=('hit10_p5', 'mean'), ch_p5=('chance10_p5', 'mean'),
                                     best_sim=('best_sim', 'mean'), ch_sim=('chance_best_sim', 'mean'),
                                     runtime=('runtime_s', 'median'))
    print(s.round(3).to_string())
    verdict = {'alive': False, 'none': True, 'centring': False}
    for n, g in d.groupby('n'):
        rows = []
        for a, b in (('hlat', 'streme'), ('hlatc', 'streme'), ('hlat', 'homer'), ('hlatc', 'hlat')):
            if {a, b} <= set(g.tool):
                (m, lo, hi), k = paired(g, a, b)
                rows.append(f'{a}-{b}: {m:+.3f} [{lo:+.3f}, {hi:+.3f}] (n={k})')
                if b == 'streme' and n in (250, 500) and m >= 0.15 and lo > 0:
                    verdict['alive'] = True
                if b == 'streme' and abs(m) > 0.05:
                    verdict['none'] = False
                if a == 'hlatc' and b == 'hlat' and n <= 1000 and m >= 0.10 and lo > 0:
                    verdict['centring'] = True
        print(f'N={n}: ' + ' | '.join(rows))
    print('VERDICTS (registered rules): accuracy thesis alive:', verdict['alive'],
          '| no accuracy gain (all |diff| <= 0.05):', verdict['none'], '| summit weighting works:', verdict['centring'])


if __name__ == '__main__':
    d = pd.read_csv(sys.argv[1], sep='\t')
    d['n'] = d.dataset.str.extract(r'__n(\d+)__').astype(float).astype('Int64') if 'n' not in d else d['n']
    report(d, 'all 24 dev factors')
    report(d[~d.dataset.str.startswith('REST__')], 'without REST (inexpressible for STREME/HOMER defaults)')
