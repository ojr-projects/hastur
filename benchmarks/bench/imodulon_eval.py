# %% imodulon_eval.py: E1 oracle test and per-call scoring (splits/protocol_E1.md)
"""
A module can carry several annotated regulators with RegulonDB truths (bench/imodulon_make.py). Each (module, TF)
pair is scored with D1's site recovery (bench/regulon_eval.py), restricted to the module's fg regions that hold at
least half of one of the TF's sites.
Oracle mode: for every pair with >= MIN_SITE_REG site-holding regions: the truth matrix, 20 column-shuffled truths and
every other E1 truth TF (those not annotated to the module). Chance threshold tau = max of the two 95th percentiles.
Evaluable pair: truth recovery >= tau + MARGIN. Evaluable module: at least one evaluable pair.
Tools mode: every call of every tool on every unit; on a module, hit = recovery >= tau on any evaluable pair of it;
chance_hit = the same rule against a fixed other evaluable module's pairs.
Usage: python bench/imodulon_eval.py oracle <root>                       -> <root>/../e1_oracle.tsv
       python bench/imodulon_eval.py tools <root> <e1_oracle.tsv> <tools> <out.tsv>
"""
import hashlib
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from regulon_eval import onehot, read_fa, site_recovery  # noqa: E402

MIN_SITE_REG, MARGIN, N_SHUF = 4, 0.15, 20


def modules(root):
    for u in sorted(os.listdir(root)):
        p = os.path.join(root, u, 'truth.json')
        if os.path.exists(p):
            j = json.load(open(p))
            if not j.get('null'):
                yield u, j


def oracle(root):
    M = dict(modules(root))
    allt = {}
    for u, j in M.items():
        for t in j['truths']:
            allt[t['tf']] = np.array(t['pwm'])
    rows = []
    for u, j in M.items():
        fg = read_fa(os.path.join(root, u, 'fg.fa'))
        own = {t['tf'] for t in j['truths']}
        for t in j['truths']:
            if len(t['fg_idx']) < MIN_SITE_REG:
                continue
            oh = [onehot(fg[k]) for k in t['fg_idx']]
            T = np.array(t['pwm'])
            rng = np.random.default_rng(int(hashlib.sha256(f'E1oracle|{u}|{t["tf"]}'.encode()).hexdigest()[:16], 16))
            rows.append(dict(module=u, tf=t['tf'], kind='truth', sites=site_recovery(T, oh, t['sites'])))
            for _ in range(N_SHUF):
                rows.append(dict(module=u, tf=t['tf'], kind='shuffled',
                                 sites=site_recovery(T[rng.permutation(len(T))], oh, t['sites'])))
            for v, P in sorted(allt.items()):
                if v not in own:
                    rows.append(dict(module=u, tf=t['tf'], kind='other_tf', other=v, sites=site_recovery(P, oh, t['sites'])))
    df = pd.DataFrame(rows)
    out = os.path.join(root, '..', 'e1_oracle.tsv')
    df.to_csv(out, sep='\t', index=False)
    tau, ev = thresholds(out)
    q = df[df.kind != 'truth'].groupby(['module', 'tf', 'kind']).sites.quantile(0.95).unstack()
    tr = df[df.kind == 'truth'].set_index(['module', 'tf']).sites
    tab = pd.concat([tr.rename('truth'), q.add_prefix('p95_')], axis=1)
    tab['tau'] = pd.Series(tau)
    tab['evaluable'] = [k in ev for k in tab.index]
    print(tab.round(2).to_string())
    print('\nsite recovery quantiles by kind:')
    print(df.groupby('kind').sites.quantile([0.5, 0.9, 0.95, 0.99]).unstack().round(3).to_string())
    em = sorted({m for m, _ in ev})
    print(f'\n{len(tab)} pairs tested (>= {MIN_SITE_REG} site-holding regions); {len(ev)} evaluable pairs; '
          f'{len(em)} evaluable modules: {" ".join(em)}')


def thresholds(oracle_tsv):
    """{(module, tf): tau} and the evaluable pairs."""
    df = pd.read_csv(oracle_tsv, sep='\t')
    tau = df[df.kind != 'truth'].groupby(['module', 'tf', 'kind']).sites.quantile(0.95).unstack().max(axis=1)
    tr = df[df.kind == 'truth'].set_index(['module', 'tf']).sites
    tau = {k: float(v) for k, v in tau.items()}
    return tau, sorted(k for k in tau if tr[k] >= tau[k] + MARGIN)


def tools(root, oracle_tsv, tool_list, out):
    from frontier import scored_calls
    tau, ev = thresholds(oracle_tsv)
    M = dict(modules(root))
    em = sorted({m for m, _ in ev})
    pairs = {m: [t for t in M[m]['truths'] if (m, t['tf']) in ev] for m in em}
    fa = {m: read_fa(os.path.join(root, m, 'fg.fa')) for m in em}
    oh = {(m, t['tf']): [onehot(fa[m][k]) for k in t['fg_idx']] for m in em for t in pairs[m]}

    def hit(pwm, m):
        best, tf = -1.0, ''
        for t in pairs[m]:
            sr = site_recovery(pwm, oh[(m, t['tf'])], t['sites'])
            if sr - tau[(m, t['tf'])] > best:
                best, tf = sr - tau[(m, t['tf'])], t['tf']
        return best >= 0, tf, best
    rows = []
    for d in sorted(os.listdir(root)):
        tj = os.path.join(root, d, 'truth.json')
        if not os.path.exists(tj):
            continue
        j = json.load(open(tj))
        for tool in tool_list:
            calls = scored_calls(os.path.join(root, d), tool)
            if calls is None:
                continue
            rt = open(os.path.join(root, d, tool, 'runtime_s')).read().strip()
            if tool.startswith('hlat'):     # native: the tool's own call flag (p_adj <= 0.05 for F2)
                nat = pd.read_csv(os.path.join(root, d, tool, 'hlat.tsv'), sep='\t').sig.astype(bool).tolist() if calls else []
            else:                           # native: E (or STREME's S) <= 0.05 for the MEME-suite tools; HOMER has none
                nat = [s >= -np.log10(0.05) - 1e-9 if tool in ('streme', 'meme', 'xstreme') else False for s, _ in calls]
            st = os.path.join(root, d, tool, 'streme.txt')
            smode = tool.startswith('streme') and os.path.exists(st) and ' S= ' in open(st).read()
            base = dict(unit=d, tool=tool, null=bool(j.get('null')), module=j.get('module'), runtime=rt,
                        n_calls=len(calls), s_mode=smode, evaluable=d in em)
            if not calls:
                rows.append(dict(base, rank=-1))
                continue
            mm = em[(em.index(d) + len(em) // 2) % len(em)] if d in em else None
            for r, (s, pwm) in enumerate(calls):
                row = dict(base, rank=r, s=s, native=bool(nat[r]))
                if d in em:
                    h, tf, margin = hit(pwm, d)
                    row.update(hit=h, hit_tf=tf if h else '', margin=margin, chance_hit=hit(pwm, mm)[0])
                rows.append(row)
    pd.DataFrame(rows).to_csv(out, sep='\t', index=False)
    print(f'{len(rows)} rows -> {out}; evaluable modules ({len(em)})')


if __name__ == '__main__':
    if sys.argv[1] == 'oracle':
        oracle(sys.argv[2])
    elif sys.argv[1] == 'tools':
        tools(sys.argv[2], sys.argv[3], sys.argv[4].split(','), sys.argv[5])
