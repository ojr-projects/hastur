# %% regulon_eval.py: D1 readouts (splits/protocol_D1.md) and the oracle test of the RegulonDB truths
"""
Readouts for a call (a PWM) on a regulon unit:
  contains / P5   the call's IUPAC against the truth IUPAC (eval/score.py rules)
  sim             pwm_sim against the truth matrix (full-length: uncovered informative columns count as uniform)
  sites           site recovery: fraction of regulon regions whose best-scoring window (both strands) overlaps an
                  annotated site of the TF by >= half the shorter of the two; matrix-construction-free
Oracle mode (no tool output needed): the truth matrix itself, 20 column-shuffled truths per TF (chance) and every
other TF's truth (cross-match), on every readout.
Usage: python bench/regulon_eval.py oracle <root>
       python bench/regulon_eval.py tools <root> <regulon_oracle.tsv> hlat_a3,streme,streme_all,homer <out.tsv>
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'eval'))
from score import as_iupac, hit_contains, hit_p5, pwm_sim  # noqa: E402

IDX = {c: i for i, c in enumerate('ACGT')}


def read_fa(p):
    out, buf = [], []
    for line in open(p):
        if line.startswith('>'):
            if buf:
                out.append(''.join(buf))
            buf = []
        else:
            buf.append(line.strip())
    if buf:
        out.append(''.join(buf))
    return out


def onehot(s):
    X = np.zeros((len(s), 4))
    for i, c in enumerate(s):
        if c in IDX:
            X[i, IDX[c]] = 1
    return X


def site_recovery(pwm, seqs_oh, sites):
    """Fraction of sequences whose best window (either strand) overlaps an annotated site by >= half the shorter."""
    M = np.asarray(pwm, float)
    M = M / M.sum(axis=1, keepdims=True)
    S = np.log2((M + 0.01) / 0.26)
    W = len(S)
    hit = 0
    for X, st in zip(seqs_oh, sites):
        n = len(X) - W + 1
        if n <= 0:
            continue
        sc = np.full(n, -np.inf)
        for P in (S, S[::-1, ::-1]):
            v = sum(X[j:j + n] @ P[j] for j in range(W))
            sc = np.maximum(sc, v)
        p = int(np.argmax(sc))
        if any(min(p + W, b) - max(p, a) >= min(W, b - a) / 2 for a, b in st):
            hit += 1
    return hit / len(seqs_oh)


def readouts(pwm, truth, seqs_oh):
    iu = as_iupac(np.asarray(pwm))
    t = truth['truth'][0]
    return dict(contains=hit_contains(iu, t), p5=hit_p5(iu, t), sim=pwm_sim(pwm, np.array(truth['pwm'])),
                sites=site_recovery(pwm, seqs_oh, truth['sites']))


def units(root):
    for u in sorted(os.listdir(root)):
        p = os.path.join(root, u, 'truth.json')
        if os.path.exists(p):
            j = json.load(open(p))
            if not j.get('null'):
                yield u, j


def oracle(root):
    import pandas as pd
    rng = np.random.default_rng(0)
    U = list(units(root))
    rows = []
    for u, j in U:
        oh = [onehot(s) for s in read_fa(os.path.join(root, u, 'fg.fa'))]
        T = np.array(j['pwm'])
        rows.append(dict(tf=u, kind='truth', **readouts(T, j, oh)))
        for k in range(20):
            rows.append(dict(tf=u, kind='shuffled', **readouts(T[rng.permutation(len(T))], j, oh)))
        for v, jv in U:
            if v != u:
                rows.append(dict(tf=u, kind='other_tf', **readouts(np.array(jv['pwm']), j, oh)))
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(root, '..', 'regulon_oracle.tsv'), sep='\t', index=False)
    print('mean readout by kind:')
    print(df.groupby('kind')[['contains', 'p5', 'sim', 'sites']].mean().round(3).to_string())
    print('\nsite recovery quantiles by kind:')
    print(df.groupby('kind').sites.quantile([0.5, 0.9, 0.95, 0.99]).unstack().round(3).to_string())
    print('\nsim quantiles by kind:')
    print(df.groupby('kind').sim.quantile([0.5, 0.9, 0.95, 0.99]).unstack().round(3).to_string())
    print('\nper TF (truth row; shuffled and other-TF 95th percentiles):')
    q = df[df.kind != 'truth'].groupby(['tf', 'kind']).sites.quantile(0.95).unstack()
    t = df[df.kind == 'truth'].set_index('tf')
    print(pd.concat([t[['contains', 'p5', 'sim', 'sites']], q.add_prefix('sites95_')], axis=1).round(2).to_string())


MARGIN = 0.15   # evaluable TF: truth site recovery >= chance threshold + MARGIN (protocol_D1.md)


def thresholds(oracle_tsv):
    """Per-TF chance threshold tau (max of shuffled and other-TF 95th percentiles) and the evaluable set."""
    import pandas as pd
    df = pd.read_csv(oracle_tsv, sep='\t')
    tau = df[df.kind != 'truth'].groupby(['tf', 'kind']).sites.quantile(0.95).unstack().max(axis=1)
    tr = df[df.kind == 'truth'].set_index('tf').sites
    return tau.to_dict(), sorted(t for t in tau.index if tr[t] >= tau[t] + MARGIN)


def tools(root, oracle_tsv, tool_list, out):
    """Every call of every tool, scored: -log10 significance (as bench/frontier.py), site recovery on its own TF,
    and on a mismatched TF (chance). Null units: scores only."""
    import pandas as pd
    from frontier import scored_calls
    tau, evaluable = thresholds(oracle_tsv)
    U = dict(units(root))
    names = sorted(U)
    oh = {u: [onehot(s) for s in read_fa(os.path.join(root, u, 'fg.fa'))] for u in names}
    rows = []
    for d in sorted(os.listdir(root)):
        if not os.path.exists(os.path.join(root, d, 'truth.json')):
            continue
        j = json.load(open(os.path.join(root, d, 'truth.json')))
        for tool in tool_list:
            calls = scored_calls(os.path.join(root, d), tool)
            if calls is None:
                continue
            rt = open(os.path.join(root, d, tool, 'runtime_s')).read().strip()
            if tool.startswith('hlat'):     # native: the Yellow Sign flag (z >= its family's zcrit), in call order
                nat = pd.read_csv(os.path.join(root, d, tool, 'hlat.tsv'), sep='\t').sig.astype(bool).tolist() if calls else []
            else:                           # native: STREME E <= 0.05; HOMER has none
                nat = [s >= -np.log10(0.05) - 1e-9 if tool == 'streme' else False for s, _ in calls]
            st = os.path.join(root, d, tool, 'streme.txt')
            smode = tool.startswith('streme') and os.path.exists(st) and ' S= ' in open(st).read()   # amendment 1
            base = dict(unit=d, tool=tool, null=bool(j.get('null')), tf=j['tf'], runtime=rt, n_calls=len(calls), s_mode=smode)
            if not calls:
                rows.append(dict(base, rank=-1))
                continue
            mm = names[(names.index(d) + len(names) // 2) % len(names)] if d in U else None   # fixed mismatched TF
            for r, (s, pwm) in enumerate(calls):
                row = dict(base, rank=r, s=s, native=bool(nat[r]))
                if d in U:
                    sr = site_recovery(pwm, oh[d], U[d]['sites'])
                    row.update(evaluable=d in evaluable, sites=sr, hit=sr >= tau[d], p5=hit_p5(as_iupac(np.asarray(pwm)), U[d]['truth'][0]),
                               sim=pwm_sim(pwm, np.array(U[d]['pwm'])),
                               chance_hit=site_recovery(pwm, oh[mm], U[mm]['sites']) >= tau[mm])
                rows.append(row)
    pd.DataFrame(rows).to_csv(out, sep='\t', index=False)
    print(f'{len(rows)} rows -> {out}; evaluable TFs ({len(evaluable)}): {" ".join(evaluable)}')


if __name__ == '__main__':
    if sys.argv[1] == 'oracle':
        oracle(sys.argv[2])
    elif sys.argv[1] == 'tools':
        tools(sys.argv[2], sys.argv[3], sys.argv[4].split(','), sys.argv[5])
