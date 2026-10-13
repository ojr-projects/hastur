# %% c3_failure_anatomy.py: why Comsequent lost on PWM-like real data (THEORY.md section 4)
"""
Reads Comsequent's frozen real-data results (read-only) for P5 ChIP, P4a HT-SELEX and P3 RBNS, and classifies every
unit where a tool's top 10 missed the truth under the phase's own rule:

  precision  a top-10 call covers >= 80% of the truth's specified positions with at most 1 true conflict, but only
             when degenerate call letters that CONTAIN the truth letter are accepted (the registered rules require the
             call's letter set to be a subset of the truth's, so W at an A, or one N inside the core, fails)
  fragment   no precision call, but a top-10 call with >= F specified positions lies wholly inside the truth with no
             conflict (F = 4 for truths of >= 8 specified positions, else 3): the motif was detected, not represented
  none       neither
A chance baseline repeats the classification against 5 letter-shuffled copies of each truth (N positions kept), because
short consistent fragments can match by chance. (Scoring against another unit's truth is not a chance baseline here:
TF and RBP families share motifs, so it measures family similarity.)
Also reported: the share of units whose rank-1 call is low-information (< 6 bits: e.g. CG, ANTNA, GC).
Run: ~/micromamba/envs/nanomotif/bin/python theory/checks/c3_failure_anatomy.py
"""
import glob
import math
import os

import numpy as np
import pandas as pd

B = '/home/oliver/code/comsequent/benchmarks'
OUT = os.path.join(os.path.dirname(__file__), '..', 'results')
IU = {'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'R': 'AG', 'Y': 'CT', 'S': 'CG', 'W': 'AT', 'K': 'GT', 'M': 'AC',
      'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'V': 'ACG', 'N': 'ACGT'}
COMP = str.maketrans('ACGTRYSWKMBDHVN', 'TGCAYRSWMKVHDBN')


def bits(p):
    return sum(math.log2(4 / len(IU[c])) for c in p if c in IU)


def classify(call, truth):
    """Best alignment class of one call against one truth: 'strict' > 'precision' > 'fragment' > None."""
    c0 = str(call).replace('x', 'N').strip('N')
    if not c0 or any(ch not in IU for ch in c0):
        return None
    spec = [j for j, ch in enumerate(truth) if ch != 'N']
    need = math.ceil(0.8 * len(spec))
    fmin = 4 if len(spec) >= 8 else 3
    best = None
    rank = {None: 0, 'fragment': 1, 'precision': 2, 'strict': 3}
    for c in (c0, c0.translate(COMP)[::-1]):
        cspec = [i for i, ch in enumerate(c) if ch != 'N']
        for off in range(-len(truth) + 1, len(c)):
            sub = sup = conf = 0
            for j in spec:
                i = off + j
                if 0 <= i < len(c) and c[i] != 'N':
                    cs, ts = set(IU[c[i]]), set(IU[truth[j]])
                    if cs <= ts:
                        sub += 1
                    elif ts <= cs:
                        sup += 1
                    elif cs & ts:
                        sup += 1
                    else:
                        conf += 1
            inside = all(0 <= i - off < len(truth) for i in cspec)
            if sub >= need and sup + conf <= 1:
                cls = 'strict'
            elif sub + sup >= need and conf <= 1:
                cls = 'precision'
            elif conf == 0 and inside and sub + sup >= fmin:
                cls = 'fragment'
            else:
                cls = None
            if rank[cls] > rank[best]:
                best = cls
    return best


def shuffle_truth(t, seed):
    r = np.random.default_rng(seed)
    letters = [c for c in t if c != 'N']
    r.shuffle(letters)
    it = iter(letters)
    return ''.join(c if c == 'N' else next(it) for c in t)


def unit_class(calls, truths):
    best = None
    rank = {None: 0, 'fragment': 1, 'precision': 2, 'strict': 3}
    for c in calls:
        for t in truths:
            k = classify(c, t)
            if rank[k] > rank[best]:
                best = k
    return best or 'none'


def load(phase):
    if phase == 'ChIP':
        calls = pd.concat(pd.read_csv(f, sep='\t') for f in glob.glob(f'{B}/phase5_negctrl/results/parts/chip_calls/*.tsv'))
        mets = pd.concat(pd.read_csv(f, sep='\t') for f in glob.glob(f'{B}/phase5_negctrl/results/parts/chip_metrics/*.tsv'))
        calls = calls.rename(columns={'factor': 'unit'})
        mets = mets.rename(columns={'factor': 'unit'})
        mets['truths'] = mets.consensus.astype(str)
    elif phase == 'SELEX':
        calls = pd.concat(pd.read_csv(f, sep='\t') for f in glob.glob(f'{B}/phase4_selex/results/parts/selex_calls/*.tsv'))
        mets = pd.concat(pd.read_csv(f, sep='\t') for f in glob.glob(f'{B}/phase4_selex/results/parts/selex_metrics/*.tsv'))
        calls = calls.rename(columns={'lib': 'unit'})
        mets = mets.rename(columns={'lib': 'unit'})
    else:
        calls = pd.concat(pd.read_csv(f, sep='\t') for f in glob.glob(f'{B}/phase3_rbns/results/parts/rbns_calls/*.tsv'))
        mets = pd.concat(pd.read_csv(f, sep='\t') for f in glob.glob(f'{B}/phase3_rbns/results/parts/rbns_metrics/*.tsv'))
        calls = calls.rename(columns={'rbp': 'unit'})
        mets = mets.rename(columns={'rbp': 'unit'})
        mets['truths'] = mets.truth.astype(str)
    mets = mets.drop_duplicates(['unit', 'tool'])
    return calls, mets


rows = []
for phase in ['ChIP', 'SELEX', 'RBNS']:
    calls, mets = load(phase)
    tools = [t for t in ['cs_allinone', 'streme', 'homer'] if t in set(mets.tool)]
    truth_of = mets.drop_duplicates('unit').set_index('unit').truths.astype(str).str.split().to_dict()
    units = sorted(truth_of)
    for tool in tools:
        top = calls[(calls.tool == tool) & (calls['rank'] <= 10)].groupby('unit').pattern.apply(list).to_dict()
        r1 = calls[(calls.tool == tool) & (calls['rank'] == 1)].set_index('unit').pattern.to_dict()
        for ui, u in enumerate(units):
            if u not in top:
                continue
            got = unit_class(top[u], truth_of[u])
            chance = [unit_class(top[u], [shuffle_truth(t, ui * 10 + s) for t in truth_of[u]]) for s in range(5)]
            rows.append(dict(phase=phase, unit=u, tool=tool, cls=got,
                             chance_frag_or_better=np.mean([c != 'none' for c in chance]),
                             **{f'chance_{k}': np.mean([c == k for c in chance]) for k in ['strict', 'precision', 'fragment']},
                             rank1_lowinfo=bits(str(r1.get(u, 'N'))) < 6))
df = pd.DataFrame(rows)
df.to_csv(os.path.join(OUT, 'c3_failure_anatomy.tsv'), sep='\t', index=False)

# %% Summary
order = ['strict', 'precision', 'fragment', 'none']
for phase in ['ChIP', 'SELEX', 'RBNS']:
    d = df[df.phase == phase]
    print(f'\n== {phase} ({d.unit.nunique()} units): share of units by best top-10 class '
          f'(ch_* = same calls vs letter-shuffled truths)')
    tab = d.groupby('tool').cls.value_counts(normalize=True).unstack().reindex(columns=order).fillna(0)
    for k in ['strict', 'precision', 'fragment']:
        tab[f'ch_{k}'] = d.groupby('tool')[f'chance_{k}'].mean()
    tab['rank1_lowinfo'] = d.groupby('tool').rank1_lowinfo.mean()
    print(tab.round(2).to_string())
