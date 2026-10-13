# %% null_eval.py: false significant calls per run on motif-free pairs (bench/null_make.py)
"""
Each tool's own significance criterion: HASTUR sig = 1 (the Yellow Sign); Comsequent all-in-one sig = 1;
STREME E-value <= 0.05; HOMER p <= 1e-12 (HOMER marks weaker motifs as possible false positives).
Usage: python bench/null_eval.py <root> <out.tsv>
"""
import json
import math
import os
import re
import sys

import pandas as pd


def tools_in(d):
    """Tool output dirs in a unit dir, discovered rather than listed (new arms were silently skipped, 2026-10-11)."""
    return sorted(t for t in os.listdir(d) if os.path.isdir(os.path.join(d, t)) and
                  (t.startswith('hlat') or t.startswith('streme') or t in ('homer', 'comse')))


def n_sig(d, tool):
    o = os.path.join(d, tool)
    if not os.path.exists(os.path.join(o, 'runtime_s')):
        return None
    if tool.startswith('hlat'):
        p = os.path.join(o, 'hlat.tsv')
        return int((pd.read_csv(p, sep='\t').sig == 1).sum()) if os.path.exists(p) else 0
    if tool == 'comse':
        p = os.path.join(o, 'cs.aio.tsv.all')
        return int((pd.read_csv(p, sep='\t').sig == 1).sum()) if os.path.exists(p) else 0
    if tool.startswith('streme'):
        p = os.path.join(o, 'streme.txt')
        if not os.path.exists(p):
            return 0
        return sum(float(e) <= 0.05 for e in re.findall(r'letter-probability matrix:.*?(?:E|S)= ([0-9.eE+-]+)', open(p).read()))
    if tool == 'homer':
        p = os.path.join(o, 'homerMotifs.all.motifs')
        if not os.path.exists(p):
            return 0
        lps = [float(l.split('\t')[3]) for l in open(p) if l.startswith('>')]
        return sum(lp <= math.log(1e-12) for lp in lps)
    raise ValueError(tool)


if __name__ == '__main__':
    root, out = sys.argv[1], sys.argv[2]
    rows = []
    for ds in sorted(os.listdir(root)):
        d = os.path.join(root, ds)
        if not os.path.exists(os.path.join(d, 'truth.json')):
            continue
        kind = json.load(open(os.path.join(d, 'truth.json')))['kind']
        for tool in tools_in(d):   # every tool output dir present (hlat*, streme*, homer, comse)
            n = n_sig(d, tool)
            if n is not None:
                rows.append(dict(dataset=ds, kind=kind, tool=tool, n_sig=n))
    df = pd.DataFrame(rows)
    df.to_csv(out, sep='\t', index=False)
    if len(df):
        print(df.pivot_table(index='kind', columns='tool', values='n_sig', aggfunc='mean').round(2).to_string())
        print(df.assign(any=df.n_sig > 0).pivot_table(index='kind', columns='tool', values='any', aggfunc='mean').round(2).to_string())
