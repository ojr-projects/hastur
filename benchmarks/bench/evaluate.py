# %% evaluate.py: score every tool's ranked calls in every dataset with eval/score.py
"""
Usage: python bench/evaluate.py <root> <out.tsv> [truth_source]      env TRUTH_MAP=<json factor -> {truth, pwm}> overrides
  truth_source: 'json' (default; <dataset>/truth.json) or a splits TSV path with columns factor, truth (ChIP), in which
  case the dataset directory name is the factor.
One row per dataset x tool: first ranks and top-1/top-10 hits under each rule (p5, contains, pwm_sim), best PWM
similarity, the shuffled-truth chance of each, n_calls, runtime.
"""
import json
import os
import re
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'eval'))
from score import pwm_to_iupac, score_unit  # noqa: E402

TOPK = 10


def tools_in(d):
    """Tool output dirs in a unit dir, discovered rather than listed (new arms were silently skipped, 2026-10-11)."""
    return sorted(t for t in os.listdir(d) if os.path.isdir(os.path.join(d, t)) and
                  (t.startswith('hlat') or t.startswith('streme') or t in ('homer', 'comse')))


def read_meme(path):
    if not os.path.exists(path):
        return []
    txt = open(path).read()
    mats = re.findall(r'letter-probability matrix:[^\n]*\n((?:\s*[0-9.eE+-]+\s+[0-9.eE+-]+\s+[0-9.eE+-]+\s+[0-9.eE+-]+\s*\n)+)', txt)
    return [np.array([[float(v) for v in r.split()] for r in m.strip().split('\n')]) for m in mats]


def read_homer(path):
    if not os.path.exists(path):
        return []
    out, rows = [], None
    for line in open(path):
        if line.startswith('>'):
            if rows:
                out.append(np.array(rows))
            rows = []
        elif line.strip() and rows is not None:
            rows.append([float(v) for v in line.split()])
    if rows:
        out.append(np.array(rows))
    return out


def read_comse(path):
    if not os.path.exists(path):
        return []
    a = pd.read_csv(path, sep='\t', keep_default_na=False)
    return [str(p) for p in a.pattern]


def calls_for(d, tool):
    o = os.path.join(d, tool)
    if tool.startswith('hlat'):
        return read_meme(os.path.join(o, 'hlat.meme'))
    if tool in ('streme', 'streme_all'):
        return read_meme(os.path.join(o, 'streme.txt'))
    if tool == 'homer':
        return read_homer(os.path.join(o, 'homerMotifs.all.motifs'))
    if tool == 'comse':
        return read_comse(os.path.join(o, 'cs.aio.tsv.all'))
    raise ValueError(tool)


if __name__ == '__main__':
    root, out = sys.argv[1], sys.argv[2]
    tsrc = sys.argv[3] if len(sys.argv) > 3 else 'json'
    truths = {}
    tmap = json.load(open(os.environ['TRUTH_MAP'])) if os.environ.get('TRUTH_MAP') else None   # factor -> {truth, pwm}
    if tsrc != 'json' and tmap is None:
        sp = pd.read_csv(tsrc, sep='\t')
        truths = {r.factor: [r.truth] for r in sp.itertuples()}
    rows = []
    for ds in sorted(os.listdir(root)):
        d = os.path.join(root, ds)
        if not os.path.exists(os.path.join(d, 'fg.fa')):
            continue
        if tsrc == 'json' and tmap is None and not os.path.exists(os.path.join(d, 'truth.json')):
            print(f'skipping {ds}: incomplete unit (no truth.json)', file=sys.stderr)
            continue
        if tmap is not None:      # matrix truths by factor (dataset name prefix before '__'): oracle-valid readout
            f = ds.split('__')[0]
            meta = {'dataset': ds, 'truth': [tmap[f]['truth']], 'pwms': [tmap[f]['pwm']]} if f in tmap else {'truth': None}
        else:
            meta = json.load(open(os.path.join(d, 'truth.json'))) if tsrc == 'json' else {'dataset': ds, 'truth': truths.get(ds)}
        if not meta.get('truth'):
            continue
        for tool in tools_in(d):   # every tool output dir present (hlat*, streme*, homer, comse)
            rt = os.path.join(d, tool, 'runtime_s')
            if not os.path.exists(rt):
                continue
            calls = calls_for(d, tool)
            vj = os.path.join(d, tool, 'hlat.json')
            version = json.load(open(vj)).get('version', 'v3') if tool.startswith('hlat') and os.path.exists(vj) else ''
            tp = meta.get('pwms') or ([meta['pwm']] if meta.get('pwm') else None)
            r = score_unit(calls, meta['truth'], topk=TOPK, truth_pwms=tp)
            top3 = [c if isinstance(c, str) else pwm_to_iupac(c) for c in calls[:3]]
            r.update({k: v for k, v in meta.items() if k not in ('pwm', 'pwms', 'truth')})
            r.update(dataset=ds, tool=tool, truth=' '.join(meta['truth']), n_calls=len(calls), top3=' '.join(top3), version=version,
                     runtime_s=float(open(rt).read()))
            rows.append(r)
    df = pd.DataFrame(rows)
    df.to_csv(out, sep='\t', index=False)
    if len(df):
        cell = [c for c in ('motif', 'n', 'pi', 'background') if c in df] or ['dataset']
        s = df.groupby(cell + ['tool']).agg(n=('dataset', 'size'), contains10=('hit10_contains', 'mean'),
                                           contains1=('hit1_contains', 'mean'), sim10=('hit10_pwm_sim', 'mean'),
                                           p5_10=('hit10_p5', 'mean'), best_sim=('best_sim', 'mean'),
                                           ch_contains=('chance10_contains', 'mean'), ch_sim=('chance10_pwm_sim', 'mean'),
                                           runtime=('runtime_s', 'median'))
        pd.set_option('display.width', 250)
        print(s.round(2).to_string())
        sim1 = df.groupby(cell + ['tool']).hit1_pwm_sim.mean()
        print('\n== overall (mean over cells): top-1 / top-10 by PWM similarity to the true matrix, and by contains')
        tot = df.groupby('tool').agg(sim1=('hit1_pwm_sim', 'mean'), sim10=('hit10_pwm_sim', 'mean'),
                                     contains1=('hit1_contains', 'mean'), contains10=('hit10_contains', 'mean'),
                                     best_sim=('best_sim', 'mean'), chance_sim10=('chance10_pwm_sim', 'mean'),
                                     runtime=('runtime_s', 'median'))
        print(tot.round(3).to_string())
