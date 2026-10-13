# %% sig_hits.py: sensitivity at the native threshold (LN7). For each unit and HASTUR arm: does any SIGNIFICANT call
# match the truth (P5 / contains on the matrix-truth IUPAC; PWM similarity >= 0.85 to the true matrix)?
# Usage: python bench/sig_hits.py <truth_map.json> <out.tsv> <root> [<root> ...] --tools hlat_a3,hlat_p5
import json
import os
import re
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'eval'))
from score import as_iupac, hit_contains, hit_p5, pwm_sim  # noqa: E402

tm = json.load(open(sys.argv[1]))
out = sys.argv[2]
i = sys.argv.index('--tools')
roots, tools = sys.argv[3:i], sys.argv[i + 1].split(',')
rows = []
for root in roots:
    for ds in sorted(os.listdir(root)):
        f = ds.split('__')[0]
        if f not in tm:
            continue
        T, t = np.array(tm[f]['pwm']), tm[f]['truth']
        n = int(ds.split('__n')[1].split('__')[0]) if '__n' in ds else 5000
        for tool in tools:
            tsv, meme = os.path.join(root, ds, tool, 'hlat.tsv'), os.path.join(root, ds, tool, 'hlat.meme')
            if not os.path.exists(tsv):
                continue
            tab = pd.read_csv(tsv, sep='\t')
            mats = re.findall(r'letter-probability matrix:[^\n]*\n((?:\s*[0-9.eE+-]+\s+[0-9.eE+-]+\s+[0-9.eE+-]+\s+[0-9.eE+-]+\s*\n)+)', open(meme).read())
            pw = [np.array([[float(v) for v in r.split()] for r in x.strip().split('\n')]) for x in mats]
            sig = [m for m, g in zip(pw, tab.sig) if g == 1] if len(tab) else []
            rows.append(dict(unit=ds, n=n, tool=tool, n_sig=len(sig),
                             sig_p5=any(hit_p5(as_iupac(m), t) for m in sig),
                             sig_contains=any(hit_contains(as_iupac(m), t) for m in sig),
                             sig_sim85=any(pwm_sim(m, T) >= 0.85 for m in sig)))
d = pd.DataFrame(rows)
d.to_csv(out, sep='\t', index=False)
print(d.groupby(['n', 'tool'])[['n_sig', 'sig_p5', 'sig_contains', 'sig_sim85']].mean().round(3).to_string())
