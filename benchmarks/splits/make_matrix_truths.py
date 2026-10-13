# %% make_matrix_truths.py: matrix truths for ChIP dev and a fresh test split (test2), with the oracle test
"""
Mechanical rules (fixed 2026-10-11, before any test2 data are built or any tool is run on them):
  dev      for each of the 24 ChIP dev factors: the first HOMER vertebrate known motif whose name (before '('), split on
           '-', '|', ':', has a token equal to the factor (case-insensitive, punctuation removed), preferring sources
           from the factor's cell line; truth =
           the matrix, consensus = eval/score.pwm_to_iupac(matrix).
  test2    every HOMER vertebrate known motif whose source is a ChIP-seq in GM12878, HEK293, HeLa(-S3), MCF-7 or A549,
           whose TF is in neither the dev nor the test split; one per TF (first in file order). ENCODE biosample term:
           GM12878, HEK293, HeLa-S3, MCF-7, A549. Factors without ENCODE peaks are reported missing, never replaced.
Writes splits/chip_dev_matrix_truth.json, splits/chip_test2_split.tsv, splits/chip_test2_matrix_truth.json and prints
the oracle test (each truth's own matrix scored as a call must pass contains, P5 and PWM similarity >= 0.85).
"""
import json
import os
import re
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'eval'))
from score import hit_contains, hit_p5, pwm_sim, pwm_to_iupac  # noqa: E402

LIB = os.path.expanduser('~/micromamba/envs/homer/share/homer/data/knownTFs/vertebrates/known.motifs')
entries, name, rows = [], None, []
for line in open(LIB):
    if line.startswith('>'):
        if name:
            entries.append((name, np.array(rows)))
        name, rows = line.split('\t')[1], []
    elif line.strip():
        rows.append([float(v) for v in line.split()])
entries.append((name, np.array(rows)))
norm = lambda s: re.sub(r'[^A-Z0-9]', '', s.upper())
sp = pd.read_csv(os.path.join(HERE, 'chip_split.tsv'), sep='\t')

# dev
dev = {}
for r in sp[sp.split == 'dev'].itertuples():
    cands = [(n, m) for n, m in entries if norm(r.factor) in {norm(t) for t in re.split(r'[-|:]', n.split('(')[0])}]
    cands.sort(key=lambda t: 0 if r.cell.upper().replace('-', '') in t[0].upper().replace('-', '') else 1)
    if cands:
        n, m = cands[0]
        dev[r.factor] = dict(source=n, truth=pwm_to_iupac(m), pwm=m.tolist())
# test2
cellmap = {'GM12878': 'GM12878', 'HEK293': 'HEK293', 'HELAS3': 'HeLa-S3', 'HELA': 'HeLa-S3', 'MCF7': 'MCF-7', 'A549': 'A549'}
used = {norm(f) for f in sp.factor}
t2, seen = [], set()
for n, m in entries:
    src = n.split('/')[1] if '/' in n else ''
    if 'ChIP' not in src:
        continue
    cell = next((v for k, v in cellmap.items() if src.upper().replace('-', '').startswith(k)), None)
    if cell is None:
        continue
    tf = n.split('(')[0].split('/')[0]
    if ':' in tf or ',' in n.split('/')[0].split(')')[-1]:      # skip composite / condensed variants
        continue
    if norm(tf) in used or norm(tf) in seen:
        continue
    seen.add(norm(tf))
    # amendment 1 (2026-10-11, before any test2 tool run): the ENCODE target label is the ChIP target named in the
    # HOMER source string (<cell>-<TARGET>-ChIP-Seq), tags such as .GFP/.eGFP/.Flag/.V5 removed; HOMER motif names
    # (e.g. AP-2gamma for TFAP2C) are not ENCODE labels
    parts = src.split('-')
    target = parts[1].split('.')[0].upper() if len(parts) > 2 else tf.upper()
    t2.append(dict(factor=target, cell=cell, truth=pwm_to_iupac(m), split='test2', source=n, pwm=m.tolist()))
pd.DataFrame([{k: v for k, v in d.items() if k != 'pwm'} for d in t2]).to_csv(os.path.join(HERE, 'chip_test2_split.tsv'), sep='\t', index=False)
json.dump(dev, open(os.path.join(HERE, 'chip_dev_matrix_truth.json'), 'w'))
json.dump({d['factor']: dict(source=d['source'], truth=d['truth'], pwm=d['pwm']) for d in t2},
          open(os.path.join(HERE, 'chip_test2_matrix_truth.json'), 'w'))
for label, items in (('dev', dev), ('test2', {d['factor']: d for d in t2})):
    ok = [(hit_contains(v['truth'], v['truth']), hit_p5(v['truth'], v['truth']),
           pwm_sim(np.array(v['pwm']), np.array(v['pwm'])) >= 0.85) for v in items.values()]
    ok = np.array(ok)
    print(f'{label}: {len(items)} factors; oracle pass rates contains {ok[:, 0].mean():.2f}, p5 {ok[:, 1].mean():.2f}, '
          f'sim {ok[:, 2].mean():.2f}')
missing_dev = sorted(set(sp[sp.split == 'dev'].factor) - set(dev))
print('dev factors without a HOMER matrix:', missing_dev)
print('test2 factors:', ' '.join(f"{d['factor']}/{d['cell']}" for d in t2))
