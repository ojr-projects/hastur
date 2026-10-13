# %% make_test3.py: test3 ChIP split for Mode A at scale (MA2), mechanical, fixed before any test3 data are built
"""
Rule (2026-10-12):
  - Candidates: every HOMER vertebrate known motif whose source is a ChIP-seq, skipping composite variants. The target
    is the ChIP target named in the source string (<cell>-<TARGET>-ChIP-Seq; tags such as .GFP/.Flag/.V5 removed),
    as in test2 amendment 1. One motif per target, the first in file order.
  - Excluded: targets in the dev or test2 splits.
  - Cell: the first of K562, HepG2, HEK293, GM12878, MCF-7, A549, HeLa-S3 with released GRCh38 conservative IDR
    thresholded peaks for the target on ENCODE. A target with none is reported missing, never replaced.
  - Truth: the HOMER matrix; consensus = eval/score.pwm_to_iupac(matrix).
Writes splits/chip_test3_split.tsv and splits/chip_test3_matrix_truth.json (accession recorded per factor).
"""
import json
import os
import re
import sys
import urllib.error
import urllib.parse
import urllib.request

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'eval'))
from score import pwm_to_iupac  # noqa: E402

LIB = os.path.expanduser('~/micromamba/envs/homer/share/homer/data/knownTFs/vertebrates/known.motifs')
CELLS = ['K562', 'HepG2', 'HEK293', 'GM12878', 'MCF-7', 'A549', 'HeLa-S3']
norm = lambda s: re.sub(r'[^A-Z0-9]', '', s.upper())  # noqa: E731


def encode_peaks(target, cell):
    q = ('/search/?type=File&output_type=conservative+IDR+thresholded+peaks&assembly=GRCh38&file_format=bed'
         f'&status=released&target.label={urllib.parse.quote(target)}'
         f'&biosample_ontology.term_name={urllib.parse.quote(cell)}&limit=all&frame=object')
    req = urllib.request.Request('https://www.encodeproject.org' + q, headers={'Accept': 'application/json'})
    try:
        with urllib.request.urlopen(req, timeout=120) as r:
            g = json.load(r)['@graph']
    except urllib.error.HTTPError as e:
        if e.code == 404:
            return None
        raise
    if not g:
        return None
    return sorted(g, key=lambda x: x.get('date_created', ''))[-1]['accession']


entries, name, rows = [], None, []
for line in open(LIB):
    if line.startswith('>'):
        if name:
            entries.append((name, np.array(rows)))
        name, rows = line.split('\t')[1], []
    elif line.strip():
        rows.append([float(v) for v in line.split()])
entries.append((name, np.array(rows)))
used = {norm(f) for f in pd.read_csv(os.path.join(HERE, 'chip_split.tsv'), sep='\t').factor}
used |= {norm(f) for f in pd.read_csv(os.path.join(HERE, 'chip_test2_split.tsv'), sep='\t').factor}
out, seen, missing = [], set(), []
for n, m in entries:
    src = n.split('/')[1] if '/' in n else ''
    if 'ChIP' not in src:
        continue
    tf = n.split('(')[0].split('/')[0]
    if ':' in tf or ',' in n.split('/')[0].split(')')[-1]:
        continue
    parts = src.split('-')
    target = parts[1].split('.')[0].upper() if len(parts) > 2 else tf.upper()
    if norm(target) in used or norm(target) in seen:
        continue
    seen.add(norm(target))
    hit = None
    for cell in CELLS:
        acc = encode_peaks(target, cell)
        if acc:
            hit = (cell, acc)
            break
    if hit is None:
        missing.append(target)
        continue
    out.append(dict(factor=target, cell=hit[0], peaks=hit[1], truth=pwm_to_iupac(m), split='test3', source=n,
                    pwm=m.tolist()))
    print(target, hit, flush=True)
pd.DataFrame([{k: v for k, v in d.items() if k != 'pwm'} for d in out]).to_csv(
    os.path.join(HERE, 'chip_test3_split.tsv'), sep='\t', index=False)
json.dump({d['factor']: dict(source=d['source'], truth=d['truth'], pwm=d['pwm'], cell=d['cell'], peaks=d['peaks'])
           for d in out}, open(os.path.join(HERE, 'chip_test3_matrix_truth.json'), 'w'))
print(f'test3: {len(out)} factors; missing ENCODE peaks: {len(missing)} ({" ".join(missing)})')
