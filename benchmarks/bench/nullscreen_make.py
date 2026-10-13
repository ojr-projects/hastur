# %% nullscreen_make.py: E6 fully null screens on the E1 promoter population (splits/protocol_E6.md)
"""
A null screen mirrors the E1 screen: one random set of promoter regions per included E1 module, at that module's
size (fg = the set, bg = every other region). Nothing in a null screen is a real module, so every call is false.
Regions are the 1,942 E1/D1 promoter regions, read back from any E1 unit (its fg + bg hold every region once).
Seeds: random.Random of sha256('E6|<screen>|<module>'), disjoint from E1 ('E1|<module>') and D1 ('<TF>').
Usage: python bench/nullscreen_make.py <e1_units_root> <out_root> <n_screens>
"""
import hashlib
import os
import random
import sys

import pandas as pd


def read_fa(p):
    out, name = {}, None
    for line in open(p):
        if line.startswith('>'):
            name = line[1:].strip()
        else:
            out[name] = line.strip()
    return out


def main():
    e1, root, n = sys.argv[1], sys.argv[2], int(sys.argv[3])
    mods = pd.read_csv(os.path.join(e1, 'modules.tsv'), sep='\t')
    mods = mods[mods.included]
    u0 = os.path.join(e1, mods.safe.iloc[0])
    regions = {**read_fa(os.path.join(u0, 'fg.fa')), **read_fa(os.path.join(u0, 'bg.fa'))}
    ids = sorted(regions, key=lambda r: int(r[1:]))
    assert len(ids) == 1942, len(ids)
    for s in range(n):
        for m in mods.itertuples():
            d = os.path.join(root, f's{s:03d}__{m.safe}')
            if os.path.exists(os.path.join(d, 'bg.fa')):
                continue
            os.makedirs(d, exist_ok=True)
            rng = random.Random(int(hashlib.sha256(f'E6|{s}|{m.module}'.encode()).hexdigest()[:16], 16))
            fg = sorted(rng.sample(range(len(ids)), m.n_regions))
            fs = set(fg)
            with open(os.path.join(d, 'fg.fa'), 'w') as fh:
                for j in fg:
                    fh.write(f'>{ids[j]}\n{regions[ids[j]]}\n')
            with open(os.path.join(d, 'bg.fa'), 'w') as fh:
                for j in range(len(ids)):
                    if j not in fs:
                        fh.write(f'>{ids[j]}\n{regions[ids[j]]}\n')
            with open(os.path.join(d, 'truth.json'), 'w') as fh:
                fh.write(f'{{"dataset": "{os.path.basename(d)}", "truth": [], "null": true, "kind": "null_screen", '
                         f'"screen": {s}, "module": "{m.module}"}}\n')
    print(f'{n} null screens x {len(mods)} sets -> {root}')


if __name__ == '__main__':
    main()
