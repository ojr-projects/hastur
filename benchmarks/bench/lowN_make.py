# %% lowN_make.py: registered low-N ChIP panel LN1 (splits/protocol_lowN.md)
"""
For each ChIP dev factor and N in {250, 500, 1000}, replicates 0-2: N foreground peaks drawn uniformly without
replacement from the 5,000 (summit-centred), 10 N background sequences from the 10x pool, both with
numpy default_rng(1000 N + r + factor index). Writes <root>/<factor>__n<N>__r<r>/{fg.fa, bg.fa, truth.json}.
Usage: python bench/lowN_make.py <chip_root> <out_root> [dev|test|test2]   (factor index = position within that split)
       env CHIP_SPLIT=splits/chip_test2_split.tsv for test2 (a factor with fewer peaks than N uses all of them)
"""
import json
import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
NS, REPS, BGX = tuple(int(x) for x in os.environ.get('LOWN_NS', '250,500,1000').split(',')), 3, 10   # LN6: LOWN_NS=25,50,100


def read_fa(p):
    out = []
    for line in open(p):
        if not line.startswith('>'):
            out.append(line.strip())
    return out


if __name__ == '__main__':
    chip, root = sys.argv[1], sys.argv[2]
    which = sys.argv[3] if len(sys.argv) > 3 else 'dev'
    import numpy as np
    sp = pd.read_csv(os.environ.get('CHIP_SPLIT', os.path.join(HERE, '..', 'splits', 'chip_split.tsv')), sep='\t')   # test2: chip_test2_split.tsv
    dev = sp[sp.split == which].reset_index(drop=True)
    for fi, r_ in dev.iterrows():
        d0 = os.path.join(chip, r_.factor)
        if not os.path.exists(os.path.join(d0, 'bg.fa')):
            print('missing', r_.factor)
            continue
        fg, bg = read_fa(os.path.join(d0, 'fg.fa')), read_fa(os.path.join(d0, 'bg.fa'))
        for n in NS:
            for rep in range(REPS):
                d = os.path.join(root, f'{r_.factor}__n{n}__r{rep}')
                if os.path.exists(os.path.join(d, 'truth.json')):
                    continue
                os.makedirs(d, exist_ok=True)
                rng = np.random.default_rng(1000 * n + rep + fi)
                fsel = rng.choice(len(fg), size=min(n, len(fg)), replace=False)
                bsel = rng.choice(len(bg), size=min(BGX * n, len(bg)), replace=False)
                for name, seqs, sel in (('fg.fa', fg, fsel), ('bg.fa', bg, bsel)):
                    with open(os.path.join(d, name), 'w') as fh:
                        for i, j in enumerate(sel):
                            fh.write(f'>s{i}\n{seqs[j]}\n')
                json.dump(dict(dataset=os.path.basename(d), factor=r_.factor, n=n, rep=rep, truth=[r_.truth]),
                          open(os.path.join(d, 'truth.json'), 'w'))
        print(r_.factor, flush=True)
