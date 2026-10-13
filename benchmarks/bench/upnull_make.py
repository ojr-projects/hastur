# %% upnull_make.py: small-N motif-free panel from bacterial upstream regions (LN6 dev; D1-independent genomes)
"""
Regions: [gene start - 250, gene start + 50] (strand-aware) around every annotated gene start, overlapping windows
merged (the D1 construction, with gene starts in place of TSSs). Runs: fg = random regions of size n, bg = every
other region; n in SIZES, REPS each, per genome. truth.json: null=True, kind=<genome>_n<n>.
Usage: python bench/upnull_make.py <root> [seed_offset]
"""
import json
import os
import sys

import numpy as np

GEN = '/home/oliver/code/comsequent/data/genomes'
GENOMES = {'pputida': 'Pseudomonas_putida_strain_KT2440', 'cglut': 'Corynebacterium_glutamicum_ATCC_13032',
           'saureus': 'Staphylococcus_aureus_subsp._aureus_MW2'}
SIZES, REPS, UP, DOWN = (10, 15, 20, 30, 50, 100, 200), 10, 250, 50


def fasta(path):
    seqs, name, buf = {}, None, []
    for line in open(path):
        if line.startswith('>'):
            if name:
                seqs[name] = ''.join(buf).upper()
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name:
        seqs[name] = ''.join(buf).upper()
    return seqs


def regions(gdir):
    d = os.path.join(GEN, gdir, 'ncbi_dataset', 'data')
    acc = [x for x in os.listdir(d) if x.startswith('GC')][0]
    base = os.path.join(d, acc)
    S = fasta(os.path.join(base, [f for f in os.listdir(base) if f.endswith('.fna')][0]))
    win = {c: [] for c in S}
    for line in open(os.path.join(base, 'genomic.gff')):
        f = line.rstrip('\n').split('\t')
        if line.startswith('#') or len(f) < 9 or f[2] != 'gene' or f[0] not in S:
            continue
        a, b = int(f[3]) - 1, int(f[4])
        win[f[0]].append((max(0, a - UP), a + DOWN) if f[6] == '+' else (max(0, b - DOWN), min(len(S[f[0]]), b + UP)))
    out = []
    for c, ws in win.items():
        cur = None
        for a, b in sorted(ws):
            if cur and a <= cur[1]:
                cur[1] = max(cur[1], b)
            else:
                if cur:
                    out.append(S[c][cur[0]:cur[1]])
                cur = [a, b]
        if cur:
            out.append(S[c][cur[0]:cur[1]])
    return [x for x in out if len(x) >= 100 and set(x) <= set('ACGT')]


if __name__ == '__main__':
    root = sys.argv[1]
    off = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    for g, gdir in GENOMES.items():
        R = regions(gdir)
        print(g, len(R), 'regions, median', int(np.median([len(x) for x in R])), 'bp', flush=True)
        for n in SIZES:
            for rep in range(REPS):
                d = os.path.join(root, f'{g}_n{n}__r{rep}')
                if os.path.exists(os.path.join(d, 'truth.json')):
                    continue
                rng = np.random.default_rng([660000, n, rep, off, 97 * len(g)])   # tuple seed: offsets never collide (2026-10-11 fix; earlier panels used an additive seed whose consecutive offsets shared 9/10 replicates)
                fs = set(rng.choice(len(R), size=n, replace=False).tolist())
                os.makedirs(d, exist_ok=True)
                with open(os.path.join(d, 'fg.fa'), 'w') as fh:
                    fh.writelines(f'>r{j}\n{R[j]}\n' for j in sorted(fs))
                with open(os.path.join(d, 'bg.fa'), 'w') as fh:
                    fh.writelines(f'>r{j}\n{R[j]}\n' for j in range(len(R)) if j not in fs)
                json.dump(dict(dataset=os.path.basename(d), kind=f'{g}_n{n}', null=True, truth=[]),
                          open(os.path.join(d, 'truth.json'), 'w'))
