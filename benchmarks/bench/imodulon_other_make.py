# %% imodulon_other_make.py: iModulon screen units for an organism without TSS maps (E1b quick look, HASTUR first)
"""
Promoter regions: operon heads. Genes (gff3 'gene') on one strand form an operon while the gap to the next gene is
<= OPGAP bp. Each operon's region is [start - UP, start + DOWN) around its first gene's start (strand-aware, forward-
strand sequence); overlapping regions merged.
Module membership: the module's imodulon_size genes with the largest |M| (the published thresholds are not
distributed; size-matching reproduces a threshold-defined membership).
Units: <root>/<module>/ (fg = the module's operon regions, bg = every other region) and 20 nulls per module (random
regions of the same count). truth.json carries the annotated regulator only (no site truth: quick look).
Usage: python bench/imodulon_other_make.py <supp.xlsx> <genome.fasta> <genome.gff3> <root> [tag]
"""
import hashlib
import json
import os
import random
import re
import sys

import numpy as np
import pandas as pd

UP, DOWN, OPGAP, MIN_REG, NULL_PER = 300, 50, 40, 5, 20


def safe(name):
    return re.sub(r'[^A-Za-z0-9_.-]+', '_', str(name)).strip('_')


def main():
    xlsx, fa, gff, root = sys.argv[1:5]
    tag = sys.argv[5] if len(sys.argv) > 5 else 'E1b'
    seqs, name = {}, None
    for line in open(fa):
        if line.startswith('>'):
            name = line[1:].split()[0]
            seqs[name] = []
        else:
            seqs[name].append(line.strip().upper())
    seqs = {k: ''.join(v) for k, v in seqs.items()}
    genes = []
    for line in open(gff):
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'gene' or f[0] not in seqs:
            continue
        a = dict(kv.split('=', 1) for kv in f[8].split(';') if '=' in kv)
        lt = a.get('locus_tag') or a.get('ID', '').replace('gene-', '')
        genes.append((f[0], int(f[3]) - 1, int(f[4]), f[6], lt))
    genes.sort()
    # operons
    op_of, heads = {}, []
    for c in seqs:
        for s in '+-':
            g = [x for x in genes if x[0] == c and x[3] == s]
            if s == '-':
                g = g[::-1]
            cur = None
            for x in g:
                gap = (x[1] - cur[2]) if (cur and s == '+') else ((cur[1] - x[2]) if cur else None)
                if cur is None or gap > OPGAP:
                    heads.append(x)
                op_of[x[4]] = len(heads) - 1
                cur = x
    win = []
    for c, a, b, s, lt in heads:
        lo, hi = (a - UP, a + DOWN) if s == '+' else (b - DOWN, b + UP)
        win.append((c, max(0, lo), min(len(seqs[c]), hi)))
    order = sorted(range(len(win)), key=lambda i: win[i])
    regions, head_reg = [], {}
    for i in order:
        c, a, b = win[i]
        if regions and regions[-1][0] == c and a <= regions[-1][2]:
            regions[-1][2] = max(regions[-1][2], b)
        else:
            regions.append([c, a, b])
        head_reg[i] = len(regions) - 1
    rseq = [seqs[c][a:b] for c, a, b in regions]
    x = pd.ExcelFile(xlsx)
    M = pd.read_excel(x, '6-M', index_col=0)
    imt = pd.read_excel(x, '4-iM table', index_col=0)
    if len(M.columns) != len(imt):
        raise ValueError('M and iM table disagree')
    M.columns = list(imt.index)
    os.makedirs(root, exist_ok=True)
    meta = []
    for m in imt.index:
        size = int(imt.loc[m, 'imodulon_size']) if pd.notna(imt.loc[m, 'imodulon_size']) else 0
        top = M[m].abs().sort_values(ascending=False).index[:size]
        pos = sorted({head_reg[op_of[g]] for g in top if g in op_of})
        reg = imt.loc[m, 'regulator']
        meta.append(dict(module=m, safe=safe(m), regulator=reg if isinstance(reg, str) else '', n_genes=size,
                         n_regions=len(pos), category=imt.loc[m, 'Category'] if 'Category' in imt.columns else '',
                         included=len(pos) >= MIN_REG))
        if len(pos) < MIN_REG:
            continue
        units = [(safe(m), pos, dict(module=m, regulator=meta[-1]['regulator'], truth=[]))]
        rng = random.Random(int(hashlib.sha256(f'{tag}|{m}'.encode()).hexdigest()[:16], 16))
        for k in range(NULL_PER):
            units.append((f'null_{safe(m)}__r{k}', sorted(rng.sample(range(len(regions)), len(pos))),
                          dict(truth=[], null=True, kind=f'{tag}_null', module=m)))
        for u, fg_idx, t in units:
            d = os.path.join(root, u)
            os.makedirs(d, exist_ok=True)
            fs = set(fg_idx)
            with open(os.path.join(d, 'fg.fa'), 'w') as fh:
                for j in fg_idx:
                    fh.write(f'>r{j}\n{rseq[j]}\n')
            with open(os.path.join(d, 'bg.fa'), 'w') as fh:
                for j in range(len(regions)):
                    if j not in fs:
                        fh.write(f'>r{j}\n{rseq[j]}\n')
            json.dump(dict(dataset=u, **t), open(os.path.join(d, 'truth.json'), 'w'))
    df = pd.DataFrame(meta)
    df.to_csv(os.path.join(root, 'modules.tsv'), sep='\t', index=False)
    gc = np.mean([(s.count('G') + s.count('C')) / len(s) for s in rseq])
    print(f'{len(genes)} genes, {len(heads)} operons, {len(regions)} promoter regions (GC {gc:.2f}); '
          f'{len(imt)} iModulons, {int(df.included.sum())} with >= {MIN_REG} regions')


if __name__ == '__main__':
    main()
