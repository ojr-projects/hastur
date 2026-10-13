# %% promoter_nulls_make.py: random gene-set null units on a genome's promoters (E5 quick look, HASTUR first)
"""
Promoters: UP bp directly upstream of every protein-coding gene's start (strand-aware, forward-strand sequence),
clipped at the replicon ends, overlapping windows merged. Units: <root>/<tag>_n<N>__r<k>/ with fg = N random
promoters and bg = every other promoter, REPS per size. Seeds: random.Random(sha256('<tag>|<N>|<k>')).
Usage: python bench/promoter_nulls_make.py <genome.fna> <genomic.gff> <root> <tag> [UP=500] [sizes=25,50,100,200] [reps=20]
"""
import hashlib
import json
import os
import random
import sys


def main():
    fna, gff, root, tag = sys.argv[1:5]
    up = int(sys.argv[5]) if len(sys.argv) > 5 else 500
    sizes = [int(x) for x in (sys.argv[6] if len(sys.argv) > 6 else '25,50,100,200').split(',')]
    reps = int(sys.argv[7]) if len(sys.argv) > 7 else 20
    seqs, name = {}, None
    for line in open(fna):
        if line.startswith('>'):
            name = line[1:].split()[0]
            seqs[name] = []
        else:
            seqs[name].append(line.strip().upper())
    seqs = {k: ''.join(v) for k, v in seqs.items()}
    win = []
    for line in open(gff):
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'gene' or 'gene_biotype=protein_coding' not in f[8] or f[0] not in seqs:
            continue
        a, b = int(f[3]) - 1, int(f[4])
        lo, hi = (a - up, a) if f[6] == '+' else (b, b + up)
        lo, hi = max(0, lo), min(len(seqs[f[0]]), hi)
        if hi - lo >= 100:
            win.append((f[0], lo, hi))
    win.sort()
    regions = []
    for c, a, b in win:
        if regions and regions[-1][0] == c and a <= regions[-1][2]:
            regions[-1][2] = max(regions[-1][2], b)
        else:
            regions.append([c, a, b])
    rs = [seqs[c][a:b] for c, a, b in regions]
    for n in sizes:
        for k in range(reps):
            d = os.path.join(root, f'{tag}_n{n}__r{k}')
            if os.path.exists(os.path.join(d, 'truth.json')):
                continue
            os.makedirs(d, exist_ok=True)
            rng = random.Random(int(hashlib.sha256(f'{tag}|{n}|{k}'.encode()).hexdigest()[:16], 16))
            fg = set(rng.sample(range(len(rs)), n))
            with open(os.path.join(d, 'fg.fa'), 'w') as fh, open(os.path.join(d, 'fg.bed'), 'w') as fb:
                for j in sorted(fg):
                    fh.write(f'>p{j}\n{rs[j]}\n')
                    fb.write('{}\t{}\t{}\tp{}\n'.format(*regions[j], j))
            with open(os.path.join(d, 'bg.fa'), 'w') as fh:
                for j in range(len(rs)):
                    if j not in fg:
                        fh.write(f'>p{j}\n{rs[j]}\n')
            json.dump(dict(dataset=os.path.basename(d), truth=[], null=True, kind=f'{tag}_n{n}'),
                      open(os.path.join(d, 'truth.json'), 'w'))
    gc = sum(s.count('G') + s.count('C') for s in rs) / sum(len(s) for s in rs)
    print(f'{len(win)} genes -> {len(regions)} promoter regions, {sum(len(s) for s in rs) / 1e6:.2f} Mb, GC {gc:.3f}; '
          f'{len(sizes) * reps} null units')


if __name__ == '__main__':
    main()
