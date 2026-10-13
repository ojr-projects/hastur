# %% regulon_make.py: D1 E. coli regulon screen units from RegulonDB (splits/protocol_D1.md)
"""
Promoter regions: [TSS - 250, TSS + 50] (strand-aware) around every RegulonDB TSS, overlapping windows merged into one
region (forward-strand sequence; every tool searches both strands). A TF's regulon = the regions containing at least
half of one of its annotated binding sites (RISet, riType tf-*). Units:
  <root>/<TF>/            fg = its regulon regions, bg = every other region; truth.json: truth matrix from its sites
                          at their modal length (strand-oriented, +0.25 pseudocount), IUPAC of it, and the site
                          intervals in fg coordinates (positional readout)
  <root>/null_<TF>__r<k>/ fg = random regions, same count as the regulon; bg = every other region (k < NULL_PER)
Included TFs: >= MIN_REG regions and >= MIN_SITES sites at the modal length.
Usage: python bench/regulon_make.py <regulondb_raw_dir> <root>
"""
import collections
import hashlib
import json
import os
import random
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'eval'))
from score import pwm_to_iupac  # noqa: E402

GENOME = ('/home/oliver/code/comsequent/data/genomes/Escherichia_coli_str_K-12_substr._MG1655/ncbi_dataset/data/'
          'GCF_000005845.2/GCF_000005845.2_ASM584v2_genomic.fna')
UP, DOWN = 250, 50
MIN_REG, MIN_SITES, NULL_PER = 10, 8, 10
COMP = str.maketrans('ACGT', 'TGCA')


def gff(path):
    for line in open(path):
        if line.startswith('#') or not line.strip():
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9:
            continue
        a = dict(kv.split('=', 1) for kv in f[8].replace('; ', ';').split(';') if '=' in kv)
        yield f, a


def main():
    raw, root = sys.argv[1], sys.argv[2]
    g = ''.join(l.strip() for l in open(GENOME) if not l.startswith('>')).upper()
    tss = [(int(f[3]), f[6]) for f, a in gff(os.path.join(raw, 'PromoterSet.gff3')) if f[2] == 'transcription_start_site']
    win = sorted((max(0, p - 1 - UP), min(len(g), p + DOWN)) if s == '+' else (max(0, p - 1 - DOWN), min(len(g), p + UP))
                 for p, s in tss)
    regions = []
    for a, b in win:
        if regions and a <= regions[-1][1]:
            regions[-1][1] = max(regions[-1][1], b)
        else:
            regions.append([a, b])
    starts = np.array([a for a, b in regions])
    sites = collections.defaultdict(list)
    for f, a in gff(os.path.join(raw, 'RISet.gff3')):
        if a.get('riType', '').startswith('tf-'):
            sites[a['Name']].append((int(f[3]) - 1, int(f[4]), f[6]))
    seqs = [g[a:b] for a, b in regions]
    os.makedirs(root, exist_ok=True)
    meta = []
    for tf in sorted(sites):
        ss = sorted(set(sites[tf]))
        L = collections.Counter(b - a for a, b, s in ss).most_common(1)[0][0]
        mod = [g[a:b] if s == '+' else g[a:b].translate(COMP)[::-1] for a, b, s in ss if b - a == L]
        mod = [x for x in mod if set(x) <= set('ACGT')]
        inreg = collections.defaultdict(list)
        for a, b, s in ss:
            i = int(np.searchsorted(starts, a, side='right')) - 1
            for j in (i, i + 1):
                if 0 <= j < len(regions):
                    ra, rb = regions[j]
                    if min(b, rb) - max(a, ra) >= (b - a) / 2:
                        inreg[j].append((max(a, ra) - ra, min(b, rb) - ra))
        if len(inreg) < MIN_REG or len(mod) < MIN_SITES:
            continue
        C = np.array([[sum(x[i] == c for x in mod) for c in 'ACGT'] for i in range(L)], float) + 0.25
        P = C / C.sum(axis=1, keepdims=True)
        pos = sorted(inreg)
        name = tf.replace('/', '_').replace(' ', '_')
        units = [(name, pos, dict(truth=[pwm_to_iupac(P)], pwm=P.round(4).tolist(), tf=tf, n_sites_modal=len(mod),
                                  site_len=L, sites=[sorted(set(inreg[j])) for j in pos]))]
        rng = random.Random(int(hashlib.sha256(tf.encode()).hexdigest()[:8], 16))
        for k in range(NULL_PER):
            units.append((f'null_{name}__r{k}', sorted(rng.sample(range(len(regions)), len(pos))),
                          dict(truth=[], null=True, kind='regulon_null', tf=tf)))
        for u, fg_idx, t in units:
            d = os.path.join(root, u)
            os.makedirs(d, exist_ok=True)
            fs = set(fg_idx)
            with open(os.path.join(d, 'fg.fa'), 'w') as fh:
                for j in fg_idx:
                    fh.write(f'>r{j}\n{seqs[j]}\n')
            with open(os.path.join(d, 'bg.fa'), 'w') as fh:
                for j in range(len(regions)):
                    if j not in fs:
                        fh.write(f'>r{j}\n{seqs[j]}\n')
            json.dump(dict(dataset=u, **t), open(os.path.join(d, 'truth.json'), 'w'))
        meta.append((tf, len(pos), len(mod), L, pwm_to_iupac(P)))
    print(f'{len(regions)} promoter regions (median {int(np.median([b - a for a, b in regions]))} bp); {len(meta)} TFs')
    for m in meta:
        print('\t'.join(map(str, m)))


if __name__ == '__main__':
    main()
