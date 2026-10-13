# %% imodulon_make.py: E1 E. coli iModulon screen units (ROADMAP_v2 E1; splits/protocol_E1.md)
"""
Promoter regions: the D1 construction (bench/regulon_make.py): [TSS - 250, TSS + 50] around every RegulonDB TSS,
overlapping windows merged (1,942 regions). Forward-strand sequence; every tool searches both strands.
Module -> regions: an iModulon's member genes (|M| > its PRECISE-1K threshold) -> every RegulonDB transcription unit
(TUSet) that contains the gene on its strand -> the region holding that TU's TSS. Genes in no TU with a mapped TSS
are dropped (coverage reported).
Units:
  <root>/<module>/             fg = the module's regions, bg = every other region
  <root>/null_<module>__r<k>/  fg = random regions, same count; bg = every other region (k < NULL_PER)
truth.json of a module: the module's annotated regulator tokens (iModulon `regulator`, split on + and /). For each
token that is a RegulonDB TF with >= MIN_SITES sites at its modal length: the truth matrix from all of its sites
(modal length, strand-oriented, +0.25 pseudocount, as D1) and the site intervals in fg coordinates, for the fg regions
that hold at least half of one of its sites.
Null seeds: random.Random of sha256('E1|<module>') (disjoint from D1's sha256('<TF>') seeds).
Included modules: >= MIN_REG fg regions.
Usage: python bench/imodulon_make.py <imodulon_raw_dir> <regulondb_raw_dir> <root>
"""
import collections
import gzip
import hashlib
import json
import os
import random
import re
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'eval'))
sys.path.insert(0, HERE)
from regulon_make import COMP, DOWN, GENOME, UP, gff  # noqa: E402
from score import pwm_to_iupac  # noqa: E402

MIN_REG, MIN_SITES, NULL_PER = 5, 8, 20


def safe(name):
    return re.sub(r'[^A-Za-z0-9_.-]+', '_', name).strip('_')


def main():
    raw_im, raw_rdb, root = sys.argv[1], sys.argv[2], sys.argv[3]
    g = ''.join(l.strip() for l in open(GENOME) if not l.startswith('>')).upper()
    # promoter regions (identical to D1)
    tss = [(int(f[3]), f[6]) for f, a in gff(os.path.join(raw_rdb, 'PromoterSet.gff3')) if f[2] == 'transcription_start_site']
    win = sorted((max(0, p - 1 - UP), min(len(g), p + DOWN)) if s == '+' else (max(0, p - 1 - DOWN), min(len(g), p + UP))
                 for p, s in tss)
    regions = []
    for a, b in win:
        if regions and a <= regions[-1][1]:
            regions[-1][1] = max(regions[-1][1], b)
        else:
            regions.append([a, b])
    starts = np.array([a for a, b in regions])

    def region_of(pos0):
        i = int(np.searchsorted(starts, pos0, side='right')) - 1
        return i if i >= 0 and regions[i][0] <= pos0 < regions[i][1] else None

    # transcription units: (strand, start0, end, region of TSS)
    tus = []
    for f, a in gff(os.path.join(raw_rdb, 'TUSet.gff3')):
        if f[2] != 'transcript' or not a.get('TSS', '').strip().isdigit() or not (f[3].isdigit() and f[4].isdigit()):
            continue
        r = region_of(int(a['TSS']) - 1)
        if r is not None:
            tus.append((f[6], int(f[3]) - 1, int(f[4]), r))
    # genes and iModulon membership
    gi = pd.read_csv(os.path.join(raw_im, 'gene_info.csv'), index_col=0)
    d = json.load(gzip.open(os.path.join(raw_im, 'precise1k.json.gz')))
    thr = d['thresholds']
    M = pd.read_csv(os.path.join(raw_im, 'M.csv'), index_col=0)
    imt = pd.read_csv(os.path.join(raw_im, 'imodulon_table.csv'), index_col=0)
    # M.csv names one column differently from the table and thresholds (Superoxide vs SoxS); the orders agree, so
    # columns are aligned by position, and any other disagreement stops the build
    bad = [(a, b) for a, b in zip(M.columns, imt.index) if a != b]
    if len(M.columns) != len(imt) or bad != [('Superoxide', 'SoxS')]:
        raise ValueError(f'M.csv and imodulon_table.csv disagree: {bad}')
    M.columns = list(imt.index)

    def gene_regions(b):
        if b not in gi.index:
            return set()
        r = gi.loc[b]
        a0, e, s = int(r.start) - 1, int(r.end), r.strand
        return {reg for st, ta, tb, reg in tus if st == s and ta <= a0 and e <= tb}
    # RegulonDB TF sites (as D1)
    sites = collections.defaultdict(list)
    for f, a in gff(os.path.join(raw_rdb, 'RISet.gff3')):
        if a.get('riType', '').startswith('tf-'):
            sites[a['Name']].append((int(f[3]) - 1, int(f[4]), f[6]))
    tf_by_lower = {t.lower(): t for t in sites}
    truths = {}
    for tf, ss in sites.items():
        ss = sorted(set(ss))
        L = collections.Counter(b - a for a, b, s in ss).most_common(1)[0][0]
        mod = [g[a:b] if s == '+' else g[a:b].translate(COMP)[::-1] for a, b, s in ss if b - a == L]
        mod = [x for x in mod if set(x) <= set('ACGT')]
        if len(mod) < MIN_SITES:
            continue
        C = np.array([[sum(x[i] == c for x in mod) for c in 'ACGT'] for i in range(L)], float) + 0.25
        truths[tf] = (C / C.sum(axis=1, keepdims=True), L, len(mod), ss)

    seqs = [g[a:b] for a, b in regions]
    os.makedirs(root, exist_ok=True)
    meta = []
    for mname in imt.index:
        if mname not in M.columns:
            raise KeyError(mname)
        genes = M.index[M[mname].abs() > thr[mname]]
        gr = {b: gene_regions(b) for b in genes}
        pos = sorted(set().union(*gr.values())) if gr else []
        cov = sum(1 for v in gr.values() if v) / max(1, len(gr))
        reg = imt.loc[mname, 'regulator']
        toks = [t for t in re.split(r'[+/]', reg)] if isinstance(reg, str) else []
        rdb = [tf_by_lower[t.lower()] for t in toks if t.lower() in tf_by_lower and tf_by_lower[t.lower()] in truths]
        meta.append(dict(module=mname, safe=safe(mname), category=imt.loc[mname, 'enrichment_category'],
                         regulator=reg if isinstance(reg, str) else '', n_genes=len(genes), gene_cov=round(cov, 3),
                         n_regions=len(pos), truth_tfs=','.join(rdb), included=len(pos) >= MIN_REG))
        if len(pos) < MIN_REG:
            continue
        T = []
        for tf in rdb:
            P, L, nmod, ss = truths[tf]
            inreg = collections.defaultdict(list)
            for a, b, s in ss:
                i = int(np.searchsorted(starts, a, side='right')) - 1
                for j in (i, i + 1):
                    if 0 <= j < len(regions):
                        ra, rb = regions[j]
                        if min(b, rb) - max(a, ra) >= (b - a) / 2:
                            inreg[j].append((max(a, ra) - ra, min(b, rb) - ra))
            idx = [k for k, j in enumerate(pos) if j in inreg]
            T.append(dict(tf=tf, pwm=P.round(4).tolist(), iupac=pwm_to_iupac(P), site_len=L, n_sites_modal=nmod,
                          fg_idx=idx, sites=[sorted(set(inreg[pos[k]])) for k in idx]))
        name = safe(mname)
        units = [(name, pos, dict(module=mname, category=meta[-1]['category'], regulator=meta[-1]['regulator'],
                                  truths=T, truth=[t['iupac'] for t in T]))]
        rng = random.Random(int(hashlib.sha256(f'E1|{mname}'.encode()).hexdigest()[:16], 16))
        for k in range(NULL_PER):
            units.append((f'null_{name}__r{k}', sorted(rng.sample(range(len(regions)), len(pos))),
                          dict(truth=[], null=True, kind='imodulon_null', module=mname)))
        for u, fg_idx, t in units:
            dd = os.path.join(root, u)
            os.makedirs(dd, exist_ok=True)
            fs = set(fg_idx)
            with open(os.path.join(dd, 'fg.fa'), 'w') as fh:
                for j in fg_idx:
                    fh.write(f'>r{j}\n{seqs[j]}\n')
            with open(os.path.join(dd, 'bg.fa'), 'w') as fh:
                for j in range(len(regions)):
                    if j not in fs:
                        fh.write(f'>r{j}\n{seqs[j]}\n')
            json.dump(dict(dataset=u, **t), open(os.path.join(dd, 'truth.json'), 'w'))
    df = pd.DataFrame(meta)
    df.to_csv(os.path.join(root, 'modules.tsv'), sep='\t', index=False)
    print(f'{len(regions)} promoter regions; {len(tus)} TUs with a mapped TSS; {len(imt)} iModulons, '
          f'{int(df.included.sum())} included (>= {MIN_REG} regions); '
          f'{int((df.included & (df.truth_tfs != "")).sum())} included with a RegulonDB truth')
    print(f'median gene coverage {df.gene_cov.median():.2f}; median regions {df[df.included].n_regions.median():.0f}')


if __name__ == '__main__':
    main()
