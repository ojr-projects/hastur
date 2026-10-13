# %% e2_make.py: E2 quick look, TF knockout responders in PRECISE-1K (ROADMAP_v2 E2; HASTUR first)
"""
Contrasts: every QC-passed PRECISE-1K strain 'del_<gene>' whose gene is a RegulonDB TF with a truth matrix (E1
construction, >= 8 sites at modal length), against the samples of the same project with the same base media, carbon
source and electron acceptor whose strain carries no deletion. Needs >= 2 samples per side.
Responders: genes with |mean log2 TPM difference| >= 1, the 100 largest, mapped to D1/E1 promoter regions (genes ->
RegulonDB TUs -> TSS regions). Unit: fg = responder regions, bg = every other region; >= 5 regions.
Nulls: 20 random region sets per contrast at the same size. truth.json: the TF's truth matrix and its site intervals
in the fg regions that hold a site (for site recovery), and the share of fg regions with a RegulonDB site of the TF
(the direct share).
Usage: python bench/e2_make.py <imodulon_raw_dir> <regulondb_raw_dir> <e1_units_root> <root>
"""
import collections
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

MIN_REG, NULL_PER, TOP, LFC = 5, 20, 100, 1.0


def main():
    raw_im, raw_rdb, e1, root = sys.argv[1:5]
    g = ''.join(l.strip() for l in open(GENOME) if not l.startswith('>')).upper()
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

    def region_of(p0):
        i = int(np.searchsorted(starts, p0, side='right')) - 1
        return i if i >= 0 and regions[i][0] <= p0 < regions[i][1] else None
    tus = []
    for f, a in gff(os.path.join(raw_rdb, 'TUSet.gff3')):
        if f[2] == 'transcript' and a.get('TSS', '').strip().isdigit() and f[3].isdigit() and f[4].isdigit():
            r = region_of(int(a['TSS']) - 1)
            if r is not None:
                tus.append((f[6], int(f[3]) - 1, int(f[4]), r))
    gi = pd.read_csv(os.path.join(raw_im, 'gene_info.csv'), index_col=0)

    def gene_regions(b):
        if b not in gi.index:
            return set()
        r = gi.loc[b]
        a0, e, s = int(r.start) - 1, int(r.end), r.strand
        return {reg for st, ta, tb, reg in tus if st == s and ta <= a0 and e <= tb}
    sites = collections.defaultdict(list)
    for f, a in gff(os.path.join(raw_rdb, 'RISet.gff3')):
        if a.get('riType', '').startswith('tf-'):
            sites[a['Name']].append((int(f[3]) - 1, int(f[4]), f[6]))
    truths = {}
    for tf, ss in sites.items():
        ss = sorted(set(ss))
        L = collections.Counter(b - a for a, b, s in ss).most_common(1)[0][0]
        mod = [g[a:b] if s == '+' else g[a:b].translate(COMP)[::-1] for a, b, s in ss if b - a == L]
        mod = [x for x in mod if set(x) <= set('ACGT')]
        if len(mod) >= 8:
            C = np.array([[sum(x[i] == c for x in mod) for c in 'ACGT'] for i in range(L)], float) + 0.25
            truths[tf] = (C / C.sum(axis=1, keepdims=True), L, ss)
    tf_by = {t.lower(): t for t in truths}
    # gene name -> TF name via gene_info gene_name
    md = pd.read_csv(os.path.join(raw_im, 'metadata_qc.csv'), index_col=0)
    X = pd.read_csv(os.path.join(raw_im, 'log_tpm_qc.csv'), index_col=0)
    md = md[md.index.isin(X.columns)]
    key = ['project', 'Base Media', 'Carbon Source (g/L)', 'Electron Acceptor']
    md['_key'] = md[key].apply(lambda r: '|'.join(map(str, r)), axis=1)
    md['_del'] = md['Strain Description'].astype(str).str.findall(r'del_([A-Za-z0-9]+)')
    seqs = [g[a:b] for a, b in regions]
    os.makedirs(root, exist_ok=True)
    meta = []
    for (k, dels), grp in md.groupby(['_key', md['_del'].apply(tuple)]):
        if len(dels) != 1 or len(grp) < 2:
            continue
        gene = dels[0]
        tf = tf_by.get(gene.lower())
        if tf is None:
            continue
        ctrl = md[(md._key == k) & (md._del.apply(len) == 0)]
        if len(ctrl) < 2:
            continue
        d = X[grp.index].mean(axis=1) - X[ctrl.index].mean(axis=1)
        resp = d[d.abs() >= LFC].abs().sort_values(ascending=False).index[:TOP]
        pos = sorted(set().union(*[gene_regions(b) for b in resp])) if len(resp) else []
        name = re.sub(r'[^A-Za-z0-9_.-]+', '_', f'{tf}__{grp.project.iloc[0]}__{len(meta)}')
        P, L, ss = truths[tf]
        inreg = collections.defaultdict(list)
        for a, b, s in ss:
            i = int(np.searchsorted(starts, a, side='right')) - 1
            for j in (i, i + 1):
                if 0 <= j < len(regions):
                    ra, rb = regions[j]
                    if min(b, rb) - max(a, ra) >= (b - a) / 2:
                        inreg[j].append((max(a, ra) - ra, min(b, rb) - ra))
        idx = [n for n, j in enumerate(pos) if j in inreg]
        meta.append(dict(unit=name, tf=tf, project=grp.project.iloc[0], n_ko=len(grp), n_ctrl=len(ctrl),
                         n_resp=len(resp), n_regions=len(pos), direct_share=round(len(idx) / max(1, len(pos)), 3),
                         included=len(pos) >= MIN_REG))
        if len(pos) < MIN_REG:
            continue
        t = dict(tf=tf, truths=[dict(tf=tf, pwm=P.round(4).tolist(), iupac=pwm_to_iupac(P), site_len=L, fg_idx=idx,
                                     sites=[sorted(set(inreg[pos[n]])) for n in idx])], truth=[pwm_to_iupac(P)],
                 direct_share=meta[-1]['direct_share'], module=name)
        units = [(name, pos, t)]
        rng = random.Random(int(hashlib.sha256(f'E2|{name}'.encode()).hexdigest()[:16], 16))
        for kk in range(NULL_PER):
            units.append((f'null_{name}__r{kk}', sorted(rng.sample(range(len(regions)), len(pos))),
                          dict(truth=[], null=True, kind='e2_null', module=name)))
        for u, fg_idx, tt in units:
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
            json.dump(dict(dataset=u, **tt), open(os.path.join(dd, 'truth.json'), 'w'))
    df = pd.DataFrame(meta)
    df.to_csv(os.path.join(root, 'contrasts.tsv'), sep='\t', index=False)
    print(df.to_string(index=False))


if __name__ == '__main__':
    main()
