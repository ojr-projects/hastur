# %% frontier.py: FR1 sensitivity vs false calls per run, thresholds swept (splits/protocol_FR1.md)
"""
Every reported call gets one significance score s on a common scale (-log10): HASTUR: one-sided normal p of its
detection z; STREME: E-value; HOMER: reported p-value. A positive unit is a hit at threshold t if any call with
s >= t matches the truth (contains, P5, or PWM similarity >= 0.85 against the true matrix); a null run's false calls at
t are its calls with s >= t. Output: curves, and sensitivity at matched false-call rates with bootstrap CIs.
Usage: python bench/frontier.py --pos <root>[:truth_map.json] ... --null <root> ... --tools hlat_a3,streme_all,homer
"""
import argparse
import json
import math
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy.stats import norm

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'eval'))
from score import as_iupac, as_pwm, hit_contains, hit_p5, pwm_sim  # noqa: E402

RATES = (0.05, 0.1, 0.5)


def scored_calls(d, tool):
    """List of (s, call) for every reported call, s = -log10 significance."""
    o = os.path.join(d, tool)
    if not os.path.exists(os.path.join(o, 'runtime_s')):
        return None
    if tool.startswith('hlat'):
        t = os.path.join(o, 'hlat.tsv')
        m = os.path.join(o, 'hlat.meme')
        if not os.path.exists(t):
            return []
        tab = pd.read_csv(t, sep='\t')
        if tab.empty:
            return []
        txt = open(m).read()
        mats = re.findall(r'letter-probability matrix:[^\n]*\n((?:\s*[0-9.eE+-]+\s+[0-9.eE+-]+\s+[0-9.eE+-]+\s+[0-9.eE+-]+\s*\n)+)', txt)
        pw = [np.array([[float(v) for v in r.split()] for r in x.strip().split('\n')]) for x in mats]
        if 'p_adj' in tab.columns:      # family-wise adjusted p (LN7 permutation calibration): one calibrated scale
            return [(-math.log10(max(float(q), 1e-300)), p) for q, p in zip(tab.p_adj, pw)]
        return [(-norm.logsf(z) / math.log(10), p) for z, p in zip(tab.z, pw)]
    if tool.startswith('streme') or tool in ('meme', 'xstreme'):
        p = os.path.join(o, {'meme': 'meme.txt', 'xstreme': 'combined.meme'}.get(tool, 'streme.txt'))
        if not os.path.exists(p):
            return []
        txt = open(p).read()
        # STREME reports E= (held-out E-value) or, when too few sequences for a hold-out set, S= (training p-value);
        # --thresh applies to whichever it reports, so both are read (protocol_D1_amend1.md)
        blocks = re.findall(r'letter-probability matrix:[^\n]*(?:E|S)= ([0-9.eE+-]+)[^\n]*\n((?:\s*[0-9.eE+-]+\s+[0-9.eE+-]+\s+[0-9.eE+-]+\s+[0-9.eE+-]+\s*\n)+)', txt)
        out = []
        for e, x in blocks:
            mat = np.array([[float(v) for v in r.split()] for r in x.strip().split('\n')])
            out.append((-math.log10(max(float(e), 1e-300)), mat))
        return out
    if tool == 'homer':
        p = os.path.join(o, 'homerMotifs.all.motifs')
        if not os.path.exists(p):
            return []
        out, rows, s = [], None, None
        for line in open(p):
            if line.startswith('>'):
                if rows:
                    out.append((s, np.array(rows)))
                f = line.rstrip('\n').split('\t')
                s, rows = -float(f[3]) / math.log(10), []
            elif line.strip():
                rows.append([float(v) for v in line.split()])
        if rows:
            out.append((s, np.array(rows)))
        return out
    raise ValueError(tool)


def call_hits(calls, truth, tpwm):
    out = []
    for s, c in calls:
        iu = as_iupac(c)
        out.append((s, hit_contains(iu, truth), hit_p5(iu, truth),
                    (pwm_sim(as_pwm(c), tpwm) if tpwm is not None else pwm_sim(as_pwm(c), truth)) >= 0.85))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--pos', nargs='+', required=True)
    ap.add_argument('--null', nargs='+', required=True)
    ap.add_argument('--tools', default='hlat_a3,streme_all,homer')
    ap.add_argument('--out', default=os.path.expanduser('~/runs/frontier'))
    ap.add_argument('--exclude-kinds', default='', help='comma list of substrings: null kinds left out (e.g. intergenic)')
    a = ap.parse_args()
    tools = a.tools.split(',')
    pos = []
    for spec in a.pos:
        root, tmap = (spec.split(':') + [None])[:2]
        tm = json.load(open(tmap)) if tmap else None
        for ds in sorted(os.listdir(root)):
            d = os.path.join(root, ds)
            if not os.path.exists(os.path.join(d, 'fg.fa')):
                continue
            if tm is not None:
                f = ds.split('__')[0]
                if f not in tm:
                    continue
                truths = [(tm[f]['truth'], np.array(tm[f]['pwm']))]
            else:
                if not os.path.exists(os.path.join(d, 'truth.json')):
                    continue
                j = json.load(open(os.path.join(d, 'truth.json')))
                pw = j.get('pwms') or ([j['pwm']] if j.get('pwm') else [None] * len(j['truth']))
                truths = [(t, None if p is None else np.array(p)) for t, p in zip(j['truth'], pw)]
            for tool in tools:
                calls = scored_calls(d, tool)
                if calls is None:
                    continue
                hits = [call_hits(calls, t, p) for t, p in truths]
                merged = [(c[0], any(h[i][1] for h in hits), any(h[i][2] for h in hits), any(h[i][3] for h in hits))
                          for i, c in enumerate(hits[0])] if hits else []
                tw = None
                if truths[0][1] is not None:
                    T = truths[0][1]
                    tw = int(((2 + (T * np.log2(np.clip(T, 1e-12, 1))).sum(axis=1)) >= 0.5).sum())
                pos.append(dict(unit=os.path.basename(root) + '/' + ds, tool=tool, calls=merged, tw=tw))
    nul = []
    for root in a.null:
        for ds in sorted(os.listdir(root)):
            d = os.path.join(root, ds)
            if not os.path.exists(os.path.join(d, 'truth.json')):
                continue
            kind = json.load(open(os.path.join(d, 'truth.json'))).get('kind', '?')
            if any(x and x in kind for x in a.exclude_kinds.split(',')):
                continue
            for tool in tools:
                calls = scored_calls(d, tool)
                if calls is not None:
                    nul.append(dict(run=ds, kind=kind, tool=tool, scores=[c[0] for c in calls]))
    os.makedirs(a.out, exist_ok=True)
    for p in pos:
        p['set'] = 'selex' if 'selex' in p['unit'].split('/')[0] else 'chip'
    for n in nul:
        n['family'] = 'human' if n['kind'].startswith('human') else ('sim' if n['kind'].split('_')[0] in ('iid', 'markov', 'cpgshift') else 'bacterial')
    json.dump(dict(pos=pos, nul=nul), open(os.path.join(a.out, 'calls.json'), 'w'), default=float)
    rng = np.random.default_rng(0)
    print(f'positives: {len({p["unit"] for p in pos})} units ({ {s: len({p["unit"] for p in pos if p["set"] == s}) for s in ("chip", "selex")} }); '
          f'nulls: {len({n["run"] for n in nul})} runs ({ {f: len({n["run"] for n in nul if n["family"] == f}) for f in ("sim", "human", "bacterial")} })')

    def threshold(N, r):
        """Smallest t whose mean false calls per null run <= r."""
        sc = [np.sort(n['scores']) for n in N]
        allv = np.concatenate([x for x in sc if len(x)] + [[]])
        grid = np.unique(np.concatenate([np.nextafter(allv, np.inf), [-np.inf, np.inf]]))   # just above each null score
        fc = np.array([np.mean([len(x) - np.searchsorted(x, t, 'left') for x in sc]) for t in grid])
        return grid[np.flatnonzero(fc <= r)[0]]

    def hits(P, t, k=1):
        """unit -> hit (any call with s >= t matching under rule k: 1 contains, 2 P5, 3 sim)."""
        return {p['unit']: int(any(c[0] >= t and c[k] for c in p['calls'])) for p in P}

    def boot(x):
        x = np.asarray(x, float)
        m = rng.choice(x, size=(10000, len(x))).mean(axis=1)
        return np.quantile(m, 0.025), np.quantile(m, 0.975)

    fams = {'real': ('human', 'bacterial'), 'human': ('human',), 'bacterial': ('bacterial',), 'sim': ('sim',),
            'all': ('sim', 'human', 'bacterial')}
    rows = []
    for fam, members in fams.items():
        for pset in ('chip', 'selex', 'all'):
            print(f'\n== nulls: {fam} | positives: {pset}')
            H = {}
            for tool in tools:
                P = [p for p in pos if p['tool'] == tool and (pset == 'all' or p['set'] == pset)]
                N = [n for n in nul if n['tool'] == tool and n['family'] in members]
                if not P or not N:
                    continue
                line = [f'{tool:12s}']
                for r in RATES:
                    t = threshold(N, r)
                    h = hits(P, t, 2)
                    H[(tool, r)] = h
                    lo, hi = boot(list(h.values()))
                    hs = hits(P, t, 3)
                    rows.append(dict(nulls=fam, pos=pset, tool=tool, rate=r, t=t, sens_p5=np.mean(list(h.values())), lo=lo, hi=hi,
                                     sens_sim=np.mean(list(hs.values())), sens_contains=np.mean(list(hits(P, t, 1).values()))))
                    line.append(f'@{r}: P5 {np.mean(list(h.values())):.3f} [{lo:.3f},{hi:.3f}] sim {np.mean(list(hs.values())):.3f} (t={t:.2f})')
                print(' | '.join(line))
            hl = [t for t in tools if t.startswith('hlat')]
            st = [t for t in tools if t.startswith('streme')]
            for h_ in hl:
                for s_ in st:
                    for r in RATES:
                        if (h_, r) in H and (s_, r) in H:
                            u = sorted(set(H[(h_, r)]) & set(H[(s_, r)]))
                            x = [H[(h_, r)][k] - H[(s_, r)][k] for k in u]
                            lo, hi = boot(x)
                            tag = ''
                            if fam == 'real' and pset == 'all' and r == 0.1:
                                tag = ('  <- BETTER FRONTIER' if np.mean(x) >= 0.05 and lo > 0 else
                                       '  <- OPERATING POINT ONLY' if abs(np.mean(x)) <= 0.05 else '  <- report as is')
                            print(f'  {h_} - {s_} @ {r}: {np.mean(x):+.3f} [{lo:+.3f},{hi:+.3f}] (n={len(u)}){tag}')
    pd.DataFrame(rows).to_csv(os.path.join(a.out, 'frontier.tsv'), sep='\t', index=False)
    print('\n== STREME calibration: share of its null-run calls with E (or S) <= 0.05, and such calls per run')
    for tool in [t for t in tools if t.startswith('streme')]:
        for fam in ('sim', 'human', 'bacterial'):
            N = [n for n in nul if n['tool'] == tool and n['family'] == fam]
            if N:
                allc = sum(len(n['scores']) for n in N)
                sig = sum(sum(s >= -math.log10(0.05) for s in n['scores']) for n in N)
                print(f'  {tool} {fam}: {sig}/{allc} calls ({sig / max(allc, 1):.3f}); {sig / len(N):.2f} per run')
    print('\nreported calls per null run (no threshold), by kind:')
    nk = pd.DataFrame([dict(tool=n['tool'], kind=n['kind'], n=len(n['scores'])) for n in nul])
    print(nk.pivot_table(index='kind', columns='tool', values='n', aggfunc='mean').round(2).to_string())
    print('\n== sensitivity by truth length (real nulls, @0.1, P5), informative truth columns <= 10 vs > 10')
    for tool in tools:
        N = [n for n in nul if n['tool'] == tool and n['family'] in ('human', 'bacterial')]
        P = [p for p in pos if p['tool'] == tool]
        if not N or not P:
            continue
        t = threshold(N, 0.1)
        h = hits(P, t, 2)
        for lab, sel in (('short', lambda w: w <= 10), ('long', lambda w: w > 10)):
            u = [p['unit'] for p in P if p.get('tw') is not None and sel(p['tw'])]
            if u:
                print(f'  {tool:12s} {lab}: {np.mean([h[k] for k in u]):.3f} (n={len(u)})')

if __name__ == '__main__':
    main()
