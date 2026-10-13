# %% score.py: the evaluator (E0). Match rules, shuffled-truth chance baseline, PWM similarity
"""
One module, used by every phase. Calls are either IUPAC strings or PWMs (L x 4 arrays over ACGT); a truth is an IUPAC
string. Rules (THEORY.md section 3 and docs/briefs in Comsequent):

  p5        >= 80% of the truth's specified positions covered by call letters that are subsets of the truth letters,
            <= 1 conflict, any offset, either strand (Comsequent P5; high chance floor for long calls)
  contains  every specified truth position covered by a subset letter, <= 2 further specified call positions
            (Comsequent P1 lenient; chance floor ~1%)            <- the PRIMARY rule for this project
  pwm_sim   PWM similarity: mean over the truth's informative columns (IC >= 1 bit) of 1 - total-variation distance
            to the aligned call column (uniform where the call does not cover), best offset/strand; a hit is
            pwm_sim >= SIM_HIT

Every scorecard reports, per unit and tool: first rank of a hit under each rule, and the same against letter-shuffled
truths (chance). Representation of a call as a PWM: IUPAC letters -> uniform over their bases.
Self-test: python eval/score.py
"""
import math

import numpy as np

IU = {'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'R': 'AG', 'Y': 'CT', 'S': 'CG', 'W': 'AT', 'K': 'GT', 'M': 'AC',
      'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'V': 'ACG', 'N': 'ACGT'}
COMP = str.maketrans('ACGTRYSWKMBDHVN', 'TGCAYRSWMKVHDBN')
SIM_HIT = 0.85
N_SHUFFLE = 10


def clean(call):
    return str(call).upper().replace('X', 'N').strip('N')


def iupac_to_pwm(s):
    m = np.array([[1.0 if b in IU[c] else 0.0 for b in 'ACGT'] for c in s])
    return m / m.sum(axis=1, keepdims=True)


def pwm_to_iupac(m, single=0.6, double=0.8):
    """Consensus: a base with p >= single; else a two-base code with combined p >= double; else N."""
    rev = {v: k for k, v in IU.items()}
    out = []
    for col in m:
        o = np.argsort(-col)
        if col[o[0]] >= single:
            out.append('ACGT'[o[0]])
        elif col[o[0]] + col[o[1]] >= double:
            out.append(rev[''.join(sorted('ACGT'[o[0]] + 'ACGT'[o[1]], key='ACGT'.index))])
        else:
            out.append('N')
    return ''.join(out).strip('N') or 'N'


def _strands(c):
    return (c, c.translate(COMP)[::-1])


def hit_p5(call, truth, frac=0.8):
    c0 = clean(call)
    spec = [j for j, ch in enumerate(truth) if ch != 'N']
    if not c0 or not spec or any(ch not in IU for ch in c0):
        return False
    need = math.ceil(frac * len(spec))
    for c in _strands(c0):
        for off in range(-len(truth) + 1, len(c)):
            match = conflict = 0
            for j in spec:
                i = off + j
                if 0 <= i < len(c) and c[i] != 'N':
                    if set(IU[c[i]]) <= set(IU[truth[j]]):
                        match += 1
                    else:
                        conflict += 1
            if match >= need and conflict <= 1:
                return True
    return False


def hit_contains(call, truth, extra=2):
    c0 = clean(call)
    if not c0 or any(ch not in IU for ch in c0):
        return False
    for c in _strands(c0):
        for off in range(-len(truth) + 1, len(c)):
            ok, used = True, set()
            for j, t in enumerate(truth):
                if t == 'N':
                    continue
                i = off + j
                if i < 0 or i >= len(c) or c[i] == 'N' or not set(IU[c[i]]) <= set(IU[t]):
                    ok = False
                    break
                used.add(i)
            if ok and sum(1 for i, ch in enumerate(c) if ch != 'N' and i not in used) <= extra:
                return True
    return False


def pwm_sim(call_pwm, truth, min_ic=1.0):
    """Mean over informative truth columns of 1 - TV distance to the aligned call column (uniform if uncovered).
    `truth` is an IUPAC string (converted to a 0/1 matrix) or, better, the true PWM (columns with IC >= 0.5 bit count)."""
    if isinstance(truth, str):
        T = iupac_to_pwm(truth)
        ic = np.array([math.log2(4 / len(IU[c])) for c in truth])
    else:
        T = np.asarray(truth, float)
        ic = 2 + (T * np.log2(np.clip(T, 1e-12, 1))).sum(axis=1)
        min_ic = 0.5
    cols = np.flatnonzero(ic >= min_ic)
    if len(cols) == 0:
        return 0.0
    M = np.asarray(call_pwm, float)
    M = M / M.sum(axis=1, keepdims=True)
    U = np.full(4, 0.25)
    best = 0.0
    for m in (M, M[::-1, ::-1]):
        W = len(m)
        for off in range(-W + 1, len(truth)):
            s = 0.0
            for j in cols:
                i = j - off
                col = m[i] if 0 <= i < W else U
                s += 1 - 0.5 * np.abs(col - T[j]).sum()
            best = max(best, s / len(cols))
    return best


def as_pwm(call):
    return iupac_to_pwm(clean(call) or 'N') if isinstance(call, str) else np.asarray(call, float)


def as_iupac(call):
    return clean(call) if isinstance(call, str) else pwm_to_iupac(np.asarray(call, float))


def first_ranks(calls, truth, topk=10, truth_pwm=None):
    """First 1-based rank of a hit under each rule (inf if none in the top k). PWM similarity uses the true matrix
    when given (a correct soft PWM should not be penalised for differing from a 0/1 consensus)."""
    ref = truth if truth_pwm is None else truth_pwm
    out = {}
    sims = [pwm_sim(as_pwm(c), ref) for c in calls[:topk]]
    for name, f in (('p5', lambda c: hit_p5(as_iupac(c), truth)), ('contains', lambda c: hit_contains(as_iupac(c), truth))):
        out[name] = next((i + 1 for i, c in enumerate(calls[:topk]) if f(c)), math.inf)
    out['pwm_sim'] = next((i + 1 for i, v in enumerate(sims) if v >= SIM_HIT), math.inf)
    out['best_sim'] = max(sims, default=0.0)
    return out


def shuffled(truth, rng):
    letters = [c for c in truth if c != 'N']
    rng.shuffle(letters)
    it = iter(letters)
    return ''.join(c if c == 'N' else next(it) for c in truth)


def shuffled_pwm(m, rng):
    return np.asarray(m)[rng.permutation(len(m))][:, rng.permutation(4)]


def score_unit(calls, truths, topk=10, n_shuffle=N_SHUFFLE, seed=0, truth_pwms=None):
    """Real and chance (shuffled-truth) readouts for one unit; truths: list of acceptable truth strings; truth_pwms:
    optional matching list of true matrices (chance then shuffles their columns and letters)."""
    rng = np.random.default_rng(seed)
    tp = truth_pwms or [None] * len(truths)
    real = {}
    for t, m in zip(truths, tp):
        for k, v in first_ranks(calls, t, topk, m).items():
            real[k] = max(real.get(k, -math.inf), v) if k == 'best_sim' else min(real.get(k, math.inf), v)
    ch = {k: [] for k in real}
    for _ in range(n_shuffle):
        for t, m in zip(truths, tp):
            s = shuffled(t, rng)
            if s == t:
                continue
            r = first_ranks(calls, s, topk, None if m is None else shuffled_pwm(m, rng))
            for k in ch:
                ch[k].append(r[k] <= topk if k != 'best_sim' else r[k])
    out = {}
    for k, v in real.items():
        if k == 'best_sim':
            out['best_sim'] = v
            out['chance_best_sim'] = float(np.mean(ch[k])) if ch[k] else np.nan
        else:
            out[f'rank_{k}'] = v
            out[f'hit10_{k}'] = v <= topk
            out[f'hit1_{k}'] = v == 1
            out[f'chance10_{k}'] = float(np.mean(ch[k])) if ch[k] else np.nan
    return out


# %% Self-test (oracle test, DEV_INSIGHTS 3.5)
if __name__ == '__main__':
    truth = 'GGAGCTGTCCATGGTGCTGA'
    member, near = 'TCAGCACCATGGACAGCTCC', 'GGAGCTGTCCATGGTGCTGA'[:10]
    pwm = iupac_to_pwm(truth) * 0.85 + 0.0375
    checks = {
        'truth p5/contains/sim': (hit_p5(truth, truth), hit_contains(truth, truth), pwm_sim(iupac_to_pwm(truth), truth)),
        'reverse complement': (hit_p5(member, truth), hit_contains(member, truth), pwm_sim(iupac_to_pwm(member), truth)),
        'soft PWM form': (hit_p5(pwm_to_iupac(pwm), truth), hit_contains(pwm_to_iupac(pwm), truth), pwm_sim(pwm, truth)),
        'half-length fragment': (hit_p5(near, truth), hit_contains(near, truth), pwm_sim(iupac_to_pwm(near), truth)),
        'IUPAC truth GRCGYC vs member GACGTC': (hit_p5('GACGTC', 'GRCGYC'), hit_contains('GACGTC', 'GRCGYC'),
                                                 pwm_sim(iupac_to_pwm('GACGTC'), 'GRCGYC')),
    }
    for k, v in checks.items():
        print(f'{k:40s} p5 {v[0]!s:5s} contains {v[1]!s:5s} pwm_sim {v[2]:.2f}')
    assert all(checks['truth p5/contains/sim'][:2]) and checks['truth p5/contains/sim'][2] == 1
    assert all(checks['reverse complement'][:2]) and checks['soft PWM form'][1] and checks['soft PWM form'][2] > 0.85
    assert not checks['half-length fragment'][1] and checks['half-length fragment'][2] < SIM_HIT
    assert checks['IUPAC truth GRCGYC vs member GACGTC'][1]
    calls = ['CG', 'TCAGSACCWYGGACAGNKCC', 'GNNCNNNGNNCNG']
    print('REST all-in-one top 3:', score_unit(calls, [truth]))
    print('self-test ok')
