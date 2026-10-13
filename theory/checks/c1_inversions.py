# %% c1_inversions.py: the three inversions of the same lattice counts (THEORY.md section 2)
"""
Claims checked (aligned windows, product background):
  T2  H_S(c_S) = pi * prod_{j in S} (p_j(c_j)/q(c_j) - 1)  for a planted product (PWM) motif, at every order,
      while I (Moebius of log f) and theta (log-linear cell-table interaction) have no product form.
  T2b H_S = 0 whenever S contains a position where the site distribution equals the background.
  T3  Under the null, the estimated H coordinates of nested patterns are uncorrelated (Hoeffding orthogonality),
      whereas raw counts, I and theta of nested patterns are strongly correlated.
  T4  The order-2 H table over (position, letter) pairs is rank one for a PWM motif and rank > 1 for an
      epistatic pair (b, sigma(b)).
Run: ~/micromamba/envs/nanomotif/bin/python theory/checks/c1_inversions.py
"""
import itertools
import os

import numpy as np
import pandas as pd

OUT = os.path.join(os.path.dirname(__file__), '..', 'results')
rng = np.random.default_rng(1)
W = 8                       # window width; motif positions 0..5, positions 6, 7 irrelevant
Q = np.array([0.3, 0.2, 0.2, 0.3])   # product background, AT-rich so q is not uniform


def column(pc, dom):
    c = np.full(4, (1 - pc) / 3)
    c[dom] = pc
    return c


def make_motif(kind):
    cons = np.array([2, 0, 1, 3, 2, 0, 0, 0])            # consensus letters (positions 6, 7: background)
    cols = [Q.copy() for _ in range(W)]
    if kind == 'crisp':
        for j in range(6):
            cols[j] = np.eye(4)[cons[j]]
    elif kind == 'pwm':                                  # mixed-strength columns, one two-letter column
        for j, pc in zip(range(6), [0.95, 0.86, 0.72, 0.60, 0.86, 0.72]):
            cols[j] = column(pc, cons[j])
        cols[3] = np.array([0.45, 0.05, 0.45, 0.05])     # R-like column; consensus letter A
        cons[3] = 0
    return cons, np.array(cols)


def sample(n, cols):
    u = rng.random((n, W, 1))
    cdf = np.cumsum(cols, axis=1)[None, :, :]
    return (u > cdf).sum(axis=2).astype(np.uint8)


def sample_mix(n, pi, cols):
    x = sample(n, np.array([Q] * W))
    k = rng.random(n) < pi
    x[k] = sample(int(k.sum()), cols)
    return x


def subsets(pos, kmax):
    for k in range(1, kmax + 1):
        yield from itertools.combinations(pos, k)


def stats_for(x, cons, S):
    """I, theta and H of pattern (S, cons[S]) from the window sample x."""
    n = len(x)
    match = (x[:, list(S)] == cons[list(S)])            # n x |S|
    f = {}
    for k in range(len(S) + 1):
        for T in itertools.combinations(range(len(S)), k):
            f[T] = match[:, list(T)].all(axis=1).sum() if T else n
    I = sum((-1) ** (len(S) - len(T)) * np.log(max(f[T], 0.5)) for T in f)
    qT = lambda T: np.prod([Q[cons[S[i]]] for i in T]) if T else 1.0
    H = sum((-1) ** (len(S) - len(T)) * f[T] / (n * qT(T)) for T in f)
    key = (match * (1 << np.arange(len(S)))).sum(axis=1)
    cells = np.bincount(key, minlength=1 << len(S)) + 0.5
    pop = np.array([bin(c).count('1') for c in range(1 << len(S))])
    theta = float(((-1.0) ** (len(S) - pop) * np.log(cells)).sum())
    return I, theta, H, f[tuple(range(len(S)))]


# %% T2 / T2b: planted motifs, N = 400k aligned windows, pi = 0.05
rows = []
PI = 0.05
for kind in ['crisp', 'pwm']:
    cons, cols = make_motif(kind)
    x = sample_mix(400_000, PI, cols)
    delta = np.array([cols[j][cons[j]] / Q[cons[j]] - 1 for j in range(W)])
    for S in subsets(range(W), 5):
        I, th, H, cnt = stats_for(x, cons, S)
        pred = PI * np.prod(delta[list(S)])
        rows.append(dict(kind=kind, S=''.join(map(str, S)), order=len(S),
                         irrelevant=any(j >= 6 for j in S), I=I, theta=th, H=H, H_pred=pred, count=cnt))
t2 = pd.DataFrame(rows)
t2.to_csv(os.path.join(OUT, 'c1_inversions.tsv'), sep='\t', index=False)

for kind in ['crisp', 'pwm']:
    d = t2[(t2.kind == kind) & ~t2.irrelevant]
    rel = (d.H - d.H_pred).abs() / d.H_pred.abs().clip(lower=1e-9)
    print(f'[T2 {kind}] relevant patterns: {len(d)}; median |H - pred|/|pred| = {rel.median():.3f}; '
          f'90th pct {rel.quantile(.9):.3f}')
    print(d.groupby('order')[['I', 'theta', 'H', 'H_pred']].mean().round(3).to_string())
    di = t2[(t2.kind == kind) & t2.irrelevant]
    print(f'[T2b {kind}] patterns with an irrelevant position: mean |H| = {di.H.abs().mean():.4f} vs '
          f'relevant mean |H| = {d.H.abs().mean():.3f}; mean |theta| {di.theta.abs().mean():.3f} vs {d.theta.abs().mean():.3f}')

# Product test: H_S / prod_{j in S} H_j / pi^(|S|-1) should be 1 if H is rank one; same ratio for theta and I.
for kind in ['crisp', 'pwm']:
    d = t2[(t2.kind == kind) & ~t2.irrelevant].set_index('S')
    out = []
    for S, r in d.iterrows():
        if len(S) < 2:
            continue
        for col in ['H', 'theta', 'I']:
            singles = np.prod([d.loc[c, col] for c in S])
            out.append((col, len(S), r[col] / (singles / PI ** (len(S) - 1)) if col == 'H' else r[col] / singles))
    o = pd.DataFrame(out, columns=['stat', 'order', 'ratio'])
    print(f'[T2 product test, {kind}] median ratio (1 = product form):')
    print(o.groupby(['stat', 'order']).ratio.median().unstack().round(3).to_string())

# %% T3: null orthogonality across replicate samples (bg only), nested patterns
cons, _ = make_motif('pwm')
chains = [(0,), (0, 1), (0, 1, 2), (0, 1, 2, 4)]
reps = []
for r in range(300):
    x = sample(20_000, np.array([Q] * W))
    reps.append({S: stats_for(x, cons, S) for S in chains})
names = ['I', 'theta', 'H', 'count']
print('[T3] null correlation between nested patterns across 300 replicates (N = 20k windows each):')
res = []
for a, b in itertools.combinations(chains, 2):
    for i, nm in enumerate(names):
        va = np.array([rp[a][i] for rp in reps])
        vb = np.array([rp[b][i] for rp in reps])
        res.append(dict(pair=f"{''.join(map(str, a))}-{''.join(map(str, b))}", stat=nm, corr=np.corrcoef(va, vb)[0, 1]))
t3 = pd.DataFrame(res).pivot(index='pair', columns='stat', values='corr')
print(t3.round(3).to_string())
t3.to_csv(os.path.join(OUT, 'c1_null_correlation.tsv'), sep='\t')


# %% T4: rank of the order-2 H table (positions x letters), PWM vs epistatic pair
def h2_table(x, pos):
    """H_ij(a, b) for i != j in pos, all letters: f_ij(a,b)/(N q_a q_b) - f_i(a)/(N q_a) - f_j(b)/(N q_b) + 1."""
    n = len(x)
    oh = np.stack([(x[:, p, None] == np.arange(4)) for p in pos], axis=1).astype(float)   # n x P x 4
    m1 = oh.mean(axis=0) / Q                                                                # P x 4
    m2 = np.einsum('nia,njb->iajb', oh, oh) / n / (Q[None, :, None, None] * Q[None, None, None, :])
    H = m2 - m1[:, :, None, None] - m1[None, None, :, :] + 1
    return H                                                                                # P x 4 x P x 4


def rank1_share(H):
    """Fraction of the off-diagonal-block energy captured by the best rank-1 symmetric fit v v^T (ALS)."""
    P = H.shape[0]
    M = H.reshape(P * 4, P * 4)
    mask = np.kron(1 - np.eye(P), np.ones((4, 4))).astype(bool)
    v = np.linalg.eigh(np.where(mask, M, 0))[1][:, -1] * 0.1
    lam = 1.0
    for _ in range(300):    # alternating fit of lam * v v^T on the masked entries
        Mf = np.where(mask, M, lam * np.outer(v, v))
        w, V = np.linalg.eigh(Mf)
        i = np.argmax(np.abs(w))
        lam, v = w[i], V[:, i]
    fit = lam * np.outer(v, v)
    return 1 - ((M - fit)[mask] ** 2).sum() / (M[mask] ** 2).sum(), np.sort(np.abs(np.linalg.eigvalsh(Mf)))[::-1][:4]


cons, cols = make_motif('pwm')
x_pwm = sample_mix(400_000, 0.05, cols)
x_null = sample(400_000, np.array([Q] * W))
sigma = np.array([3, 2, 0, 1])           # epistatic pair at positions (0, 5): b5 = sigma(b0), flat margins
x_epi = sample(400_000, np.array([Q] * W))
k = rng.random(len(x_epi)) < 0.05
b0 = rng.choice(4, size=int(k.sum()), p=Q)
x_epi[k, 0] = b0
x_epi[k, 5] = sigma[b0]
for name, xx, pos in [('pwm (6 informative positions)', x_pwm, range(6)), ('epistatic pair (0,5)', x_epi, [0, 5]),
                      ('null', x_null, range(6))]:
    share, ev = rank1_share(h2_table(xx, list(pos)))
    print(f'[T4] {name}: rank-1 share of order-2 H energy = {share:.3f}; top |eigenvalues| {np.round(ev, 3)}')
