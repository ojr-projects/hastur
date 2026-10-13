# %% c7_spanning_shapes.py: translation pooling and spanning shapes (THEORY.md sections 5.1-5.2)
"""
Unaligned sequences with a planted PWM of L = 6 informative positions. For every spacing d, the order-2 H table over
letter pairs (translation-pooled contrast, fg vs bg) is predicted to be
    H_d(a, b) = pi_w * sum_u delta_u(a) delta_{u+d}(b)      (sum over placements u with u, u+d inside the motif)
so it has rank <= L - d: rank one at the spanning spacing d = L - 1. With sites on both strands the reverse-complement
motif adds its own placements (rank <= 2 at d = L - 1 for a non-palindromic motif).
Readouts per spacing: correlation of observed with predicted H_d; share of energy in the top 1 and top 2 singular values.
Run: ~/micromamba/envs/nanomotif/bin/python theory/checks/c7_spanning_shapes.py
"""
import numpy as np

rng = np.random.default_rng(7)
LSEQ, N, PI = 100, 60_000, 0.5
Q = np.full(4, 0.25)
cons = np.array([2, 0, 3, 1, 2, 3])            # GATCGT (not a palindrome)
pwm = np.full((6, 4), 0.06)
pwm[np.arange(6), cons] = 0.82
L = len(pwm)


def sim(strands):
    X = rng.integers(0, 4, (N, LSEQ))
    cdf = np.cumsum(pwm, axis=1)
    for i in np.flatnonzero(rng.random(N) < PI):
        site = (rng.random((L, 1)) > cdf).sum(axis=1)
        if strands == 2 and rng.random() < 0.5:
            site = 3 - site[::-1]
        t = rng.integers(0, LSEQ - L + 1)
        X[i, t:t + L] = site
    return X


def pair_H(fg, bg, d):
    def tab(X):
        a, b = X[:, :-d].ravel(), X[:, d:].ravel()
        return np.bincount(a * 4 + b, minlength=16).reshape(4, 4) / len(a)
    Rf = tab(fg) / tab(bg)
    single = lambda X: np.bincount(X.ravel(), minlength=4) / X.size
    r1 = single(fg) / single(bg)
    return Rf - r1[:, None] - r1[None, :] + 1


delta = pwm / Q - 1
for strands in (1, 2):
    fg, bg = sim(strands), rng.integers(0, 4, (N, LSEQ))
    pi_w = PI / (LSEQ - 1)          # sites per pair position, approximately
    print(f'== sites on {"one strand" if strands == 1 else "both strands"}')
    print(' d  rank bound  corr(obs, pred)  top-1 share  top-2 share  |H|_F')
    for d in range(1, L + 2):
        H = pair_H(fg, bg, d)
        pred = np.zeros((4, 4))
        motifs = [delta] if strands == 1 else [delta, pwm[::-1, ::-1] / Q - 1]   # forward and reverse-complement sites
        w = 1.0 / len(motifs)
        for dl in motifs:
            for u in range(L - d):
                pred += w * pi_w * np.outer(dl[u], dl[u + d])
        s = np.linalg.svd(H, compute_uv=False)
        share = (s ** 2).cumsum() / (s ** 2).sum()
        corr = np.corrcoef(H.ravel(), pred.ravel())[0, 1] if np.abs(pred).sum() > 0 else float('nan')
        bound = max(L - d, 0) * strands
        print(f'{d:2d}  {bound:10d}  {corr:15.3f}  {share[0]:11.3f}  {share[1]:11.3f}  {np.sqrt((H ** 2).sum()):.4f}')
