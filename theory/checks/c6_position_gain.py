# %% c6_position_gain.py: data value of a summit-relative position axis (THEORY.md Proposition 9)
"""
Dense-regime gain from weighting windows by their position r relative to the peak summit: the optimal weight is the
density ratio w = f/g, worth (1 + chi2(f || g)) times the data; a hard central window is the best 0/1 weight.
Sites are Gaussian around the summit (sd = summit precision); windows are uniform over a 200-bp peak.
Run: ~/micromamba/envs/nanomotif/bin/python theory/checks/c6_position_gain.py
"""
import numpy as np
from scipy.stats import norm

PEAK = 200.0
r = np.linspace(-PEAK / 2, PEAK / 2, 20001)
dr = r[1] - r[0]
g = np.full_like(r, 1 / PEAK)
print('summit sd (bp)  optimal weighting (x data)  best hard window (x data)  half-width (bp)')
for sd in [10, 15, 20, 30, 50]:
    f = norm.pdf(r, 0, sd)
    f /= f.sum() * dr
    soft = (f ** 2 / g).sum() * dr
    hw = np.arange(2, 100)
    hard = (norm.cdf(hw / sd) - norm.cdf(-hw / sd)) ** 2 / (2 * hw / PEAK)
    print(f'{sd:14d}  {soft:26.2f}  {hard.max():25.2f}  {hw[np.argmax(hard)]:15d}')
