# Registered protocol LN2: closing the low-N representation gap, then confirming STREME parity (2026-10-10)

**Context.** LN1 (dev): HASTUR − STREME located top-10 = −0.14 [−0.25, −0.04] at 250 peaks, +0.08 / +0.06 (n.s.) at
500 / 1,000. Diagnosis (on LN1 dev units, therefore development): in 12 of 13 STREME-only units HASTUR's call is at
rank 1–3 with PWM similarity 0.75–0.93, and fails the located rule on consensus precision (degenerate letters at exact
truth positions, or one edge position trimmed). PWM similarity at 250 peaks: 0.844 vs STREME 0.853.

## Part A: development (dev units only)

**Variants** (one factor at a time, then combined):
- **A1 site-based final matrix:** after the representation is chosen, recompute it from the called sites (each fg
  sequence's best window above the held-out-chosen threshold), as MEME and STREME report matrices.
- **A2 trimming threshold:** edge-column G test at p = 0.01 (G ≥ 11.34) instead of p = 0.001.

**Dev panel:** the LN1 units (24 dev factors × N 250/500/1,000 × 3) plus the N = 5,000 dev inputs; HASTUR only
(comparators cached).

**Acceptance** (located top-10, paired against HASTUR v6 plain): a variant is kept if it improves N = 250 by ≥ 0.05 and
loses ≤ 0.03 at every other N; P5 top-10 and PWM similarity are reported and must not drop by more than 0.03.
Regression checks: the sim panel (no drop > 0.05 in top-1) and the motif-free panel (false `sig` calls per run must
not rise; these variants do not touch detection).
**Stopping rule:** at most 4 variant rounds; stop at the first round whose best change is < 0.03 at N = 250.

## Part B: confirmation (locked test factors, run once, configuration frozen and hashed before the run)

**Units:** the 34 ChIP **test** factors (`splits/chip_split.tsv`), built with `data/prep_chip.py test`; N = 250, 500,
1,000 (3 random subsets, LN1's rule) and 5,000. Tools: frozen HASTUR, STREME, HOMER (LN1 settings).

**STREME parity criterion (all must hold):** at every N, HASTUR − STREME located top-10 ≥ −0.05 (paired mean) with the
95% bootstrap lower bound ≥ −0.10; the same on the P5 rule; and PWM similarity not lower by more than 0.02.
**Speed:** median runtime ≥ 3× lower than STREME at every N.
**Calibration:** on the motif-free panel, false `sig` calls per run under the CpG shift ≤ 1/5 of STREME's.
**Outcome statement:** "matches STREME's accuracy, faster and better calibrated" only if all three hold; otherwise
report which failed.
