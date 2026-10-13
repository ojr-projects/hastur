# Registered protocol LN4: exact small-count calibration (2026-10-10, before any LN4 run)

**Problem.** False `sig` calls per run on motif-free pairs: iid 0.0, same-chain Markov 0.2, CpG shift 0.2 (target
≤ 0.05). Residual false calls are low-count cells where the normal approximation's tail is too light.
**Change (E1, flag `--exact`).** For lattice cells with F + B ≤ 100, z is computed from the exact conditional test:
given T = F + B, F ~ Binomial(T, Ef / (Ef + Eb)) under the null (Ef, Eb the double-centred expectations), mid-p,
converted to a signed z. Larger cells keep the current statistic. Genomic control and Bonferroni are unchanged.

**Development units:** the existing motif-free panel (`bench/null_make.py`, 30 pairs) for calibration; accuracy
regression on SELEX dev and ChIP dev low-N (located / P5 top-10 must not drop by more than 0.02 against the same
configuration without `--exact`).
**Confirmation (calibration "test"):** a fresh motif-free panel generated with new seeds (`null_make.py` with seed
offset 900000; 10 pairs per kind), run once with the frozen configuration. Pass if false `sig` calls per run ≤ 0.05
for each kind (95% upper bound reported).
