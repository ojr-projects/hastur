# Registered protocol LN5: matrix accuracy (registered 2026-10-11, before any arm run)

**Why.** Part B (post hoc) and the low-N matrix readout show HASTUR behind STREME on stating motifs (P5, PWM
similarity) while level on locating them. The planned diagnosis was "full-length / width". The step-0 anatomy on
low-N dev (`bench/repr_anatomy.py`, results in ROADMAP §5) refuted it.
- **Width:** HASTUR's best calls cover the truth's informative columns as well as STREME's (0.95–0.96 vs
  0.94–0.95), and they are wider.
- **Where the gap is:** all of it is in the covered columns (column similarity 0.895–0.910 vs 0.921–0.933).
- **What the column error is not:** softening toward uniform or toward the peaks' composition, ±1 phase blur, or a
  smaller site count (HASTUR's site count is 93% of STREME's).
- **Contamination:** contamination by an earlier related motif (JUND TRE/CRE) explains 29% of the gap.
- **What helps:** one more hard-EM site-refit round on HASTUR's output matrix (best window per sequence at the
  enrichment-optimal threshold) raises similarity 0.877 → 0.893 (N = 250) and 0.885 → 0.897 (N = 500), against
  STREME's 0.895 and 0.908. Further rounds add nothing.

## Arms
The base is frozen v8 (`--sites --gtrim 0.01`); prototype v11. Calls, ranking and significance are unchanged by
construction, since keep/dedupe decisions are made on the pre-refit matrix (LN2 A3).

| arm | tool key | flags added | tests |
|---|---|---|---|
| L1 | `hlat_l1` | `--refit 2` | a second site-refit round |
| L2 | `hlat_l2` | `--plainz` | site threshold from plain fg-vs-bg enrichment instead of the double-centred contrast |
| L12 | `hlat_l12` | `--refit 2 --plainz` | both |

## Dev data and readouts
- **ChIP dev under the matrix readout** (`splits/chip_dev_matrix_truth.json`): low-N panel (N = 250, 500, 1,000)
  and N = 5,000.
- **SELEX dev** (JASPAR matrices). This is the guard against fitting ChIP truths: site matrices hurt SELEX in LN2.
- **Per unit:** best PWM similarity among the top 10, P5 top-10 and located top-10 (`bench/evaluate.py`).
- **Sanity:** the motif-free panel's significant calls per run must be identical to v8.

## Keep rule (each arm against frozen v8 `hlat_a3`, paired by unit, bootstrap 95% CIs)
An arm is kept only if all of these hold:
1. **ChIP similarity:** best PWM similarity on ChIP dev (all N pooled) improves by ≥ +0.005, with the CI
   excluding 0.
2. **ChIP P5 and located:** P5 top-10 and located top-10 do not drop by more than 0.02 at any N.
3. **SELEX:** best PWM similarity does not drop by more than 0.01, and P5 top-10 does not drop by more than 0.02.
4. **Nulls:** significant calls per motif-free run are unchanged.

Among kept arms, the larger ChIP similarity gain wins, and the simpler arm wins a tie within 0.003.

The kept arm is frozen as **v11** (flags and git hash recorded in `splits/freeze_LN5.md`) and goes forward to the
FR1 confirmation on test2 and the SELEX test split. It replaces frozen v8 there only if kept here. Otherwise v8
stays the confirmation configuration.
