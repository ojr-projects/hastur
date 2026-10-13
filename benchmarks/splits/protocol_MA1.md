# Registered protocol MA1: PROMETHEUS Mode A, matched backgrounds (registered 2026-10-12, before any run on these units)

**Question (docs/PROMETHEUS_PLAN.md step 3).** Do genome-matched backgrounds, built from the hg38 covariate index,
remove F2's false calls on composition-shifted (GC-selected) nulls without costing ChIP sensitivity? Does the gain
belong to the background, so that STREME gets it too, or to HASTUR?

**Dev evidence (disclosed; not reused here).**
- **Nulls:** on the existing GC-selected real nulls (realnull_conf and realnull_dev3, 80 runs), false calls per run
  fell from 3.58 / 2.72 (N = 2,000 / 500) to 0.00 / 0.05 with the matched background. Exchangeable nulls stayed at
  ≤ 0.10 → 0.00.
- **ChIP:** dev ChIP native sensitivity was unchanged (17/19 for F2 and both arms).
- **Bucket design:** informed by step 2b. Unmasked false calls on GC-selected nulls sat in SINE copies; masked ones sat
  in GC-rich unique sequence.

## Units (`infra/box/run_modeA_round.py`)
- **Fresh nulls:** `bench/realnull_make.py` with `SEED_OFFSET=9000`, at human exch and GC N = 500 and 2,000,
  20 replicates each (80 runs).
- **Freshness:** measured before the round is read, as the share of fg windows (chrom, start) shared with any
  earlier human panel (realnull, realnull_conf, realnull_dev3). It must be < 1%.
- **ChIP:** the 27 held-out test2 sets (`chip5k_test2`; truth `splits/chip_test2_matrix_truth.json`). They were not
  used in the dev look.
- **Matched backgrounds (`rust/prometheus/tools/match_sets.py`, unchanged from the dev look):**
  - K = 10 index tiles per fg region, drawn within the region's bucket;
  - buckets: GC in 0.05 bins, CpG o/e (4 bins), soft-masked share (3 bins), SINE at the centre when repeat-rich, and
    TSS distance (4 bins, ChIP only);
  - coarsening rule as in the script;
  - tiles overlapping fg regions or ENCODE peaks are never drawn.

## Arms
| arm | statistic | background | permutation |
|---|---|---|---|
| `hlat_f2r` | F2 | original | global |
| `hlat_mA1` | F2 | matched | global |
| `hlat_mA2` | F2 | matched | stratified (`--strata`) |
| `streme` | STREME default | original | — |
| `streme_mA` | STREME default | matched | — |

## Readouts
1. False calls per null run at the native threshold, per null kind and arm, with Poisson 95% intervals.
2. Test2 ChIP native sensitivity per arm (hit rules as FR1: contains, P5 or PWM similarity ≥ 0.85), with paired
   differences against `hlat_f2r`.
3. The STREME pair: false calls per run and sensitivity, original vs matched background.
4. Matching diagnostics: coarsening events, regions short of K, strata counts.

## Decision rule (for adopting Mode A into PROMETHEUS)
- **Adopt `hlat_mA1` (or `hlat_mA2`)** only if all three hold:
  - GC-selected false calls ≤ 0.1 per run (both N);
  - exchangeable false calls ≤ 0.1 per run;
  - test2 ChIP sensitivity within 0.03 of `hlat_f2r`.
- **Prefer `mA2` over `mA1`** only if it is strictly better on one of these with no loss on the others.
- **Ownership of the gain:**
  - If `streme_mA` also reaches ≤ 0.1 false calls per run on GC-selected nulls, the gain is reported as the
    background's ("matched backgrounds fix composition shift for both tools").
  - Otherwise it is reported as HASTUR's.
