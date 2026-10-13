# Registered protocol LN6: calibration at small N (registered 2026-10-11, before any arm run)

**Why.** D1 (`results/d1_report.txt`) showed 1.5 native false calls per run on random E. coli promoter sets of 10–206
regions. Gapped-lattice false calls scale with smallness: 2.29 per run at ≤ 15 regions, 0.23 at 61–250. The soft
L-mer table adds 0.2–0.3 per run, with z up to 67.

**Mechanism.** Both tables use a Wald z for a log ratio, with the variance taken from the *observed* counts. That is
badly anti-conservative when the expected fg count is small:
- one observed 8-mer against 0.09 expected gives z ≈ 4.2, where the exact Poisson tail gives z ≈ 1.4;
- three observed against 0.2 expected gives z ≈ 5.4 vs 3.1.

Genomic control cannot correct this, because λ is estimated only from cells with expected counts above 5, while
candidate seeds come from all cells. LN4's exact test (E1) applies only when F + B ≤ 100, which never holds with a
large background. No earlier calibration panel had N < 500.

## Arms (prototype v12; the base is frozen v8 `--sites --gtrim 0.01`)
| arm | tool key | flags added | change |
|---|---|---|---|
| C1 | `hlat_x1` | `--exact-all` | exact conditional binomial (mid-p) test for every gapped cell |
| C2 | `hlat_x2` | `--lmer-valid` | soft L-mer seeds only from cells with expected fg and bg counts > 5 (the domain λ already uses) |
| C12 | `hlat_x12` | both | |

## Dev data (independent of D1: no E. coli)
- **Small-N null panel** (`bench/upnull_make.py`): random sets of upstream regions ([gene start − 250, + 50],
  merged; median 300 bp) from P. putida (GC-rich), C. glutamicum and S. aureus (AT-rich). Sizes 10, 15, 20, 30, 50,
  100, 200, with 10 replicates each: 210 runs. bg = every other region.
- **Simulated nulls** (`null_v2`, N = 2,000).
- **Positives:**
  - ChIP dev small-N panel (`LOWN_NS=25,50,100 bench/lowN_make.py`, 216 units);
  - low-N panel (250, 500, 1,000; 216 units);
  - N = 5,000 (24 units);
  - all under the matrix readout.

## Keep rule (against frozen v8 `hlat_a3`)
An arm is kept only if all of these hold:
1. **Small-N nulls overall:** native false calls per run ≤ 0.1, pooled over all 210 runs.
2. **Small-N nulls by size:** ≤ 0.2 in every size bin (≤ 15, 16–30, 31–60, 61–200).
3. **Simulated nulls:** false calls per run not above v8's.
4. **ChIP at N ≥ 250:** P5 top-10 and located top-10 do not drop by more than 0.03 at any N.
5. **Small-N ChIP (25, 50, 100):** the sensitivity loss is reported, not gated. A calibrated tool is expected to call
   less at N = 25.

Among kept arms, the simplest wins (C1 or C2 alone over C12 when both pass). If none passes, the closest arm is
reported with its false-call profile by size, and small-N calibration stays an open problem. D1 is not re-run on
any arm until it is frozen. The re-run is then a separate, registered D1 confirmation on fresh random-promoter
nulls.
