# Registered protocol E6: screen-level FDR on fully null screens (registered 2026-10-12, before any run)

**Question (ROADMAP_v2 E6).** If each run is calibrated, does the calibration carry over to a whole screen? Can a
screen-level FDR be controlled across runs, and does that hold for HASTUR's p_adj and for STREME's E-values?

## Data (`bench/nullscreen_make.py`)
- **Null screens:** 100 of them. Each mirrors the E1 screen, with one random set of promoter regions per included E1
  module (94), at that module's size. The population is the same 1,942 regions, with bg = every other region.
- **No real modules:** every call in a null screen is false.
- **Seeds:** sha256 of `E6|<screen>|<module>`, disjoint from E1 and D1. These screens are fresh relative to both.

## Tools (workstation, same machine for both; `bench/run_tools.py`)
- **`hlat_f2r`:** frozen F2 (Rust), with `--threads 2`.
- **`streme`:** default settings, `--nmotifs 10`.
- **Native thresholds:** as in E1, p_adj ≤ 0.05 for F2, and E (or S) ≤ 0.05 for STREME.
- **HOMER** may be added later as filler. If it is, it is reported descriptively only.

## Readouts (`bench/nullscreen_report.py`)
1. **False calls per null screen at the native threshold.** Report the mean, the distribution, and P(screen has ≥ 1
   false call). Compare with the E1-derived prediction: E1 false calls per null run × 94.
2. **Screen-level FDR control.**
   - In each screen, take each run's minimum significance, converted to a p-value: F2 p_adj, or STREME min(1, E or
     S) used as a p-value.
   - Apply Benjamini–Hochberg across the screen's 94 runs at q = 0.05 and q = 0.10.
   - Under the full null, FDR equals P(any discovery). The realised FDR is the share of null screens with ≥ 1
     discovery, with its 95% binomial (Clopper–Pearson) interval.
3. **Screen-level FWER (Bonferroni over 94 runs)** at 0.05: the share of screens with any discovery.
4. **Runtime per screen,** per tool.

## What may be claimed
- **"Screen-level FDR controlled at q"** for a tool, only if the share of null screens with any BH discovery is ≤ q,
  or its 95% interval includes q.
- **"Per-run calibration predicts screen behaviour"** for HASTUR, only if the observed mean of false calls per screen
  lies within the 95% interval of the E1-derived prediction (Poisson interval on the E1 count, scaled to 94 runs).
- **E1 screen FDR.** When E1 is read, its estimated screen FDR uses the E6 null-screen mean as the false-call
  estimate, alongside E1's own nulls.
- Otherwise results are reported as they stand.
