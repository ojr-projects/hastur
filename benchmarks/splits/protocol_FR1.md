# Registered protocol FR1: the matched-threshold frontier (registered 2026-10-11, before any FR1 run)

**Question.** Does HASTUR have a better sensitivity / false-call frontier than STREME (and HOMER), or only a
different operating point? A tightened STREME threshold might reach HASTUR's false-call rate; the question is how
much sensitivity each tool keeps there. This is the decisive comparison for the "calibrated, fast, high-throughput"
positioning (ROADMAP §6).

## Scores swept
| tool | score per reported motif | run settings |
|---|---|---|
| HASTUR | detection z (ranking statistic; Yellow Sign = family Bonferroni) | the configuration under test |
| STREME | held-out E-value | `--thresh 1 --nmotifs 25` (report everything; thresholds applied afterwards) |
| HOMER | reported log p-value | `-len 8,10,12`, all reported motifs |

## Positives (sensitivity)
Development: ChIP dev under the matrix readout (`splits/chip_dev_matrix_truth.json`; low-N panel and N = 5,000) and
SELEX dev (JASPAR matrices). A unit is a hit at threshold t if any call with score passing t matches the truth
(reported under `contains`, P5 and PWM similarity ≥ 0.85). Confirmation: ChIP test2 and the SELEX test split, run once
with frozen settings.

## Nulls (false calls per run), each with 20 runs per kind
1. **Simulated** (`bench/null_make.py`, stationary burn-in): iid, same-chain Markov, CpG shift.
2. **Real human genome:** (a) exchangeable: two disjoint random sets of genomic windows (201 bp, no N), N = 500 and
   2,000 foreground, 10× background; (b) composition shift: foreground = random windows with GC in [0.55, 0.65],
   background = random genome-wide windows.
3. **Real bacterial genomes** (AT-rich and GC-rich): for C. jejuni (AT-rich), E. coli and B. cenocepacia (GC-rich):
   exchangeable splits of random intergenic windows, and of random coding windows (100 bp).
Any reported motif passing t in a null run is a false call.

## Readout
For each tool, the curve of sensitivity against mean false calls per null run, as t varies. Headline: sensitivity
at matched false-call rates of 0.05, 0.1 and 0.5 per run (separately for each null family, and pooled), with bootstrap
CIs over units.
**Better frontier** if, on real nulls and real positives, HASTUR's sensitivity at 0.1 false calls per run exceeds
STREME's by ≥ 0.05 with the CI excluding 0. **Operating point only** if within ±0.05. Otherwise report as is.
Also reported: STREME's own calibration (fraction of its false calls with E ≤ 0.05), replicate reproducibility across
the low-N subsamples, and accuracy by truth-motif length.
