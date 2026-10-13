# Registered protocol LN6 round 2 (registered 2026-10-11, before any round-2 arm run on the decision panel)

**Round 1** (`results/ln6_report.txt`): no arm was kept.
- **C12** (exact test for every gapped cell, plus an expected-count > 5 cutoff on L-mer cells) cut small-N false
  calls from 1.75 to 0.086 per run. It failed on the 31–60 bin (0.27) and on located −0.053 at ChIP N = 250.
- **The residual false calls came from clumping:** a tandem repeat gave four calls in one run. Occurrence counts
  treat hits from one locus as independent evidence.
- **The cutoff switches the L-mer table off** below about 1,000 sequences.

**A rejected design.** A quasi-Poisson exact test on the soft L-mer counts (`--lmer-qp`) was rejected before
registration. Its null dispersion κ ≈ 0.35 makes one real occurrence count as about 2.9, which is anti-conservative
in the tail; it gave marginal false calls in a sanity run.

## Arms (prototype v13; the base is frozen v8 `--sites --gtrim 0.01`; each also has `--exact-all`)
| arm | tool key | flags added | change |
|---|---|---|---|
| R1 | `hlat_y1` | `--exact-all --seqcount` | both tables count sequences holding a cell, not occurrences. Null expectation Σ_s 1 − exp(−2 W_s p); exact test |
| R2 | `hlat_y2` | `--exact-all --seqcount --lmer-ball` | R1, plus L-mer cells tested by exact radius-1 Hamming-ball counts (integer; pools one-mismatch variants) instead of the soft Wald z |
| R3 | `hlat_y3` | `--exact-all --lmer-ball` | the ball test without presence counting (isolates the two changes) |

## Decision data
- **A fresh small-N null panel** (`bench/upnull_make.py <root> 1`): same three genomes, sizes and replicate counts
  as round 1 (210 runs), new random sets. The round-1 panel was used for diagnosis, so it does not decide.
- **The same simulated nulls and ChIP dev panels** as round 1 (N = 25 to 5,000, matrix readout).

## Keep rule
The keep rule is identical to round 1's (`splits/protocol_LN6.md`):
- small-N nulls ≤ 0.1 pooled and ≤ 0.2 in every size bin;
- simulated nulls not above v8;
- P5 and located top-10 within −0.03 of v8 at every ChIP N ≥ 250;
- small-N ChIP reported, not gated.

The simplest passing arm wins (R1, then R3, then R2). The winner is combined with LN5's L2 at the freeze
(FR1 amendment 2).
