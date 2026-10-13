# Registered protocol: low-N ChIP panel (LN1), registered 2026-10-10 before any low-N run

**Question.** Does HASTUR gain accuracy over STREME (and HOMER) when data are limited, and does summit weighting
(THEORY Prop. 5) add to it? At N = 5,000 peaks the dev result is parity (THEORY §11b).

**Units.** The 24 ChIP **dev** factors (`splits/chip_split.tsv`; test factors stay locked). Inputs: `prep_chip.py`
outputs (top 5,000 conservative-IDR summits ±100 bp; GC-matched genomic background).

**Downsampling.** For N ∈ {250, 500, 1000}: 3 replicates r ∈ {0, 1, 2}; foreground = N peaks drawn uniformly without
replacement from the 5,000 (numpy default_rng(1000·N + r + factor index)); background = 10·N sequences drawn the same
way from the 10× background pool. Sequences keep their summit-centred coordinates (position 100 = summit).

**Arms (configurations fixed here).**
| arm | command |
|---|---|
| HASTUR plain | `proto/hlat.py` VERSION v6, defaults (k 4–5, w 13, L 8, 12 seeds per family, Markov order 1), no position weights |
| HASTUR centred | the same with `--center 40` (summit-relative window weights, Prop. 5; fixed prior, not fitted) |
| STREME | 5.5.9, `--dna --nmotifs 10` (defaults otherwise: widths 8–15) |
| HOMER | `findMotifs.pl fasta -noknown -len 8,10,12 -p 2` |

HASTUR also runs at N = 5,000 (the existing dev inputs) in both arms; comparators there are the existing runs.

**Readouts** (`bench/evaluate.py`, truth = the split's consensus): located (`contains`) top-10 and top-1, P5 rule
top-10, best PWM similarity, each with its shuffled-truth chance; per tool; paired by unit (factor × N × replicate).
Expressibility: truths longer than a tool's maximum width (REST, 20 bp, for STREME and HOMER) are flagged and the
headline is also reported without them.

**Decision rules (breakpoints), on located top-10 recall, paired over the 24 × 3 units per N:**
- **Accuracy thesis alive** if, at N = 250 or 500, the better HASTUR arm beats STREME by ≥ 0.15 (paired mean
  difference, 95% bootstrap CI excluding 0).
- **No accuracy gain** if the difference is within ±0.05 at every N: the claim becomes speed + calibration only.
- **Summit weighting works** if HASTUR centred beats HASTUR plain by ≥ 0.10 at any N ≤ 1,000 (same CI rule).
- In between: report as is; no further low-N variants until a new protocol is registered.

**Calibration readout (alongside).** The motif-free panel (`bench/null_make.py`): false `sig` calls per run per tool.
