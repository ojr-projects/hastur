# HASTUR

**H**oeffding-**A**ligned **S**parse-**T**ensor **U**nified **R**ank-one motif discovery.

HASTUR is an exhaustive, calibrated de novo motif finder for a foreground set against a background set.
- **Counting:** it counts every sparse pattern exactly: the gapped k-of-w lattice (k = 4–5, w = 13) and the 8-mer
  table, on both strands.
- **Ranking:** cells are ranked in composition-free (double-centred) coordinates.
- **Testing:** the top cells are tested with exact statistics. Their thresholds are set by permutation max-T, so a call
  carries a family-wise adjusted p.
- **Representation:** each call is reported as the matrix that held-out data prefers (exact frame or EM PWM), followed
  by a site refit.

The theory is in `theory/THEORY.md`.

This repository holds the frozen **F2** configuration, in Rust. Development is finished; this is the version to use and
cite.

## Build and run

```bash
cargo build --release                      # .cargo/config.toml sets target-cpu=native (PEXT counting)
target/release/hastur fg.fa bg.fa -o out/hastur --threads 8 > out/hastur.json
```

**Options:**
- `-o PREFIX`: writes `PREFIX.tsv` (one row per call) and `PREFIX.meme`. Without `-o`, only detection runs.
- `--perm R`: the number of label permutations for max-T calibration (default 20).
- `--seed S`
- `--threads T` (0 = all cores)
- `--selftest`: checks the NumPy-compatible RNG.

**Outputs:**
- `PREFIX.tsv` columns: `rank consensus z sig z_heldout bits width representation source seed p_adj`. A call is
  significant when `sig = 1` (p_adj ≤ 0.05).
- stdout: a summary JSON with permutation maxima, thresholds, all seeds and timings.

## What the evidence says

The full record is in `benchmarks/`; the dated history is in `docs/history/`.

| claim | evidence | strength |
|---|---|---|
| Calibrated on exchangeable real sequence | 0.00–0.15 false calls per motif-free run on fresh human and bacterial nulls | registered confirmation (FR1) |
| Calibrated under composition shift | simulated CpG shift: 0.05 false calls per run vs STREME 3.45 | registered confirmation |
| Calibrated on small sets | random *E. coli* promoter sets (N = 10–206): 0.08 vs STREME 7.65 | dev + fresh panel |
| Better frontier with composition-shifted nulls | +0.24 [+0.20, +0.28] sensitivity at 0.1 false calls per run vs STREME | registered confirmation (FR1) |
| Level with STREME on exchangeable data | +0.01 overall | post-hoc decomposition |
| Long motifs (> 10 informative columns) | 0.31 vs 0.09 at matched precision | confirmation, secondary |
| Fast | 0.8–7 s per full run on 12 cores; 15–40× STREME per core | measured |
| Real screens (E-series) | calibrated everywhere (≈ 0.04 false calls per run), but no sensitivity advantage on small sets; SEA owns known-motif tests | quick looks, 2026-10-12 |

**Positioning:** a calibrated, very fast motif finder for screen-scale work (many small or composition-shifted
contrasts), where false calls from other tools pile up. Its strength is trust at scale. It does not find more motifs in
a single clean ChIP experiment.

**Known weak spots:**
- GC-selected human windows with a random background. Supply a composition-matched `bg.fa` instead; see below.
- Sets of about 25 sequences are calibrated but conservative.

## Retired: matched genome backgrounds (Mode A)

An hg38 matched-background builder (Python, formerly `tools/matched_background/`) was retired on 2026-10-13 and removed
from the repository. It is not part of HASTUR.
- **What it showed** (`benchmarks/splits/protocol_MA1.md`): on GC-selected human null sets, a background matched on
  GC, CpG o/e and repeat share cut false calls per run from 3.70 / 2.70 to 0.05 / 0.05 (N = 2,000 / 500).
- **Why it was retired:** the failure it fixes is narrow, and on real ChIP it gained nothing (23/27 → 22/27 on test2;
  40/58 in every arm on test3).
- **If your foreground is composition-selected,** build a GC/CpG-matched `bg.fa` with any tool you trust. The engine
  only reads the two FASTA files.

## Layout

| path | what |
|---|---|
| `src/` | the engine: counting, statistics and max-T in `main.rs`; representation in `represent.rs`; NumPy-compatible RNG in `nprng.rs` |
| `benchmarks/` | registered protocols, splits and truths, the harness, result reports (see its README) |
| `theory/` | THEORY.md and the numerical checks behind it (c1–c8) |
| `docs/` | development insights, the deep-model literature review, and `history/` (roadmaps, the PROMETHEUS plan, the Comsequent background) |

Naming, from *The King in Yellow*: the calibrated significance flag is the **Yellow Sign**. The companion genome
visualiser is **Hassering**.
