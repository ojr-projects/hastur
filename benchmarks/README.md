# Benchmarks

This is the record of every benchmark and comparison HASTUR went through. The protocols were registered (hashed in
`splits/SHA256SUMS`) before their outputs were read, and each one names its own decision rule.

## Confirmations and comparisons (read these first)

| what | protocol | result |
|---|---|---|
| Frozen configuration F2 (= this engine) | `splits/freeze_v18.md` | `results/ln7e_report.txt`, `results/frontier_ln7e.txt` (dev evidence) |
| **FR1: sensitivity vs false calls per run against STREME, held-out ChIP + SELEX** | `splits/protocol_FR1*.md` | `results/frontier_conf_f2*.txt`, `results/frontier_masked_f2*.txt` |
| Small-N calibration (fresh panel) | `splits/protocol_lowN.md`, `splits/protocol_FR1_amend*.md` | `results/smallN_conf_f2.txt`, `results/smallN_conf5_f2.txt`, `results/lowN_eval.tsv` |
| Reproducibility across subsamples | `splits/protocol_D1*.md` | `results/repro_conf_f2.txt`, `results/d1_report*.txt` |
| Repeat-masked secondary | `splits/secondary_masked.md` | `results/frontier_masked_f2*.txt` |
| *E. coli* regulon oracle | — | `results/regulon_oracle.txt` |
| Speed (Rust vs the Python prototype) | — | `results/speed_*.txt`; scripts `bench_full.sh`, `bench_local.sh` (these call the retired Python prototype) |

## Real-data quick looks (E-series, 2026-10-12)

Each quick look ran HASTUR first; comparators were added only where a win looked possible.

| experiment | data | HASTUR | comparator | verdict |
|---|---|---|---|---|
| E1 | *E. coli* iModulons, 32 evaluable modules (`splits/protocol_E1*.md`) | 8/32 native, 0.044 false calls/run | STREME 30/32 native but 9.6 false calls/run; level with HASTUR at matched rates | calibrated, no frontier win |
| E1, known-motif arm | same (`results/e1_known_motif.tsv`) | own known-motif test 16/32 | SEA 29/32 at 0.059 false calls/run | SEA owns known motifs |
| E1b | *P. putida* | 0.032 false calls/run; 12/48 modules called, FNR box among them | not run | calibration only |
| E6 | 100 null screens (`splits/protocol_E6.md`) | per-run calibration predicts screen counts (4.7 per screen) | not run | p_adj is anti-conservative in the deep tail: BH across runs fails (46% of null screens) |
| E5 | *P. falciparum*, 13% GC nulls | 0.04 false calls/run | STREME calibrated at N ≥ 50, fails only at N = 25 | composition does not separate the tools; small N does |
| E2 | 23 TF knockouts | own motif in 1/14 scoreable contrasts | not run | no |

## Matched backgrounds (M-series)

- **MA1** (`splits/protocol_MA1.md`):
  - GC-selected human nulls: false calls per run fell from 3.70 / 2.70 to 0.05 / 0.05;
  - exchangeable nulls stayed ≤ 0.10;
  - held-out ChIP sensitivity went from 23/27 to 22/27.
  - Not adopted as the default; the builder was kept as an option, then retired and removed on 2026-10-13.
- **MA2:** 94-factor test3 ChIP, stopped at 58 factors. Sensitivity was 40/58 in every arm.
- **Status:** retired. The matched-background code is no longer in this repository.

## Development rounds (LN-series)

The `results/ln*_report.txt` and `results/frontier_ln7*.txt` files, with `splits/protocol_LN*.md`, are the dev cycles
that produced F2. The keep-rule scripts for these rounds were not carried over; the reports are the record.

## Harness

| path | what |
|---|---|
| `bench/run_tools.py` | runs each tool on `<unit>/fg.fa` + `bg.fa`; the HASTUR arm is `hlat_f2r`, set by `HASTUR_BIN` (default `~/.cache/hastur_target/release/hastur`) |
| `bench/*_make.py` | build the panels: simulation, nulls, real nulls, small-N, regulons, iModulons, E2/E6 |
| `bench/evaluate.py`, `null_eval.py`, `frontier.py`, `repro.py`, `regulon_*`, `imodulon_*` | scoring and reports |
| `eval/score.py` | match rules, shuffled-truth chance baseline, PWM similarity |
| `data/prep_chip.py`, `prep_selex.py` | real-data preparation (ENCODE ChIP, HT-SELEX) |

**Notes on the harness:**
- `run_tools.py`'s `hlat`, `hlat_a*`, … arms called the retired Python prototype. Only `hlat_f2r` and the comparators
  (`streme`, `meme`, `xstreme`, `sea`, `homer`, `comse`) still run.
- Output files keep the historical name `hlat.tsv` / `hlat.meme`.
- Raw run outputs are not in git. On the workstation they are in `~/runs/` and `~/data/{e1,e2,e5,e6,imodulon_ppu}`; the
  rest is on the EC2 box's persistent volume.
