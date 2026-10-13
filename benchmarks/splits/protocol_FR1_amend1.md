# FR1 amendment 1 (2026-10-11, before any confirmation run)

1. **All tools on one machine.** FR1 dev runs every tool on the box (`infra/box/run_fr1_box.sh`), not the
   workstation, so that runtimes are comparable.
2. **The STREME parser reads `E=` or `S=`.** This is D1 amendment 1. No FR1 dev unit has `S=` output.
3. **Confirmation positives:**
   - ChIP test2 at N = 5,000, for the 27 factors whose inputs built. The other 26 have no conservative IDR peaks
     in their cell line.
   - A test2 low-N panel: N = 250, 500, 1,000 × 3 replicates (`bench/lowN_make.py`, `CHIP_SPLIT`). ZNF416 has
     780 peaks, so its N = 1,000 units use all 780.
   - The SELEX test split (186 TFs).
   - Matrix readout for test2 (`splits/chip_test2_matrix_truth.json`).
4. **Confirmation nulls are fresh replicates,** so no null run is shared with the dev frontier:
   - simulated: `bench/null_make.py` seed offset 910000;
   - real genomes: `bench/realnull_make.py` with `SEED_OFFSET=5000`;
   - the same kinds and counts as dev (20 runs each).
5. **HASTUR configuration:** frozen v8 (`hlat_a3`), unless LN5 keeps an arm, in which case that arm (frozen v11) is
   the configuration confirmed, and v8 is also reported. Comparator runs do not depend on this choice. No
   confirmation output is read until the LN5 decision is recorded.
