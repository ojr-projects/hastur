# FR1 amendment 3 (2026-10-11, before any confirmation output is read)

1. **Disclosure.** Frozen v8 failed FR1 on dev: −0.314 at 0.1 false calls per run on real nulls
   (`results/frontier_dev.txt`). Development paused, then resumed as an exploratory line (LN7, ROADMAP §5). The
   configuration confirmed here is **F2** (`splits/freeze_v18.md`).
   - **F2's dev numbers are optimistic.** It was developed on dev panels that include the FR1 dev panel. Its
     calibration gates used other dev panels.
   - **v8 is also reported** on the confirmation data.
2. **Confirmation data, as amendments 1 and 2 define it:**
   - positives: ChIP test2 (27 factors), its low-N panel and N = 5,000 under the matrix readout
     (`splits/chip_test2_matrix_truth.json`), and the SELEX test split (JASPAR matrices);
   - nulls: fresh replicates (`null_conf`, seed offset 910000; `realnull_conf`, SEED_OFFSET 5000).
3. **Intergenic nulls are excluded from the decision cell, for every tool,** because of a construction flaw:
   100-bp windows were sampled with overlap from small intergenic spaces, so fg windows share sequence (found on
   the dev panel: 376 of 500 C. jejuni windows). A second readout with them included is reported alongside.
4. **Comparators:**
   - STREME (`streme_all`) was run on the confirmation data earlier and is read only now.
   - HOMER was dropped after the FR1 preliminary (user decision, to save compute).
5. **Scores per call:**
   - F2: −log10 of the family-wise adjusted p;
   - v8: −log10 of the one-sided normal p of its z;
   - STREME: −log10 E (or S).
6. **The verdict rules are unchanged** (`protocol_FR1.md`). The decision is F2 − STREME sensitivity at 0.1 false
   calls per run, real nulls, all positives:
   - **better frontier** if ≥ +0.05 with the CI above 0;
   - **operating point only** if within ±0.05;
   - otherwise reported as it stands.
