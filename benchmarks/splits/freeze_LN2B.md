# Freeze record for LN2 Part B (confirmation on locked ChIP test factors), 2026-10-10

Development stopped by the registered stopping rule after round 2 (THEORY §11b, ROADMAP log).

**Frozen HASTUR:** `proto/hlat.py` VERSION v8, run as tool key `hlat_a3` = `--sites --gtrim 0.01` with all other defaults
(k 4–5, w 13, L 8, 12 seeds per family, Markov order 1, no centring, min 6 bits).
**Comparators:** STREME 5.5.9 `--dna --nmotifs 10`; HOMER `findMotifs.pl fasta -noknown -len 8,10,12 -p 2` (as in LN1).
**Units:** the ChIP test factors in `splits/chip_split.tsv` with ENCODE conservative-IDR peaks: 32 of 34
(RARA and RFX3 have none in K562; reported missing, not replaced). N = 250, 500, 1,000 (3 subsets each,
`bench/lowN_make.py ... test`, LN1's sampling rule with the factor index taken within the test list) and N = 5,000.
**Readout and criteria:** as registered in `splits/protocol_LN2.md` Part B (STREME parity on located top-10 and P5
top-10 at every N: paired mean ≥ −0.05 with bootstrap lower bound ≥ −0.10; PWM similarity not lower by > 0.02;
runtime ≥ 3× lower than STREME; calibration from the motif-free panel: under the CpG shift HASTUR 0.2 false calls per
run vs STREME 3.8, already measured on the frozen configuration).

**File hashes (sha256):**
```
b3eb6a36cb3e679c0a0c34943cd7fef1ea291573844a2ddf78e9183ad496c49a  proto/hlat.py
dd3fa9bbb17a48961d3d56c2b12121dcd738b7f0d99cd434364ceb3ef0a67092  bench/run_tools.py
6b24a83f18a9bc4d404782f9d148e2aa656b5daaa3d2d4c776cfed282e211736  bench/evaluate.py
319eb29389ef53269969fee40d014cbf248dd1eefcf413fc81aad22fe6dc8964  eval/score.py
d61af86055f40253adcc7e21b546f11658d1b8d2fb7a189ffa015f33ef2c2f6a  bench/lowN_make.py
944f30b8b8f7de536c06b66f87da8f46ba2b15b7a395c44a0d51f47d45300e9b  data/prep_chip.py
13fe7e84c731a0c3b7bd5a1b4b155f62da36c78a3f0f683b7f9542733aee23e3  bench/lowN_report.py
```
