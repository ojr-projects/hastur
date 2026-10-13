# ROADMAP: master status for HASTUR (repo idea_coins)

> **Superseded for forward planning by `ROADMAP_v2.md` (2026-10-12).** This file remains the dated development log
> (§5) up to and including the FR1 confirmation, the masked comparison and the Rust port.


> **Read this first after a restart or compaction.** This file supersedes every other doc's "status" section. Theory
> and design: `theory/THEORY.md`. Process rules: `docs/DEV_INSIGHTS.md`.

## State at 2026-10-10 (evening)

- Day 0 theory is done and pushed (`theory/THEORY.md`, checks c1–c8).
- Eval corrections were handed to the Comsequent agent (`../comsequent/docs/briefs/comsequent_eval_corrections.md`).
- Splits registered (§2). Box set up (§4). Now building: the evaluator (E0), the Python prototype, the dev panel, then
  E1/E2.

## 1. Operating rules

- Commit and push to `main` myself; add paths explicitly; `*.pem` stays ignored.
- Heavy compute runs on the box when one is up; the box is currently lost (§5), so work runs locally. The local machine is busy with another agent until ~2026-10-11 05:00; keep local runs light.
- No real-data CRISPR/PAM work (user request).
- Develop on **dev** units only. **Test** units stay locked until a freeze. Any post hoc rule change makes a phase
  development data: label it.
- Every scorer reports a shuffled-truth chance baseline (THEORY §3), and **every truth source passes an oracle test
  (its own matrix scored as a call) before it is registered**; a readout whose oracle fails yields no verdict.
- **Simulation is for debugging only** (THEORY §11b): the sim panel inflated HASTUR's margin by ~0.3 against real ChIP
  dev. Decisions use real dev units; claims use locked test units only.
- **Register the configuration** (version tag, flags, background multiple, comparator flags) with a hash before any
  test run; one output list per tool; harness self-check that every tool can express every truth.
- Stopping rule: measure replicate noise on the dev panel first; stop iterating when two consecutive variants move the
  headline by less than that noise, or only trade wins between classes.

## 2. Splits (registered 2026-10-10, before any run of this project; hashes in `splits/SHA256SUMS`)

| phase | dev | test (locked) | rule |
|---|---|---|---|
| ChIP (ENCODE, top 5k summits ±100 bp, GC-matched genomic bg) | Comsequent's 24 P5 factors (all already seen) | 34 new factors: every HOMER-library motif from ENCODE-era K562/HepG2 TF ChIP not among the dev 24 (truth = HOMER consensus); factors without ENCODE peaks are reported missing, never replaced | `splits/chip_split.tsv` |
| HT-SELEX | Comsequent dev (186 TFs) | Comsequent test (186 TFs; never scored by either project) | `splits/selex_split.tsv` |
| RBNS | Comsequent dev (56 RBPs) | Comsequent test (55) | `splits/rbns_split.tsv` |
| Methylation | Comsequent's 10 dev strains | Comsequent's 25 held-out strains (its own registered run) | Comsequent `HELDOUT_STRAINS.tsv` |
| Simulation | dev panel (§3) | none: sims are for mechanisms only | — |

## 3. Plan (THEORY §11)

| step | what | status |
|---|---|---|
| E0 | Evaluator (`eval/score.py`): p5 / contains / PWM-similarity rules + shuffled-truth chance; oracle self-test passes. Representation is scored against the TRUE matrix where one exists (a 0/1 consensus penalises correct soft PWMs) | done |
| P0 | Prototype `proto/hlat.py` | v2 running (see log) |
| D0 | Sim dev panel (`bench/sim_make.py`, 11 cells × 5 reps, bg 5×) with STREME, HOMER, Comsequent v1.0 all-in-one | running (5 reps) |
| LN1 | **Registered low-N ChIP panel** (`splits/protocol_lowN.md`, hash in SHA256SUMS): 24 dev factors × N 250/500/1000 × 3 reps; HASTUR plain and summit-centred (v6) vs STREME and HOMER; registered breakpoints decide whether an accuracy claim exists | **done: accuracy thesis not alive; summit weighting no effect** (THEORY §11b LN1) |
| LN2 | Registered (`splits/protocol_LN2.md`): Part A dev variants for low-N representation; Part B frozen confirmation on the 34 locked test factors | Part A round 1 done (A1+A2 kept); round 2 (A3) running; test inputs building |
| E1 | Representation layer on ChIP dev (24) and SELEX dev subset | ChIP dev done at N = 5,000 (parity with STREME, §11b); SELEX running |
| E2 | Double-centred ranking on SELEX dev (high drift) and ChIP dev (CpG) | todo |
| E3 | Full-length seeds at small N | todo |
| E4 | c2 robustness with matrices refitted on all peaks | todo |
| Go | Port the winning components into a Go engine | after E1–E3 |

## 3b. Positioning and demonstrations (adopted 2026-10-11, from an outside review)

**Positioning:** *calibrated, fast motif discovery for high-throughput and automated use*: "the tool you can run
10,000 times and trust the output". It is not a STREME replacement: STREME is stronger on full-length motifs
(THEORY §11b), and a head-to-head on raw accuracy is the wrong fight. Speed matters because it multiplies with
calibration: at 3.8 false calls per run, STREME's output always needs a person; at ~0.1, output can be taken at face
value, and false calls stop accumulating across thousands of runs.

**Decisive check first: FR1** (`splits/protocol_FR1.md`): sensitivity vs false calls per run for every tool, with
thresholds swept, on real positives and real nulls. Better frontier → methods contribution; same frontier →
"good defaults + speed" (weaker). Known so far: STREME's 3.8 false calls per run under a composition shift all have
E ≤ 0.05, so its held-out test calls them significant. That is miscalibration for "is this a motif" (honestly:
STREME's test is correct about enrichment; a composition shift is enrichment, so the claim is that HASTUR separates
motif structure from composition).

**Workflows the calibration enables, and how each will be demonstrated:**

| # | workflow | demonstration (data) | status |
|---|---|---|---|
| D1 | **Discovery on everything**: per regulon / operon set, per genome across a pangenome | E. coli regulon screen: every RegulonDB regulon's upstream regions vs other operons, plus thousands of random operon sets of matched sizes as nulls; true vs false discoveries per tool. Then the same regulons' orthologs across the 21 local bacterial genomes (pangenome) | to build (genomes and GFFs local) |
| D2 | **Empty results that mean something** (screening, QC) | ENCODE ChIP experiments with and without audit flags: does a calibrated "nothing found" track failed experiments? | to build |
| D3 | **Many contrasts without multiplying false calls** | differential discovery across many conditions (e.g. per-cell-type accessible regions vs the rest), false calls counted across all contrasts | to build |
| D4 | **Motif calls as features** | calibrated motif presence/absence matrix across all ENCODE TF ChIP sets (a motif atlas): cluster it, recover cofactor modules; compare STREME's matrix | to build |
| D5 | **Interactive and automated use** | wall-clock per run on an idle machine (LN3 amendment); parameter sweeps (window, background) in seconds; one pipeline with no human checkpoint | timing in FR1 |
| D6 | **Run-level FDR across a screen** (new) | because each run is calibrated, a global FDR over thousands of runs can be computed from per-run significance; show it holds on the D1 nulls | new |
| D7 | **Auditing deep-model motifs** (new) | check TF-MoDISco motifs from a public ChromBPNet/BPNet model against calibrated enrichment in the same peaks: which "motifs" survive? | new |
| D8 | **Replicate reproducibility as an output** (new) | same motif across subsamples / biological replicates; per-motif reproducibility score (low-N subsamples already exist) | in FR1 |
| D9 | **Allele- or variant-centred contrasts** (new) | windows around caQTL/ASB alleles (bound vs unbound allele) as many small matched contrasts | later |
| D10 | **Core-then-refine pipeline for long motifs** | HASTUR finds cores fast; MEME/STREME restricted to the windows around its sites recovers the full frame; report accuracy by motif length | after FR1 |

**Stress tests required before any claim:** real unbound genomic sequence with realistic composition (human and
AT-/GC-rich bacterial genomes), null-only runs many times over (the false-call distribution), replicate
reproducibility, and accuracy reported by motif length (the long-motif gap shown, not averaged away).

**Order of work:** (1) FR1 on dev (needs STREME/HOMER reruns with permissive thresholds and the real nulls); (2) LN5
matrix accuracy (the step-0 anatomy showed the gap is column accuracy, not width; registered arms: extra site-refit
round, plain site threshold); (3) D1 regulon
screen as the headline demonstration; (4) FR1 confirmation on test2 / SELEX test with frozen settings; (5) D2–D4, D6,
D7 as time allows.

## 4. Compute

- Box: c7a.16xlarge spot (64 cores, 123 GB), host in memory/`BOX` env (changes on restart), setup
  `infra/box/box_setup.sh`, repo synced to `~/idea_coins` (rsync, no `.git`, no `*.pem`).
- Box jobs run in tmux sessions named after the job; results are pulled back into `results/` (small tables only).
- From 2026-10-11 the box has a persistent 500 GB EBS volume at `/data` that
  survives spot termination; envs, repo copy, data, runs and queue state live there. See `infra/box/BOX.md`.

## 5. Log

- 2026-10-11 morning: **New spot box with a persistent EBS volume at `/data`** (`infra/box/BOX.md`). FR1 moved
  entirely to the box: HASTUR, permissive STREME and HOMER on every null and dev positive, one machine, one code
  version (`infra/box/run_fr1_box.sh`).
- 2026-10-12: **Masked comparison (secondary; `splits/secondary_masked.md`; window identity verified, 0 of 2.37M
  mismatching; 45% of bases masked).**
  - **Registered-style cell:** F2 − STREME @0.1 = +0.266 [+0.224, +0.310], about the same as unmasked (+0.241).
  - **Exchangeable nulls only:** +0.033 [−0.006, +0.073], parity. ChIP −0.019 (unmasked −0.059), SELEX +0.101.
  - **GC-selected human nulls stay miscalibrated for both tools** after masking (F2 4.0–8.0, STREME 6.5–8.0 per run).
    Removing the Alus exposes other GC-rich structure (CpG islands) that an order-1 Markov model does not capture.
  - **Conclusion:** parity with STREME on exchangeable sequence, now including ChIP once repeats are masked. HASTUR
    keeps a large frontier advantage whenever composition-shifted backgrounds are in the null.
  - Files: `results/frontier_masked_f2*.txt`.
- 2026-10-12: **Full Rust HASTUR (end to end), outputs identical to Python F2** (`rust/hastur_core`; validated
  byte-identical TSV and MEME on 8 units: regulon, ChIP 250 / 500 / 1,000 / 5,000, SELEX, human and bacterial nulls).
  - **Ported for identity:** numpy's PCG64 and SeedSequence, `choice`, `permutation`, quantile `_lerp`, `linspace`
    and pairwise summation.
  - **Full run, 12-core workstation, 24 threads:** 0.8–1.8 s (≤ 1,000 peaks, regulons) and 8.2 s (5,000 + 50k bg).
    Python F2 takes 9–40 s; STREME takes 124–2,150 s single-threaded.
  - **Single thread:** 5.5–70 s. Representation dominates (about 70%: held-out z scores all four sets, including two
    synthetic copies).
  - **Next speed steps change the numerics** (prefix-sum background, smaller synthetic sets for the frame/EM
    choice), so they would be validated on accuracy, not identity.
- 2026-10-12: **Speed floor: Rust detection core** (`rust/hastur_core`; `results/speed_rust_vs_python.txt`).
  - **Validated identical to Python's `find_seeds`** (24/24 seeds, z equal to 3 d.p., 4 units). A rigorous pruned
    exact maximum gives identical permutation maxima with far fewer tail evaluations.
  - **Detection with 20 permutations, 12-core workstation, 24 threads:** 0.37 s (20 regions), 0.41 s (CRP regulon),
    0.37 s (250 peaks), 0.63 s (1,000), 2.15 s (5,000 + 50k bg). Python F2 takes 6–18 s with 24 threads; STREME takes
    75–2,150 s single-threaded.
  - **Single thread:** 1.8–21 s, about 3–5× numba thread for thread.
  - **Floor:** counting runs at about 2.3e9 increments/s per core, near the random-increment limit. Tiling might
    give about 2×.
  - **Not ported:** representation (EM, site refit; 2–19 s in Python).
- 2026-10-12 ~01:00: **F2's false calls on GC-selected human nulls are Alu elements** (`CTACTAAAAATACAAAAAATTAGCCGG`,
  `TTTGTATTTTTAGTAGAGA`, `CCCAGCTACTCGGGAGGC`, all Alu consensus segments).
  - **Cause:** Alus concentrate in GC-rich isochores, so GC-selected windows really are Alu-enriched; any motif
    finder reports them (STREME's 5–9 per run very likely the same). This is not a composition-model failure, so a
    higher-order Markov chain will not fix it.
  - **Remedy: repeat masking** (hg38 is soft-masked; HOMER masks by default). All panels were upper-cased at build
    time.
  - **Next step if resumed:** build masked human panels and add a lowercase-as-N option; re-read the GC-selected
    calibration and ChIP sensitivity with both tools masked. The registered FR1 advantage partly reflects how each
    tool handles unmasked repeats in GC-shifted backgrounds.
- 2026-10-12 ~00:50: **Decomposition of the FR1 confirmation (post hoc, secondary; labelled as such).**
  - **On exchangeable real nulls only** (human random windows + bacterial coding; GC-selected human excluded):
    **F2 − STREME @0.1 = +0.010 [−0.031, +0.052], parity** (`results/frontier_conf_f2_exchangeable_posthoc.txt`).
    ChIP −0.059 (STREME better), SELEX +0.101 (F2 better).
  - **So the registered +0.241 is driven by the GC-selected nulls,** where STREME is badly miscalibrated (5–9 false
    calls per run) and must be thresholded very strictly. F2 is miscalibrated there too (2.6–3.6), but less.
  - **Accurate claim for F2:**
    - calibration at STREME's level on exchangeable real sequence, and much more robust to composition shift;
    - sensitivity at parity with STREME on exchangeable data (better on SELEX, worse on ChIP);
    - clearly better when composition-shifted backgrounds are part of the null;
    - about 3–4× faster.
  - **Also reported (FR1):**
    - replicate reproducibility on test2 at parity (significant calls reproduced in both other subsamples:
      F2 0.52 / 0.60 / 0.67 vs STREME 0.53 / 0.59 / 0.55; `results/repro_conf_f2.txt`);
    - by motif length: long truths 0.31 vs 0.09 at @0.1 on the registered nulls.
  - **LN7 cycle f (union of rankings):** no gain at N ≤ 250, +0.017 at 500, and more calls. Dropped
    (`results/ln7f_union_report.txt`). At N = 25 the limit is the de novo search's multiplicity; only a small family
    (known-motif testing) would help.
- 2026-10-12 ~00:15: **Small-N calibration confirmed for frozen F2 on a verified-fresh panel** (tuple seed,
  offset 5000; fg overlap with earlier panels about 2%, chance level; `results/smallN_conf5_f2.txt`).
  - 0.052 false calls per run pooled; by set size ≤ 15 0.050, 16–30 0.050, 31–60 0.167, 61–200 0.000. Passes the
    LN6 gates.
  - **Runs with any false call: 7 of 210 (0.033)**, below the family-wise target of 0.05. The 31–60 bin is one
    S. aureus run with five correlated calls on a palindromic repeat family. Max-T controls *any* false call per
    run; when a run fails, a repeat family yields several calls.
- 2026-10-11 ~23:55: **Seed-collision flaw found; the small-N "confirmation" above is withdrawn.**
  - **The flaw:** `bench/upnull_make.py` seeded with `… + rep + offset`, so consecutive offsets share 9 of 10
    replicates. The "fresh" panel (offset 3) was 90% the dev3 panel, with identical false calls at shifted replicate
    indices.
  - **Same flaw in `bench/null_make.py`:** `seed0 + 1000·rep` makes offsets 10,000 apart share 10 of 20 replicates.
    The FR1 confirmation's simulated nulls (910000) half-overlap the LN6/LN7 dev panel null_v3 (920000).
  - **Unaffected:**
    - `realnull_make` (offsets 5,000 apart against 20 replicates), so the FR1 decision cell (real-genome nulls) is
      clean;
    - the test positives;
    - D1 (hash seeds).
  - **The FR1 verdict stands. Caveats:** the confirmation's simulated-null rows are half dev data. The LN6 rounds'
    "fresh" small-N decision panels overlapped earlier ones (rounds 2 and 3); their conclusions (no arm kept) do
    not change.
  - **Fix:** tuple seeds (`default_rng([…, offset, rep, …])`), plus a freshness check (fg sequences shared with
    earlier panels). A genuinely fresh small-N panel is running.
- 2026-10-11 ~23:40: **Small-N calibration confirmed for frozen F2** on a never-seen panel (seed offset 3; 210
  runs; `results/smallN_conf_f2.txt`).
  - 0.038 false calls per run pooled (v8 1.72), passing the LN6 gates.
  - By set size: ≤ 15 0.017, 16–30 0.017, 31–60 **0.167**, 61–200 0.017.
  - The 31–60 bin sits at about 0.17 in every panel, so it is systematic. Next: diagnose it.
- 2026-10-11 ~23:15: **FR1 CONFIRMATION (registered; frozen F2 = freeze v18; fresh nulls; unseen test positives;
  `results/frontier_conf_f2.txt`): BETTER FRONTIER.**
  - **Decision cell** (real nulls, intergenic excluded per amendment 3, all 478 positives, @0.1): **F2 − STREME =
    +0.241 [+0.201, +0.280].** With intergenic included: +0.211 [+0.172, +0.251]. v8: −0.031 (operating point only).
  - **By set:** ChIP test2 +0.237 (0.30 vs 0.06), SELEX test +0.245 (0.68 vs 0.44). Human nulls +0.215, bacterial
    +0.010 (parity), simulated +0.525. Long truths: 0.31 vs 0.09.
  - **Native false calls per run, fresh nulls:**
    - exchangeable real: F2 0.00–0.15, STREME 0.00, v8 0.5–1.8;
    - CpG shift: F2 0.05, STREME 3.45;
    - GC-selected human: F2 2.65–3.60, STREME 5.2–9.25.
  - **Caveats:**
    - STREME is the only comparator (HOMER dropped).
    - On exchangeable real nulls, STREME's native operating point is slightly more conservative than F2's (0.00 vs
      about 0.06 per run).
    - Both tools are miscalibrated on GC-selected sets.
    - F2's dev numbers were optimistic (dev +0.190; confirmation +0.241 on different data).
- 2026-10-11 ~22:00: **LN7 cycle e (dev): genomic control on the tested statistics was the SELEX ceiling.**
  - **Cause:** selection shifts the whole word spectrum between SELEX cycles, so the median-based λ was large and
    divided true motifs from z ≈ 37 to z ≈ 4. Under permutation max-T, λ is redundant for exchangeable nulls; the
    thresholds stay near Bonferroni in SELEX.
  - **P7 (no GC on tests) and F2 (P7 + L2) pass all gates.** Simulated CpG shift is 0.05, at the bound.
  - **FR1 dev decision cell:** P5 +0.130, P7 +0.167, **F2 +0.190 [+0.142, +0.237]**. SELEX goes from −0.12 to
    +0.19, bacterial nulls from −0.15 to +0.02, human nulls +0.23.
  - **Native-threshold sensitivity of F2 equals v8's at N = 250–5,000:** calibration no longer costs sensitivity
    there.
  - **Costs:** GC-selected human sets (not gated) worsen, 2.8–3.55 false calls per run (STREME 4.9–9.2). D1:
    0.079 false calls per run, 10/21 regulons.
- 2026-10-11 ~22:15: **Freeze v18 = F2** (`splits/freeze_v18.md`; hash checked on the box) **and FR1 amendment 3**
  (disclosure; intergenic nulls out of the decision cell; HOMER dropped), both committed before the confirmation
  outputs were read. **FR1 confirmation launched** (`infra/box/run_fr1conf_hastur.sh hlat_f2`).
- 2026-10-11 ~20:30: **LN7 cycle d (dev): background-only permutations (P6) pass all gates but do not help.**
  - FR1 cell +0.092 vs P5's +0.130. SELEX is unchanged (0.34 vs 0.35), so the pool-contamination diagnosis of the
    SELEX ceiling was wrong; the cause is open.
  - D1 is slightly better (11/21 regulons, 0.062 false calls per run).
  - **P5 is chosen.** Freeze sequence (FR1 amendment 2): P5 + LN5's L2 (`--plainz`) as tool `hlat_f1`, a dev sanity
    check of the combination, then `splits/freeze_v17.md` and FR1 amendment 3, then the confirmation
    (`infra/box/run_fr1conf_hastur.sh`).
- 2026-10-11 ~19:40: **LN7 cycle c (dev): P5 clears the FR1 bar on the dev panels** (`results/frontier_ln7c.txt`;
  HASTUR-only rerun against the stored STREME outputs; calls scored by the family-wise adjusted p).
  - **FR1 decision cell** (real nulls, all positives, @0.1): **P5 − STREME = +0.130 [+0.067, +0.190]** (v8:
    −0.314). Registered bar: ≥ +0.05 with the CI above 0.
  - **By positives and nulls:**
    - ChIP +0.41, SELEX −0.12;
    - human nulls +0.15, bacterial nulls −0.15, simulated +0.36;
    - long truths 0.42 vs STREME 0.21.
  - **Weakness:** SELEX sensitivity tops out around 0.41 even at loose thresholds (v8 0.59). Diagnosis: the
    permutation pool mixes fg into the null, about 9% of the signal per permuted set, which inflates thresholds in
    strong-signal data.
  - **D1 rerun with P5** (`results/d1_report_ln7c.txt`): 0.076 false calls per random-promoter run (v8 1.50,
    STREME 7.65), 10/21 regulons at the native threshold (STREME 15/21), screen FDR 0.15 vs 0.996. @0.1 false calls
    per run: 0.57 vs 0.48. Runtime 7,446 s vs STREME 27,529 s.
  - **Dev only.** Confirmation = the FR1 confirmation sets (unread STREME runs exist) with a frozen configuration.
    Cycle d tests background-only permutations (`--perm-bg`) first.
- 2026-10-11 ~18:45: **LN7 cycle b (dev; `results/ln7b_report.txt`): P5 (rank by v8, test exact, permutation
  max-T) passes all four gates.**
  - **Calibration:** small-N 0.038 per run, real-genome 0.03, simulated 0.00 (CpG shift included). The ChIP motif
    list is identical to v8's.
  - **Sensitivity at the native threshold** (significant calls matching the truth, P5 rule):
    - v8 / P5: 0.32 / 0.04 (N = 25), 0.53 / 0.25 (50), 0.49 / 0.44 (100), 0.75 / 0.61 (250), 0.83 / 0.75 (500),
      0.79 / 0.77 (1,000), 0.84 / 0.84 (5,000).
    - Calibration is free from N = 1,000. It costs 0.08–0.14 at 250–500, and N = 25 calls almost nothing.
    - P2 (tests its own ranking) is better at small and mid N (0.07 / 0.33 / 0.49 / 0.67).
  - **Next:** a family-wise adjusted p per call (the Gumbel tail; one calibrated scale across runs) for the frontier.
    Cycle c reads the FR1 dev frontier and D1, HASTUR-only against the existing STREME outputs.
- 2026-10-11 ~18:00: **LN7 cycle a (dev; `results/ln7a_report.txt`).**
  - **P2 (R3 statistic + permutation max-T) passes all three calibration gates for the first time:**
    - small-N nulls 0.048 per run (v8 1.72);
    - real-genome nulls 0.03 (v8 0.68);
    - simulated nulls 0.00, CpG shift included;
    - GC-selected human sets 0.50 at N = 500 (v8 3.0).
  - **P2 fails ChIP gate 4** (located −0.05 to −0.07 at 250–500, best similarity about −0.02). R3 lost the same
    without permutation, so the cause is the exact ball test changing which L-mer seeds are chosen.
  - **P3 (+ presence counts)** fails CpG shift (5.15), as predicted.
  - **Speed (outputs identical):** about 2.3× from closed-form order-1 expectations, batched exact tests, a subset
    cache, numba revcomp and a max-only permutation pass. Dev roots now run concurrently (`infra/box/run_arms.sh`).
  - **Cycle b:** P5 = rank by v8's statistics, test the exact statistics against the permutation thresholds
    (`--rank-v8`). Max-T thresholds bound every cell's tested statistic, so FWER holds whatever chose the cell. The
    motif list is v8's (verified), so gate 4 is near-automatic. Sensitivity at the native threshold is read with
    `bench/sig_hits.py`, and the decisive check will be a HASTUR-only FR1 frontier rerun.
- 2026-10-11 evening: **Resumed (user): HASTUR-only exploratory dev, box only, short cycles (LN7).**
  - **New ground: calibrate the maximum empirically, not each cell analytically.** Westfall–Young max-T: about 20
    label permutations of the pooled sequences, each family's maximum z recorded, a Gumbel fitted by moments, and
    the family-wise 2.5% point used as the threshold instead of Bonferroni (`--perm R`, v15). Exact under
    exchangeability whatever the clumping, repeats or N. Cheap because counts are additive: only the permuted fg is
    counted.
  - **Composition shift still has to come from the statistic.** Permutation does not cover non-exchangeable nulls.
  - **Sanity runs:** the known false calls vanish (C. jejuni z = 44; CTCF's low-complexity call).
    - v8's soft L-mer table has an absurd null tail on bacterial and small sets (thresholds 55–86). R3's exact
      statistics stay near Bonferroni there.
    - Human sets need thresholds of 6.3–7.9; repeats are the residual. CTCF at N = 50 is lost (gapped threshold 9.0).
  - **Cycle a:** P2 (R3 statistic + perm) and P3 (+ presence counts) on the four gates, using the round-3 panels as
    dev. A pass gets a fresh-panel confirmation.
- 2026-10-11 ~17:00: **LN6 round 3 (the last round) fails → development PAUSED** (user rule;
  `results/ln6r3_report.txt`). R4 = R2 + presence expectations from synthetic Markov replicates, on fresh panels.
  - **Gate 1, small-N nulls: pass,** 0.014 per run (v8 1.72).
  - **Gate 2, real-genome nulls: fail,** 0.35 per run (v8 0.68). human_exch N = 2,000 is *worse* than v8: 1.10 vs
    0.70.
  - **Gate 3, simulated nulls: fail,** 0.10–0.25 per run vs v8 0.00. The CpG-shift fix itself worked (4.45 → 0.10).
  - **Gate 4, ChIP: fail,** located −0.158 at N = 5,000 (n = 19).
  - **Diagnosis:** the synthetic-replicate expectation carries its own counting noise, and the exact test treats it
    as fixed, which is anti-conservative at moderate counts. A fix would propagate that uncertainty (e.g. an exact
    test on real vs synthetic counts with both sets as data). That would be a fourth round, so not run.

**Status at pause (2026-10-11):**
- **Accuracy:** HASTUR v8 locates motifs about as well as STREME. It states them slightly worse; LN5's L2 narrows
  the gap.
- **Speed:** 4–14× faster than STREME today, with 10–100× more available from a native rewrite.
- **Calibration:** **not achieved on real sequence.** v8 makes 0.45–1.05 false calls per run on exchangeable
  real-genome nulls, against STREME's 0.00. Every repair so far fixes one regime (small N, repeats, the bacterial
  L-mer table, composition shift) and breaks another.
- **Most promising next step, if resumed:** a held-out significance test of the few preselected candidates
  (STREME's design), which needs tail accuracy around 0.05 instead of 1e-8. HASTUR detection would then be used only
  for speed. Unread FR1 confirmation runs and the built SELEX test split remain available.
- 2026-10-11 ~16:15: **FR1 dev, final (frozen v8 vs permissive STREME, no HOMER; `results/frontier_dev.txt`):**
  - **Decision cell:** real nulls, all positives (190 ChIP + 212 SELEX), @0.1: **−0.314 [−0.359, −0.269]** → report
    as is. Same as the preliminary.
  - **Missing units:** two ChIP 5k STREME runs were killed twice (HNF1A and one other) and are missing, not misses.
  - **Scope:** this is v8's verdict. The LN6 configuration, if one is kept, gets its own HASTUR-only rerun on the
    same panels.
  - **Also done:** the FR1 confirmation runs (v8 + STREME on test2, SELEX test and fresh nulls; not read, per
    amendment 1) and the SELEX test-split build (212 libraries).
- 2026-10-11 ~15:00: **LN6 round 2 result (registered keep rule; `results/ln6r2_report.txt`): no arm kept.**
  - **R2 (presence + ball) passes both real-data null gates:** 0.05 per run on the real-genome nulls, 0.019 on the
    small-N panel.
  - **It fails the simulated CpG-shift null: 4.45 false calls per run vs v8's 0.05,** about STREME's level (R1 the
    same; R3, the ball test alone, 0.00). The likely cause: presence expectations from 1 − exp(−λ) assume Poisson
    occurrences, while under a CpG-rich chain CG words clump. Expected presence then differs between fg and bg, and
    the double-centring no longer cancels. This loses HASTUR's composition-shift robustness, its clearest advantage
    over STREME.
  - **It also fails located at N = 5,000:** −0.105 [−0.263, 0.000], n = 19. At N = 100, P5 rises +0.25; at
    250–1,000 it is within noise.
  - **Reading:** calibration on real sequence moved a lot, but the passing configuration trades away composition
    robustness. User's rule: pause unless calibration moves. It moved, but no arm is keepable as registered. Decision
    with the user.
- 2026-10-11 ~14:40: **LN6 round 2, null gates (decision panels): R2 passes both. The ChIP sensitivity guard is
  still running.**
  - **Exchangeable real-genome nulls** (100 runs), native false calls per run: v8 0.71, R1 (presence) 0.31, R3 (ball)
    0.25, **R2 (presence + ball) 0.05**, STREME 0.00.
    - Presence counting fixes the human nulls (repeats); the ball test fixes the bacterial ones (soft L-mer table).
    - Under GC composition shift, R2 makes 2.0–2.4 calls per run against STREME's 4.9–9.2. Both are miscalibrated
      there.
  - **Fresh small-N panel** (210 runs, 10–200 regions): v8 1.76, R1 0.22, R3 0.11, **R2 0.019**, with at most 0.033
    in any size bin. v8 makes 3.0 at ≤ 15 regions.
  - Tables in `results/ln6r2_realnull.tsv` and `results/ln6r2_upnull.tsv`.
- 2026-10-11 (user decision): **LN6 round 2 is the stopping point for calibration.**
  - **Pause:** if it shows no movement on calibration, development pauses. That means native false calls not
    ≤ 0.1 per run on the exchangeable real-genome nulls (amendment 1) and on the fresh small-N panel, with ChIP
    sensitivity within 0.03 of v8.
  - **Continue:** if it passes, rerun the FR1 frontier with HASTUR only for the new configuration, then decide.
  - **Learnings** recorded in `docs/DEV_INSIGHTS.md` (§3.7 and §8, building on Comsequent's).
- 2026-10-11 ~13:30: **FR1 preliminary (HASTUR v8 vs permissive STREME, no HOMER and no ChIP 5k;
  `results/frontier_prelim.txt`): fails.**
  - **Decision cell:** real nulls, all positives, @0.1 false calls per run: **−0.31 [−0.36, −0.26]**.
  - **By null family:**
    - simulated: +0.50 (STREME badly miscalibrated under composition shift);
    - human: −0.02 (ChIP +0.22, SELEX −0.21);
    - bacterial: **−0.72**.
  - **Native false calls per run on exchangeable real-genome nulls: HASTUR 0.45–1.05, STREME 0.00.** So v8's
    calibration claim held only on simulated Markov nulls. The cause is the one LN6 found: real k-mer counts are
    clumpy, so occurrence-count Poisson tests are anti-conservative for the clumpiest cells, and median-based
    genomic control cannot fix a tail. STREME counts sequences.
  - **Null-construction flaw:** the bacterial intergenic nulls sampled overlapping windows from small intergenic
    spaces (376 of 500 C. jejuni fg windows share a 20-mer). Rebuild by tiling. The coding and human nulls are
    unaffected and fail too.
  - **Action (user: save HOMER compute if this fails):** HOMER dropped from FR1 dev and confirmation. FR1 dev
    completes with HASTUR and STREME only. LN6 round 2 is now also gated on the exchangeable real-genome nulls
    (amendment 1).
- 2026-10-11: **Box incident (no data lost).** The LN5 job's audit, running on `lowN` at the same time as FR1,
  deleted an in-progress FR1 HOMER output dir (it had no `runtime_s` yet). That FR1 run then crashed on writing its
  marker, and the queue restarted FR1 at 09:55, resuming from cached outputs. Fix: the audit now leaves dirs without
  `runtime_s` alone. `run_tools` reruns them anyway, and they may belong to another job. HOMER on SELEX units takes
  up to ~70 min each, so FR1 dev runs into the evening.
- 2026-10-11: **LN6 round 2 registered** (`splits/protocol_LN6_round2.md`): per-sequence presence counts
  (`--seqcount`) and an exact radius-1 Hamming-ball test for L-mer cells (`--lmer-ball`). It is decided on a
  fresh-seed small-N null panel, because round 1's panel was used for diagnosis. A quasi-Poisson soft test was
  rejected before registration: its null dispersion makes one occurrence count as about 2.9, which is
  anti-conservative in the tail.
- 2026-10-11: **LN6 round 1 result (registered keep rule; `results/ln6_report.txt`): no arm kept.**
  - **Small-N null panel** (three non-E. coli genomes, 210 runs), false calls per run pooled: **v8 1.75**,
    C1 0.31, C2 1.54, **C12 0.086**. C12 passes the pooled bound but fails the 31–60 bin (0.27 > 0.2).
  - **ChIP:** C2 and C12 also fail on located −0.053 at N = 250. C2's expected-count > 5 cutoff switches the soft
    L-mer table off below roughly 1,000 sequences, and with it the motifs that table found at N = 100–250.
  - **Simulated nulls:** C1 and C12 lower false calls (0.067 → 0.033 per run).
  - **Residual C12 false calls** are marginal (z 5.2–6.2), low-complexity patterns, plus one tandem repeat
    (`TGSCAYSR`×4) giving four calls in one run. That is clumping: occurrences from one locus treated as
    independent evidence.
  - **Round 2** (to register): count gapped cells per sequence (presence, as STREME and HOMER do) with the exact
    test; replace C2's cutoff with a quasi-Poisson exact test on the soft L-mer counts.
- 2026-10-11: **LN5 result (registered keep rule; `results/ln5_report.txt`): L2 kept** (`--plainz`, the site
  threshold from plain fg-vs-bg enrichment).
  - **L2:** ChIP best PWM similarity +0.008 [+0.005, +0.010] pooled; P5 and located do not fall at any N; SELEX
    flat. Significant calls are identical to v8 in every null run.
  - **L1 fails** (P5 −0.035 at N = 250, SELEX P5 −0.033).
  - **L12 fails** on located −0.035 at N = 250, the consensus readout the D1 oracle work showed to be noisy. Its
    gain is the largest, +0.014 [+0.011, +0.017], which leaves it −0.002 to −0.006 behind permissive STREME on best
    similarity. It is on record as the round-2 candidate.
  - **Two evaluator bugs found on the way:** `bench/evaluate.py` and `bench/null_eval.py` used hard-coded tool
    lists, so new arms were silently skipped, and `null_eval` would have crashed on `streme_all`. Both now discover
    tools from the unit directories. Earlier LN reports listed their arms explicitly and are unaffected.
- 2026-10-11: **D1 result (E. coli regulon screen, registered readout; `results/d1_report.txt`).**
  - **Native points:** HASTUR (frozen v8) and STREME both find 15 of 21 evaluable regulons (0.71). False calls per
    random-promoter null run are **1.50 for HASTUR** vs 7.65 for STREME. Estimated screen FDR: 0.54 vs 0.996.
  - **Frontier:**
    - at 0.1 false calls per run: **HASTUR 0.00** vs STREME 0.48 (HOMER 0.52), because HASTUR's null scores reach
      −log10 p ≈ 970;
    - at 0.5: HASTUR 0.57 vs 0.48.
  - **Cost:** HASTUR 5,592 s for all 374 runs, against STREME 27,529 s (default) / 64,180 s (permissive) and HOMER
    41,090 s.
  - **Verdict under the registered rules:** "calibrated screen" fails (1.5 > 0.1), "better frontier" fails, and
    "matches STREME at the native point" passes (equal sensitivity, 5× fewer false calls).
  - **Diagnosis:** HASTUR is not calibrated at small N. Gapped-lattice false calls scale with smallness: 2.29 per
    run at ≤ 15 promoters, 1.27 at 16–30, 0.88 at 31–60, 0.23 at 61–250. That is the normal approximation for
    cells with few fg counts against a large bg; LN4's exact test applies only when F + B ≤ 100. The soft L-mer
    table adds a size-independent 0.2–0.3 per run, including every z > 20 call (up to z ≈ 67). Small sets were
    never in the calibration panels (minimum N = 250 in the ChIP work, 500 in the nulls).
- 2026-10-11: **D1 registered** (`splits/protocol_D1.md`). E. coli regulon screen from RegulonDB v12: 34 TFs and
  340 random promoter sets of matched size.
  - **Oracle test:** the IUPAC and PWM-similarity readouts fail it; site recovery passes, and is evaluable for 21 TFs.
  - **Amendment 1 (pre-run):** STREME reports `S=` instead of `E=` when it has no hold-out set. The parser read only
    `E=` and would have scored those runs as having no calls. No earlier panel was affected.
- 2026-10-11: **LN5 step 0 (anatomy, low-N dev): the long-motif/width hypothesis is refuted.**
  - **Width:** HASTUR covers the truth's informative columns as well as STREME, with wider calls.
  - **Where the gap is:** the whole PWM-similarity gap is in the covered columns.
  - **Not the cause:** softening toward uniform or toward the peaks' composition, ±1 phase blur, or site count
    (93% of STREME's).
  - **Contamination** by an earlier related motif (JUND TRE/CRE hybrid) explains 29% of the gap.
  - **What helps:** one more hard-EM site-refit round on HASTUR's output closes over half of the gap: 0.877 → 0.893
    at N = 250 and 0.885 → 0.897 at N = 500, against STREME's 0.895 and 0.908. Further rounds are flat.
  - **Next:** LN5 registered with arms L1 (second refit round), L2 (plain fg/bg site threshold) and L12. Selection
    is untouched by construction (verified).
  - **Erasure (local dev check, not a registered arm):** refitting later motifs with the earlier motifs' sites
    masked, STREME-style, changes best similarity by −0.0003 on top of L1 (30 of 171 units move, in both
    directions). Dropped. The JUND hybrid comes from weak half-site variants, not from rank-1 sites. The same check
    puts L1 at +0.006 over v8 (low-N pooled).

- 2026-10-11 early: **LN2 Part B (locked ChIP test) readout failed its oracle test** (truth strings over-specified:
  own matrix passes `contains` 6%, P5 30%) → registered verdict uninformative. Post hoc against the true HOMER
  matrices: HASTUR behind STREME on P5 by 0.09–0.19 (chance-corrected 0.11–0.17) and PWM similarity by 0.02–0.03,
  located −0.08 to +0.03 (n.s.); 4–14× faster. Dev parity did not transfer (dev truths were short cores). THEORY §11b.
  Rule added: **every new truth source gets the oracle test before registration.**

- 2026-10-10 late: **LN4 calibration (dev).** The exact test (E1) cut Markov-null false calls 0.2 → 0.1 per run;
  CpG shift stayed at 0.2. Remaining false calls were ordinary cells (z 5.4–6.2), not low-count. Diagnosis: the
  null generator started every sequence from a uniform base (non-stationary start). Regenerated with a stationary
  burn-in, the CpG-shift nulls give **0.0 false calls per run** (v8 and E1 alike). The residual was a simulator
  artefact. `null_make.py` now burns in; the fresh-seed LN4 confirmation panel uses the fixed generator.
  Calibration readout to date: HASTUR 0.0 under composition shift (stationary nulls) vs STREME 3.8, HOMER 29,
  Comsequent 53 (earlier panel; comparators to be rerun on the fixed panel for a like-for-like table).

- 2026-10-10 ~20:15 UTC: **the EC2 box became unreachable** (connection closed by remote host, then port 22 timeouts;
  no spot-notice line before the last util tick at 20:05). Lost with it unless the volume survives a *stop*: the built
  ChIP/SELEX inputs and the in-progress LN2 Part B and LN3 Part A runs. **Neither Part B nor LN3 had produced any
  evaluated result**, so rerunning them from scratch is clean; partially pulled Part B tool outputs in
  `results/box/runs/lowN_test` are left unread. All earlier result tables were pulled (last pull 20:08).
  Continuing locally (user freed the machine): local `bench` env, ChIP inputs rebuilt with the same seeded builder,
  heavy jobs capped at ~6 parallel (15 GB RAM). Disk: vhdx 618.8 GB, 163 GB used inside, C: 130 GB free.

- 2026-10-10: theory + checks (c1–c8) committed; splits registered; box set up.
- 2026-10-10 night: **LN2 round 2** (A3, v8: selection on the pre-refit matrix): KEEP; null calls identical to plain
  (0.0 / 0.2 / 0.2); located vs plain +0.10 / +0.01 / 0.00 / +0.17. Stopping rule → **frozen** (`splits/freeze_LN2B.md`,
  hashes verified on the box); **Part B launched** on 32 of 34 locked ChIP test factors (RARA, RFX3: no ENCODE peaks).
- 2026-10-10 night: **SELEX dev** (213 libraries; HOMER partial): HASTUR plain (v6) located top-10 0.25 / top-1 0.23 /
  P5 0.72 vs STREME 0.24 / 0.13 / 0.71, HOMER (109) 0.30 / 0.12 / 0.70, Comsequent 0.07 / 0.05 / 0.23; 16 s vs 168 s.
  **The ChIP-developed site matrix hurts SELEX** (frozen v8: 0.22 / 0.19 / 0.64): partly de-duplication (fixed in
  v8), partly the site refit itself (short reads, palindromic dimers). Lesson: A1 should be a held-out-chosen
  representation candidate (THEORY §6.2), not a blanket refit; HASTUR also outputs few motifs (median 2 vs STREME 10).
  Next development round (to register): sites as a third candidate in the held-out choice; more distinct calls.
- 2026-10-10 night: **LN2 round 1** (dev, vs v6 plain, located top-10): A1 site-based matrix +0.125 at 250 peaks but
  −0.056 at 500 (rejected); A2 trimming p 0.01 no effect (rejected); **A1+A2 kept**: +0.097 / −0.028 / 0.000 / +0.167
  at 250 / 500 / 1,000 / 5,000 peaks; vs STREME −0.04 [−0.14, +0.06], +0.06, +0.06, +0.21. Sim panel top-1 0.78 → 0.86.
  One extra null call (Markov 0.2 → 0.3/run) traced to de-duplication on the refitted matrix → round 2 (A3, v8):
  selection decisions on the pre-refit matrix.
- 2026-10-10 night: **LN1 verdict** (registered rules): HASTUR − STREME located top-10 = −0.14 [−0.25, −0.04] at 250
  peaks, +0.08 / +0.06 (n.s.) at 500 / 1,000; ahead of HOMER from 500; 3–6× faster than STREME; summit weighting no
  effect (sites sit in the sparse regime at the conjunction optimum: Prop. 5 corrected).
- 2026-10-10 evening: focus moves to limited data (where any accuracy gain must show, and runs are cheap). LN1
  registered before any run; HASTUR v6 adds summit-relative weights. Real-data genomic-control λ is 2.8–6.5 on a ChIP
  unit (sims ~1.0): real fg/bg contrasts carry structure an order-1 chain does not absorb; recorded as a calibration
  readout, not tuned mid-protocol.
- 2026-10-10 evening: Comsequent's final reflections reviewed (THEORY §11b). First real data, ChIP dev (24 factors,
  10× bg): HASTUR v3 P5-rule top-10 0.79 vs STREME 0.83 / HOMER 0.74; contains top-1 0.46 = STREME; best PWM
  similarity 0.887 vs 0.861 / 0.825; 84 s vs ~1,000 s. R1 (methylation) downgraded, R6 (gapped/bipartite) withdrawn.
- 2026-10-10 evening, prototype iterations on the sim dev panel (3 reps; top-1 = PWM similarity ≥ 0.85 to the true
  matrix at rank 1):
  - v0 (lattice seeds, EM, cross-fitted plain contrast): top-1 0.64 vs STREME 0.58, Comsequent 0.48, HOMER 0.33; failed
    rank 1 under CpG drift (composition PWMs) and lost the crisp gapped motif (EM smears near-palindromes).
  - Motif-level double centring needs synthetic sequences from each set's own Markov-2 chain (hit-rate contrasts
    against a set-specific background are not composition-invariant).
  - Full cross-fit (seeds on one half) was calibrated but lost too much power at small N; adopted instead: detect and
    rank by the double-centred lattice z on all data (exact family Bonferroni, genomic-control λ), held-out halves only
    choose the representation (frame vs EM). λ ≈ 1.0–1.1 on the panel.
  - Fragment seeds of a kept motif are absorbed (explaining away); EM may refine, not abandon, its seed (drift guard).
  - Working name adopted (user, 2026-10-10): **HASTUR**, Hoeffding-Aligned Sparse-Tensor Unified Rank-one.
  - v1: top-1 **0.70** (STREME 0.58, Comsequent 0.48, HOMER 0.33), median 23 s; CpG cell 1.0 at rank 1; crisp gapped =
    Comsequent. Weak: long motifs (REST, L12 low) → v2 adds frame widening (MEME-style width search on frame counts).
  - v2 (5 reps, 55 datasets): top-1 0.78. GATC missed: the lattice enumerated only k = 5 shapes, so a 4-mer was
    diluted 4-fold → v3 enumerates every order 4..k; double-centring order m = 1 by default (order 2 partly absorbs
    short motifs into their trinucleotide margins). Null pair: λ ≈ 1.0, no sig calls, top z 4.6.
  - **v3: top-1 0.80, top-10 0.82, best similarity 0.89** vs STREME 0.55 / 0.73 / 0.84, Comsequent 0.49 / 0.55 / 0.80,
    HOMER 0.31 / 0.35 / 0.77; median 26 s vs 70–110 s. Behind only at N = 500, π = 0.1 (0.4 vs STREME 0.8).
