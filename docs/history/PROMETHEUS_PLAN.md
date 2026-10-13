# PROMETHEUS: plan for the experimental successor to HASTUR F2

*2026-10-12. Scrutiny of `docs/PROPOSED_EXPANSIONS.md` against the measured state of HASTUR F2 (ROADMAP_v2 §1).
Crate: `rust/prometheus`.*

**Name.** Large models take fire HASTUR lacks: parameter sharing across sites, genome-informed context, and priors
learned from other experiments. PROMETHEUS steals that fire for a calibrated, algorithmic tool. In
PROPOSED_EXPANSIONS terms, the extra information comes from the **genome** (the null) and from **other experiments**
(the alternative), never from re-comparing a set with itself (§1 sufficiency).

**Ground rules.**
- `rust/hastur_core` stays frozen as F2, the confirmed reference.
- `rust/prometheus` starts as a byte-identical copy (verified 2026-10-12). Every new behaviour sits behind a flag,
  and with all flags off the output must stay identical to F2 (tested on every build).
- Each component is adopted only through a registered round with the ROADMAP_v2 §2 rules, and on fresh panels.

---

## 1. Verdicts, ranked by expected boost × fit with the theory and current state

| rank | proposal | measured weak spot it targets | fit with theory/state | expected boost | build cost | verdict |
|---|---|---|---|---|---|---|
| 1 | **§3 Mode A: matched-set background + stratified permutation** | GC-selected human nulls: F2 2.6–3.6 false calls per run, 4–8 masked (STREME 5–9) | exact extension of LN7. Max-T needs exchangeability only *within* strata, and counts stay additive, so the permutation trick still works | calibration on composition-shifted nulls: the remaining calibration failure | small: a matched-set builder plus a strata-aware permutation loop | **build first** |
| 2 | **§6.2 weighted family (prior from a motif bank)** | small N: N = 25 calls almost nothing; D1 10/21 regulons at native vs STREME's 15/21 | the same max-T with per-cell weights fixed in advance. FWER holds whatever the prior; only power depends on it | large at small N, which is exactly ROADMAP_v2 E1/E2 territory | small for bacteria (RegulonDB/RegPrecise bank exists); medium for human (needs an atlas) | **build second**; bacterial bank first |
| 3 | **§3 Mode B: tabled expectations from a genome index** | speed: representation dominates (held-out z scores 50k bg + 50k synthetic at 5k); background counting | expectations become a lookup; matched genome draws replace bg counting | about 3–5× at 5k peaks; screen and atlas scale (E3, E12) | medium: the genome index | **build with the index**, after 1 |
| 4 | §2 importance-weighted counts by context class | none measured yet; a power claim | additive-model optimal weighting; calibration preserved by permutation | 1.0–1.6× effective data, *only if* within-set context heterogeneity is large (σ ≥ 0.5) | needs the genome index plus coordinates for every fg region | **measure §2.6 first**; build only if median predicted gain ≥ 1.3 |
| 5 | §5 flank contrast against genome copies of the core | ChIP motif statement: P5 −0.06 unmasked, −0.02 masked vs STREME | small exact families; representation only, so calls are unchanged | modest; the gap is small once repeats are masked | medium; composition confounding needs two controls | after 1–3, if the P5 gap persists on E-series data |
| 6 | §6.1 stratified pooling across experiments | weak shared cofactors; site counts for representation | Mantel–Haenszel over per-experiment residuals | situational | medium | with the atlas (E12) |
| 7 | §4 multi-radius ball seeds | — | a small extension of the radius-1 ball | small (doc's own verdict) | small | later |
| 8 | §7 grammar tables | — | needs Mode A strata and repeat classes | new capability, prior art (SpaMo) | medium | later |
| 9 | §7b count-space fitting | speed at atlas scale | engineering | moderate; EM is about 20% of runtime at 5k, held-out scoring about 60% | medium | after Mode B |

**Why this order.** It follows the evidence:
- **Mode A fixes the one calibration failure we can measure.** Composition shifts beyond order-1 Markov are why
  GC-selected nulls fail, and repeat masking did not fix them.
- **The weighted family attacks the measured small-N limit,** the main remaining sensitivity gap in the bacterial
  screens that ROADMAP_v2 leads with.
- **Mode B is the main remaining speed lever.**
- **§2's power gain is a prediction** whose size depends on a quantity nobody has measured.
- **§5 targets a gap that masking already shrank to noise.**

---

## 2. Feasibility: data, memory, speed

**Lattice sizes.** Gapped cells (k = 4: 220 shapes × 256 cells; k = 5: 495 × 1,024) total 563,200. L-mer cells are
65,536 at L = 8, 1.05 × 10⁶ at L = 10 and 1.7 × 10⁷ at L = 12.

| item | size or cost | where | notes |
|---|---|---|---|
| hg38, soft-masked | 0.98 GB gz; 3.1 GB as 2-bit codes in RAM | box (`/data`) | already on the box (`~/data/chip_masked/_downloads`) |
| RepeatMasker families (UCSC `rmsk.txt.gz`) | ~150 MB | box | for repeat class beyond soft-mask yes/no |
| TSS annotation (GENCODE) | ~50 MB | box | distance-to-promoter stratum |
| Context-class track (repeat class × composition class) | 1 byte/position = 3.1 GB, or computed on the fly from flanks | box | §2.2 rule 6: classes exclude the motif footprint plus a margin |
| Genome index: per-class gapped tables | 563k cells × 64 classes × 4 B = 144 MB per table; leave-one-chromosome-out ×24 = 3.5 GB | box | L-mer tables: L ≤ 10 (64 classes: 270 MB); L = 12 only unstratified (67 MB) |
| Genome index build | ~4 × 10¹² increments (2.9 × 10⁹ positions × 2 strands × 715 shapes). At the measured 2.3 × 10⁹/s per core: ~30 core-minutes, **~1–3 min on the 64-core box** if tiled | box | must be cache-tiled (genome chunks × shape groups). Streaming the 17 GB window-code array once per shape would be ~25 TB of memory traffic |
| Matched-set index | ~15 M tiled 201-bp windows × 4 covariates (GC, CpG o/e, repeat fraction, TSS distance), bucketed | box, ~0.5 GB | K = 10 matches per region by bucket sampling: ms per set |
| Mode A run cost | same as F2 (bg = matched windows, 10×); stratified permutation uses the same subtraction trick | anywhere | no slowdown |
| Mode B run cost | counting: fg only (bg tabled). Permutations: fg-sized matched draws. Representation: K = 2–5 matched bg instead of 10× bg + 10× synthetic | anywhere | expected ~3–5× faster at 5k peaks; small sets unchanged (already ~1 s) |
| Bacterial genomes | 2–8 Mb each; every table trivial | local or box | the weighted-family prior comes from curated banks (RegulonDB, RegPrecise), no atlas run needed |
| Human motif atlas for the §6.2 prior (E12) | ~2,000 ENCODE TF ChIP sets × ~7 s at 5k peaks | box: ~2 h | peak files small; leave-factor-out priors are a sum minus the factor's own row |

**WSL.** Nothing in the genome-scale items is built locally (see the disk notes: the vhdx is at its high-water
mark). The workstation is for code, unit tests and small-set benchmarks only.

**Coordinates.** §2, §3 and §5 need genomic coordinates for every fg region. The current builders write sequences
only. `data/prep_chip.py` and `bench/realnull_make.py` are deterministic, so a BED is added alongside each FASTA
(same windows) before any of this is measured.

---

## 3. Component designs (as they will sit in `rust/prometheus`)

**C1. Matched sets and stratified permutation (Mode A; flag `--strata FILE`).**
- **Input:** a strata file mapping every sequence (fg and bg) to a stratum id. Built by `prometheus match` from fg
  coordinates and the matched-set index.
- **Permutation:** each permuted fg draws, *within each stratum*, as many sequences as the real fg holds in it.
  Background counts are still pool − fg, so the cost is unchanged.
- **Statistic:** unchanged (F2). Only the null distribution changes.
- **Registered covariates (§3.2):** GC bin, CpG o/e bin, repeat fraction bin, and distance to TSS for ChIP. Fixed
  before any run; never tuned.

**C2. Weighted family (flag `--prior BANK`).**
- **Weights:** fixed in advance from a motif bank. A cell or ball whose consensus matches a bank motif gets weight
  w > 1, while a reserved share (1 − ρ) of α stays with de novo cells. The weights are normalised so that the
  weighted family-wise error stays at α.
- **Max-T:** taken over p_c / w_c.
- **Safeguards (§6.2):** the query's own factor is held out of the bank; each prior weight is capped; calls report
  whether they used prior weight.
- **Banks:** bacterial ones first (RegulonDB, RegPrecise), then a human atlas.

**C3. Genome index (subcommand `prometheus index`).**
- **Contents:** per-assembly, leave-one-chromosome-out class tables (gapped shapes, L ≤ 10), plus the matched-set
  window index and covariate buckets.
- **Build:** cache-tiled counting on the box.
- **Format:** flat little-endian u32 arrays with a JSON header (assembly, shapes, classes, sha256).

**C4. Tabled expectations (Mode B; flag `--genome-index DIR`).**
- **Expectations:** E(w) = Σ_r n_r(fg) · g_r(w) for the double-centred contrast. The background side becomes the
  genome's stratum-weighted rates.
- **Permutations:** fg-sized matched genome draws.
- **Representation:** a small matched background (K = 2–5) and its synthetic copy.

**C5 (conditional on §2.6).** Context-class weighted counts. These are class-stratified fg histograms, giving
S(w) = Σ_c a_c (n_c − E_c) with cross-fitted a_c ∝ s_c / q_c, permutation-calibrated.

**C6 (later).** Flank contrast against genome copies of a core (§5). It scans the 2-bit genome for the core's copies
(a few seconds per core on the box) and applies the two required controls.

---

## 4. Order of work, with gates

| step | what | gate to continue | compute |
|---|---|---|---|
| 0 | `prometheus` = F2 baseline (done); identity test in the build script | outputs identical to `hastur_core` | — |
| 1 | BED outputs for the ChIP and real-null builders; matched-set index for hg38 (box) | the rebuilt windows are identical to the existing units | box, minutes |
| 2 | **No-build measurements** (PROPOSED §9.1): (a) §2.6 context heterogeneity → predicted gain on dev ChIP; (b) copy-number and composition profile of F2's false calls and permutation maxima on GC-selected nulls; (c) share of true dev sites in high-copy contexts | (a) decides C5; (b) predicts what Mode A must capture | box, about an hour |
| 3 | **C1 Mode A as a registered round:** arms = matched sets; + stratified permutation; + copy-number collapse. Decided on fresh GC-selected and exchangeable nulls (tuple seeds, freshness check), plus ChIP dev sensitivity | GC-selected false calls ≤ 0.1 per run; exchangeable unchanged; ChIP sensitivity within 0.03 | box, ~1 day |
| 4 | **C2 weighted family with the bacterial bank:** registered on D1-style regulon sets (fresh random-promoter nulls) | native regulon recovery up, false calls per run ≤ 0.1 | local or box, hours |
| 5 | C3 + C4: genome index and Mode B | outputs within registered tolerance of Mode A; ≥ 2× faster at 5k | box |
| 6 | C5 if step 2a passed; C6 and §6.1 against the E-series data | each by its own registered rule | box |

**Hand-off to ROADMAP_v2.** Steps 3–4 feed E1 (iModulon), E2 (perturbations), E5 (extreme-composition genomes),
E7 (ATAC contrasts with matched backgrounds) and E8 (low-input). Step 5 is the speed basis for E3 (pangenome) and E12
(atlas).

---

## 5. Risks recorded now

- **Mode A changes the question.** A matched null says "beyond what the covariates explain". Promoter-generic motifs
  disappear from promoter-rich sets by design. Reports state the strata, and the unmatched F2 run stays available.
- **Strata too fine.** With K < 10 per stratum the permutation cannot resolve the threshold. The builder refuses,
  coarsens the strata, and logs it.
- **Prior leakage and usual-suspect bias** (§6.2): leave-factor-out banks, capped weights, and calls flagged when they
  used prior weight.
- **Identity drift.** Every new flag defaults off. CI compares outputs against `hastur_core` on the eight validation
  units.

---

## 6. Status 2026-10-12: development paused

**M-series (Mode A, matched genome backgrounds).**

- **MA1** (`splits/protocol_MA1.md`; fresh nulls with freshness < 0.03%; held-out test2 ChIP):
  - GC-selected false calls per run fell from 3.70 / 2.70 to 0.05 / 0.05 (N = 2,000 / 500).
  - Exchangeable nulls stayed at ≤ 0.10.
  - ChIP sensitivity went from 23/27 to 22/27 (lost E2F1, a CpG-promoter factor). That is −0.037 against a gate of 0.03,
    so **not adopted as default**.
  - Stratified permutation gave no benefit over the plain matched background.
  - The STREME-on-matched-background arm was stopped by the user, so it remains open whether the gain is the
    background's.
- **MA2** (test3 ChIP, 94 factors; stopped by the user at 58 completed): sensitivity was identical, 40/58 in all three
  arms. Mode A gained GATA3 and lost KLF1, and made about 0.6 fewer native calls per factor.
- **Verdict:** Mode A removes composition-shift false calls (a real but narrow failure). On real ChIP with GC-matched
  backgrounds the effect is small and costs some promoter-bound factors. It is worth keeping as an explicit option,
  not as a step change.

**Other components.**
- **C2, cell-level bank prior:** rejected on dev (the prior was too broad).
- **C5, context weighting:** failed the §2.6 rule (median predicted gain 1.23 < 1.3).
- **Known-motif mode:** dominated by SEA.
- **C3/C4 (Mode B), C6 and §6.1:** not started. Their expected gains are speed only, or small.

**Literature check (`docs/DNA_deep_models_motif_lit_review.md`).** The review's large documented gains are in site
modelling (occupancy scoring, higher-order context plus flanks, multiple modes, spacing kernels), measured as
site-classification AUC or AvRec. Existing tools (BaMM, Dimont, ProBound) already own that space. HASTUR's measured
limit is detection power at small N, and the review's §7.6 information floor explains it rather than removing it. No
item was judged a large expected gain, so development is paused pending a new angle (user decision).

**Infrastructure left in place.**
- `rust/prometheus`: the `--strata` and `--prior` flags, each verified identical to F2 when off.
- `tools/matched_index.py`, `match_sets.py`, `step2_measure.py`.
- The hg38 index on the box's EBS volume (`/data`).
