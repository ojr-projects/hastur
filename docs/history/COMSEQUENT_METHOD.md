# Comsequent: the method, what works, what doesn't, and why

> **Jumping-off point for a new project.** This is a self-contained digest of the Comsequent method as of its freeze
> (`v1.0-freeze`, 2026-10-09) and the real-data results to 2026-10-10. It covers the theory, the components and their
> verdicts, the evidence, the limits, and the ideas left untested.
> The source repository is `../comsequent` (Go engine in `lattice/` and the repo root; benchmarks in `benchmarks/`).
> Its key documents are listed in §11. For process learnings, see `DEV_INSIGHTS.md` in this folder.

---

## 1. Summary

Comsequent **counts every sparse DNA pattern exactly** (up to k specified bases at any offsets within a window of
width w) and scores them with **lattice statistics**: interaction terms obtained by inclusion–exclusion (Möbius
inversion) over each pattern's sub-patterns. In **contrast mode** (foreground vs a matched background) it is a
calibrated, fast, exhaustive motif finder.

**Bottom line from the benchmarks:**
- It **wins** where the specificity is **exact and discrete** (a crisp, possibly gapped or multi-part pattern) **and** the
  contrast is **near-binary and matched by construction**. The proven case is bacterial methylation motifs from nanopore
  data, especially at low coverage and for bipartite Type I motifs.
- It is uniquely **exact on gapped multi-part motifs**, **calibrated** where dyad tools are not, and **20–100× faster**
  than STREME/HOMER at large N.
- It **loses** wherever the signal is **PWM-like**: long, low-information, additive sites (eukaryotic TF ChIP, HT-SELEX,
  RBNS, bacterial σ70 promoters). There it should seed or complement PWM tools, not replace them.

---

## 2. Data structure and engine

- **Pattern:** up to k specified bases at chosen offsets in a window of width w. The *shape* is the offset set; the
  *pattern* is shape + letters. Typical settings: k5/w13, k6/w13–15, k7/w15 for long bipartite motifs.
- **Lattice:** patterns ordered by sub-pattern inclusion. Every sub-pattern of a counted pattern is also counted.
- **Implicit indexing:** layer (l, w), then shape rank × 4^l + 2-bit packed bases. There is no explicit object per
  pattern (k7/w15 = 58.5M patterns in ~1 GB instead of ~4 GB).
- **Counting:** a direct CPU shape enumeration over packed 2-bit windows, block-wise (2^26 positions per block, ~0.6 GB),
  both strands by default, with non-ACGT bases ending a run. It replaced a per-pattern CUDA kernel that could not scale.
  - E. coli genome (4.6 Mb, both strands): 5/13 in 4 s, 6/15 in 21 s, 7/15 in 46 s (4.4 GB peak).
  - k = 7 contrast runs peak around 4.75 GB after memory work (8.1 GB before).
- **Möbius stage:** precomputed per-shape sub-pattern tables, about 13× faster than naive.
- **Commands:** `comse discover` (exhaustive statistics), `comse extend` (upward walk, refinement, all-in-one output),
  `comse anchored` (windows centred on a base), `comse query` (targeted counts for a motif list), `comse mmseed`
  (mismatch-neighbourhood seeds, prototype). The legacy `calc/apply/run` sequence-scoring pipeline still exists but is
  not the useful path (§8).

---

## 3. Theory

### 3.1 Two inversions from the same counts

Let f(T) be the count of windows matching every specified position of T.

- **Marginal-log (Möbius) interaction:** I(S) = Σ_{T⊆S} (−1)^{|S|−|T|} log f(T). This is the original method's
  statistic.
- **Log-linear interaction θ:** recover the 2^|S| match/mismatch **cell table** n(x) from the down-set by
  inclusion–exclusion, then θ_S = Σ_x (−1)^{|S|−|x|} log n(x).

**Closed form for a planted whole motif** (i.i.d. background, per-position match probability b, motif in a fraction π of windows):
- f(T)/N = (1−π)·b^{|T|} + π, so I_k = Δ^k h(0) with h(j) = log((1−π)b^j + π). h is hinge-shaped, so its
  high-order finite differences **alternate in sign and stay small, and can be negative**.
- Only the all-match cell is perturbed, so **θ_k = log(1 + π / ((1−π)·b^k)), which grows monotonically with order.**

| π (b = 0.25) | I₂ | I₃ | I₄ | I₅ | θ₂ | θ₃ | θ₄ | θ₅ |
|---|---|---|---|---|---|---|---|---|
| 0.01 | 0.08 | 0.16 | 0.03 | −0.28 | 0.15 | 0.50 | 1.28 | 2.43 |
| 0.05 | 0.28 | 0.16 | −0.27 | 0.18 | 0.61 | 1.47 | 2.67 | 4.01 |
| 0.30 | 0.42 | −0.20 | 0.05 | 0.05 | 2.06 | 3.35 | 4.71 | 6.09 |

**Consequence: raw Möbius I cannot see whole enriched motifs.** Measured: I_con on GCNNNANNNNGC was −0.007 against
θ_con +0.55, and strict recall for I-based arms was about 0 for every class. **Use θ, not I, for enrichment.**

**Depletion is the opposite case.** Removing most full-pattern occurrences barely moves the lower marginals, so all
the signal lands in the top-order term, and **both inversions see it cleanly** (k = 5, 95% depletion: I₅ = −2.46,
θ₅ = −3.49). In simulation, Möbius I_con, θ_con and ENR reached recall@1 = 1.0 at 50% depletion, against 0.0 for the
Markov contrast and 0.6 for STREME. (On real genomes the Möbius avoidance idea still lost to Markov expectations; see §7.)

### 3.2 Attribution: θ peaks at the true pattern, z peaks at low order

7-position motif in a 9-wide window, π = 0.05:

| pattern | θ | z |
|---|---|---|
| 2-sub-pattern | 0.61 | 82 |
| 5-sub-pattern | ~4.0 | ~43 |
| 4 motif positions + 1 irrelevant | 0.06 | 0.7 |
| **full 7-position motif** | **6.03** | 11 |
| 7 + 1 irrelevant | −1.35 | −1.1 |

- Every k-subset of the motif has the same θ_k; adding an irrelevant position sends θ back to about 0. **The effect
  size identifies the true pattern.**
- SE(θ) ≈ sqrt(Σ_x 1/n(x)) is dominated by the smallest cells, which shrink with k. **Ranking by significance therefore
  favours fragments.** This is the **fragment effect** in statistical form.
- **What worked in practice:** rank by a sensitive statistic (ENR z) and **filter by θ-maximality**: keep a call only if
  no immediate sub-pattern has a θ significantly higher (1.96 SE) and no one-position extension raises θ significantly.
  This cut redundancy about 5× at equal recall. On a 6-position fragment test it called 0 of 49 fragments.
- **What failed:** ranking directly by a θ effect-size lower bound. The winner's curse put junk high-order patterns on top.

### 3.3 Variance and calibration

- **Exact delta-method variance of I for every pattern** in O(l·2^l) via the cell form:
  Var = Σ_W cell(W)·A(W)² − 1, with A the subset sum of μ_T/p_T. The earlier diagonal approximation **overstated SE about
  9× at l = 7** (null λ ≈ 0.01), which crushed every k = 7 z. After the fix, λ = 0.96–1.02 at l = 3–7 on i.i.d. nulls.
- θ uses the Woolf SE.
- **Background-expectation variance must be included** in any statistic that compares foreground counts with a background
  model ("mobbg": fg count vs the background's Möbius expectation). Without it, false calls on same-source splits were
  2,522; with it, 0.
- **Score (Pearson) z, not Wald z,** for log odds ratios. Under separation (near-zero cells, common in methylation
  windows) the Wald z collapses: GATC ranked about 10,000th.
- **Optional `-calib` table:** λ per (mode, statistic, l, gap class), built from motif-free splits; z is divided by √λ.
- **Multiple testing at k = 7 (~49M tests):** the all-in-one `sig` column is a two-sided Bonferroni bound over discover's
  tests, or a call table's own BH ≤ 0.05. This gives about 0.3 sig rows per motif-free run.

### 3.4 Contrast statistics and robustness to composition

| statistic | definition | behaviour |
|---|---|---|
| **ENR** | log odds of occurrence, fg vs bg (score z) | most sensitive detector; **floods under composition mismatch** (~157k BH calls on a motif-free GC 0.5 vs 0.3 pair) |
| **mobbg** | fg count vs background Möbius expectation (+ background variance) | sensitive; inherits the background model's misspecification |
| **MKV / mkvbg** | count vs order-(l−2) contiguous Markov expectation from the same counts | the natural baseline; anti-conservative on structured backgrounds |
| **I_con, θ_con** | I_fg − I_bg, θ_fg − θ_bg | **immune to composition shifts** (0 false calls under the GC shift, λ 1.03), because interaction terms cancel main-effect differences. Still leaks under coding vs non-coding (λ 1.7) |
| **robust ENR** | ENR gated on θ-contrast z ≥ 3.29 with the same sign | ENR's power with θ's robustness (~190 false calls vs 157k; planted recall 0.59 vs 0.51) |
| **factorial dd** (`-i2/-bg2`) | (fg1 − bg1) − (fg2 − bg2) for ENR, Möbius or θ | cancels a confounder shared by two contrasts; removed coding leakage in a test (ENR 20,000-call cap → dd_theta 6 calls, motif first) |

**Theory result (I-projection):** the minimum-discrimination-information projection keeps the background's top-order θ
while matching the foreground's lower margins. Its likelihood-ratio test is therefore asymptotically **θ_fg = θ_bg**,
the θ-contrast arm itself. No separate I-projection arm is needed.

**Single-set mode (no control) is not usable on real genomes.** Markov and codon structure inflates low orders
(λ up to 1,700 at l = 2 in coding sequence). The top calls are composition.

### 3.5 Why exact patterns lose to PWMs: the dilution argument

Every motif finder models the joint distribution of a window as a hierarchical log-linear expansion:
log p(x) = c + Σ f_i + Σ g_ij + Σ h_ijk + …

| | PWM tools (MEME, STREME, HOMER) | Comsequent |
|---|---|---|
| terms kept | main effects only, all 4 letters per position | interactions of all orders ≤ k, but only for one letter combination's match indicators |
| unit of evidence | a site *scored softly* (variants still count) | a site *matched exactly* |
| site model | latent mixture, EM or discriminative alignment | none: marginal counts over all windows |
| search | heuristic, complete over letters | exhaustive and exact over offsets × letters, bounded by k |

- A moderate-IC PWM has many site variants. The exact consensus at 5 positions matches only Π p(consensus) of sites
  (e.g. 0.6⁵ ≈ 8%). Each exact pattern carries a small slice of the signal; a PWM sums over all variants.
  This is **conjunctive-exact vs additive-soft**, the same mechanism as the fragment effect.
- **Where Comsequent wins instead:** position dependencies (invisible to first-order models), gapped/bipartite structure
  (exhaustive offsets, no seed failures, no EM local optima), and robust contrasts (interaction terms cancel composition).
- **What stays hard:** long, low-information additive motifs (~0.7 bits/column over 12+ columns). No single sparse or
  degenerate pattern carries enough signal; only aggregated scoring sees them.

### 3.6 Unifying moves (how the lattice relates to other model families)

- **A. Counts are sufficient statistics.** Exact sparse-pattern counts are the marginal tables of every hierarchical
  log-linear model of order ≤ k within w. A PWM is the order-1 model; with the alignment fixed (anchored mode, or sites
  anchored on a core) a PWM *is* a set of lattice marginals.
- **B. The lattice selects terms.** θ-contrast and θ-maximality say which higher-order terms are needed beyond the additive
  model: "PWM main effects + significant interactions", fitted by IPF or discriminatively (prior art: MaxEntScan,
  TFFMs, BEEML-style models).
- **C. Seeds → soft refinement** (how STREME and MEME work). Comsequent is a richer seed generator: gapped, bipartite,
  exhaustive, calibrated.
- **D. Degeneracy without recounting.** An IUPAC pattern's count is the sum of its members, so a degenerate lattice
  (10 symbols per position) is derivable from existing counts.
- **Mismatch identities.** The number of windows matching ≥ k−1 of P's positions is Σᵢ n(P∖i) − (k−1)·n(P). The full
  "number of matching positions" distribution follows from P's 2^k cell table. A soft-scoring ladder therefore costs
  no recount.
- **Extensions that worked were functions of counts the lattice already holds** (e.g. the PWM from the walk's frame
  histograms). Bolting on external models did not.

---

## 4. Components and verdicts

| Component | What it does | Verdict |
|---|---|---|
| Exact CPU sparse counter, implicit index | count all patterns ≤ k in w | **Works.** Fast and exact; the foundation |
| Contrast ENR (score z) | detection | **Works.** The most sensitive detector; needs a matched background |
| θ statistic + θ-maximal filter (v3 symmetric rule) | attribution, redundancy removal | **Works** |
| Exact delta-method SEs, background variance, score z | calibration | **Works** (essential) |
| Robust ENR; θ_con; factorial dd | robustness to composition / confounders | **Works** as safeguards; does not fix a badly mismatched control |
| Upward walk (`extend`): single-base extensions at z > 4.5 on conditional log-OR contrast, then IUPAC prune/merge by excess share | grows patterns beyond k; degenerate letters | **Works** (EcoKI, GCNNNANNNNGC, GRCGYC exact). Guards needed: pattern-repeat stop, MinShare prune, pair-share guard, no merges at self-complementary positions |
| Cross-fitted walk (accept only if significant in both halves) | fewer chance extensions | **Works**, modest (25 vs 23 of 40 truth motifs) |
| Seed joins (`-join N`, lattice join of overlapping seeds) | bipartite / mutated bipartite | **Works** |
| U1 core PWM from walk frame histograms | PWM output | **Works**: PWM class strict recall 0 → 0.8–1.0 |
| U2 cross-fitted discriminative PWM refinement (+ ZOOPS) | PWM quality | **Works at N ≥ 2k** (matched STREME at 10k: 1.00 vs 0.97); **fails in small samples** (locks onto seed letters) |
| All-in-one ranked output (u14): pool walk patterns, refined PWMs, θ calls, ENR calls; rank by z; promote one refinement within √2 z ahead of its parent; reserve slots 4/7/10 for BH-significant θ-maximal calls; `sig` column | one output | **Works**; frozen. Core grid strict @1 0.33 → 0.41 |
| Anchored lattice (`comse anchored`) | positional patterns around a modified base (`GxTC`) | **Works**; the winning real-data mode. Also the right command for any positional assay (PAM: 0.08 → 0.50 vs `discover`) |
| Explaining-away rounds (methylation): find motif, remove its windows, repeat; **10 rounds, ≥ 3 specified bases per round** | multiple MTases per genome | **Works**; biggest single gain (pilot strict 0.61 → 0.83; bipartite 0.11 → 0.86). A 5-round cap and 1–2-base composition rounds were the defect |
| Residual cross-fitted walk (methylation) | each round's seed walks on its residual windows | **Works** |
| Soft labels / min-coverage 2 (U11) | use low-coverage sites | **Works at 5×** (0.74 vs Nanomotif 0.47; bipartite 0.58 vs 0.09); soft weights alone diluted by single-read sites |
| Pooled pair-table epistasis test (U7) | interaction motifs | **Partial.** Detects epi_pair at π·N ≥ ~3,000 and the triple at large N; an earlier "inroad" counted both strands and was anti-conservative (retracted) |
| Ranking by θ effect-size bound | attribution | **Rejected** (winner's curse) |
| Raw Möbius I / I_con as ranking statistic | enrichment | **Rejected** (§3.1) |
| Single-set discovery on real genomes | label-free | **Rejected** (composition dominates) |
| Möbius expectation for RM-site avoidance (Tier 0b, 17 genomes) | depletion on real genomes | **Negative**: Markov beat Möbius (AUROC 0.69–0.71 vs 0.63–0.65); likely because avoided motifs' sub-patterns are also depleted, so the Möbius expectation absorbs the deficit |
| u10–u12, u15, u20 ordering guards | rank tweaks | **Trade-offs only** (±0.02); plateau |
| u17 softened refinement, u18/u19 + `mmseed` mismatch seeds | small-sample PWMs | **Negative**: detection improved, refinement did not |
| Parsimony pruning; one-pass refit; fg-only refinement | various | **Negative / reverted** |
| Empirical-Bayes shrinkage per order (E4) | low-count stability | **Deferred**: τ² ≈ 0 because motifs are a tiny fraction of patterns; needs a two-group mixture prior |
| Seeding HOMER with Comsequent calls | hybrid workflow | **Withdrawn**: promoters 0.44 → 0.56, never beat HOMER alone (0.75) |
| Original sequence-scoring pipeline (promoter detection) | label-free promoter scoring | **Not supported**: weak or inverted on genomic backgrounds (AUC 0.26–0.57); with a Markov expectation it only matched simple baselines (~0.55) |

---

## 5. The frozen pipeline (`v1.0-freeze`)

**General contrast discovery:**
1. `comse discover -i fg.fa -bg bg.fa -k 5 -w 13` (k7/w15 for long bipartite motifs) → call tables incl. `calls_theta`,
   `calls_theta_tmax`, `calls_enr*`.
2. `comse extend -join 30 -theta <calls_theta> -thetamax <calls_theta_tmax>` → `<o>.all`, the all-in-one ranked list.
   **Only `sig = 1` rows are discoveries.**

**Methylation (the winning application):** dorado SUP v5 with 6mA and 4mC/5mC models → alignment → modkit pileup →
foreground windows around modified sites, background around unmodified sites of the same base (thresholds matched to
Nanomotif: Nvalid ≥ 5, unmethylated < 0.3) → `comse anchored` with θ-maximal calls → explaining-away rounds
(`P1_EXPLAIN_ROUNDS=10`, `P1_EXPLAIN_MIN_L=4`) → residual cross-fitted walk.
Pipeline scripts: `benchmarks/phase1_methylation/11_meth_pipeline.py`, `15_meth_downsampling.py`, `17_meth_heldout.py`.

**Positional assays (e.g. PAM depletion):** use `comse anchored`, not `discover` (which pools a pattern over all offsets).

---

## 6. Evidence: simulations

Strict recall = exact motif (IUPAC members count) at rank 1 or in the top 10. "Best comparator" = per-dataset best of
STREME/HOMER/MEME (an oracle choice, conservative against Comsequent).

**Core grid** (1,050 datasets, Markov-2 E. coli background, π 0.01–1, N 500/2k/10k), all-in-one, @10:

| class | Comsequent | best comparator |
|---|---|---|
| gapped GCNNNANNNNGC | 0.46 / 0.66 / 0.80 | 0.22 / 0.20 / 0.00 |
| bipartite EcoKI | 0.58 / 0.80 / 0.80 | 0.38 / 0.44 / 0.62 |
| contiguous GATC | 0.44 / 0.60 / 0.72 | 0.26 / 0.38 / 0.50 |
| degenerate GRCGYC | 0.42 / 0.66 / 0.78 | 0.36 / 0.52 / 0.72 |
| PWM 8 bp mid IC | 0.42 / 0.64 / 0.80 | 0.60 / 0.70 / 0.80 |
| epistatic pair / triple | ≤ 0.22 / ≤ 0.12 | ~0 / 0 |
| **all** | **0.48** (@1 0.41) | 0.32 (@1 0.31) |

The artefact-free grid (flanks resampled under the background): all 0.68 vs 0.39 @10.
Sweeps (2,720 datasets, 13 backgrounds): @1 0.454.

**Gapped-capable panel** (RSAT dyads + matrices at monad 3 and 2, Weeder 2, GADEM, XXmotif, GLAM2; per-dataset best of 8), strict @10:
gapped **+0.40 / +0.73 / +1.00**; bipartite +0.13 / 0 / 0; contiguous and degenerate ≈ 0; **PWM −0.40 / −0.13 / 0**.
The panel detects the GC…GC pair but never the internal base.

**Calibration** (motif-free splits): Comsequent 0 false motifs (ENR) or ~0.3 sig rows (all-in-one); RSAT dyads ~21,
XXmotif ~5, HOMER 1.7–2.7, STREME 0–0.3. Under composition mismatch every enrichment method floods; θ/I_con resist.

**Epistasis:** a three-way motif needs about 10⁴ planted sites (N = 50k / π 1: 1.0); every comparator scores 0.
Exactly pure higher-order structure (all lower margins pinned) is intrinsically near the detection floor: at k = 4,
p(all match) only moves 0.0039 → 0.0116 even at a ×1000 tilt.

**Speed** (median per dataset, N 500 / 2k / 10k): all-in-one 8 / 30 / 130 s; STREME 2 / 8 / 78 s; HOMER 40 / 80 / 214 s.
Contrast ENR alone: 4 s at 10k, 22 s at 50k (STREME 35 min, HOMER 70 min). The walk with joins dominates all-in-one time.

---

## 7. Evidence: real data

| Phase | Data | Result | Verdict |
|---|---|---|---|
| **P1 methylation** (dev, 10 strains, Nanomotif truth) | PRJEB74343 nanopore, 100× cap, downsampled 5–40× | strict **0.69 vs 0.66**, lenient **0.79 vs 0.71**; bipartite 0.70 vs 0.56; **5×: 0.70 vs 0.47**; full coverage Nanomotif ahead on strict (0.74 vs 0.65), level on lenient | **Win** (tuning on these strains disclosed) |
| P1b REBASE truth | 47 MTase sites, same 10 strains | strict **0.64 vs 0.46**, lenient 0.79 vs 0.47; 5× 0.58 vs 0.30 | **Win** (species-level truth; not proven independent) |
| P1 held-out | 25 new strains, 102 motifs, registered | pending | predicted win |
| P6 promoters | σ70 boxes, 16 bacteria | −10 box top 10: 0.44 vs STREME 0.69 / HOMER 0.75; −35 box ~0 vs 0.44 / 0.69 | Loss (PWM-like) |
| P5 ChIP (negative control) | 24 ENCODE TFs | top 10: 0.67 vs STREME 0.92 (HOMER 0.58); misses long/degenerate (REST, SP1, ZNF143, NRF1, …) | Loss (predicted) |
| P3 RBNS | 53 dev RBPs | top 10: 0.79 vs STREME 0.89 / HOMER 0.87 | Loss |
| P4a HT-SELEX | 162–249 dev libraries | top 10: 0.35 vs 0.72 / 0.72; top calls are GC composition (CCCCC, GGNGGG) from cycle drift | Loss |
| P4b CAP-SELEX | 83 dev pairs | composite readout **0 for every tool**; half-sites 0.06 vs 0.39 / 0.33 | Loss / null |
| P7 PAM depletion | HT-PAMDA, 12 samples | anchored 0.50 vs HOMER 0.75 (registered rule) | Parity-to-loss (graded depletion, ~2 informative positions) |
| Tier 0b RM avoidance | 17 genomes, REBASE motifs | Markov expectation beats Möbius | Negative for the Möbius idea |

### Domain-of-applicability rule (derived from the above)

Comsequent is strongest when **all** of these hold:
1. **The specificity is discrete and exact** (crisp pattern, not a PWM);
2. **The contrast is matched by construction** (same genome, same pool; differs only in the signal);
3. **The contrast is near-binary** (sites essentially present vs absent, not 20–50% graded depletion);
4. **The motif has structure to exploit** (several parts, a gap, or real length; not ~2 informative positions).

**Pre-check, computed before any discovery run:** split the fg/bg divergence (KL at k = 3) into **drift** (the reference
motif masked) and **signal** (total − drift).

| phase | drift KL₃ | signal / drift | outcome |
|---|---|---|---|
| methylation | 0.022 | **3.5** | win |
| HT-SELEX | **0.104** | 0.39 | loss (high drift) |
| ChIP | 0.008 | 0.21 | loss (low signal) |
| RBNS | 0.018 | 0.11 | loss (low signal) |
| PAM | 0.005 | 0.11 | parity-to-loss (low signal) |

This separates two loss mechanisms: **high drift** and **low signal**. The hypothesis is that signal/drift > ~1
predicts a win. It is fitted to five phases and still needs confirmation (tools: `drift_check.py`,
`composition_check.py` in `eval/`).

---

## 8. Limitations

| Limit | Evidence | Why it is structural |
|---|---|---|
| Low-site / low-IC PWMs | N 500, π ≤ 0.1: 0 vs STREME 0.4–0.8, MEME 0.67–1.0; all PWM-like real phases lost | Exact patterns dilute soft signal (§3.5); refinement on cross-fitted halves fails in small samples; fixing it means MEME-style EM, i.e. rebuilding existing tools |
| Composition-driven ranking | HT-SELEX GC calls; RBNS low-information calls (ANNA, RRNNR); θ-quota slots admitting 2–3-base calls | When fg and bg differ in composition, the frozen ranking puts the shift first. A candidate fix exists (§10) |
| Background mismatch in general | every enrichment arm floods | No contrast method can fix a mismatched control; θ/robust ENR/factorial only mitigate |
| Short motif vs its own extensions | GATC at rank 1 trades against degenerate/gapped cells under every rule tried | Evidence genuinely ambiguous at moderate z |
| Epistasis needs large N | ~10⁴ sites | Interaction tests need data |
| Very low-information gapped motifs | mutated CANNNNNNTG: 0 for every tool | ~80 sites against ~700 chance matches |
| k ceiling | k ≤ 7 exhaustive; longer via the walk only | 4^k × shapes; k7/w15 is already 58.5M patterns |
| Single-set mode | composition dominates | Needs a control |
| Positional signals in `discover` | PAM 0.08 with `discover` | `discover` pools over offsets; use `anchored` |
| Both-strand palindromes | counted twice | Must be stated in methods |
| Runtime vs Nanomotif | 2–4× slower per methylation run (~400 s vs 160 s) | Walk and rounds dominate |
| Adoption | at parity, STREME/HOMER are entrenched and faster at small N | The case for use has to be the winning domain, not generic motif finding |

---

## 9. Implementation pitfalls already hit (check these in any re-implementation)

- Angle wrap with a negative `math.Mod`; layers assumed contiguous ("hollow" models dropped layers); off-by-one z
  bins; NaN z passing every filter; a cache keyed too loosely (real counts reused for scrambled runs); non-deterministic
  map order. All in the original code.
- Diagonal SE approximation above l = 6 (9× overstatement).
- Wald z under separation.
- Background-model variance omitted (anti-conservative mobbg/mkvbg).
- Counting both strands in a pooled test (anti-conservative).
- Walk cycling (an extension and a prune undoing each other forever); strand mixing at self-complementary positions.
- Explaining-away rounds consumed by 1–2-base composition calls (one removed 21k windows); a hard round cap exhausted
  by multi-MTase genomes (each strand of a non-palindromic motif costs a round).
- Input de-duplication on by default in the legacy CLI (`-dr 0.8`); it must be off for benchmarks.
- Legacy GPU counter encoded non-ACGT as A, and its w = 2 counts included a padding base.
- Simulator: planted N positions drawn uniformly into a Markov background create a context artefact; plant them from the
  background chain.

---

## 10. Open directions (not done; candidate starting points)

**Method:**
- **Composition-routing fix (v1.1 candidate):** measure fg/bg divergence; above a threshold rank by robust ENR or θ
  rather than plain ENR, and require minimum information content for reserved slots. Develop on dev splits (RBNS,
  HT-SELEX dev), then confirm on the locked test splits.
- **U7 full shape-table log-linear model:** fit main effects + selected interactions per shape (fg vs bg) and choose the
  order by held-out likelihood. Replaces the all-in-one heuristics with likelihood-based model selection ("PWM +
  dependencies" with exhaustive, calibrated term selection).
- **U3 exhaustive degenerate lattice** from summed counts (R/Y/S/W/K/M at up to 2 positions) with no recount.
- **U8 mismatch ladder** as a calibrated soft detector (the identities in §3.6); detection helped, refinement didn't, so
  pair it with a better refiner.
- **U10 spacer merge** for variable-linker bipartite motifs (needs a variable-spacer simulation class).
- **U12 extra count axes:** sequence-level presence (ZOOPS) for repeat robustness; position bins for centrality.
- **U13 L1-logistic over θ-maximal patterns:** an interpretable gkm-SVM.
- **R-b background re-estimation** with discovered sites masked (a motif partly explains itself away in its own Möbius
  expectation); might revive Möbius vs Markov for depletion.
- **Continuous-label methylation GLM:** binomial outcome per site, Möbius-structured predictors; the expected gain is at 5–10×.
- **Order-wise shrinkage with a two-group (null + signal) mixture prior.**
- **IUPAC merging as a Breslow–Day homogeneity test** across concrete instances.

**Applications in the winning class (registered as next tests, data availability to verify):**
- Held-out methylation strains (25 strains, registered; ready to run).
- Natural Cas ortholog PAMs with long crisp motifs (Gasiunas 2020, 79 orthologs) rather than relaxed engineered variants.
- Phosphorothioate (Dnd) modification motifs (anchored lattice directly applicable).
- Base-editor / deaminase context (anchor on edited base vs unedited of the same type).
- Transposon/integrase target sites; nuclease off-target context.
- Methylome-based strain tracking; metagenome MGE–host linkage (needs paired nanopore + Hi-C).

**Strategic reading:** the defensible product is not a general motif finder. It is a **calibrated exact-pattern engine
for enzyme-specificity problems with matched controls**: methyltransferases, restriction systems, PAMs, modification
and editing contexts, with a pre-check that tells users in advance whether their data is in the domain.

---

## 11. Source documents (in `../comsequent`)

| Document | Content |
|---|---|
| `docs/ENGINE.md` | Command reference and which statistic to use |
| `docs/UNIFIED_MOTIF_MODEL.md` | Dilution theory, unifying moves, U1–U13 and R-a–R-d with an update log |
| `docs/briefs/comsequent_theory_addendum.md` | I vs θ closed forms, attribution, contrast/calibration theory, phase-specific theory |
| `docs/briefs/comsequent_differential_brief.md` | Differential (fg vs bg) Möbius theory: contrast inversion, product lattice, I-projection, factorial lattice |
| `benchmarks/phase0_sim/DECISION.md` | Phase 0 go/no-go, claims and non-claims |
| `docs/SIM_DEV_FREEZE.md` | Frozen method, competitor tables, known limits, held-out protocol |
| `docs/REAL_DATA_PLAN.md` | Domain rule, prediction registry, drift pre-check, all real-data phases |
| `docs/COLLABORATION_STORY.md` | Chronology, metric tracks, mistakes |
| `docs/working/NOTES.md` | Legacy CLI/output facts (original sequence-scoring pipeline) |
| `docs/briefs/comsequent_methods.pdf` | Original methods write-up (code is authoritative where they differ) |
