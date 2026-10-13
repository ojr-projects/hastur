# PROPOSED_EXPANSIONS: genome-scale backgrounds and exhaustive comparisons for HASTUR

> **What this is.** Proposals for two new operating modes and several components. None of them has been built or
> tested. Each idea went through one serious round of self-doubt, and only what survived is written up here. Some
> survived only in narrower form, and the narrowing is recorded with each one. §8 lists what was dropped and why.
>
> **Working assumptions.** HASTUR runs about 1000× faster than STREME at roughly the same calibration and accuracy.
> The current dev numbers in the repo are deliberately not used, because development is ongoing.
>
> **Labels** follow THEORY.md: *Derivation* (argued here), *Prediction* (to be tested), *Prior art* (known elsewhere).
> Started 2026-10-11.

---

## 0. Summary

| # | proposal | what it buys | main limit found under doubt | first measurement | priority |
|---|---|---|---|---|---|
| 2 | Context-class background with importance-weighted counts | power: gain mean(q)·mean(1/q) ≥ 1 | gain is zero where contexts are homogeneous; needs genomic coordinates | heterogeneity of q across dev peak positions (§2.6) | 1 |
| 3 | Matched-set background + stratified permutation; tabled expectations | calibration under shifts the strata capture; faster runs | exact only for the covariates matched; the null defines the question | GC-selected / repeat-rich null sets, stratified vs global | 1 |
| 5 | Flank contrast conditioned on the core, against its genome copies | flank and core×flank specificity, calibrated | composition confounding between bound and genomic copies | held-out site likelihood with vs without flank terms | 2 |
| 6 | Borrowing across experiments (stratified pooling, empirical-Bayes family) | the only real source of extra samples; small-N detection | leakage, usual-suspect bias, pooling across unlike backgrounds | small-N panels with a factor-held-out prior | 2 |
| 4 | Hamming kernels as mode products; multi-radius density-ratio seeds | exhaustive all-vs-all at O(d·L·4^L) | adds no information; modest step beyond the radius-1 ball | multi-radius vs radius-1 on degenerate-motif sets | 3 |
| 7 | Exhaustive spacing/orientation tables for motif pairs | grammar as calibrated tests | prior art (SpaMo); self-overlap and repeats confound | known composite elements recovered, nulls clean | 3 |
| 7b | Fitting representations in count space | representation cost independent of N | engineering, not new statistics | equal matrices to sequence-space EM, faster | 4 |

---

## 1. Framing: where extra information can come from

Deep models are often said to benefit from "more sampling of motif structure". On a single dataset, though, a CNN
sees the same peaks HASTUR sees. What it does better is **share parameters**:
- every degenerate near-instance contributes a little to a filter;
- filters are fitted jointly, so motifs explain each other away;
- many instances of one core teach it flank context;
- pretraining or multi-task heads bring in other data.

**Sufficiency limits what within-set comparison can add (Derivation).** THEORY Prop. 1 makes the L-mer lattice
sufficient for the position-pooled window distribution. So any statistic that treats windows symmetrically, by
content alone, is a function of the counts. That includes all-vs-all Hamming comparisons of windows within the set.
Comparing a set with itself gives better *estimators* (weighting, smoothing, low-rank fits) but no new samples.

**The genuinely extra data comes from two places:**
- **The genome**, which informs the **null**: what a word's rate is in a given context, and what copies of a core
  look like when they are not bound.
- **Other experiments**, which inform the **alternative**: which motifs exist, and with what flanks.

The proposals below are organised around that split.

---

## 2. Context-class background with importance-weighted counts (operating mode)

### 2.1 Claim (Derivation)

HASTUR's model is additive: P_fg = (1 − π)Q + πP. If the background varies by position, Q_t, then the likelihood
ratio of one occurrence of word w at t scales with P(w)/q_t(w). Here q_t(w) is the null rate of w at t. An occurrence
where w is locally rare is strong evidence; one where w is locally common (a repeat copy, a w-rich isochore) is weak.

Sites are rare (πP ≪ q), so the null is approximately Poisson. In that regime, the locally most powerful statistic
for one cell is the weighted count

  S(w) = Σ_t a_t (x_t(w) − q_t(w)), with a_t = 1/q_t(w).

Its efficiency relative to the unweighted count (a_t = 1, which is what a global null does) is

  **gain = mean_t(q) · mean_t(1/q) ≥ 1** (Jensen), which equals e^{σ²} when log q_t is normal with spread σ.

| σ | 0.3 | 0.5 | 0.7 | 1.0 |
|---|---|---|---|---|
| gain in effective data | 1.09× | 1.28× | 1.63× | 2.72× |

This is a different gain from enlarging the background. A bigger background only removes background counting noise,
which is about 1/m of the variance at an m× background.

### 2.2 Doubts, and what survived

1. **"A word's genome-wide Hamming neighbourhood is a per-position background."** *Not as stated.* The count of words
   within n mismatches of w depends only on w, not on t. Summed over positions it becomes a global (unmatched)
   background. A background is position-local only when the expectation depends on some **context c_t** of the
   position.
2. **"Neighbours of the flanking context give a genome language model's knowledge."** *Narrowed.* A 30 bp flank has
   genome-wide neighbours within a few mismatches almost only in **repetitive** sequence. In unique sequence it has
   none, and the expectation falls back to the global one. The useful context is therefore two-part:
   - **repeat family / copy number**, which Hamming neighbourhoods measure well;
   - **composition class** (GC/isochore, CpG o/e), a low-dimensional summary, not a neighbourhood.

   Contexts are defined as (repeat class) × (composition class). The "k-NN genome language model" framing is dropped
   (§8).
3. **Heavy weights.** 1/q explodes where q is small. Mitigations: shrink each class table toward the global table,
   cap the weights, and calibrate only by permutation (no analytic tail). The power claim in §2.1 is asymptotic.
   Measure it, don't assume it.
4. **The additive alternative may be wrong.** If sites are placed *because of* context (a GC-loving factor in
   GC-rich classes), pure 1/q weighting hurts. The general optimum is a_c ∝ s_c / q_c, with s_c the site density in
   class c. s_c is unknown, so either cross-fit it, or re-learn it inside every permutation (§2.4).
5. **Positional windows absorb signal.** A ±1 kb local window treats a homotypic site cluster as "w is locally
   common" and down-weights the signal. Context classes must come from **genome-wide** tables, never from a positional
   window around t.
6. **The motif leaks into its own context.** Flanks next to a site predict the site: one half of a palindrome
   predicts the other half, and extended motifs predict their core. Contexts must exclude the motif footprint plus a
   margin, e.g. [t−40, t−10] and [t+L+10, t+L+40]. Composition classes use the wider window minus the footprint.
7. **Inputs.** Class assignment needs genomic coordinates (BED, or windows mapped back). SELEX/RBNS libraries have
   none, and their contexts are homogeneous anyway (σ ≈ 0), so the mode is meant for genomic region sets.

### 2.3 Design

- **Genome index, built once per assembly:** per-class word tables q_c(w) for the gapped k4–5/w13 shapes and
  contiguous L ≤ 12, plus per-position class tracks. The counting is about 4·10¹² increments, minutes on a 64-core
  box. Build leave-one-chromosome-out versions, so fg regions never inform their own null.
- **Counting:** class-stratified fg histograms (one per shape × class), then S(w) = Σ_c a_c(w) · (n_c(w) − E_c(w)).
  The number of increments is unchanged. Memory is about shapes × 4^k × classes, e.g. 715 × 1024 × 64 ≈ 47M counts.
- **Calibration:** the weights come from the genome, not from labels, so any permutation scheme stays valid under its
  own exchangeability assumption. Use sequence-level permutation within sequence-level strata (§3), with hundreds of
  permutations and empirical tails.
- **Two readouts per call:**
  - *enrichment beyond context* (weighted);
  - *enrichment explained by context* (unweighted minus weighted), e.g. "carried by an LTR family".

  Real repeat-borne sites (ERV/LTR-borne binding) are reported, not discarded.

### 2.4 Adaptive weighting inside the permutation loop

The class weights, number of classes, shrinkage and cap are all data-adaptive. Recomputing the whole adaptive
procedure on each permuted labelling, and taking the max-T over its output, keeps family-wise control under
exchangeability. Speed is what makes this affordable, and STREME cannot do it.

**Caveat:** adaptive weights learned on the same fg still overfit in terms of *power* (not validity). Prefer
cross-fitted weights, with the loop as the validity guarantee.

### 2.5 Predictions

- **P2a:** on human ChIP peak sets, gain > 1.3 for the truth motif's core words. Smaller on promoter-only sets, which
  have homogeneous composition.
- **P2b:** at matched false calls per run, weighted detection beats unweighted on human ChIP. No change on SELEX
  (control).

### 2.6 First measurement (before any build)

For dev ChIP sets:
- Assign each fg position a context class from a quick genome table (one chromosome is enough to start).
- For the truth motif words, compute mean(q)·mean(1/q) across peak positions.
- Also compute the efficiency with weights a_c ∝ s_c/q_c, where s_c comes from where the truth's best sites actually
  fall. This shows whether the 1/q direction helps or hurts.

**Rule:** build only if the median predicted gain across factors is ≥ 1.3.

---

## 3. Matched-set backgrounds, stratified permutation and tabled expectations

### 3.1 Claim

Any per-position expectation that averages genome positions sharing a context equals, in expectation, the count in a
background sampled per fg region from genome positions with that context. So a "local background" can be delivered
as **matched sets**: K genome windows per fg region, matched on registered covariates.

Permuting labels *within* each matched set (fg window + its K matches) is exact whenever windows are exchangeable
within matched sets. That is far weaker than global exchangeability. Under global permutation, a non-exchangeable
null (e.g. GC-selected sets) can only be handled by the statistic. Stratified permutation also handles every shift
that the strata capture.

### 3.2 Doubts, and what survived

- **Exactness is conditional.** Confounding left within a stratum passes straight through. Continuous covariates
  need binning, and strata that are too fine leave too few matches. K ≥ 10 is needed for the permutation to resolve
  the threshold. *Survives as stated.*
- **The null defines the question.** An unmatched genome null calls promoter-generic motifs (GC box, ETS, NRF1,
  CCAAT) in any promoter-rich set. Those are real motifs, but not the factor's. For ChIP, the covariates must include
  distance to promoter and accessibility where available. Covariates are **registered, not tuned**.
- **Repeats carry real binding sites.** Collapse rather than discard: with copy-number collapse, occurrences whose
  ±15 bp contexts lie within d mismatches count as one locus, which keeps integer counts and exact tests. Report
  repeat-derived calls as a flagged class.

### 3.3 Two forms

- **Mode A (matched sets):** the background is the K matched windows per region, written out as a FASTA. The
  statistic is unchanged; only the permutation becomes stratified. This is the minimum viable experiment.
- **Mode B (tabled expectations):** expectations come from stratum tables, E(w) = Σ_r n_r(fg) · g_r(w). The background
  is never counted, so runs get faster. Max-T comes from matched genome draws, each costing one fg-sized count.

### 3.4 Smaller hypotheses that survived

- **Presence counts with genome-empirical expectations.** Presence (per-sequence) counting resists clumping but needs
  expectations that reflect clumping. Earlier attempts used Poisson or noisy synthetic expectations. Matched genome
  windows of the same length give empirical presence rates with negligible noise. *Prediction, to test.*
- **Long-word expectations.** hg38 holds about 300 copies of the average 12-mer, so real expectations exist for long
  L-mers. Expect a small effect on detection (the optimal-conjunction law: long exact words rarely match soft sites).
  Any gain is in representation. *Low priority.*

---

## 4. Hamming kernels as mode products (enabling fact) and multi-radius seeds

### 4.1 Fact (Derivation)

All-vs-all comparison of windows with a Hamming kernel is a quadratic form in the counts:

  Σ_{i≠j} K(x_i, x_j) = nᵀKn − Σ_i K(x_i, x_i).

The kernel factorises over positions, K = ⊗_j (I + z(J − I)), with the power of z marking the distance. So K·n is a
sequence of mode products: the same operation as `mode_smooth` and `ball()` in `proto/hlat.py`. All radii up to d cost
O(d·L·4^L), independent of N. **Exhaustive all-vs-all comparison is free given the lattice.**

### 4.2 What survived

- **Multi-radius density-ratio seeds.**
  - Compute fg ball mass against background (or genome-q) ball mass at radii 0…d over all of L-mer space.
  - Seeds are the local maxima. The radius at which a mode's ratio peaks estimates the motif's degeneracy.
  - This is a modest extension of the existing radius-1 ball. The extra radii enlarge the family, so radius choice
    goes inside the permutation loop.
- **Omnibus test, narrowed.** (f − b)ᵀK(f − b) is a kernel two-sample (MMD) statistic.
  - On raw counts it detects composition shift, which HASTUR must not report. It is only usable on double-centred
    residuals, calibrated by permutation.
  - Its only advantage over max-T is power against **diffuse** signals that no cell or ball captures, such as
    low-complexity RNA preferences.
  - *Low priority.*

---

## 5. Flank contrast conditioned on the core, against its genome copies

### 5.1 Claim

For a detected core c (say an 8-mer), the genome holds about 10⁵ copies, mostly unbound. Contrast the flank letters
at each offset across fg occurrences of c with the flanks of c's genome copies. The core is identical by
construction, so only bound-ness differs. Extensions:
- **core variant × flank tables:** ball members × flank letters, giving core–flank dependencies;
- **half-site × spacer × half-site tables** for dimers.

The families are small, so the tests are powerful, exact and calibrated. These dependencies are where deep models are
known to beat PWMs.

### 5.2 Doubts, and what survived

- **Composition confounding (serious).** fg copies sit in peaks (open, GC-richer, promoter-proximal), while genome
  copies sit anywhere. Flank differences may then reflect regions, not specificity. Required:
  - weight the genome copies to the fg's context-class distribution (§2);
  - double-centre each flank column on the copy's own wider-window composition;
  - add a second, within-fg control: flanks of c in high-signal vs low-signal peaks.

  A flank effect counts only if both controls agree.
- **Dilution.** Some fg occurrences of c are not the bound site (a second copy in the window, or indirect binding).
  This makes the test conservative, not invalid.
- **Cofactors.** An adjacent partner motif appears as a "flank preference". That is still a finding (a composite
  element), but it is reported as such and checked against §7.
- **Genome copies bound elsewhere.** Copies bound in other cell types make the test slightly conservative.
- **Prior art.** Flank effects on binding are known from in-vitro and genomic studies. What would be new is a
  calibrated, exhaustive tool component. Claim that, not the phenomenon.
- **Scope.** This is for genomic sets. For SELEX, the analogue is flanks within the library itself, not genome copies.

### 5.3 First measurement

On ChIP dev factors, compare the held-out site log-likelihood (or held-out peak-vs-background discrimination) of the
core PWM with and without the significant flank and core×flank terms.

**Rule:** keep only if the held-out gain beats a matched number of random extra parameters.

---

## 6. Borrowing across experiments

This is the one route that brings **new samples** of motif structure, the true analogue of multi-task learning and
pretrained filters. It is also the one that speed makes practical: every public ChIP, SELEX and RBNS set can be run.

### 6.1 Stratified pooling for shared motifs

Experiments that share a motif (replicates, paralogs, cofactor partners) can pool evidence. Doubt: **raw count
addition is wrong** when backgrounds differ (cell type, peak caller, composition). Use a stratified combination
instead:
- a Mantel–Haenszel-style sum of per-experiment observed-minus-expected values, each experiment against its own
  background;
- calibration by permutation within each experiment.

**Use:** more sites per motif for representation (flanks, §5), and detection of weak shared cofactors. Pooled and
per-experiment calls are reported separately, because co-binding can blur specificity.

### 6.2 Empirical-Bayes family for small N

A motif prior built from an atlas of calls across experiments supports a **weighted family**: more α for motifs
already seen elsewhere, plus a reserved share of α for de novo calls. This is the main lever at small N, where the
limit is the multiplicity of de novo search and a better null cannot help. Doubts and safeguards:
- **Leakage:** the prior must be built with the query's factor (and its replicates) held out.
- **Usual-suspect bias:** CTCF, AP-1 and ETS dominate atlases and would gain α everywhere. Cap each prior weight and
  report which calls used prior weight.
- **Validity:** weights must be fixed before the query is seen. Under weighted Bonferroni or weighted max-T,
  family-wise control holds whatever the prior's quality, so only power depends on it.

### 6.3 Narrowed: the experiment × motif matrix

Factorising z-tensors across experiments (about 10⁶ cells each) would be dominated by composition and pipeline
artefacts. Factorise the **calibrated call matrix** (the D4 atlas) instead, for cofactor modules.

---

## 7. Grammar: exhaustive spacing and orientation tables

For each pair of motifs in the bank, tabulate co-occurrences by spacing and orientation in fg against matched
expectations (§3), with max-T over pairs × spacings × orientations.

- **Prior art:** SpaMo and similar tools. The contribution would be exhaustiveness, calibration and speed, not the
  idea.
- **Confounds found under doubt:**
  - self-overlap: a palindrome or overlapping motifs matched by two motif models at fixed offsets;
  - repeats: an Alu carries several motifs at fixed spacings;
  - peak-width selection, which gives broad co-occurrence enrichment, not spacing-specific enrichment.

  Required: exclude overlapping matches, use repeat classes from §3, and test against a smooth spacing baseline
  rather than a flat one.
- **First measurement:** known composite elements recovered from ChIP dev; zero spacing calls on nulls.

### 7b. Fitting representations in count space

A ZOOPS/TCM-style EM, or a rank-one fit, over the **weighted table of distinct L-mers**, rather than over sequences.
Its cost depends on the number of distinct words, not on N. This matters for genome-scale and atlas-scale use. It is
engineering, not new statistics, and it drops the one-occurrence-per-sequence constraint (TCM-like). Accept it only if
its matrices equal sequence-space EM's on dev.

---

## 8. Dropped or narrowed after doubt

| original claim | verdict | why |
|---|---|---|
| A word's genome-wide n-mismatch neighbourhood count is a per-position background | **reframed** (§2.2.1) | it is a property of the word, not the position; position-locality needs a context |
| Genome-wide k-NN over flanks recovers a genome language model's knowledge | **dropped** | flank neighbours exist almost only in repeats; in unique sequence the context reduces to composition. Kept as repeat class × composition class |
| Pairs from different sequences give a clump-robust statistic | **dropped (wrong)** | for one cell, Σ_{s≠s'} n_s n_s' still scales with the clump size n_s; only presence (n_s capped at 1) is clump-robust |
| Method of moments / tensor decomposition replaces EM, with guarantees | **dropped** | site rates per window are tiny (sample complexity grows polynomially in 1/π); a motif shorter than L forms a sum of shifted rank-one terms in the L-mer tensor; and EM from lattice seeds already reaches oracle-start EM in most conditions, so local optima are not the bottleneck. Narrowed to count-space fitting (§7b) |
| A genome-sized background gives a large precision gain | **narrowed** | background noise is about 1/m of the variance at an m× background; the gain is small except for small backgrounds |
| An MMD omnibus test as a screening statistic | **narrowed** (§4.2) | max-T already gives a calibrated "nothing found"; the omnibus only helps diffuse signals, and only on double-centred residuals |
| All-vs-all within a set gives deep-model-style "extra sampling" | **dropped** | sufficiency (§1): within-set comparisons add estimators, not samples; extra samples come from the genome (null) and other experiments (alternative) |
| Gains quoted from current dev thresholds | **dropped** | dev is ongoing; this doc uses derivations and registered first measurements only |

---

## 9. Suggested order

1. **Measurements that need no build:**
   - §2.6 (context heterogeneity → predicted gain);
   - copy-number profile of the contexts behind null false calls and permutation maxima;
   - fraction of true dev sites in high-copy contexts (what collapsing would cost).
2. **Genome index** (class tables and tracks, leave-one-chromosome-out). It serves §2, §3, §5 and §7.
3. **Mode A** (§3) as a registered round: arms = matched sets, + stratified permutation, + copy-number collapse.
4. **Weighted counting** (§2), if §2.6 clears its rule.
5. **Flank contrast** (§5) and **stratified pooling with a held-out prior** (§6), the two routes to extra samples.
6. **Multi-radius seeds** (§4), **grammar** (§7) and **count-space fitting** (§7b), as time allows.

Each step is either a linear map or a low-rank fit on additive counts, or a draw from the genome index. Every
data-adaptive choice can be recomputed inside the permutation loop, so calibration is preserved by construction.
