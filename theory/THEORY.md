# The lattice and the PWM: a unified theory of motif detection, and the design it implies

> **What this is.** The day-0 mechanism document for a new motif tool, written in the order `docs/DEV_INSIGHTS.md` asks for:
> mechanism first, then predictions, then the plan. It reconciles Comsequent's exact sparse-pattern lattice with the PWM
> tools it lost to (STREME, HOMER, MEME), explains what each family captures and discards, and derives the design of a
> tool that keeps the lattice's wins and closes its losses.
>
> **How claims are labelled.** *Theorem/Proposition*: proved here (proofs in Appendix A, or standard). *Checked*: verified
> numerically by a script in `checks/` (outputs in `results/`, index in Appendix B). *Analogy*: a result from another
> field used to guide design, not proved for motifs. *Prediction*: to be tested; listed in §10.4.
>
> Started 2026-10-10. Sources: `../docs/COMSEQUENT_METHOD.md`, `../docs/DEV_INSIGHTS.md`, `../../comsequent/docs/*`
> (UNIFIED_MOTIF_MODEL, theory addendum, differential brief) and Comsequent's frozen real-data results (read-only).

---

## 0. The answer in one page

1. **One object, three coordinate systems.** Comsequent's counts are the translation-invariant moment tensors of the
   one-hot sequence: the sufficient statistics of every statistic of degree ≤ k. The Möbius transform of their *log*
   gives I (Comsequent's original statistic). The log-linear cell-table transform gives θ (the frozen method's statistic). A third transform,
   the Möbius inversion of the *linear* enrichment ratio, gives **H**: the Hoeffding (ANOVA) decomposition of the
   foreground/background likelihood ratio (§1).
2. **A PWM is a rank-one object in H.** For a PWM motif in a fraction π of windows, H_S(b) = π ∏_{j∈S} δ_j(b_j), with
   δ_j = p_j/q_j − 1, exactly, at every order and for every letter combination (Theorem 2; checked to 1–2%). So the
   whole lattice of a PWM motif is fixed by 3 numbers per position. Comsequent tests the coordinates of this
   tensor one at a time; a PWM tool fits its rank-one factor. **That is the reconciliation:** the two families are a
   coordinate test and a factor fit on the same tensor. Dependencies (epistasis) are what the rank-one factor fails to
   explain (rank-one share 0.998 for a PWM, 0.45 for an epistatic pair).
3. **θ is the right attribution statistic only for crisp motifs.** On a soft (PWM) motif, θ peaks at order 3 and turns
   negative by order 5, while H keeps growing (checked). θ-maximality and the walk's "extend while θ grows" therefore
   stop early on soft motifs. This is the mechanism behind "the walk prunes TATAAT to ANNNTA".
4. **The dilution gap in detection is real but small on real TF motifs.** With oracle knowledge, the best exact 7-position
   cell needs a median **1.54×** the data of the optimal PWM scan, across 472 HOMER vertebrate motifs at a ChIP-like site
   rate (IQR 1.29–1.91; 2.0× at k = 5), and **1.37×** across 1,293 JASPAR matrices. Real motifs are "crisp core plus
   soft flanks" (median 6 of 12 columns with p ≥ 0.8). The gap is large only at k = 5 with sparse sites, and for long
   low-information motifs (§2). The
   **optimal-conjunction law** says where an exact cell sits on the Neyman–Pearson frontier: planted matches outnumber
   background matches by u\* = ln(p²/q)/ln(1/p), and conjunction helps only while p² > q (§2.2).
5. **So the real-data losses were mostly not detection.** In Comsequent's ENCODE ChIP misses the motif was usually
   found but represented imprecisely: REST as a full 20-bp degenerate call at rank 2, ZNF143 the same, MAFK two
   positions short, NRF1 and CREB1 as fragments. Meanwhile a low-information composition call (CG) took rank 1 in 54% of
   runs. The RBNS readout has a chance floor of ~0.7–0.8 for tools that emit long calls (a letter-shuffled truth is
   "found" 69–80% of the time), so that loss is not interpretable as detection (§3).
6. **Composition is a coordinate problem with a clean fix.** A *double-centred* contrast (each set against its own
   order-m Markov expectation, then fg against bg) removes everything an order-m chain can express. ENR makes
   140–5,561 false calls on motif-free pairs with a CpG/GC shift; the double-centred contrast makes 0–0.6. With a 10×
   background it keeps the motif in the top 10 in 100% of runs, at rank 1 in 80–100% (checked, §4).
7. **Position has three lattice-native forms.** Spanning shapes (a shape as long as the motif has one placement, so it
   aligns the motif), anchors (known shift) and position axes (non-uniform shift). This parallels the
   multireference-alignment hierarchy: uniform shifts need third-order invariants, non-uniform shifts second-order,
   known shifts first-order (§5, analogy). Exactly: a position axis is worth 1 + χ² of the site-position distribution,
   i.e. **2.8–5.6× data for ChIP summits with sd 10–20 bp**, more than any letter-model change (Prop. 5).
8. **Estimation, not detection, is what U2 got wrong.** Hard thresholds, cross-fitted halves and exact starts bias and
   waste data. The theory says to detect with the exhaustive calibrated statistic, then estimate with all-data EM from
   a lattice start (Le Cam one-step). Checked against STREME: lattice seeds + EM reach the oracle-start EM in 9 of 11
   conditions and beat STREME's top motif in 67% of datasets (lose 13%). Where they fail, at small N, the k = 5 seed is
   missing, and full-length 8-mer seeds recover most of the gap (0.59 → 0.83, oracle 0.85). STREME's best-of-3 still
   leads in the two hardest cells (§6.3).
9. **Why the field went deep.** Deep sequence models fit exactly the unified model class implicitly (§7): a bank of
   rank-one lattice weightings (filters = PWMs), a mixture nonlinearity, latent alignment (pooling), sparse interactions
   (deeper layers) and quantitative labels. They won on *representation* at large N, not because algorithmic
   discovery is exhausted. They gave up calibration, small-N efficiency, speed and stability. Room remains for an
   algorithm that fits low-rank + sparse structure directly from exact sufficient statistics, calibrated, in
   seconds (§9).
10. **The design.** Count once, exactly, against a large matched background. Rank in composition-free (double-centred)
    coordinates. Detect with three matched detectors: a cell test (crisp), a rank-one/soft test (PWM) and a pooled test
    (diffuse interaction). Detect at full length (contiguous tables to L = 12, spanning shapes) to remove the k-ceiling.
    Represent with all-data EM from lattice starts, chosen by the detection statistic. Test residual dependencies with
    θ. Explain away by deflation. Use quantitative labels when they exist (§10).

---

## 1. One object, three coordinate systems

### 1.1 Notation

- Alphabet Σ = {A, C, G, T}; one-hot X_j(b) = 1[x_j = b] for position j of a window of width w; both strands counted.
- A **pattern** (S, b): offsets S = {o_1 < … < o_k} ⊆ {0, …, w−1} (the *shape*) and letters b ∈ Σ^k.
- **Lattice count** n(S, b) = Σ_t ∏_i X_{t+o_i}(b_i): occurrences of the gapped word at every position t.
- Foreground and background window distributions P_fg and Q. The site model: a fraction π of windows (the
  *per-window site rate*) carry a site drawn from the site distribution P. Typical π: 10⁻² (SELEX), 10⁻³ (ChIP, 200 bp
  peaks), 10⁻⁴ (sparse; Comsequent's N = 500 sims are ~5·10⁻⁴).
- **Marginal enrichment ratio** R_T(b_T) = P_fg(x_T = b_T) / Q(x_T = b_T), estimated by [n_fg(T,b)/N_fg] / [n_bg(T,b)/N_bg]
  or against a model expectation.

### 1.2 The lattice is the moment tensor

**Proposition 1.** The table {n(S, ·) : |S| = k} is the order-k translation-invariant moment tensor of the one-hot
sequence. For |S| ≤ k these are all the marginal tables of the (position-pooled) window distribution up to order k.
Every statistic that is a polynomial of degree ≤ k in the window one-hot is therefore a linear function of lattice
counts. In particular, the lattice holds the sufficient statistics of every hierarchical log-linear model of order ≤ k
(Comsequent "move A").

Two consequences. First, a "degree-k" statistic is the natural notion of what the lattice can see; this connects to the
low-degree theory of detection (§2.6). Second, any per-position re-weighting of letters (IUPAC, mismatches, PWM weights)
is a linear map of the counts (§1.7), so **soft letters never need a recount**.

### 1.3 Three inversions of the same counts

| | definition | scale | what it measures |
|---|---|---|---|
| **I** (Comsequent original) | I(S) = Σ_{T⊆S} (−1)^{\|S∖T\|} log f(T) | log of marginals | interaction information of match indicators |
| **θ** (frozen method) | θ(S) = Σ_{x∈{0,1}^S} (−1)^{\|S\|−\|x\|} log n_S(x) on the match/mismatch cell table | log-linear | Amari θ-coordinate: interaction in the log-linear model |
| **H** (new here) | H(S) = Σ_{T⊆S} (−1)^{\|S∖T\|} R_T | **linear** | the Hoeffding/ANOVA component of the likelihood ratio |

All three are Möbius transforms over the sub-pattern lattice; they differ in the scale on which they invert. Under
the null all three are first-order equivalent. Checked (c1, T3): across 300 null replicates, nested patterns have
correlations within ±0.1 for I, θ and H alike, against 0.09–0.49 for raw counts. Null orthogonality is therefore *not*
what distinguishes H. What distinguishes it is Theorem 2.

### 1.4 H is the orthogonal decomposition of the likelihood ratio

**Theorem 1 (Hoeffding).** Let Q = ⊗_j q_j (product background) and ℓ(x) = P_fg(x)/Q(x). Then ℓ(x) = Σ_{S} H_S(x_S) with
H_∅ = 1, and E_Q[H_S | x_T] = 0 whenever S ⊄ T. Hence the H_S are L²(Q)-orthogonal and
Var_Q(ℓ) = χ²(P_fg‖Q) = Σ_{S≠∅} E_Q[H_S²].
*Proof.* For product Q, E_Q[ℓ | x_T] = P_fg(x_T)/Q(x_T) = R_T(x_T). The Hoeffding (Efron–Stein) components of any
function are the Möbius inversion of its conditional expectations over subsets. ∎

So the lattice, read in H coordinates, decomposes the evidence for "foreground ≠ background" into non-overlapping
pieces, one per pattern. Evidence from different sub-patterns can be *added* without double counting.

### 1.5 A PWM is rank one in H

**Theorem 2.** For windows aligned to the site (§5 treats unknown positions), with a product background Q, if
P_fg = (1 − π)Q + πP with a product (PWM) site distribution P = ⊗_j p_j, then for every S ≠ ∅ and every letter
combination b,

  **H_S(b) = π ∏_{j∈S} δ_j(b_j),  δ_j(b) = p_j(b)/q_j(b) − 1.**

*Proof.* R_T(b_T) = (1 − π) + π ∏_{j∈T}(1 + δ_j(b_j)); the Möbius sum of ∏_{j∈T}(1 + δ_j) over T ⊆ S is ∏_{j∈S} δ_j, and
the constant cancels for S ≠ ∅. ∎

Consequences:
- **Rank one at every order.** The entire lattice of a PWM motif (every shape, every order, every letter combination)
  is fixed by 3 free numbers per position. Checked (c1, T2): H matches π∏δ with median relative error 0.6% (crisp)
  and 1.8% (PWM). The order-2 (position × letter) H table has a rank-one share of **0.998** for a PWM, against 0.45 for
  an epistatic pair and 0.23 for null noise (T4).
- **Exact attribution.** If p_j = q_j at a position (off-motif), δ_j ≡ 0 and H vanishes for *every* pattern containing
  that position (checked: mean |H| 0.018 with an irrelevant position vs 0.58–3.1 without).
- **Several motifs** at fixed alignment: H = Σ_m π_m ⊗_j δ^{(m)}_j, a CP tensor decomposition with one term per motif.
- **Dependencies are the residual.** A non-product site distribution leaves structure the best rank-one fit cannot
  explain. That residual is exactly what Comsequent's interaction statistics are for.
- **Crisp limit.** δ_j(c_j) = 1/q_j − 1 (= 3 for uniform q): H grows ×3 per position, like θ, and the full motif has
  the largest |H|. Crisp motifs are rank one as well; the PWM is the general case and the exact pattern its corner.

**This is the reconciliation in one sentence:** the lattice tensor of a motif is rank one in linear (Hoeffding)
coordinates. Comsequent tests its coordinates one at a time; PWM tools fit its rank-one factor by EM or discriminative
search; the unified tool should do both (§7).

### 1.6 θ in the weak and strong regimes, and why θ-maximality fails on soft motifs

**Proposition 2.** Expanding log((1 − π)bg + π·site) to first order in π·site/bg gives, in the weak regime,
θ_S ≈ π ∏_{j∈S} [p_j(c_j)/q_j(c_j) − (1 − p_j(c_j))/(1 − q_j(c_j))]. This is product form, like H. In the strong regime
(sites dominate the all-match cell) the logarithm saturates. For crisp motifs θ_S = log(1 + π/((1 − π)q^k)) is still
increasing (Comsequent's closed form). For soft motifs the saturation is uneven across cells and high-order θ becomes
non-monotone.

Checked (c1): for a 6-position PWM (π = 0.05), the mean θ by order 2/3/4/5 is 0.31 / 0.48 / 0.35 / **−0.24**, while H is
0.24 / 0.47 / 0.88 / 1.48 (predicted 0.24 / 0.48 / 0.88 / 1.47). The θ product test fails (median ratio 0.27–0.36 vs 1.0).

**Consequence.** θ-maximality ("keep S if no sub-pattern has larger θ and no extension raises it") is a sound attribution
rule for crisp motifs and an unsound one for soft motifs: it stops at order ~3 and prunes informative soft positions.
This is the documented failure where "the walk prunes even a TATAAT seed to ANNNTA" on promoters. For soft motifs,
attribution should come from H (which grows with every informative position and vanishes with every uninformative one)
or from the rank-one fit.

### 1.7 Soft letters are linear maps of the lattice: the soft lattice transform

For a shape table n_S (4^k cells) and a 4×4 kernel K_j per position, define the mode product

  **S_S(c) = Σ_x n_S(x) ∏_j K_j(x_j, c_j) = (n_S ×_1 K_1 ×_2 ⋯ ×_k K_k)(c).**

It scores *every* consensus c at once. Cost k·4^{k+1} per shape (≈ 0.8 Gflop for a full 12-mer table). Variances
follow from the elementwise-squared kernels on the same counts. Special cases:

| kernel K(x, c) | what S computes | existing name |
|---|---|---|
| 1[x = c] | exact count | the lattice |
| 1[x ∈ B(c)] | IUPAC letter sets | degenerate lattice (U3) |
| 1 + t(1[x ≠ c]) | generating function of Hamming distance (coefficient of t^d) | mismatch ladder (U8), Weeder/HOMER oligos |
| p_c(x)/q(x) for a column "dominant c" | LR of a homogeneous PWM, every consensus at once | soft lattice (§6, check 4) |
| 1[x = c]/q(c) − 1 | H coordinates (elementary symmetric functions of centred matches) | Hoeffding lattice |

So the "positional letter weights" a PWM uses are a kernel applied to counts the lattice already holds. The
soft/exact distinction is a choice of K, not of data structure.

---

## 2. Detection: what fraction of the evidence does a lattice cell carry?

### 2.1 Every window statistic is a set; the Neyman–Pearson frontier

A count statistic is defined by a window set A (or a weighting). With site mass P(A) = P(x ∈ A) and background mass
Q(A), a count contrast over N windows has (Poisson, known background)

  **z²/N ≈ π² P(A)² / (π P(A) + Q(A)).**

This is the dense regime z² ≈ Nπ²P²/Q when background matches dominate, and the sparse regime z² ≈ NπP when
planted matches dominate. By the Neyman–Pearson lemma, at each Q(A) the largest P(A) is attained by a superlevel set
of the likelihood ratio p(x)/q(x), i.e. by the PWM score threshold with the true PWM. Every exact cell, IUPAC
pattern or truncated PWM is a point under this frontier. Its **data multiple** is r_NP/r_A, where r = max over the
class of the right-hand side. This generalises U6 (which compared two oracle detectors) to nested classes.

### 2.2 The optimal-conjunction law

**Proposition 3.** For a homogeneous motif (consensus probability p per position, background q), an exact conjunction of
k positions has r(k) = π² p^{2k}/(π p^k + q^k). An interior optimum exists iff **p² > q**. At the optimum, planted
matches outnumber background matches by

  **u\* = π P/Q = ln(p²/q) / ln(1/p),  so  k\* = [ln(1/π) + ln u\*] / ln(p/q).**

*Proof.* Set d ln r/dk = 0 with u = πP/Q: 2 ln p = (u ln p + ln q)/(u + 1) ⇒ u = ln(p²/q)/ln(1/p). ∎

Checked numerically (q = 0.25):

| p (consensus) | u\* | best k, π = 10⁻² | best k, π = 10⁻³ | best k, π = 10⁻⁴ | sites captured at best k, π = 10⁻³ |
|---|---|---|---|---|---|
| 0.95 (high IC) | 25 | 6 | 8 | 9 | 0.66 |
| 0.86 (Comsequent "mid") | 7.2 | 5 | 7 | 9 | 0.35 |
| 0.72 ("low") | 2.2 | 5 | 7 | 9 | 0.10 |
| 0.60 | 0.7 | 5 | 7 | 10 | 0.03 |

Reading: a conjunction buys specificity (Q falls by q per position) at a geometric price in sensitivity (P falls by p).
The best cell buys just enough specificity for planted matches to outnumber chance matches u\*-fold, then stops. Only
p^{k\*} of the sites are inside it. For crisp motifs that is most sites; for soft motifs it is a few percent, and
the rest of the evidence is left on the table.

**Per-position criterion.** In the dense regime, adding position j multiplies z by p_j(c)/√q_j(c). It helps iff
p > √q (0.5 for uniform q, which for "dominant + equal rest" columns means IC > 0.21 bits). A two-letter column
(p ≈ 0.45 + 0.45) loses as a single letter (0.9) and wins as its IUPAC letter (0.9/√0.5 = 1.27). In H coordinates the
factor is (p − q)/√(q(1 − q)) (> 1 iff p > 0.683): the threshold for a pure interaction to strengthen with order.

### 2.3 Decomposing the gap: letter, weighting, truncation

The classes EXACT_k ⊂ IUPAC_k ⊂ SOFT_k ⊂ SOFT_L (NP frontier) are nested. Here SOFT_k is the LR threshold on the k most
informative positions, i.e. the rank-one weighting of one shape table. For each motif the three factors multiply
exactly to the total gap (their medians across motifs need not):

- **letter** = r_IUPAC_k / r_EXACT_k: one letter per position vs letter sets;
- **weighting** = r_SOFT_k / r_IUPAC_k: hard sets vs LR weights;
- **truncation** = r_SOFT_L / r_SOFT_k: k positions vs the whole motif.

### 2.4 Real motifs: the gap is 1.3–2×, mostly truncation and letters

Checked (c2): 472 HOMER vertebrate motifs (median length 12, median 12.6 bits), human-like background, aligned, oracle
weights. The table gives the median data multiple vs the oracle PWM scan, with the share of motifs above 2×:

| per-window site rate π | EXACT_5 | EXACT_7 | IUPAC_7 | SOFT_7 | k = 7 split: letter / weighting / truncation |
|---|---|---|---|---|---|
| 10⁻² (SELEX-like) | 1.29 (8%) | 1.28 (8%) | 1.11 | 1.05 | 1.11 / 1.05 / 1.05 |
| 10⁻³ (ChIP-like) | 2.00 (49%) | **1.54 (21%)** | 1.26 | 1.20 | 1.16 / 1.03 / 1.20 |
| 10⁻⁴ (sparse) | 5.58 (98%) | 2.01 (51%) | 1.69 | 1.66 | 1.06 / 1.01 / 1.66 |

- **The motifs Comsequent missed on ChIP are not hard for exact detection**: REST 1.2×, SP1 1.2×, ZNF143 1.4×,
  NRF1 1.3× (EXACT_7, π = 10⁻³). p53 (2.9×), Sox21 (2.6×) and ATF1 (2.1×) are among the hardest of the named set.
- **Crisp references** come out at 1.00 (GATC, GANTC, EcoKI, GCNNNANNNNGC). CCWGG comes out at 1.07–1.66, because its
  W is a two-letter column; IUPAC removes that.
- The gap shrinks with mean column χ² (Spearman −0.66) and grows with the share of soft columns (+0.37).
- **Why it is small: real motifs are crisp cores with soft flanks.** The median motif has 6 of 12 columns with
  p_max ≥ 0.8, and 69% of motifs have ≥ 5 such columns. 34% of all columns have p_max ≥ 0.9; 35% have p_max < 0.6,
  of which 30% are two-letter. An exact conjunction on the core (§2.2) captures most of the detection evidence. What
  it misses is mainly the *representation* of the flanks.
- **Where the gap is large:** k = 5 at sparse rates (5.6×, nearly all of it truncation); Comsequent's simulated soft
  classes (every column p = 0.72 or 0.86, no crisp core); long low-IC motifs (U6: up to 20×).
- **Robust to the matrix source.** The same analysis on 1,293 JASPAR vertebrate matrices (count-based, many from
  HT-SELEX; `C2_LIB=jaspar`) gives a *smaller* gap: EXACT_7 median **1.37×** at π = 10⁻³ (IQR 1.19–1.69; > 2× for 13%)
  and 1.69× at 10⁻⁴; EXACT_5 1.94× and 5.86×. Split at k = 7, π = 10⁻³: letter 1.08, weighting 1.01, truncation 1.19.
- **Caveats.** Any fitted matrix may be sharper than the true in-vivo site distribution (matrices are learned from the
  best sites). The analysis is aligned (no search cost), with a known background (no control noise) and one strand.
  Day-1 task: refit matrices on all peaks rather than top peaks, and repeat.

### 2.5 Control noise: equal-size backgrounds cost 2× data

For a contrast with N_bg background windows, the variance of the count difference carries a factor (1 + N_fg/N_bg).
An equal-size control therefore needs 2× the data of a known expectation; a 10× control needs 1.1×. Comsequent's
simulation benchmarks gave every tool an equal-size control (`shared_sim.make_dataset`). Since counting is cheap, **the new tool should default to a
large matched background**. Check 5 shows where this matters most: words depleted in the background (CpG words in a
CpG-depleted genome), whose background count sets the variance.

### 2.6 The local view: efficiency is a correlation, and the lattice is the low-degree space

In the small-π limit, the noncentrality of any window statistic g is Nπ² Cov_Q(g, Λ)²/Var_Q(g), with Λ = P/Q.
By Cauchy–Schwarz:

  **efficiency(g) = Corr_Q(g, Λ)²,**

the squared background correlation between the statistic and the site likelihood ratio. For a product site model
χ²(P‖Q) = ∏_j(1 + χ²_j) − 1 = Σ_{m≥1} e_m(χ²_1, …, χ²_L), where e_m is the elementary symmetric polynomial. The
**order-m share of the local detection information is e_m / (∏(1 + χ²_j) − 1)**, spread over C(L, m) position subsets.
A degree-≤k statistic captures at most Σ_{m≤k} e_m. This is the *low-degree likelihood-ratio norm* of Hopkins and of
Kunisky–Wein–Bandeira; the lattice of order k is exactly the degree-k space of that theory.

- In sliding windows, order 1 is composition (no positional information), so motif information starts at order 2.
- The order-m shares peak near m ≈ Lχ²/(1 + χ²). **The per-position constant is χ²_j = 1**: for a "dominant letter"
  column that is p_c ≈ 0.68, the same threshold at which an H coordinate starts to strengthen with order (§2.2).
- **Very soft columns** (χ²_j < 1) put their mass at moderate orders. For χ²_j = 0.5 and L = 12, orders 2–4 over all 12
  positions carry 58% of the information, while all orders of the best 7-position subset carry 12.5%. **Such motifs
  want a wide, shallow lattice** (orders 2–4 over spans of 20–30).
- **Columns with χ²_j > 1** (most real columns: p_max ≥ 0.8 gives χ²_j ≥ 1.5) put their mass at high orders. With
  χ²_j = 1.18 and L = 12, orders 2–4 carry 12% and the best 7-subset 2%; with χ²_j = 1.98 and L = 8, the best 7-subset
  carries 34%. **These want full-length tables**, i.e. all orders over all positions (L-mer tables, spanning shapes).
  Crisp motifs are the extreme: a narrow, deep lattice (Comsequent's k5–7/w13–15) holds all their information.
- Caveat: the local view ignores the heavy tail of Λ, so in the sparse regime §2.1 governs. Both views agree on the
  design: cover the motif's full length, and pool letters where columns are soft.

### 2.7 Verdict on the dilution argument

Conjunctive-exact vs additive-soft is the right mechanism, but on real TF motifs at realistic site rates it costs
detection only 1.3–2× at k = 7. That cannot explain losses like 0.67 vs 0.92 (ChIP) or 0.35 vs 0.72 (SELEX). The next
section shows what did.

---

## 3. Anatomy of the real-data losses

Checked (c3): every top-10 list in Comsequent's frozen P5 ChIP (24 TFs), P4a HT-SELEX (503 libraries) and P3 RBNS
(104 units) results, for Comsequent's all-in-one, STREME and HOMER. Each unit is classified by its best call:
**strict** (the phase's own rule: ≥ 80% of the truth's specified positions, ≤ 1 conflict, call letters ⊆ truth
letters); **precision** (the same, but accepting degenerate call letters that *contain* the truth letter); **fragment**
(≥ 3–4 specified positions inside the truth, no conflict); or **none**. The chance column gives the same classification
against letter-shuffled truths, 5 per unit.

| phase | tool | strict | precision | fragment | none | chance strict | chance precision | chance fragment | rank-1 call < 6 bits |
|---|---|---|---|---|---|---|---|---|---|
| ChIP | Comsequent | 0.67 | **0.12** | 0.08 | 0.12 | 0.00 | 0.03 | 0.12 | **0.54** |
| ChIP | STREME | 0.92 | 0 | 0.04 | 0.04 | 0.06 | 0 | 0.02 | 0 |
| ChIP | HOMER | 0.58 | 0 | 0.33 | 0.08 | 0.03 | 0 | 0.05 | 0 |
| SELEX | Comsequent | 0.37 | 0.08 | 0.39 | 0.17 | 0.03 | 0.03 | 0.32 | 0.05 |
| SELEX | STREME | 0.70 | 0.05 | 0.14 | 0.11 | **0.13** | 0.03 | 0.02 | 0 |
| SELEX | HOMER | 0.71 | 0.05 | 0.17 | 0.07 | **0.14** | 0.03 | 0.03 | 0 |
| RBNS | Comsequent | 0.60 | 0.25 | 0.11 | 0.05 | 0.30 | 0.28 | 0.18 | 0.21 |
| RBNS | STREME | 0.89 | 0 | 0.01 | 0.10 | **0.80** | 0 | 0.01 | 0.06 |
| RBNS | HOMER | 0.90 | 0 | 0 | 0.10 | **0.69** | 0 | 0 | 0 |

What it shows:

- **ChIP: representation and ranking, not detection.** Of Comsequent's 8 misses:
  - *Found, but imprecise* (3): REST at rank 2 as TCAGSACCWYGGACAGNKCC, the full 20-bp reverse complement of the truth. It
    failed on letters S/W/Y/K, where the truth is exact, and one N. STREME and HOMER returned 10–13-bp fragments and
    *also* missed REST. ZNF143 was found at rank 2 as a full-length degenerate call; MAFK at rank 1 was two positions short.
  - *Fragments* (2): NRF1 (GCGCA, CNNNNGCGC, GNNNNNGCGC) and CREB1 (TNACGT, ACGT, CGTCA).
  - *Swamped by composition* (1): SP1.
  - *Found a related factor's motif* (2): ATF3 → AP-1 (TGASTCA), RFX5 → NFY (CCAAT); STREME found the canonical one at rank 2.
  - And in 54% of runs Comsequent's rank-1 call is a low-information composition pattern (CG); STREME/HOMER never do this,
    partly because STREME's minimum width is 8.
- **SELEX: a real gap, smaller than reported.** STREME/HOMER's strict hits include a 13–14% chance floor (long calls
  cover many offsets of a 6-mer truth). Net of chance: STREME ≈ 0.57, Comsequent ≈ 0.34 + 0.05 precision. Comsequent's
  fragment rate (0.39) is close to its chance rate (0.32), so most SELEX "fragments" are not evidence of detection.
  The documented cause is composition from cycle drift (GC calls). This is §4's problem, not a dilution problem.
- **RBNS: the readout cannot separate the tools.** With a 5-mer truth and a 4-of-5 rule, STREME's top 10 "hits" a
  letter-shuffled 5-mer 80% of the time (HOMER 69%), against 89–90% for the real truth. The reported 0.79 vs 0.89 is
  therefore not interpretable as a detection difference. (Caveat: the rule here is P5's, applied uniformly. It
  reproduces STREME/HOMER's reported RBNS rates, 0.89/0.90 vs 0.89/0.87, but gives Comsequent 0.60 vs its reported 0.79,
  so the phase rule differs in details. The chance floor for long calls does not depend on them.)

**Conclusion.** On PWM ground, the deficit decomposes into: (i) **representation** (degenerate letters, short calls,
fragments: what a soft motif model fixes), (ii) **ranking** (composition at rank 1; fragments occupying slots: what
composition-free coordinates and informativeness fix), (iii) **evaluation** (rules with high chance floors for long
calls), and (iv) a residual **detection** gap that the theory puts at 1.3–2× for real motifs and larger only in sparse,
soft, or k = 5 regimes. The new tool must fix (i) and (ii); the evaluator must fix (iii) before any comparison (§11).

---

## 4. Composition, drift and the ranking problem

### 4.1 Where composition lives

By Amari's hierarchical decomposition, KL(P_fg‖Q) splits into orthogonal per-order contributions. Composition and
cycle drift live in the low contiguous orders: base composition (order 1), dinucleotides including CpG (order 2),
codon periodicity (contiguous order 3 with phase). A motif is coherent structure at higher orders: rank one in H
across many orders. ENR ranks by total enrichment, so when drift is large, composition words outrank the motif.
θ-contrast cancels order-1 differences but not dinucleotide ones. The factorial contrast cancels a confounder only
when a second contrast shares it.

### 4.2 The double-centred Markov contrast

**Proposition 4.** Let E^{(m)}_f(w) be the expected count of word w under the order-m Markov chain fitted to the
foreground's *own* (m+1)-mer counts, and E^{(m)}_b(w) likewise for the background. Define

  **D(w) = log(O_f(w)/E^{(m)}_f(w)) − log(O_b(w)/E^{(m)}_b(w)),  Var D ≈ 1/O_f + 1/O_b.**

If neither set has a motif, E[D(w)] ≈ 0 for every word, however different the two chains are: composition cancels
*within* each set before the sets are compared. A motif adds excess beyond its own (m+1)-mer margins in the
foreground; its contribution to those margins is small (a rough estimate for check 5's setting: ~2% of one
trinucleotide's count at π = 0.1). The
cost is variance, dominated by the smaller count. Words depleted in the background (CpG words in a CpG-depleted
genome) need a larger background (§2.5). The same construction with a phase-aware chain removes codon structure.

### 4.3 Checked

c5: CpG-island-like foreground (GC ≈ 0.6, CpG not depleted) vs genome-like background (GC ≈ 0.41, CpG depleted),
planted PWM (L = 8, mid IC), 5 replicates, all 8-mers:

| N, π, background | statistic | truth top 1 | truth top 10 | top-10 composition share | false calls (z > 5) on a motif-free pair |
|---|---|---|---|---|---|
| 2,000, 0.1, 1× | ENR | 1.0 | 1.0 | 0.90 | **149** |
| | double-centred, m = 1–3 | 0.4 | 0.6 | 0.00 | 0 |
| 2,000, 0.1, **10×** | ENR | 0.8 | 1.0 | 0.88 | 146 |
| | double-centred, m = 1–3 | **1.0** | **1.0** | 0.00 | 0.6 |
| 2,000, 0.3, 1× | ENR | 1.0 | 1.0 | 0.46 | 151 |
| | double-centred | 0.2–0.4 | 0.8 | 0.00 | 0 |
| 2,000, 0.3, 10× | ENR | 1.0 | 1.0 | 0.40 | 140 |
| | double-centred | 0.8 | 1.0 | 0.00 | 0.4–0.6 |
| 10,000, 0.05, 1× | ENR | 0.0 | 0.2 (median rank 134) | 0.98 | **5,561** |
| | double-centred | 0.8 | 1.0 | 0.00–0.02 | 0.2 |

The two replicates that failed at 1× had a CpG inside the motif (CTGCATCG, GTCGCAAA): their background counts were tiny,
so the variance swamped the z. At 10× background both reach rank 1. ENR's flood grows with N (149 → 5,561), because
drift z grows like √N. This is the "high drift" loss mechanism of Comsequent's signal/drift pre-check, now with a
statistic that removes it.

---

## 5. Position: alignment, anchors and the multireference-alignment view

### 5.1 Lattice counts are translation-invariant autocorrelations

n(S, b) = Σ_t ∏_i X_{t+o_i}(b_i) is the order-|S| autocorrelation of the one-hot signal at lags S. For a motif at
uniformly random positions, the centred order-m excess is π Σ_u ∏_i δ_{u+o_i}(b_i), summed over placements u. That is
the order-m autocorrelation of δ: a *sum of shifted rank-one terms*, not a single rank-one term. Translation pooling is
what makes unaligned motif discovery harder than fitting Theorem 2 directly.

### 5.2 Spanning shapes are alignments

A shape whose span equals the motif's informative extent has exactly one placement, so Theorem 2's rank-one structure
is restored for it. **Exhaustive shape enumeration contains self-aligning coordinates: the span carries the positional
information.** This is why Comsequent's high-order gapped cells win on bipartite motifs (one spanning shape aligns both
halves at once, with no seed and no local optimum), and why a full-length contiguous table (L-mer, §6) aligns a
contiguous PWM.

Checked (c7; unaligned 100-bp sequences, a 6-position PWM, pair tables at every spacing d):
- **Prediction matches everywhere.** The observed H_d correlates with π Σ_u δ_u δ_{u+d}ᵀ at 0.98–0.998 for every
  spacing inside the motif.
- **One strand of sites.** The spanning spacing (d = 5) is rank one (top singular share **0.982**). Shorter spacings,
  which hold several placements, have shares of 0.52–0.81.
- **Both strands.** The spanning table is rank two (top-2 share **0.997**): forward plus reverse complement (§5.4).
- **Beyond the span** (d ≥ 6), |H| falls to the noise floor (0.005–0.008, against 0.028–0.083 inside).

### 5.3 The moment hierarchy (analogy: multireference alignment)

In multireference alignment, a signal is observed under random shifts at noise level σ. With uniform shifts,
third-order invariants (the bispectrum) are needed, and sample complexity grows like σ⁶ at low SNR (Perry et al.
2019). With a non-uniform (aperiodic) shift distribution, second moments suffice (σ⁴; Abbe et al. 2019). With known
shifts, first moments suffice (σ²). The motif analogue:

| shift information | lattice form | lowest useful order | applications |
|---|---|---|---|
| known shift | **anchored windows** (modified base, edit site, core occurrence) | 1 (frame histograms: U1) | methylation, base editing, PAM (Comsequent's winning mode) |
| non-uniform shift | **position axis** (summit-relative bins; U12b) | 2 (centrality as a positional order-1/2 moment) | ChIP/ATAC summits, CentriMo-like evidence |
| uniform shift | **spanning shapes / full-length tables** | ≥ 3 coherent orders | SELEX random regions, RBNS |

This is an analogy: MRA theory concerns continuous signals in Gaussian noise, and motif discovery is a sparse
categorical mixture. The ChIP row, though, can be made exact:

**Proposition 5 (positional information is worth 1 + χ² of the site positions).** In the dense regime, weight each
window by w(r), where r is its position relative to the summit, f(r) is the site-position density and g(r) the window
density (uniform). The noncentrality becomes Nπ²(P²/Q)·(∫wf)²/∫w²g. By Cauchy–Schwarz this is at most
(1 + χ²(f‖g)) times the unweighted value, with equality at w = f/g. *Positional information multiplies the effective
data like one more motif column.*

For Gaussian summit errors in 200-bp peaks, the optimal weighting gives **5.6× / 3.8× / 2.8× / 1.9× / 1.2×** the data at
sd 10 / 15 / 20 / 30 / 50 bp; the best hard central window gives 5.0 / 3.4 / 2.5 / 1.7 / 1.0×. At realistic summit
precision this exceeds every letter-model refinement in §2.4 (1.2–1.5×). **For ChIP, a summit-relative position axis
is likely the largest single lever**, and it is cheap: one extra count axis (U12b), or position-binned tables.

### 5.4 Reverse complement

The symmetry group is dihedral (shift plus reverse-complement). Counting both strands symmetrises. A non-palindromic
motif is rank two in H (δ and its reverse complement); a palindrome is rank one. Dihedral MRA results
(reflection-invariant bispectrum; Bendory et al. 2024) carry over in the same analogical sense.

### 5.5 What "the positional information a PWM uses" is, and the lattice equivalents

| PWM ingredient | what it is | lattice-native equivalent |
|---|---|---|
| per-position letter weights | soft letters, all 4 per column | kernels on shape tables (§1.7); rank-one H factor (§1.5) |
| a contiguous frame of all L positions | evidence summed over the whole site | full-length tables (L ≤ 12), spanning shapes, extension beyond k by anchoring |
| latent site positions | alignment of each site (EM posteriors) | spanning shapes (implicit), anchors (explicit), position axes |

A PWM tool gets the third by EM, with local optima and seed dependence. The lattice gets it by exhaustive
enumeration, without them. The new tool should take letters from the PWM side and alignment from the lattice side.

---

## 6. Estimation: from a detected cell to a represented motif

### 6.1 Why Comsequent's refinement (U2) failed in small samples

U2 refined a seed PWM by thresholding scores, re-estimating columns from fg hits minus scaled bg hits, and
cross-fitting on halves. The theory finds four defects:
1. **Truncation bias.** Columns estimated from windows selected by a hard score threshold are biased toward the seed's
   letters ("locks onto seed letters": 0.99-pure columns from 73 of 375 promoters).
2. **Halved data.** Cross-fitting doubles the variance of the estimate; with ~50 sites (N = 500, π = 0.1) that leaves 25.
3. **Exact starts.** The initial threshold admits only near-exact matches, so step 1 is self-reinforcing (a softened
   start raised sites from 73 to 139).
4. **Seeds from a k = 5 lattice.** At small N the motif-consistent k = 5 cell is often outside the top 10, while a
   full-length 8-mer table ranks it first (check 4, §6.3).

### 6.2 What the theory prescribes

- **EM on all data is the maximum-likelihood estimator** of the ZOOPS mixture (MEME's model) and is asymptotically
  efficient. Soft posteriors avoid truncation bias.
- **Le Cam one-step.** From any √n-consistent start, one Newton/EM step gives an efficient estimator. With a good
  start, a few EM iterations reach the MLE basin; this is the "spectral initialisation + EM" pattern of the
  tensor-method literature (Anandkumar et al. 2014). The lattice supplies the start. Anchored frame counts around seed
  occurrences (U1) are consistent for the flank columns. Full-length or spanning-shape tables give the core columns
  (aligned by §5.2).
- **Decouple detection from estimation.** A motif's existence and rank come from the exhaustive, calibrated detection
  statistic, which has a valid p-value by design. EM supplies only the representation, using all the data.
  Choosing among candidate EM fits by their in-sample likelihood brings back a winner's curse; choose by the detection
  statistic. If a test of the *refined* representation is needed, use a held-out split for that test only.

### 6.3 Checked

c4: iid background, 100-bp sequences, equal-size control, 5 replicates. Similarity is the mean over truth columns of
1 − total-variation distance, at the best offset and strand (null similarity of a flat matrix in brackets). Ranks are
of the first motif-consistent call.

| motif | N | π | reps | k5 seed in top 10 | soft 8-mer seed in top 10 | LAT+EM (k5 seeds, by gain) | LAT8+EM (top 8-mer seed) | ORACLE+EM | STREME top 1 | STREME best of 3 | flat |
|---|---|---|---|---|---|---|---|---|---|---|---|
| L8 mid | 500 | 0.1 | 5 | 0.4 | 1.0 | 0.59 | **0.83** | 0.85 | 0.77 | 0.89 | 0.39 |
| L8 mid | 500 | 0.3 | 5 | 1.0 | 1.0 | 0.93 | 0.93 | 0.93 | 0.90 | 0.90 | 0.39 |
| L8 mid | 2,000 | 0.05 | 5 | 0.4 | 1.0 | 0.54 | **0.83** | 0.86 | 0.77 | 0.93 | 0.39 |
| L8 mid | 2,000 | 0.1 | 4 | 1.0 | 1.0 | 0.91 | 0.91 | 0.92 | 0.84 | 0.92 | 0.39 |
| L12 low | 2,000 | 0.1 | 5 | 1.0 | 1.0 | **0.90** | 0.80 | 0.91 | 0.75 | 0.83 | 0.53 |
| L12 low | 2,000 | 0.3 | 5 | 1.0 | 1.0 | 0.96 | 0.96 | 0.96 | 0.89 | 0.89 | 0.53 |
| NRF1 (L 12) | 2,000 | 0.1 | 5 | 1.0 | 1.0 | 0.97 | 0.93 | 0.96 | 0.95 | 0.95 | 0.38 |
| REST (L 20) | 2,000 | 0.1 | 5 | 1.0 | 1.0 | **0.98** | 0.77 | 0.98 | 0.85 | 0.85 | 0.41 |
| CTCF (L 20) | 2,000 | 0.1 | 5 | 1.0 | 1.0 | 0.96 | 0.89 | 0.96 | 0.89 | 0.89 | 0.52 |
| Sox2 (L 12) | 2,000 | 0.1 | 5 | 1.0 | n/a | 0.89 | 0.83 | 0.89 | 0.84 | 0.85 | 0.57 |
| p53 (L 20) | 2,000 | 0.1 | 5 | 1.0 | 1.0 | 0.91 | 0.79 | 0.91 | 0.88 | 0.88 | 0.64 |

(One replicate of L8 mid / 2,000 / 0.1 was lost to a STREME output-decoding crash, since fixed. "n/a": Sox2 has
columns with no letter at p ≥ 0.4, so no 8-mer can count as motif-consistent.)

Reading:
- **All-data EM from lattice seeds reaches the oracle wherever the seed pool contains the motif.** With k5 seeds
  (chosen by in-sample gain), the result is within 0.01 of EM started from the true PWM in 9 of 11 conditions. Against
  STREME's top motif it wins in 67% of datasets and loses in 13% (mean +0.02). This is the representation step that
  U2 lacked: REST at 0.98 vs STREME's 0.85.
- **Where it fails, the k-ceiling is the cause.** In the two small-sample cells (N = 500 / π 0.1 and N = 2,000 / π
  0.05) the k5 seed is in the top 10 in only 40% of replicates, while the soft full-length 8-mer seed is in the top 10
  in 100%. Starting EM from the top 8-mer seed lifts these cells from 0.59 / 0.54 to **0.83 / 0.83** (oracle 0.85 /
  0.86; STREME top 1 0.77 / 0.77).
- **Neither seed pool wins everywhere.** 8-mer seeds give a frame too narrow for 20-bp motifs (REST 0.77 vs 0.98), and
  in-sample gain is the better chooser at N = 2,000 but the worse one at N = 500. The design therefore takes seeds from
  both pools, keeps only detection-significant ones, and chooses among them by *held-out* likelihood (§6.2, last
  point). That rule is not yet tested.
- **STREME's best-of-3 is still ahead in the two hardest cells** (0.89 / 0.93), above even the oracle-started EM (0.85
  / 0.86). So this bare ZOOPS EM is less well regularised at small N than STREME's estimator: headroom for the
  estimation layer (MEME-style Dirichlet priors, tuned pseudocounts).
- **Cost.** In numpy, EM from 10 seeds takes 60–70 s per dataset, against STREME's 21 s and the lattice's 2.4 s. The
  Go port and seed de-duplication are needed for the speed claim of §10.2.

---

## 7. The unified model: low-rank plus sparse on the H-lattice

### 7.1 Statement

  **ℓ(x) − 1 = Σ_m π_m [∏_j (1 + δ_j^{(m)}(x_j)) − 1]  +  Σ_{S∈𝒮} H^{sparse}_S(x_S)**

The first sum is the low-rank part: one rank-one term per motif and strand, at all its placements, covering all orders.
The second is the sparse part: a few interaction coordinates that no product term explains (epistasis, pure
higher-order structure).

| family | what it fits | how |
|---|---|---|
| Comsequent | sparse coordinates (cells), one at a time | exhaustive enumeration, calibrated tests |
| PWM tools (MEME, STREME, HOMER) | one rank-one term at a time, greedy masking | EM / discriminative search from seeds |
| FastMotif (Colombo & Vlassis 2015) | several rank-one terms, ungapped, single set | method of moments, spectral decomposition |
| dinucleotide / TFFM / BaMM / MaxEnt | rank one + fixed low-order corrections | EM / IPF |
| ProBound / BEEML | rank one + a link function (biophysics) | likelihood over all reads |
| gkm-SVM | a full-rank weighting of gapped k-mer coordinates (mismatch kernel = a soft lattice) | convex (SVM) |
| CNNs (BPNet, ChromBPNet, …) | a bank of rank-one terms (filters) + nonlinear sparse composition + position | SGD end-to-end |
| **the new tool** | low-rank + sparse, from exact counts | exhaustive support, EM representation, calibrated tests |

**What each PWM tool optimises, in these terms.**
- **MEME** maximises the ZOOPS/OOPS mixture likelihood over a width-W frame: the MLE of one rank-one term, with soft
  site posteriors (efficient once in the right basin; seeded from substrings, so basin-finding is the hard part).
- **STREME** chooses PWM weights and a score threshold that maximise the enrichment of sequences with a hit, then tests
  on a hold-out. It searches for the best *point on the Neyman–Pearson frontier* (§2.1) of the contrast, seeded from
  exact k-mers. That is the lattice's cell test, made soft and full-length.
- **HOMER** seeds from oligos with mismatches (a Hamming kernel on full-length tables, §1.7) and then optimises a
  thresholded PWM for hypergeometric enrichment: again an NP-frontier point, from mismatch seeds.
- **Comsequent's walk and joins** grow a support one position or one join at a time: a greedy, combinatorial way of
  finding the rank-one factor's support with hard letters. The rank-one fit finds the same support in one step, with
  soft letters, from the same counts.

### 7.2 Matched detectors

For one shape's H table (signal S plus approximately Gaussian noise of variance ≈ 1/N per whitened cell), the optimal
test depends on the structure of S:
- **One cell** (crisp motif): the max-norm test, i.e. the best cell's z. This is Comsequent's ENR/θ, and the right
  test here.
- **Rank one** (PWM): the spectral test, i.e. the top singular value / tensor eigenvalue. It pools all 3k
  letter-position directions, with a BBP-type threshold: detectable once the signal exceeds the noise operator norm
  (≈ √(d/N)).
- **Diffuse** (spread interaction): the Frobenius/χ² test (U7's pooled G²).

Comsequent's own results fit this picture. ENR wins on crisp motifs; θ beats the pooled test on epi_pair; neither is
matched to a four-cell, rank-three structure, which is why neither dominated. The tool should run all three per shape
and Bonferroni over the three (cost ≈ 0.3 in z).

**How much can the rank-one detector add?** Within one shape it can recover at most the letter and weighting factors
of §2.3: on real motifs at k = 7 and π = 10⁻³ the medians are 1.16× (letter) and 1.03× (weighting), i.e. EXACT_7 1.54×
vs SOFT_7 1.20× of the oracle. Its main value is therefore representation and attribution (the δ it estimates *is*
the PWM restricted to the shape), not raw power. The large detection gains come from covering the full length (§2.6) and from removing control noise (§2.5).

### 7.3 Several motifs: deflation as well as masking

Explaining-away by window removal was Comsequent's single biggest real-data gain (0.61 → 0.83) and should stay. Theorem 2
adds a second operation: subtract a fitted motif's predicted contribution, π∏δ at every placement, from every lattice
table, then re-test. This removes the motif without removing other motifs' windows, and needs no recount.

### 7.4 Quantitative labels

With a per-sequence label y (peak strength, methylated fraction, SELEX round, concentration), the efficient
score-test statistic for a pattern feature g is Σ_i (y_i − ȳ) g(x_i). That is a label-weighted lattice count, and
Comsequent already has the machinery (`CountAnchoredWeighted`). Dichotomising a continuous label into fg/bg throws
efficiency away (a median split of a Gaussian covariate keeps 2/π ≈ 64%). Dimont, the overall winner of the
2025 GRECO-BIT/Codebook benchmark, learns discriminatively from per-sequence signal weights, in the same spirit.

### 7.5 Biophysical grounding of the domain rule

With an additive binding energy E(x) = Σ_j ε_j(x_j) and occupancy φ(x) = 1/(1 + e^{E(x) − μ}) (Berg & von Hippel 1987;
BEEML; ProBound):
- **Unsaturated** (E ≫ μ): bound sites ∝ Q(x)e^{−E(x)}, which is product form. **TF sites are rank one in H; the PWM is exact.**
- **Saturated** (high concentration, strong sites): bound sites ∝ Q(x)σ(μ − E(x)), a single-index model (rank one plus a
  monotone link: "global epistasis", Otwinowski et al. 2018). The log-linear θ will show apparent interactions that are
  not specific epistasis. Fit the direction and the link, as ProBound does.
- **Stringent enzymes** (MTases, REases, crisp PAMs): φ(x) ≈ 1[x ∈ A], a single (IUPAC) cell. **The lattice cell is the
  correct model.**
- **SELEX rounds** multiply by φ(x) each cycle: sites sharpen (towards cells) while composition drifts (§4).

So Comsequent's domain-of-applicability rule (crisp, matched, near-binary, structured) has a biophysical root:
recognition stringency (step vs exponential) sets crisp vs soft.

---

## 8. Search and calibration

| family | tests M (canonical) | two-sided Bonferroni z at α = 0.05 |
|---|---|---|
| gapped k5/w13 lattice | 2.5·10⁵ | 5.20 |
| gapped k7/w15 lattice | 2.5·10⁷ | 6.00 |
| contiguous 8-mers | 3.3·10⁴ | 4.81 |
| contiguous 10-mers | 5.2·10⁵ | 5.34 |
| contiguous 12-mers | 8.4·10⁶ | 5.82 |
| per-shape tests (rank-one / pooled), k5/w13 | 495 | 3.89 (more df per test) |
| per-shape tests, k7/w15 | 3,003 | 4.31 |

- The exhaustive lattice has an exact family size. Nesting makes the effective number of tests smaller, and Tarone
  pruning removes untestable cells (LAMP), so Bonferroni is conservative but valid.
- PWM tools search a continuous family. STREME reports p-values on a 10% hold-out (valid, at a power cost); MEME's
  E-values and HOMER's corrections are approximations. Only the lattice gives exact multiplicity on all the data.
- Calibration recipe (inherited from Comsequent, and required): λ ≈ 1 on motif-free splits for every statistic;
  background variance included; score z under separation; `sig` from Bonferroni or BH over the exact family.
- **Evaluation needs the same care as inference** (§3): every scorer gets a chance baseline from shuffled truths, and
  long calls must not buy hits with offsets.

---

## 9. Why the field moved to deep models, and the room that leaves

### 9.1 The reasons

1. **The question changed.** Single-TF specificity was largely catalogued in vitro (PBM, HT-SELEX; Codebook closing the
   gaps for poorly studied TFs). The open problems moved in vivo: which sites are bound in which cell type,
   cooperative grammar, quantitative effects, variant effects. These need many motifs modelled jointly, with spacing,
   position and quantitative output. That is beyond "find the top PWM".
2. **The benchmark changed.** The field evaluates sequence-to-function prediction (peak classification, base-resolution
   profiles, variant effects). Under that metric any flexible model trained on 10⁴–10⁶ examples beats a set of PWMs.
   Motif discovery became an interpretation step (attribution → TF-MoDISco seqlets → motifs), judged qualitatively.
3. **Data volume moved the binding constraint.** Classical algorithms were built for small samples (tens to hundreds of
   sites), where search, multiplicity and local optima are the bottleneck. With 10⁴–10⁵ peaks, detecting the main
   motif is easy, and the bottleneck becomes **representation**: several motifs, dimers, spacing, flanks, low-affinity
   sites. Gradient descent on a flexible model attacks representation directly.
4. **Classical tools plateaued on representation, not search.** MEME, STREME and HOMER all output single, fixed-width,
   additive matrices, found one at a time with greedy masking. Since ~2010 their improvements have been speed and
   search heuristics. No classical tool fits many motifs, interactions and positional effects jointly. This document's
   checks agree with the field's implicit judgement: detection is within 1.3–2× of optimal even for exact counting,
   and the losses are representation and ranking.
5. **Deep models fit the unified model class implicitly.** A first convolutional layer is a bank of rank-one weightings
   of the window one-hot: each filter is a PWM score. The nonlinearity plays the mixture likelihood ratio, pooling plays
   latent alignment (ZOOPS-like), deeper layers express sparse interactions among motif detections (grammar), and
   multi-task heads use quantitative labels. That is §7's low-rank + sparse model with §7.4's labels, fitted end to end.
   The field moved to deep learning for the model class, not because algorithms ran out.
6. **Incentives.** Machine-learning talent and tooling, consortium-scale data, publication incentives favouring
   architectures, and a sense since DREAM5 that "PWMs are near the ceiling for in-vitro specificity" all pushed effort
   to prediction.

### 9.2 What the deep route gave up

- **Data efficiency.** Training needs 10⁴+ sequences; the small-sample regimes (rare cell types, bacteria, archaea,
  phage, modification and enzyme specificities, low coverage) are left to classical tools.
- **Calibration.** No p-values for discovered motifs. TF-MoDISco takes hours and "reports many redundant and false
  motifs" (Yaish et al. 2025), and motifs vary across seeds and are split across filters.
- **Speed.** GPU training and attribution per dataset.
- **Explicit controls.** Confounders are learned unless factored out by hand (ChromBPNet's Tn5 bias model, and its own
  advice to check MoDISco output for bias motifs).
- **Exhaustive gapped support.** Filters have fixed width; bipartite motifs with long spacers are learned only if the
  architecture happens to compose them.

### 9.3 The algorithmic opening, stated precisely

gkm-SVM, a gapped-k-mer (lattice) feature space with a mismatch kernel (a soft lattice, §1.7) and a learned full-rank
weighting, matches CNNs at typical training-set sizes and is more robust (CNNs pull ahead only with larger training
sets; 2020 PLOS ONE comparison). That is direct evidence that **exact sufficient statistics plus the right weighting
get most of the way.** The missing pieces are structure (low rank: interpretable motifs instead of 10⁵ weights),
calibration, exhaustive gapped support and speed.

**What "extracting motifs from a deep model" is, in these terms.** The trained network is a flexible estimate of the
log likelihood ratio log ℓ(x) for every sequence. Attribution (DeepLIFT/SHAP/ISM) spreads that estimate over positions.
TF-MoDISco cuts out high-attribution seqlets (the sites), aligns them and clusters them, and the cluster averages are
the motifs. That is a **two-stage estimator of the rank-one components of §7.1**: stage 1 needs 10⁴+ sequences to fit
ℓ; stage 2 solves the alignment problem of §5 by clustering. The lattice route estimates the same components directly
from the moment tensor, where alignment comes from spanning shapes and anchors. It skips stage 1, so it needs far less
data, and its search is exhaustive, so its errors can be calibrated.

The opening is therefore an estimator of the deep model's first two layers (a motif bank plus sparse grammar) computed
directly from exact counts: calibrated, exhaustive over gapped supports, and running in seconds. It would also work
where deep models cannot be trained. Comsequent already has the counting engine, calibration and gapped search; this
document supplies the rank-one structure, the composition-free coordinates and the estimation step.

---

## 10. Design of the new tool

### 10.1 Principles (each tied to a result)

| # | principle | from |
|---|---|---|
| P1 | Count once, exactly, against a **large matched background** (default ≥ 5–10× foreground) | §2.5, §4 |
| P2 | Rank in **composition-free coordinates**: double-centred Markov contrasts (m = 1–3; phase-aware for coding) | §4 |
| P3 | Three **matched detectors** per shape: cell (crisp), rank-one/soft (PWM), pooled (diffuse); each calibrated | §7.2, §2 |
| P4 | Detect at **full length**: contiguous tables to L = 12, spanning gapped shapes; soft kernels by mode products. Removes the k-ceiling at small N | §2.2, §2.4, §6 |
| P5 | **Represent with all-data EM** from lattice starts; choose by the detection statistic, never by in-sample fit | §6 |
| P6 | **Align** via spanning shapes, anchors and position axes (summit-relative weighting for ChIP: Prop. 5) | §5 |
| P7 | **Residual dependencies**: θ/pooled tests on the residual from the rank-one fit, so "PWM + dependencies" only where significant | §1.5, §7 |
| P8 | **Several motifs**: explaining-away rounds (masking) plus deflation in H | §7.3 |
| P9 | **Quantitative labels** via weighted lattices when available (peak scores, fractions, rounds) | §7.4 |
| P10 | Report one list: detection z and `sig`, representation chosen by held-out likelihood (cell / IUPAC / PWM / PWM + dependencies), composition calls flagged; pre-checks printed first | §3, §8 |

### 10.2 Pipeline and compute budget (Go, N = 10k sequences of 100 bp)

| stage | operation | estimated cost |
|---|---|---|
| 0 pre-checks | signal/drift (existing), crispness index (per-position factor of top seeds, §2.2), background adequacy (bg counts of top words) | seconds |
| 1 count | gapped lattice k ≤ 7, w ≤ 15 on fg and a 10× bg; contiguous L-mer tables to L = 12; optional anchor/position axes and label weights | ~1–2 min at k7 with a 10× bg (E. coli k7/w15: 46 s per 9 Mb); seconds at k5 |
| 2 transform | Markov expectations (m = 1–3), soft kernels (mode products), H | seconds (12-mer mode products ≈ 0.8 Gflop) |
| 3 detect and rank | cell / spectral / pooled tests per shape, double-centred, calibrated | seconds |
| 4 represent | ZOOPS EM from the top ~20 de-duplicated seeds, all data | ~10⁹–10¹⁰ simple ops ≈ 10–60 s in Go |
| 5 residuals and rounds | θ on residuals; deflation / masking rounds | × rounds |
| **total** | | **target ≤ 2–4 min**, against HOMER 214 s and STREME 78 s at N = 10k (Comsequent sims); ENR alone is 4 s |

The speed headroom is real: everything except EM is a function of counts, and EM is linear in the data.

### 10.3 What is new, and what is prior art

New here (as far as this survey found): the H transform as the linear Möbius inversion of lattice enrichment ratios,
with Theorem 2's rank-one structure *in a fg/bg contrast over a gapped exhaustive lattice*; the nested-class
decomposition of the dilution gap and the optimal-conjunction law; the soft lattice transform as mode products on
shape tables; the matched-detector triple per shape; double-centred Markov contrasts as the composition fix for
exhaustive lattices; and the failure anatomy of Comsequent's real-data losses.

Prior art to cite and build on: Hoeffding/Efron–Stein (decomposition); Amari (η/θ, order decomposition); Anandkumar et
al. (tensor methods for latent-variable models); **FastMotif** (Colombo & Vlassis 2015: method-of-moments PWM discovery
from ungapped windows, single set, HT-SELEX only, sensitive to secondary motifs); Weeder/HOMER (mismatch oligos); gkm-SVM
(gapped k-mer kernels); MEME/STREME (ZOOPS EM, hold-out p-values); Dimont (weighted labels); ProBound/BEEML
(biophysical links); multireference alignment (Perry, Abbe, Bandeira, Bendory, Singer et al.); low-degree method
(Hopkins; Kunisky, Wein & Bandeira); sparse mixture detection (Ingster; Donoho & Jin); motif detection limits (Keich &
Pevzner 2002).

### 10.4 Prediction registry (register before running)

R2–R4 are informed by the c3 anatomy of Comsequent's real-data results, which included dev units. Those units are
therefore development data for the new tool. These predictions must be confirmed on the locked test splits.

| # | regime / application | prediction vs STREME/HOMER | mechanism | readout |
|---|---|---|---|---|
| R1 | ~~crisp: methylation, RM sites~~ | **downgraded (§11b): not a target**; parity at best against continuous-signal specialists | — | — |
| R2 | ChIP (ENCODE P5 dev factors) | **parity on strict recall** (from representation + composition fixes); win on calibration and speed | §3 anatomy: precision/fragment/composition misses are fixable | strict and chance-corrected recall; PWM similarity |
| R3 | HT-SELEX (P4a dev) | parity at late cycles; win over Comsequent by ≥ 0.2 | double-centred ranking removes cycle drift | chance-corrected strict recall |
| R4 | RBNS | no claim until the evaluator is chance-corrected | §3 | — |
| R5 | low-site PWM sims (N = 500, π = 0.1) | close to STREME; MEME may still lead | full-length seeds + all-data EM | PWM similarity, strict recall |
| R6 | ~~gapped / bipartite~~ | **withdrawn (§11b)**: variable spacers break spanning-shape alignment; needs U10 spacer unions | — | — |
| R7 | ChIP with a summit position axis | **failed (LN1)**: no gain; Prop. 5 holds only in the dense regime (§11b, LN1) | §5.3 | — |
| R8 | weighted labels (peak strength) | gain over fg/bg dichotomy, largest at small N | §7.4 | recall at N = 500–1,000 peaks |

### 10.5 Kill criteria

- If, on ChIP dev factors, strict recall with the representation layer stays more than 0.1 below STREME at N ≥ 2k
  after the evaluator fix, the "losses were representation" thesis is wrong. Stop and re-examine detection.
- If the double-centred ranking loses more than 0.1 recall against ENR on matched-composition dev data (where ENR is
  not flooded), the variance cost outweighs the robustness; route by the drift pre-check instead.
- If all-data EM from lattice seeds is worse than STREME's representation at N ≥ 2k (check 4 and dev data), the
  estimation layer needs MEME-style multi-start, and the speed advantage shrinks.

---

## 11. Following the DEV_INSIGHTS pathway (day 0 → day 5)

**Day 0 (this document plus the following):**
- [x] Mechanism: what the lattice exploits and discards, and what the PWM does (§1–§7). Per-regime predictions (§10.4).
- [x] Evaluation corrections sent to Comsequent: `../comsequent/docs/briefs/comsequent_eval_corrections.md`.
- [ ] **Evaluator first** (DEV_INSIGHTS 3.5). Add a shuffled-truth chance baseline to every scorer (c3 has the code);
  decide in advance whether degenerate call letters that contain the truth letter count; add a representation-level
  metric (PWM similarity) next to the IUPAC rule; oracle test (truth, IUPAC member, PWM form, near miss → 1/1/1/0).
- [ ] Comparator panel per class, at defaults: STREME, HOMER, MEME; add **Dimont** (GRECO-BIT's best) for ChIP and
  **ProBound** for SELEX (hostile-review answer to "did you beat the current best?").
- [ ] Reuse Comsequent's frozen dev/test splits; touch dev only.
- [ ] Pre-checks: signal/drift (exists), crispness index (new), background adequacy (new).
- [x] Matrix-source robustness of §2.4 (JASPAR vs HOMER): done, gap smaller with JASPAR.
- [ ] Dev panel ≤ 15 min: sims (L8 mid N = 500/2k; L12 low; planted HOMER motifs; crisp and bipartite classes; a
  CpG/GC-drift pair) plus real dev units (5 ChIP, 10 SELEX, 3 methylation strains).
- [ ] Stopping rule: replicate noise of the headline metric measured first; stop after two variants within it.

**Days 1–2: real-data smoke panel of the three new layers, in order of expected gain:**
1. **E1 representation layer on P5/P4a dev.** Re-run the all-in-one's top seeds through all-data EM; score strict,
   chance-corrected and PWM similarity. Tests R2 directly and cheaply.
2. **E2 double-centred ranking on P4a dev** (high drift) and P5 dev (CpG). Tests R3.
3. **E3 full-length seeds at small N** (sims; then the low-site real units).
4. **E4 c2 robustness**: matrices refitted on all peaks (JASPAR vs HOMER is already done: the gap is smaller with JASPAR).

**Days 2–4:** develop only in the regimes where E1–E3 win; Go ports of the soft lattice, double centring and EM.
**Days 4–5:** freeze; held-out once.

---

## 11b. Revision after Comsequent's final reflections and the first real-data results (2026-10-10, evening)

Inputs: `docs/final_reflections.md` (Comsequent's post-mortem, with partial held-out methylation results) and HASTUR's
first real-data run (24 ChIP dev factors, 10× GC-matched background; HOMER 23 of 24 done, Comsequent 11 of 24).

**ChIP dev, chance-corrected (development data: these factors informed the design, §3):**

| tool | P5 rule top 10 (chance) | contains top 10 | contains top 1 | best PWM similarity (chance) | median runtime |
|---|---|---|---|---|---|
| HASTUR v3 | 0.79 (0.03) | **0.54** | **0.46** | **0.887** (0.57) | **84 s** |
| STREME | **0.83** (0.03) | 0.46 | **0.46** | 0.861 (0.59) | 988 s |
| HOMER | 0.74 (0.01) | 0.35 | 0.26 | 0.825 (0.58) | 1,024 s |

HASTUR misses (P5 rule): CREB1 (RTGACGT, one letter short), MAFK (MTGACTCAGCA, two positions short), REST (16 of 20
positions), RFX5 (NFY found instead), ZNF143 (partial). All but RFX5 are representation-length misses; none is a
failure to detect.

### What the reflections confirm (the theory predicted these)

| reflection | theory | status in HASTUR |
|---|---|---|
| §2 "locates motifs but will not commit to them" (precision gap 0.08–0.30, one-sided) | §3: losses are representation and ranking, not detection | the representation layer (frame vs EM, sample-size trimming) is the fix; on ChIP dev HASTUR's located/stated rates exceed STREME's |
| §9 `sig` saturates under any composition mismatch; the Markov-background arm is the best-calibrated | §4: composition is a coordinate problem; double-centred contrasts | HASTUR ranks and gates in double-centred coordinates. Motif-free panel (10 pairs each), false `sig` calls per run, v5: iid **0.0**, same-chain Markov **0.2**, CpG shift **0.2** (STREME 3.8 under the shift; Comsequent saturates at 20,000). Better than every alternative, but above the nominal 0.05 on structured nulls: residual false calls are low-count cells (Poisson tails); next fix is exact small-count tests. v4 was 0.0 / 0.3 / 0.6: 6 of 9 false calls were reverse-complement palindromes counted twice (now z/√2) |
| §5 two readouts measured nothing; three of five rules cannot separate tools | §3, §8: chance floors | chance baseline in `eval/score.py` from day 0 |

### What changes

1. **Simulation is demoted from evidence to debugging.** The sim panel draws sites from a product PWM with fixed
   spacing in iid or order-1 backgrounds, which is exactly the model Theorem 2 inverts. HASTUR's sim margin over STREME
   at rank 1 is +0.25 (+0.31 on cells every tool can express). On real ChIP dev it is 0.00 at rank 1 and −0.04 on the
   P5 rule. **The simulator inflated our margin by about 0.3, as Comsequent's did.** From now on the sim panel is a
   regression and debugging tool; claims come only from real units, and the headline dev metric is real dev data.
2. **The rank-one model is an approximation whose scope is empirical.** Theorem 2 is exact for product sites, but real
   sites are occupancy-weighted (saturating, §7.5), come in several modes, and carry flank dependencies. The math is
   unchanged; the claim that a motif "is" a rank-one tensor becomes "is well approximated by one where real data says
   so". The representation choice by held-out evidence (frame vs EM) is the safeguard.
3. **Variable spacing breaks spanning-shape alignment.** Real bipartite methylation motifs lost by 0.27–0.42 despite
   the simulated win. In H terms, a motif with a spacer distribution w_s is Σ_s w_s (shifted rank-one terms): it is
   rank one in *no* single shape, so §5.2's "the span is the alignment" fails. The lattice-native fix is a spacer
   union (U10: a union over shapes, like an IUPAC letter is a union over bases), not yet built. **R6 (gapped/bipartite
   win) is withdrawn** until a variable-spacer real dataset says otherwise.
4. **The crisp niche is occupied and better served by continuous signal.** Nanomotif beats the lattice on its home
   ground by using modification probabilities directly; §7.4 predicts exactly this (dichotomising a continuous label
   costs efficiency). **R1 is downgraded: methylation is not a target for HASTUR.** If revisited, it needs the
   weighted (continuous-label) lattice, not thresholded windows.
5. **The winning thesis flips.** Comsequent's thesis was "win where specificity is crisp". The evidence supports the
   opposite: the lattice machinery is competitive *on PWM ground* once representation and composition are handled,
   and its advantages there are **speed (≈10× STREME/HOMER at a 10× background), calibration under drift, and
   exhaustive gapped support**, not accuracy. The realistic success criterion is: STREME-level accuracy on locked
   ChIP and SELEX test units, at ≥ 5× lower runtime, with a `sig` column that holds on drifted motif-free pairs.
6. **Optimism correction.** Comsequent's registry: every predicted loss held, every predicted win or parity failed.
   This registry was written by the same process, so its parity and win entries should be read as upper bounds.
   R2 (ChIP parity) has dev support but the dev factors shaped the design; it stands or falls on the 34 locked factors.

### Process changes (adopted)

- **Register the configuration**, not only the protocol: HASTUR version tag, flags, k/w/L, background multiple, Markov
  order, and comparator flags, pinned in `splits/` with a hash before any test run. One output list per tool; no
  union of arms.
- **Harness self-checks before scoring:** assert each tool can express each truth (max width ≥ the truth's
  informative span; HOMER defaults stop at 12 bp, STREME at 15). Report results on all units and on the expressible
  subset. The sim panel already had this asymmetry in HASTUR's favour (REST and p53 are 20 bp).
- **Report located vs stated** (contains / P5 / PWM similarity, each with its chance rate) for every phase.
- **Same denominators:** a tool with no output on a unit scores a miss (already the case in `bench/evaluate.py`).

### Expectations, restated

| target | expectation now | why |
|---|---|---|
| ChIP (locked test, 34 factors) | parity with STREME on located/representation, slightly below on strict; ≥ 5× faster | dev table above; design was tuned on dev, so allow a drop |
| HT-SELEX (locked test) | open; parity would be a strong result | high drift; short reads; the double-centring is the untested key |
| RBNS | no claim until a readout with a low chance floor exists | §3 |
| methylation, bipartite | not targets | items 3–4 |
| novel claim | factorial interaction (CAP-SELEX: θ between two TFs' contributions), never validly tested | reflections §12; no comparator answers it |

### LN1 result: the registered low-N ChIP panel (2026-10-10, protocol hash 9578ad1c)

Located top-10 recall, 24 dev factors × 3 random subsets per N, paired differences with 95% bootstrap CIs:

| N peaks | HASTUR | HASTUR centred | STREME | HOMER | HASTUR − STREME | runtime HASTUR / STREME / HOMER |
|---|---|---|---|---|---|---|
| 250 | 0.32 | 0.29 | **0.46** | 0.36 | **−0.14 [−0.25, −0.04]** | 13 / 37 / 293 s |
| 500 | **0.54** | 0.51 | 0.46 | 0.31 | +0.08 [−0.04, +0.21] | 18 / 84 / 506 s |
| 1,000 | **0.54** | **0.54** | 0.49 | 0.19 | +0.06 [−0.06, +0.17] | 33 / 194 / 703 s |

Registered verdicts: **accuracy thesis not alive**; "no gain" rule not triggered (the 500/1,000 differences are
positive but not significant); **summit weighting does not work** (−0.03 to 0.00). Same conclusions without REST.
On the P5 rule STREME is ahead at every N (0.75–0.81 vs 0.60–0.79). HASTUR is significantly ahead of HOMER from
N = 500 (+0.24, +0.35) and 3–6× faster than STREME.

**Why Prop. 5 failed (measured).** Real sites are summit-concentrated: the median site-position sd is 27 bp, so
1 + χ² ≈ 2.4 (CEBPB 5.1, MAFK 6.9, SP1 1.1). But at the consensus level 77–99% of foreground matches are already sites:
the lattice's best cells sit in the *sparse* regime (§2.2: planted matches outnumber chance ones u\* = 7–25-fold).
Prop. 5 is a dense-regime result. In the sparse regime positional weighting adds about log(1 + χ²) ≈ 0.9 nats of
likelihood ratio per site, less than one more specified base (≈ 1.2 nats). **Correction to §5.3 and Prop. 5:
positional information is worth (1 + χ²)× data only for statistics in the dense regime (soft, low-specificity
patterns); at the conjunction optimum it is worth less than one base.** It may still help very soft motifs with no
specific cell; that is untested.

**Standing of the accuracy question.** On ChIP, HASTUR is level with STREME from 500 peaks up, behind below that,
clearly ahead of HOMER, faster than both, and better calibrated under composition shift (§11b). The registered rules
allow no further low-N variants without a new protocol.

### LN2 Part B: the locked ChIP test (2026-10-11)

**The registered readout failed its own oracle test**, so the registered verdict is **uninformative** (neither pass
nor fail). The test truths are HOMER-library consensus strings, whose low-information flank letters their own
matrices do not support: scored as a call, each truth's own matrix passes the `contains` rule 6% of the time, P5 30%,
PWM similarity ≥ 0.85 18%. Located recall was ~0 for both tools. Process failure: the new truth source was registered
without the day-0 oracle check (DEV_INSIGHTS 3.5).

**Post hoc, secondary** (same runs, truths = the evaluator's projection of the true HOMER matrices, so the oracle
passes by construction; information, not confirmation): HASTUR (frozen v8) − STREME, 31 factors:

| peaks | located top-10 | P5 top-10 | PWM similarity | speed |
|---|---|---|---|---|
| 250 | −0.07 [−0.15, +0.02] | −0.18 [−0.28, −0.09] | −0.030 | 4× |
| 500 | +0.03 [−0.07, +0.13] | −0.09 [−0.18, 0.00] | −0.020 | 6× |
| 1,000 | −0.08 [−0.16, +0.01] | −0.17 [−0.27, −0.09] | −0.030 | 7× |
| 5,000 | −0.03 [−0.16, +0.10] | −0.19 [−0.36, −0.07] | −0.019 | 14× |

(P5 chance rates: STREME 0.04–0.09, HASTUR 0.02–0.03; the chance-corrected P5 gap is 0.11–0.17.)

**Reading.** On unseen factors HASTUR is behind STREME on stating full-length motifs (P5, PWM similarity), and level
to slightly behind on locating them. The dev "parity" came partly from the dev truth form: short canonical cores
reward a crisp short call, while full ChIP-derived matrices reward a full-length PWM, which is STREME's strength and
HASTUR's remaining weakness (shorter, softer matrices). This is the dev-to-test optimism the Comsequent reflections
warned about, measured. **"Matches STREME all round" is not supported.** Supported: 4–14× faster, and calibrated
under composition shift where STREME is not.

**LN5 step 0: anatomy of the matrix gap (2026-10-11, low-N ChIP dev, matrix readout).** The reading above blamed
"shorter, softer matrices". Measured directly, both halves of that are wrong, and the cause is narrower.
- **Not shorter.** HASTUR's best call covers the truth's informative columns as well as STREME's (0.95–0.96 vs
  0.94–0.95), and it is wider (17–24 columns vs 13–14). Uncovered columns cost both tools the same (≈ 0.025).
- **Not softer.** Fitting the covered columns as (1 − λ)·truth + λ·reference gives λ ≈ +0.01 for HASTUR and −0.01
  for STREME, whether the reference is uniform or the peaks' own composition. The same holds separately for core
  columns (IC ≥ 1) and flank columns.
- **Not ±1 phase blur, and not site count.** HASTUR's matrices come from 93% as many sites as STREME's.
- **What it is:** a direction-free column error, 0.080–0.093 total-variation residual vs 0.061–0.073.
  - **About a third** sits in units whose best HASTUR call is below rank 1. JUND is the clearest: a TRE call at
    rank 1, then a CRE/half-site hybrid `TGAC n TCA`. Erasing the rank-1 motif's sites does not fix it; the blurred
    column comes from weaker half-site variants below rank 1's threshold.
  - **The rest** sits in rank-1 matrices. One more hard-EM site-refit round (best window per sequence at the
    enrichment-optimal threshold) closes over half of the gap: 0.877 → 0.893 at N = 250 and 0.885 → 0.897 at N = 500,
    against STREME's 0.895 and 0.908. Further rounds are flat, and the same refit pulls STREME's own matrices down to
    ≈ 0.898.
- **Interpretation:** HASTUR's reported matrix is one refinement step short of the site-count fixed point. Beyond
  that fixed point, STREME holds a small edge that site counting alone does not reach. A plausible reason is that
  HOMER's truth matrices and STREME's are both discriminatively optimised, while a site-count matrix is not (open).
  Registered as LN5 (`splits/protocol_LN5.md`).

### 11c. Calibration by permutation max-T (LN7; confirmed in FR1, 2026-10-11)

**The failure.** v8 controlled significance analytically. Every lattice cell got a z on a Poisson/Wald model with
Markov double-centring, the maximum over about 10⁶ cells was taken, and Bonferroni plus genomic control (λ) was
applied. That needs each cell's null tail to be right at p ≈ 1e-8. On real sequence it is not:
- **Clumping:** occurrence counts are clumpy (low complexity, repeats, codon structure).
- **Small expected counts:** the Wald z on a log ratio with observed-count variance is badly anti-conservative.
- **λ cannot fix a tail:** it is a median, so it rescales the bulk.

v8 made 0.45–1.8 false calls per run on exchangeable real-genome nulls and 1.5–3 on small sets (STREME: 0.00).
Three rounds of per-statistic repair (LN6) each fixed one regime and broke another.

**The design change.** Calibrate the *maximum* empirically. Under exchangeability of fg and bg sequences
(Westfall–Young max-T):
- Relabel the pooled sequences R = 20 times (fg size kept), recompute each lattice family, and record its maximum.
- Fit a Gumbel to the maxima by moments, and use its upper α quantile as the family threshold. That gives exact
  family-wise error control under exchangeability, whatever the clumping, repeats or N.
- **It is cheap because counts are additive:** each permutation counts only the permuted fg (bg = pool − fg). The
  cost is about 1.5–2× a plain run.

**Three consequences:**
1. **Ranking and testing separate.** A max-T threshold bounds every cell's tested statistic at once, so seeds can be
   *chosen* by any statistic (v8's Wald and soft-kernel z, which pick good motifs) and *tested* by another (exact
   conditional binomial for gapped cells, exact radius-1 Hamming-ball counts for L-mers), with FWER intact.
2. **Genomic control becomes redundant, and harmful with strong signal.** In SELEX, selection shifts the whole word
   spectrum, so λ grew large and deflated true motifs from z ≈ 37 to ≈ 4. Without λ, the permutation thresholds stay
   near Bonferroni there.
3. **Composition robustness must still come from the statistic.** Permutation does not cover non-exchangeable
   nulls. Occurrence counts with v8's normalisation (expectation = observed shape total × Markov share) keep the
   simulated CpG shift at ≤ 0.05 false calls per run. Presence counts with absolute expectations broke it (4.45).

**Result** (frozen as F2, `splits/freeze_v18.md`; FR1 confirmation on fresh nulls and 478 unseen test positives):
- **Matched-threshold sensitivity** at 0.1 false calls per run is **+0.24 over STREME** on the registered real-null
  family (ChIP +0.24, SELEX +0.25).
- **The advantage comes from composition shift.** Post hoc, on exchangeable real nulls alone, the two tools are at
  parity (+0.01 [−0.03, +0.05]; ChIP −0.06, SELEX +0.10). The registered advantage comes from the GC-selected null
  sets, where STREME's E-values are badly miscalibrated.
- **Calibration:** about 0.06 false calls per run on exchangeable real sequence, 0.05 under CpG shift (STREME 3.45).
- **Native-threshold sensitivity** equals v8's at N ≥ 250.

**Open:**
- GC-selected sets (both tools miscalibrated);
- small N (N = 25 calls almost nothing);
- a held-out-candidate design remains an alternative, not tested.

## 12. Open problems and risks

- **In-vivo sites may be softer than any fitted matrix** (HOMER or JASPAR), which would widen the detection gap of
  §2.4. E4 measures this.
- **Real backgrounds are not Markov.** Repeats and CpG islands leave structure that double centring will not remove;
  combine it with matched genomic backgrounds and the factorial contrast.
- **The MRA hierarchy is an analogy.** A sample-complexity theory for sparse categorical mixtures under shifts is open;
  so is the exact null distribution of the spectral (rank-one) per-shape test with overlapping windows (calibrate
  empirically).
- **Dimers and variable spacers** are sums of shifted rank-one terms with a spacing distribution; U10's spacer merge
  is the lattice form. Not yet theorised here.
- **EM cost** dominates the budget; MEME-style multi-start would erase the speed advantage. The bet is that lattice
  seeds make one start per motif enough (check 4 tests it).
- **Saturation (single-index) regimes** need a link function; ProBound shows how, at a cost.

---

## 13. References

- Abbe, Bendory, Leeb, Pereira, Sharon, Singer (2019). Multireference alignment is easier with an aperiodic translation distribution. *IEEE Trans Inf Theory*.
- Amari (2001). Information geometry on hierarchy of probability distributions. *IEEE Trans Inf Theory*.
- Anandkumar, Ge, Hsu, Kakade, Telgarsky (2014). Tensor decompositions for learning latent variable models. *JMLR*.
- Avsec et al. (2021). Base-resolution models of transcription-factor binding reveal soft motif syntax (BPNet). *Nat Genet*.
- Bailey & Elkan (1994). MEME. *ISMB*. Bailey (2021). STREME. *Bioinformatics*.
- Bendory et al. (2024). The reflection-invariant bispectrum: signal recovery in the dihedral model. arXiv:2408.09599.
- Berg & von Hippel (1987). Selection of DNA binding sites by regulatory proteins. *J Mol Biol*.
- Colombo & Vlassis (2015). FastMotif: spectral sequence motif discovery. *Bioinformatics* 31(16):2623.
- Donoho & Jin (2004). Higher criticism for detecting sparse heterogeneous mixtures. *Ann Stat*.
- Efron & Stein (1981). The jackknife estimate of variance. *Ann Stat*. Hoeffding (1948). *Ann Math Stat*.
- Ghandi et al. (2014). Enhanced regulatory sequence prediction using gapped k-mer features (gkm-SVM). *PLoS Comput Biol*.
- Grau et al. (2013). A general approach for discriminative de novo motif discovery from high-throughput data (Dimont). *NAR*.
- Heinz et al. (2010). HOMER. *Mol Cell*.
- Hopkins (2018). Statistical inference and the sum of squares method. PhD thesis, Cornell.
- Keich & Pevzner (2002). Subtle motifs: defining the limits of motif finding algorithms. *Bioinformatics*.
- Krützfeldt, Schubach, Kircher (2020). The impact of different negative training data on regulatory sequence predictions. *PLoS ONE*.
- Kunisky, Wein, Bandeira (2022). Notes on computational hardness of hypothesis testing: predictions using the low-degree likelihood ratio.
- Le Cam (1956); van der Vaart (1998), *Asymptotic Statistics* §5.7 (one-step estimators).
- Otwinowski, McCandlish, Plotkin (2018). Inferring the shape of global epistasis. *PNAS*.
- Pampari et al. (2024). ChromBPNet. bioRxiv 2024.12.25.630221.
- Perry, Weed, Bandeira, Rigollet, Singer (2019). The sample complexity of multireference alignment. *SIAM J Math Data Sci*.
- Rube et al. (2022). Prediction of protein–ligand binding affinity from sequencing data with interpretable machine learning (ProBound). *Nat Biotechnol*.
- Shrikumar et al. (2018). TF-MoDISco. arXiv:1811.00416. Shrikumar, Prakash, Kundaje (2019). GkmExplain. *Bioinformatics*.
- Terada et al. (2013). Statistical significance of combinatorial regulations (LAMP). *PNAS*.
- Vorontsov et al. (2025). Cross-platform motif discovery and benchmarking to explore binding specificities of poorly studied human transcription factors (GRECO-BIT / Codebook). *Commun Biol*.
- Weirauch et al. (2013). Evaluation of methods for modeling transcription factor sequence specificity (DREAM5). *Nat Biotechnol*.
- Yaish et al. (2025). F-MoDA: efficient Fourier-based motif discovery in attribution maps. *PMLR* 311 (MLCB).
- Zhao & Stormo (2011). Quantitative analysis demonstrates most transcription factors require only simple models of specificity (BEEML-PBM). *Nat Biotechnol*.

---

## Appendix A. Proofs and derivations

**A1 (Theorem 1).** For product Q and any T, E_Q[ℓ(x) | x_T] = Σ_{x_{−T}} Q(x_{−T}) P_fg(x)/Q(x) = P_fg(x_T)/Q(x_T) = R_T(x_T).
The Efron–Stein components of a square-integrable function f of independent coordinates are
f_S = Σ_{T⊆S}(−1)^{|S∖T|} E[f | x_T], and they satisfy E[f_S | x_T] = 0 for S ⊄ T. Apply this with f = ℓ.

**A2 (Theorem 2).** Given in §1.5. For several motifs, R_T is affine in the mixture weights, so H is the corresponding sum.

**A3 (θ weak regime, Prop. 2).** The fg cell probability for match mask x on S is bg(x)[1 + π(site(x)/bg(x) − 1)]
(to first order in π, with bg = (1 − π)·product). θ_bg = 0 for product bg. Then
θ_S ≈ π Σ_x (−1)^{|S|−|x|} ∏_j ρ_j(x_j) = π ∏_j (ρ_j(1) − ρ_j(0)), with ρ_j(1) = p_j/q_j and ρ_j(0) = (1 − p_j)/(1 − q_j).

**A4 (Prop. 3).** Given in §2.2. The existence condition: d ln r/dk at k = 0 has the sign of ln(p²/q) when Q ≫ πP.

**A5 (local efficiency).** For g with finite variance, E_{P_fg}[g] − E_Q[g] = π Cov_Q(g, Λ), so the noncentrality is
Nπ² Cov_Q(g, Λ)²/Var_Q(g) ≤ Nπ² Var_Q(Λ) by Cauchy–Schwarz, with equality iff g is affine in Λ. For product P,
Var_Q(Λ) = ∏(1 + χ²_j) − 1, and expanding the product by subsets gives the order decomposition Σ_m e_m.

**A6 (double centring, Prop. 4).** If a set is exactly order-m Markov, the plug-in expectation E^{(m)}(w) is the MLE of
the expected count of w, so log(O/E^{(m)}) has mean O(1/O) for every w; this holds separately in each set, whatever
the two chains are. A motif's excess in a long word w enters E^{(m)} only through the (m+1)-mer margins, which change
by O(π × motif occurrences of that (m+1)-mer / total (m+1)-mers).

**A7 (positional weighting, Prop. 5).** With window weights w(r), the weighted count has mean excess
Nπ P ∫w f and null variance N Q ∫w² g (dense regime). Its squared z is Nπ²(P²/Q)(∫w f)²/∫w² g. Writing ∫w f = ∫(w√g)(f/√g)
and applying Cauchy–Schwarz gives (∫w f)² ≤ ∫w² g · ∫f²/g, with equality at w ∝ f/g, and ∫f²/g = 1 + χ²(f‖g).

## Appendix B. Numerical checks

| script | claim | runtime | output |
|---|---|---|---|
| `checks/c1_inversions.py` | Theorems 1–2, Prop. 2, rank-one share | 4 s | `results/c1_*.tsv` |
| `checks/c2_real_motif_efficiency.py` | §2.3–2.4 decomposition on 472 HOMER motifs (`C2_LIB=jaspar`: 1,293 JASPAR) | 2–5 s | `results/c2_real_motif_efficiency[_jaspar].tsv` |
| `checks/c3_failure_anatomy.py` | §3 anatomy (reads Comsequent results read-only) | 14 s | `results/c3_failure_anatomy.tsv` |
| `checks/c4_lattice_to_pwm.py` | §6 ranking and representation vs STREME | ~35 min on 11 processes | `results/c4_lattice_to_pwm*.tsv` |
| `checks/c5_composition_ranking.py` | §4 double-centred contrast | 40 s | `results/c5_composition_ranking.tsv` |
| `checks/c6_position_gain.py` | Prop. 5 positional weighting | < 1 s | printed |
| `checks/c7_spanning_shapes.py` | §5.1–5.2 placements and spanning shapes | 20 s | printed |
| `checks/c4_summarize.py` | merges c4 parts, prints the §6.3 table | 1 s | `results/c4_lattice_to_pwm.tsv` |
| `checks/c8_p1_chance_floor.py` | P1 methylation lenient rule: chance floor 1.7% / 1.2% (the win survives) | 5 s | printed |

Run with `~/micromamba/envs/nanomotif/bin/python theory/checks/<script>` (numpy, scipy, pandas; STREME from the
`meme` env; HOMER's motif library from the `homer` env).
