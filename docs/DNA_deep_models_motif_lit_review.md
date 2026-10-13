# Deep / Large Sequence Models × DNA: How They Represent Motifs, What They Gain Over Simple Algorithms, and What That Means for a Frontier Motif Detector

**Literature review, compiled 12 Oct 2026 for Oliver Rimington, as a hand-over to the coding agent working on the HASTUR / Comsequent motif-detection tools.**

> **How to use this document (for the coding agent)**
> - §1 is the executive summary: the one-page answer.
> - §2 is the quantitative "where does the gain come from" decomposition.
> - §§3–6 cover how models ingest, store and represent DNA.
> - §7 is the theory and maths, which is the most algorithmically useful part.
> - §§8–9 cover how motifs are pulled out of trained models and what that reveals.
> - §10 is the benchmark compendium with numbers.
> - §11 lists failure modes.
> - §12 is the concrete design and evaluation roadmap for a frontier motif detector.
> - §13 lists known gaps. §14 is the bibliography.
>
> **Verification legend.** Numbers were extracted from primary sources fetched during this review unless marked:
> - **[unverified]**: from background knowledge or a secondary source; check before relying on it.
> - **[derived]**: worked out in this review, not taken from a paper.
> - **[fig-only]**: the paper reports the value only in a figure, so no number is given.
>
> Several important numbers live only in figures. They are flagged rather than guessed.

---

## Contents
1. Executive summary
2. Where large/deep models actually gain over PWMs and k-mers (decomposition with numbers)
3. Input transformation: tokenisation
4. Positional encoding, strand symmetry, long-range operators
5. Model zoo and scaling laws
6. Internal representations: what is stored where
7. Theory and maths
8. Extracting motifs from trained models (algorithms)
9. What deep-model interpretation finds that PWMs miss
10. Benchmark compendium (numbers)
11. Failure modes and pitfalls
12. Design roadmap for a frontier motif detector
13. Open gaps in the literature
14. Bibliography

---

## 1. Executive summary

**1. "Large" is not where the motif-detection gains are.**
- Five independent 2024–2026 evaluations agree on this: DART-Eval, Tang & Koo 2025, Vishniakov et al. 2024, the IBIS challenge 2025, and Civale et al. (NeurIPS 2026).
- Self-supervised genomic language models (gLMs: DNABERT-2, Nucleotide Transformer, HyenaDNA, Caduceus, Evo 2) do **not** beat small supervised CNNs, and often do not beat k-mer features, on TF-binding, motif and regulatory tasks.
- Plain k-mer features (k = 4–6) beat frozen HyenaDNA embeddings on 86% of 57 datasets, DNABERT-2 on 95%, and Evo 2-1B on 93%. Only NTv3-650M beats k-mers consistently, on 63% of datasets.
- Per-motif discrimination from gLM embeddings is near chance (40–60%).
- No gLM reached the IBIS final.
- Randomly initialised models often match pretrained ones; for example, Caduceus scored +0.114 MCC higher when randomly initialised than when pretrained on GUE.

**2. The real gains come from supervised sequence-to-function models (CNNs, sometimes with attention) trained on functional data.** In rough order of size, they add:
- (a) **Many-site occupancy / avidity**: sum or log-sum-exp pooling over many weak sites, as opposed to best-hit scoring.
- (b) **Multiple binding modes per TF**: 2–4 PWMs combined beat one.
- (c) **Flanks and higher-order (Markov-order 2–5) context** around the core.
- (d) **Soft syntax**: spacing preferences, orientation, ~10.5-bp helical phasing, directional cooperativity.
- (e) **Non-equivalence of instances of the same motif**, driven by flanks and context.
- (f) **Bias/composition factorisation** (GC content, assay bias).

Within-core dinucleotide epistasis is real but small: a median gain of +0.006 to +0.019 in DREAM5, and no gain in vivo.

**3. Quantitatively, additive motif counts capture roughly half of what a good deep model captures.**
- ProCapNet: motif-count R² = 0.25 vs model R² = 0.57.
- The lentiMPRA lasso on features reaches r ≈ 0.72 vs EnformerMPRA r ≈ 0.84.
- The remainder is positions, flanks, interactions and weak sites.

**4. On pure motif discovery, classical PWM tools remain extremely competitive.**
- In GRECO-BIT/Codebook (4,237 experiments, 394 TFs), **Dimont** produced about half the top motifs, **ProBound** won on HT-SELEX, and **ExplaiNN** (a deep model) could not compete on genomic data.
- In IBIS, the PWM track was won by **CNN-proposed motifs (ZMotif) refined by STREME**. Deep "any-architecture" models led the PWM baseline only "by a limited margin".
- **BaMM** (5th-order Bayesian Markov motifs + ±4 bp flanks) beats lower-order models on 97% of 427 ChIP-seq sets (+13.6% median AvRec) On a genome-wide CTCF scan it calls ~27× more sites than the JASPAR PWM at about 1/30th of the FDR (1.8M sites at <0.5% vs 67k at ~15%).

**5. Theory explains why.**
- A first-layer conv filter is a PWM log-odds scanner.
- exp → sum-pool → log → sigmoid **is** the thermodynamic occupancy model; ProBound is exactly a 2-layer, weight-tied, RC-aware CNN with a count likelihood.
- Infinite-width CNNs are gapped-k-mer kernels.
- Transformers trained as language models converge on variable-order Markov / context-tree-weighting estimators, built from n-gram/induction heads.
- Masked LM pretraining "sees" a motif only through the mutual information *within* the motif, so short low-IC TF motifs are almost invisible to it.
- SGD learns "staircase" (leap-1) structure, such as motif plus extension, efficiently, but pays $d^{\Theta(k)}$ for zero-marginal interactions. That is no better than exhaustive enumeration, which favours exact combinatorial search for those.

**6. Implication for HASTUR / Comsequent.** The frontier is beaten by an **exact, statistically calibrated, occupancy-scored, multi-mode, flank-aware, syntax-aware** motif model, not by a large model. The ingredients:
- **Proposals:** exact gapped/sparse k-mer enumeration (Comsequent's strength), optionally plus small-CNN proposals.
- **Interaction filter:** replace Möbius-on-log-counts (which is *pointwise co-information*) with **maximum-entropy (IPF) residuals**, to test excess beyond all sub-patterns jointly.
- **Background:** a **variable-order (CTW/BaMM back-off)** model.
- **Scoring:** **LSE/occupancy** (sum over sites), not best hit.
- **Model refinement:** **mixtures of 2–4 PSAMs** plus **higher-order flank** terms plus an explicit **pairwise spacing kernel** ω(δ, strand), as in ProBound.
- **Evaluation:** AUPRC/AvRec at realistic priors, cross-platform/cross-cell-line transfer, first-motif-is-best rate, and false-discovery calibration via held-out E-values.

---

## 2. Where the gain over simple algorithms comes from: a quantitative decomposition

| Source of gain (vs single PWM, best-hit) | Representative evidence | Size of effect |
|---|---|---|
| Sum/occupancy over many sites (avidity) | DeepBind paper: MEME-SUM vs best single MEME motif; GRECO-BIT: Dimont's lead depends on sum-occupancy scoring | +0.04 AUC (0.82 vs 0.78); Dimont's lead mostly disappears under best-hit [fig-only] |
| Multiple modes per TF (2–4 PWMs) | GRECO-BIT "Archipelago" (random forest on 2–4 PWMs, 140 TFs) | Significant auROC/auPRC gains, mostly saturating by 4 PWMs [fig-only] |
| Higher-order Markov context + flanks | BaMMmotif2, 427 ENCODE ChIP-seq | 5th vs lower order: +13.6% median AvRec, 97% of sets; ±4 bp flanks: +11.4% (5th order) vs +3% (PWM) |
| Higher-order (cross-cell-line) | BaMMmotif2 | +32.7% vs ChIPMunk |
| Genome-scale precision | BaMM CTCF hg19 scan | 1.8M sites at FDR < 0.5% vs 67k at ~15% FDR (JASPAR PWM) |
| Within-core dinucleotides (in vitro) | DREAM5 PBM | Median per-TF +0.019 (FeatureREDUCE) / +0.006 (BEEML); **no gain in vivo** |
| k-mer vs PWM (in vitro PBM) | DREAM5 | Median gap 0.014; PWMs fail for ~10% of TFs (mostly C2H2 zinc fingers) |
| Supervised CNN vs PWM on ChIP | DeepBind | 0.90 vs 0.82 (MEME-SUM) mean AUC, 137 TFs |
| Supervised CNN vs PWM on HT-SELEX | DeepBind | 0.71 vs 0.63 |
| Cross-platform transfer (SELEX → ChIP) | DeepBind | 0.85 vs 0.82: the gap shrinks out of distribution |
| Syntax + flanks + instance effects (MPRA/STARR) | ProCapNet; lentiMPRA; DeepSTARR | Motif-count R² 0.25 vs model 0.57; lasso r ≈ 0.72 vs EnformerMPRA r ≈ 0.84; DeepSTARR PCC 0.68/0.74 vs replicate ceiling 0.73/0.76 |
| Cell-type context + chromatin | ENCODE-DREAM reanalysis | Removing ATAC features −0.12 AUPR; removing motif features −0.07 |
| Variant effects | Boytsov 2022 (SNP-SELEX) | Well-chosen PWM r = 0.531, ≈ ΔPWM 0.534; deltaSVM wins on ~74% of TFs; PWM wins on 34/129 |
| Variant effects (accessibility QTLs) | DART-Eval | ChromBPNet AUROC 0.772 (caQTL) / 0.892 (dsQTL) vs ≤ 0.65 for every gLM |
| Self-supervised gLM vs k-mer features | Civale 2026 (57 datasets) | k-mers win on 86–95% of datasets vs HyenaDNA / DNABERT-2 / Evo 2-1B |

**Reading the table.** The biggest, most robust gains come from:
1. counting all the sites (occupancy);
2. allowing more than one mode;
3. modelling local context to order 3–5 with flanks;
4. syntax between motifs, which is mostly what the deep sequence-to-function models capture.

All four can be built explicitly, interpretably and exactly.

---

## 3. Input transformation: tokenisation

### 3.1 Schemes

| Scheme | Models | Vocab | Notes |
|---|---|---|---|
| One-hot / single nucleotide | All supervised CNNs (DeepSEA, Basset, BPNet, Enformer, Borzoi, AlphaGenome), HyenaDNA, Caduceus, AIDO.DNA (A/T/C/G/N + specials), NTv3 (11 tokens), Evo/Evo 2 (byte-level) | 4–~16 (Evo char tokenizer 512 [unverified]) | Single-base resolution; a SNP changes one token |
| Overlapping k-mer | DNABERT (k = 3–6), SpeciesLM (6-mer) | 4^k + 5 (k = 6: 4,101) | **Leakage** under MLM (see below) |
| Non-overlapping 6-mer | Nucleotide Transformer v1/v2, GENERator | ~4,105 [unverified] | A 1-bp shift changes every downstream token (frame problem) |
| BPE | DNABERT-2 (4,096), GROVER (601), GENA-LM (32k) | 601–32k | ~5× compression for DNABERT-2; GENA-LM median token 9 bp |
| Learned / VQ / chunking | VQDNA (VQ-VAE codebook), MxDNA (mixture of conv experts + deformable conv), dnaHNet (H-Net differentiable chunking), DNACHUNKER, DNAMotifTokenizer | learned | Learned tokenisers gain ~1–2.5 points |
| Mixed modality | gLM2 (amino acids for CDS + nucleotides for intergenic + strand tokens) | — | Metagenomic contigs |

**Overlapping k-mer leakage.**
- Token $t_i$ covers bases $i..i+k-1$. If only $t_i$ is masked, $t_{i-1}$ and $t_{i+1}$ reveal all of its bases, so the task is trivial.
- With DNABERT's contiguous span masking of $k$ tokens, the flanking unmasked tokens still reveal all but one base. The $4^k$-way softmax collapses to roughly a 4-way choice [derived; matches the DNABERT-2 argument].
- Lesson: in any self-supervised k-mer scheme, the mask must cover more than the k-mer span or the objective learns nothing about motifs.

**GROVER BPE token statistics** (Sanabria et al., *Nat Mach Intell* 2024):
- 601 tokens, mean 4.07 nt, range 1–16.
- Only 32% of possible 4-mers appear as tokens.
- No CG-containing tokens formed.
- ~97% of tokens start with A or T.

So BPE vocabularies reflect genome composition and repeats, not TF motifs. In Lindsey et al. 2025, only **1.54%** of learned BPE tokens matched JASPAR 2024 motifs.

### 3.2 Controlled tokeniser comparisons

**DNABERT-2 ablation** (same data and architecture): BPE averages 65.33 vs 60.92 for 6-mers on GUE, winning 21/28 tasks. A larger vocabulary does not help monotonically.

**MxDNA ablation** (NT-v2 100M backbone, 100k steps):

| Tokeniser | NT benchmark | Genomic Benchmarks |
|---|---|---|
| 1-mer | 75.07 | 88.56 |
| Overlapping 6-mer | 74.35 | 88.55 |
| Non-overlapping 6-mer | **67.65** | 86.83 |
| BPE-4096 | 74.96 | 87.30 |
| MxDNA (learned) | **77.52** | 88.89 |

**VQDNA** (SARS-CoV-2 variant macro-F1, linear probe / fine-tune):

| Tokeniser | Linear probe | Fine-tune |
|---|---|---|
| 6-mer | 23.5 | 55.5 |
| One-hot | 5.5 | 54.1 |
| BPE | 36.5 | 71.0 |
| VQ-HRQ | 48.9 | 74.3 |

**Lindsey et al. 2025, *Bioinformatics*** (44 tasks, mean MCC):
- Caduceus 0.674, DNABERT-2 0.673, HyenaDNA 0.643, NT 0.634, Mamba-char 0.621, DNABERT 0.617, Mamba-BPE 0.599, CNN 0.594.
- Controlled Mamba, char vs BPE:
  - splice sites: **+0.224 MCC for character tokens** (p < 1e-4)
  - promoters: +0.029 for character tokens
  - SARS-CoV-2: +0.073 for BPE

**GENERator (autoregressive):** the opposite preference. The 6-mer beats BPE (BPE recovery < 0.30 vs a 0.25 random baseline). Random 0–5 nt frame offsets during training give robustness.

**Vishniakov 2024:** HyenaDNA with character tokens reached lower pretraining loss (1.180 vs 1.215) yet scored **0.187 lower downstream MCC** than with 6-mers. Pretraining loss is not a proxy for motif usefulness.

**Synthesis for motif work:**
- **Single-nucleotide resolution wins on position-precise tasks** such as splice sites and motif boundaries.
- Non-overlapping k-mers are the worst choice.
- Tokens learned by compression (BPE) mostly capture composition and repeats, not motifs.
- A motif tool should work at base resolution and treat k-mers as *features*, not *tokens*.

---

## 4. Positional encoding, strand symmetry, long-range operators

### 4.1 Positional encodings

Attention logits are $a_{ij}=q_i^\top k_j/\sqrt d$ plus any positional bias.

| Scheme | Mechanism | DNA models |
|---|---|---|
| Learned absolute | Added embedding | DNABERT, GENA-LM-BERT |
| ALiBi | $-m_h\lvert i-j\rvert$ bias | DNABERT-2 (trained at 700 bp, runs at 10 kb) |
| RoPE | $q_i^\top R_{i-j}k_j$ | NT-v2 [unverified], AIDO.DNA, GPN-MSA, GENA-BigBird, Evo 2 attention layers |
| Relative basis functions | Exponential + central-mask + gamma bases of $i-j$ | Enformer; Borzoi uses central-mask only |
| None (order in operator) | Implicit long conv / SSM recurrence | HyenaDNA, Caduceus, Evo (mostly) |

- **Evo 2 context extension:** positional interpolation plus a **10× RoPE base frequency per context doubling**, taking context from 8,192 to ~1M tokens.
- **Why it matters for motifs:** regulatory syntax is *relative* (helical phasing ~10.5 bp; TATA ~25–30 bp upstream of the TSS) and *oriented*. Relative schemes make attention translation-equivariant, like a convolution. Absolute positional encodings must learn shift invariance from data, which costs samples [derived].
- **Gong et al. 2025:** RoPE best for periodic motifs and length extrapolation; ALiBi good for local dependencies. Depth 3 → 12 layers helps; 24 layers adds little or overfits.

### 4.2 Reverse-complement (RC) symmetry

**Definition.** $\mathrm{RC}(X)$ flips the length and channel axes when the channel order is A, C, G, T. A model is RC-equivariant if $f(\mathrm{RC}(x)) = \mathrm{RC}(f(x))$.

**Mallet & Vert 2021 (NeurIPS)** characterise *all* linear translation+RC-equivariant layers.
- The group is $\mathbb Z\rtimes\mathbb Z_2$. Each representation decomposes as $\rho(s)=P\,\mathrm{Diag}(I_a, sI_b)P^{-1}$.
- A conv kernel is equivariant iff $\kappa(-x)=\rho_{n+1}(-1)\,\kappa(x)\,\rho_n(-1)$.
- Arbitrary pointwise nonlinearities are allowed only on "regular" and "+1" channels; "−1" channels need odd functions (e.g. tanh).
- "Irrep" models with $a/(a+b)=0.75$ beat RCPS (reverse-complement parameter sharing) on TF tasks, most in low-data regimes.
- "−1" channels can carry **strand-asymmetric** features, such as oriented promoter elements.

**Zhou, Shrikumar & Kundaje 2022 (MLCB):** RC augmentation alone still gives strongly strand-discordant predictions. *Augment during training, then conjoin (average the two strands) post hoc* was best or second-best throughout.

**Caduceus (ICML 2024):**
- BiMamba plus RCPS channel-splitting (MambaDNA); Caduceus-PS is equivariant by construction.
- Results: mouse enhancers 0.793 vs HyenaDNA 0.780; H3K9ac MCC 0.622 vs NT-v2 0.562; but NT-v2 is better on splice sites (0.983 vs 0.940).
- Beats Enformer on eQTLs more than 100 kb from the TSS [fig-only].

**Other models:**
- GPN-MSA averages both strands for its variant score.
- Evo/Evo 2 are single-strand autoregressive with no RC symmetry.

### 4.3 Long-range operators

- **Hyena:** $y = D_{x_2}T_hD_{x_1}v$. The $D$ are data-dependent diagonal gates; $T_h$ is the Toeplitz matrix of a long implicit filter $h_t=\gamma_\theta(t)$ (parameters do not grow with $L$). It runs in $O(L\log L)$ via FFT, up to 160× faster than attention at 1M tokens.
- **StripedHyena (Evo):** 29 Hyena layers plus 3 attention layers.
- **StripedHyena 2 (Evo 2):** short explicit, medium regularised and long implicit Hyena layers, plus RoPE attention. 1.3× faster than an optimised Transformer at 16k context, and 3× faster at 1M.
- **Selective SSM / Mamba-2 (SSD):** $M_{ij}=C_i^\top A_i\cdots A_{j+1}B_j$. With scalar $A_t$ this becomes $M=L\circ(CB^\top)$, i.e. masked linear attention with a data-dependent decay mask.
- **What each can express** [derived, plus Zoology/Arora 2023 unverified]:
  - Time-invariant long convolutions are linear and shift-equivariant, so they suit *fixed-spacing* kernels over motif-presence channels.
  - Content-based retrieval ("the earlier occurrence of this k-mer") needs gating or softmax attention.
  - Gated convolutions need state that grows with the number of key–value pairs to match attention on associative recall.

---

## 5. Model zoo and scaling laws

### 5.1 Summary table

| Model (year) | Objective | Params | Context | Tokens / output resolution | Training data | Notable numbers |
|---|---|---|---|---|---|---|
| DeepSEA (2015) [unverified] | Supervised | small CNN | 1 kb | one-hot → 919 labels | ENCODE/Roadmap | — |
| Basset (2016) | Supervised | CNN | 600 bp | one-hot → 164 DNase | ENCODE | ~45% of 300 first-layer filters match CIS-BP |
| DeepBind (2015) | Supervised | 1-layer CNN | ~100 bp | one-hot | PBM / SELEX / ChIP | ChIP AUC 0.90 vs MEME-SUM 0.82 |
| BPNet (2021) | Supervised profile | ~small CNN | ~1 kb | base-res profile (multinomial) + counts | ChIP-nexus | 51 motifs; Nanog ~10.5–10.9 bp periodicity |
| Basenji2 (2020) [unverified] | Supervised | dilated CNN | 131 kb | 128 bp bins | ENCODE etc. | — |
| Enformer (2021) | Supervised | 246M | 196,608 bp | 128 bp bins; 5,313 human tracks | ENCODE / FANTOM | 11 transformer blocks; CAGE r ~0.85 vs 0.81 Basenji2 [unverified] |
| Borzoi (2025) | Supervised | n/s | 524 kb | 32 bp bins, RNA-seq | 7,611 human tracks | U-Net upsampling |
| ChromBPNet (2024–25) | Supervised + bias model | small | ~2 kb | base-res ATAC/DNase | — | 41–60 non-redundant motifs per cell line |
| AlphaGenome (2026, *Nature*) | Supervised + distillation | n/s | 1 Mb | 1 bp, plus 128 bp and 2,048 bp pairwise | ~5,930 tracks, human + mouse | best on 22/24 track tasks and 24/26 VEP tasks |
| DNABERT (2021) | MLM | ~110M | 512 tokens | overlapping k-mer | GRCh38 | attention motifs: 1,595/1,999 match JASPAR |
| DNABERT-2 (2024) | MLM | 117M | ALiBi | BPE 4,096 | multispecies | GUE 66.80 vs NT-2.5B 66.93 at ~21× fewer params |
| GROVER (2024) | MLM | BERT-base | 510 tokens | BPE 601 | hg19 | — |
| Nucleotide Transformer v1 (2024) | MLM | 0.5–2.5B | 6 kb | non-overlapping 6-mer | hg38 + 3,202 genomes + 850 species | — |
| NT v2 [unverified] | MLM | 50–500M | ~12 kb | 6-mer | 850 species | — |
| NTv3 (2025) | MLM + supervised post-training | 8M–650M | 1 Mb | 1 nt, U-Net | OpenGenome2 ~9T bp + ~16k tracks | only gLM that consistently beats k-mers (Civale 2026) |
| GENA-LM (2025) | MLM | 110M / 336M | 4.5–36 kb (RMT ~50 kb) | BPE 32k | T2T + multispecies | — |
| HyenaDNA (2023) | Autoregressive | 0.44–6.6M | up to 1M | 1 nt | hg38 | — |
| Caduceus (2024) | MLM | 0.47–1.9M | 131k | 1 nt, RC-equivariant | hg38 | — |
| GPN (2023) | MLM (dilated CNN) | small | ~512 bp | 1 nt | 8 Brassicales genomes | ~160 motifs from TF-MoDISco on predictions |
| GPN-MSA (2025) | MLM on alignments | ~86M | 128 bp × 89 species | 1 nt/column | 100-way vertebrate alignment | beats CADD / phyloP on ClinVar [fig-only] |
| SpeciesLM (2024) | MLM + species token | ~90M | — | overlapping 6-mer | 1,500 fungal genomes | motif recovery ≈ STREME |
| Evo (2024, *Science*) | Autoregressive | 7B | 131k | byte | OpenGenome, 300B tokens (prokaryote + phage) | 5S rRNA fitness r = 0.64; essentiality AUROC 0.86 |
| Evo 2 (2025/26) | Autoregressive | 1 / 7 / 40B | 1M | byte | OpenGenome2 ~8.84T nt; 40B trained on 9.3T tokens | BRCA1 probe AUROC 0.94; SAE features |
| gLM2 (2025) | MLM | 150M / 650M | 4,096 tokens | amino acid + nucleotide + strand | OMG metagenomes, 315B tokens | — |
| AIDO.DNA (2024) | MLM | 300M / 7B | 4 kb | 1 nt | 10.6B nt, 796 species | 300M beat 7B on variant effects |
| GENERator (2025) | Autoregressive | 1.2B / 3B | ~98 kb | 6-mer | 386B nt of eukaryotic genic sequence | training on all 2T nt gave lower loss but worse downstream results |
| dnaHNet (2026) | Autoregressive, tokenizer-free | 10M–1B | — | learned chunks | prokaryotes | learns codon periodicity unsupervised |

### 5.2 Scaling

- **Evo:** over 300 models compared. Transformer++ had worse perplexity at every compute budget; Hyena and StripedHyena scaled best. The compute-optimal data budget for 7B was ~250B tokens. Exponents are given only in figures.
- **dnaHNet (2026):** the only explicit DNA scaling exponents found. With $\mathrm{PPL}=A\,C^{-\alpha}$, $\alpha=0.06$ (dnaHNet), 0.04 (StripedHyena 2) and 0.01 (Transformer). This was on prokaryotic data up to 1B parameters.
- **Perplexity improvements do not transfer reliably to downstream tasks.**
  - AIDO.DNA 300M beat 7B on variant effects.
  - GENERator trained on more data got lower loss but did worse downstream.
  - Lower loss with character tokens went with worse MCC (Vishniakov).
- **Gap:** no published head-to-head of gLM perplexity against an order-k Markov or context-tree-weighting (CTW) baseline on the same held-out DNA was found (§13). That makes it hard to know how much of gLM likelihood is beyond k-mer statistics.

---

## 6. Internal representations: what is stored where

- **Layer-wise probing.**
  - Nucleotide Transformer v1: the best layer is task-dependent, usually intermediate and **never the last**.
  - Datta et al. 2026 (DNABERT-2, NT, HyenaDNA, GENERator-v2, Omni-DNA): frozen probes recover 95–100% of fine-tuned performance on promoters but only 60–88% on splice sites, with splice information most accessible in middle layers.
- **Mostly k-mer composition.**
  - **Civale, Bagdanov & Magi (NeurIPS 2026)** split each embedding into a part linearly predictable from k-mer counts (k = 4–6) and a residual, across 57 datasets.
  - k-mer features alone beat HyenaDNA on 86% of datasets, DNABERT-2 on 95% and Evo 2-1B on 93%.
  - For HyenaDNA and DNABERT-2 **the non-k-mer residual hurts** classification.
  - Only NTv3-650M wins consistently (63.2% of datasets; MCC 0.562), mostly on splice and promoter tasks, which is likely due to its supervised post-training.
- **Tang & Koo 2025.** CLS embeddings are only marginally better than a dinucleotide-frequency MLP. Attribution maps look like "diffuse low-level statistics". Zero-shot CAGI5 Pearson r (HepG2):

  | Model | r |
  |---|---|
  | NT | 0.125 |
  | HyenaDNA | 0.064 |
  | GPN | 0.002 |
  | One-hot CNN | 0.324 |
  | Enformer | 0.510 |

- **Decode-gLM (2025).** Preserving 5-mer frequencies alone recovered **85%** of the probe probability in a CMV-enhancer probe of a Nucleotide Transformer SAE. Only E-box ablation mattered; ablating CREB/ATF, NF-κB or CCAAT did not.
- **Geometry.**
  - GROVER token embeddings: PC1 tracks frequency (ρ 0.88), PC2 tracks GC (−0.96), PC3 tracks AG content (0.94).
  - DNABERT-2 embeddings are strongly anisotropic, and whitening helps (Cuadrat 2026).
  - Embeddings of sequences differing by one SNP have cosine similarity ≥ 0.9, up to ~0.999 for NT-500M (Vishniakov). That means **single-base motif-breaking changes barely move the embedding**.
- **Pooling.** Options are CLS, mean, max, attention or top-k. A CNN trained on full per-token embeddings beats a probe on the CLS token (Tang & Koo).
- **Interpretable features (SAEs).**
  - Evo 2 (BatchTopK SAE, layer 26 of the 7B model) has features for exon, intron, exon boundaries, prophage/CRISPR spacers, frameshift/stop and α/β structure. Promoter-motif-like activations are mentioned, but **no TF-motif benchmark** is given.
  - An ICML 2026 SAE study on Evo 2 found SAE features beat raw embeddings on 42 of 55 concept evaluations.
- **DART-Eval motif probe** (1,443 HOCOMOCO v12 motifs inserted into neutral backgrounds):
  - Using *likelihoods*, gLMs generally prefer the true motif, with high per-motif variance.
  - Using *embeddings*, per-motif accuracy is mostly **40–60% (near chance)**.
  - A FIMO motif-count baseline gave the best cell-type embedding.

**Interpretation.** Self-supervised gLM representations are dominated by k-mer composition, GC content and repeats, with gene structure on top for the largest models. Specific TF-motif identity is weakly encoded. Supervised sequence-to-function models trained on TF or accessibility data are what encode motifs cleanly.

---

## 7. Theory and maths

Notation: $X\in\{0,1\}^{L\times4}$ is a one-hot sequence, $x_t=X_{t:t+w-1}$ is the window starting at $t$, and $q$ is the background distribution.

### 7.1 A conv filter is a PWM log-odds scanner

A width-$w$ filter $W\in\mathbb R^{w\times4}$ with bias $c$ computes

$$z_t=\sum_{j=0}^{w-1}\sum_b W_{jb}X_{t+j,b}+c.$$

- With $W_{jb}=\log(p_{jb}/q_b)$, $z_t$ is exactly the log-likelihood ratio $\log P(x_t\mid\text{motif})/P(x_t\mid\text{bg})$ (Stormo).
- Under a two-state binding model it equals $-\Delta\Delta G/RT$ plus a constant (Berg & von Hippel).
- Adding a constant to every entry of a column changes only the bias. CNN filters therefore have to be normalised before they can be compared with PWMs [derived].

**First-layer pooling decides what the filters learn** (Koo & Eddy 2019, *PLoS Comput Biol*; synthetic data with 12 JASPAR motifs):

| Pool size | Filter behaviour |
|---|---|
| ≥ 10 | Filters learn whole motifs: at least 88–98% match the ground truth (worst case ≥ 97% match some JASPAR motif) |
| ≤ 4 | Filters learn distributed fragments: at most 57% match any JASPAR motif, as few as 13% match the true motifs; CNN-4 has only 3 ground-truth filters |
| In vivo, 1 kb inputs | ≤ 57% (small pool) vs ≥ 83% (large pool) match JASPAR |

Pool size matters more than filter width. A width-9 filter gets 93% with pool 25 and none with pool 4. Deep layers still encode whole motifs, which saliency maps recover in both settings.

**Exponential first-layer activations** (Koo & Ploenzke 2021, *Nat Mach Intell*):

| Measure | Exponential | ReLU |
|---|---|---|
| Relevant JASPAR match, real DNA | 0.584 | 0.011 |
| Any JASPAR match, real DNA | 0.859 | 0.259 |
| Basset filters matching | 0.617 | 0.370 |
| Motif localisation AUROC from filter scans | 0.889 | 0.391 |

A true PWM scan gives 0.884 on the localisation measure. Predictive accuracy is unchanged (AUPR 0.489 vs 0.486). The reason is that $\exp$ is the Boltzmann weight (§7.2).

### 7.2 Pooling is an occupancy model; ProBound is a shallow CNN

For one binding mode, sum over offsets and strands $x$:

$$Z(S)=\sum_x e^{-\Delta\Delta G(S_x)/RT},\qquad P_{\text{bound}}=\sigma\big(\log Z(S)+\log\alpha\big).$$

That is the pipeline **conv (energy) → exp → sum-pool → log → bias → sigmoid**, i.e. log-sum-exp pooling followed by a sigmoid. Max-pooling is the zero-temperature limit:

$$\max_t z_t\le\mathrm{LSE}(z)\le\max_t z_t+\log N_{\text{sites}}.$$

Max-pooling and best-hit scoring therefore **cannot represent avidity from many weak sites**. LSE and sum pooling can [derived]. This is the theoretical basis of the sum-occupancy gains in §2.

**ProBound** (Rube et al. 2022, *Nat Biotechnol*):

$$Z_{\text{bound},r}=\alpha_{NS,r}+\sum_a\alpha_{a,r}\sum_x\omega_a(x)\,e^{\vec\beta_a\cdot\vec X(S_x)}+\sum_{a,b}\alpha_{a:b,r}\sum_{x_1,x_2}e^{\vec\beta_a\cdot\vec X(S_{x_1})+\vec\beta_b\cdot\vec X(S_{x_2})}\,\omega_{a:b}(x_1,x_2)$$

The terms are:
- **Features $\vec X$:** mononucleotides, with optional adjacent or non-adjacent dinucleotides.
- **$\omega_a(x)$:** a positional bias.
- **$\omega_{a:b}$:** a spacing/orientation cooperativity term.
- **$\alpha=[P]/K_D$:** activity.

Selection is modelled as

$$f_{i,r}=f_{i,r-1}\,\frac{Z_{i,r}^{\rho}}{(1+Z_{i,r})^{\gamma}},$$

and counts are Poisson. Profiling out the input frequencies leaves a multinomial over rounds:

$$\log L=\sum_{r,i}k_{i,r}\log\frac{\eta_r f_{i,r}}{\sum_{r'}\eta_{r'}f_{i,r'}}.$$

Across ~2,100 datasets ProBound beat JASPAR, HOCOMOCO, DeepBind, DeepSELEX and Jolma motifs on AUPRC and $R^2$, mostly at $P<10^{-3}$. In network terms it is a **2-layer, weight-tied, RC- and shift-aware CNN with exp activation, sum pooling, a pairwise relative-offset term that plays the role of attention, and a proper count likelihood**. Earlier models in this line are BEEML and BEESEM (Zhao & Stormo; Ruan & Stormo).

### 7.3 String kernels, the representer theorem, and CNNs as kernels

**Spectrum kernel:** $K_k=\langle\phi_k(S),\phi_k(S')\rangle$, where $\phi_k$ counts k-mers.

**Mismatch kernel:** $\phi_{(k,m)}(S)_\beta=\sum_{\alpha\in S}\mathbb 1[d_H(\alpha,\beta)\le m]$.

**Gapped k-mer kernel** (Ghandi 2014). For two $l$-mers $u,v$ that agree at $a(u,v)$ positions:

$$K(S,S')=\sum_{u\in S}\sum_{v\in S'}\binom{a(u,v)}{k},\qquad \binom{a}{k}=\sum_{T\subseteq[l],\,|T|=k}\ \prod_{j\in T}\mathbb 1[u_j=v_j].$$

So the gkm kernel is an inner product over **all order-$k$ positional match indicators** in an $l$-window. A weighted sum $\sum_k c_k\binom{a}{k}$ assigns a prior variance $c_k$ to each interaction order, which is the variance-component Gaussian-process view (Zhou, McCandlish et al. 2022 [unverified]).

**Representer theorem.** A fitted kernel model has the form

$$f(S)=\sum_i\alpha_iK(S_i,S)=\sum_t g(x_t),$$

where $g$ is a fixed-order polynomial in the one-hot window. gkm-SVM is therefore a **sum-pooled single-layer filter bank**, and that is what deltaSVM and GkmExplain exploit [derived].

**Infinite-width CNNs** (Novak et al. 2019 [unverified]). One conv layer with Gaussian weights, a nonlinearity and global sum or average pooling has the kernel

$$K(S,S')=\sum_{t,t'}\kappa\big(a(x_t,x'_{t'})\big).$$

That is the same structure as the gkm kernel with a different order spectrum $\kappa$ [derived].

**What CNNs add beyond kernel methods**:
- feature learning at finite width,
- max or LSE pooling instead of sums over windows,
- depth, which composes motifs into grammar.

**Strong cheap baseline:** exact kernel ridge regression over gapped matches with learned per-order variances.

### 7.4 What transformers and SSMs compute on sequences

- **Attention as a kernel smoother.** Each output is $y_i=\sum_j\mathrm{softmax}_j(q_i^\top k_j/\sqrt d)\,v_j$, a Nadaraya–Watson estimator. Uniform attention therefore outputs the *frequency* of any per-token feature, so counting "k-mer ends here" takes one head.
- **Statistical induction heads** (Edelman et al. 2024, NeurIPS).
  - **Task:** sequences come from a random Markov chain with Dirichlet-distributed rows. The Bayes-optimal in-context rule is $\hat P_{ij}=(c_{ij}+\alpha_j)/(N_i+\sum_l\alpha_l)$.
  - **Two-layer mechanism:** the first layer copies the previous token; the second attends to positions whose predecessor matches the current token. Together they compute an in-context bigram count.
  - **Training dynamics:** learning passes through uniform, then unigram, then bigram stages. The unigram solution is an attractor that delays the bigram stage.
  - **Higher order:** trigrams need 3 heads.
- **n-gram heads** (Akyürek et al. 2024, ICML). Transformers beat RNNs and CNNs at in-context learning of regular languages by using higher-order induction ("n-gram") heads. Hard-wiring n-gram heads into a 340M-parameter language model cut perplexity by up to 6.7%.
- **Variable-order Markov chains and CTW** (Zhou, Tian & Diggavi 2024).
  - **Construction:** a $(D+2)$-layer transformer can implement context-tree weighting (CTW), the Bayes-optimal predictor for variable-order Markov sources.
  - **Empirical behaviour:** one-layer models barely beat unigram. Models with 2 or more layers track CTW, saturating at about 4 layers.
  - **What is needed:** the MLP is required to learn the blending weights.
  - **Baselines:** fixed-order PPM and Kneser–Ney do worse than CTW.
- **Causal structure learning** (Nichani, Damian & Lee 2024, ICML). The first-layer attention gradient is proportional to the **mutual information** between token pairs, so first-layer attention picks out the highest-MI position pairs.
- **Associative memory** (Bietti et al. 2023). Weight matrices behave as $W=\sum_j v_ju_j^\top$. Global bigram statistics are learned before in-context ones.
- **Expressivity limits.**
  - **Hahn 2020:** soft attention cannot robustly compute parity or other maximally sensitive functions at length, because the effect of one input is $O(1/n)$.
  - **Merrill & Sabharwal 2023 [unverified]:** log-precision transformers are in TC⁰.
  - **What motif grammars need:** "contains M", "A within $d$ bp of B on the same strand" and "≥ n copies" are regular or threshold functions. They are easily within reach [derived].
  - **The hard case:** parity-like epistasis, which is hard for the optimiser rather than for the architecture (§7.6).

**Implication [derived].** A genomic LM is, to first order, a *learned variable-order Markov/CTW estimator with long-range retrieval*. A motif finder can compute the same statistics, pairwise MI and variable-order context likelihoods, **directly and exactly**.

### 7.5 What MLM and autoregressive pretraining can and cannot learn about motifs

The population minimiser of masked cross-entropy is the true conditional:

$$\mathbb E[-\log q_\theta(x_i\mid x_{\setminus i})]=H(X_i\mid X_{\setminus i})+\mathbb E\,\mathrm{KL}\big(p(\cdot\mid x_{\setminus i})\,\|\,q_\theta\big).$$

**Planted-motif analysis [derived].** Take a Markov background with PWM instances planted at prior rate $\pi$. For position $i$ at column $j$ of a possible site at offset $s$:

$$p(x_i\mid x_{\setminus i})=\gamma_s\,p_{j,x_i}+(1-\gamma_s)\,p_B(x_i\mid\text{ctx}),\qquad \gamma_s=P(\text{site at }s\mid x_{\setminus i}).$$

Here $\gamma_s$ is the posterior occupancy computed from the other $w-1$ columns. Two consequences follow:
1. When $\gamma_s\to1$, the LM's variant score reduces to PWM column log-odds.
2. The expected log-likelihood gain per instance is at most $\sum_j\mathrm{KL}(p_{j\cdot}\|q)=\mathrm{IC}_{\text{rel}}$, discounted by uncertainty in $\gamma$.

So **an MLM sees a motif only through mutual information within the motif.** Short, low-IC TF motifs, whose columns barely predict one another, are almost invisible to it. That matches the near-chance embedding results in §6.

**Zero-shot variant score (GPN):** $s=\log q(\text{alt}\mid\text{ctx})-\log q(\text{ref}\mid\text{ctx})$, averaged over strands.

**Where density estimation does recover motifs**:
- **GPN** (CNN, 8 Brassicales genomes): TF-MoDISco on the *predicted nucleotide distributions* found ~160 promoter motifs. Correlation with phyloP is r = 0.22.
- **SpeciesLM** (1,500 fungal genomes): masked reconstruction approaches close-genus alignment. TF-MoDISco on the reconstructions recovered about as many known motifs as STREME, and reconstruction fidelity predicted ChIP-exo binding better than phastCons.

Species conditioning and compact genomes help.

**Takeaway.** Use LMs as *background or conservation priors*, not as motif discovery engines. A motif is then the **residual information beyond the LM**.

### 7.6 Identifiability and sample complexity

**Planted (l,d) motifs** (Pevzner & Sze 2000; Buhler & Tompa 2002). Take $t$ uniform sequences of length $n$.
- Probability that a random $l$-mer lies within distance $d$ of a pattern: $p=\sum_{i\le d}\binom li(3/4)^i(1/4)^{l-i}$.
- Probability per sequence: $P=1-(1-p)^{n-l+1}$.
- Expected number of spurious patterns present in all $t$ sequences: $\mathbb E=4^lP^t$.

**Spurious counts for $t=20$, $n=600$** [derived; recomputed for this review]:

| (l,d) | log₁₀ E[spurious] | Verdict |
|---|---|---|
| (8,1) | −9.1 | easy |
| (10,2) | −7.2 | easy |
| (12,3) | −6.5 | solvable |
| (14,4) | −6.4 | solvable |
| (15,4) "challenge" | −14.7 | solvable in principle |
| (9,2) | +0.2 | **at the wall** |
| (15,5) | +0.45 | **at the wall** |
| (11,3) | +0.67 | **at the wall** |

When $\mathbb E\gtrsim1$, no algorithm can distinguish the planted motif from chance. **A tool should compute this information floor for each dataset's $(t,n)$ and refuse to report motifs below it.**

**Random projection** (Buhler & Tompa). Hash $l$-mers on $k$ random positions. Planted instances collide when the $k$ chosen positions avoid all mismatches, which happens with probability $\binom{l-d}{k}/\binom lk$ per instance.

**Detectability of PWMs.**
- **Schneider's rule:** $R_{\text{sequence}}\approx\log_2(G/\gamma)$ bits, where $G$ is genome size and $\gamma$ the number of sites.
- **Background hit rate:** roughly $2n\cdot2^{-\mathrm{IC}}$ false hits per sequence [derived].

**Significance practice in current tools.**
- **STREME:** Fisher or binomial test of positives vs controls. E-value = $p$ × number of motifs tested, evaluated on held-out data. It is accurate to within 2× and slightly conservative, with ≤ 5% type I error at 0.05.
- **MEME:** E-value of the log-likelihood ratio.
- **CentriMo:** binomial test of central enrichment.

**Size of the gapped search space.** There are $\binom{l-1}{k-1}4^k$ patterns with $k$ fixed positions in a span $l$ (first position anchored). Nested patterns are highly dependent, so Bonferroni over this space is over-conservative. Held-out or permutation E-values are tighter.

**How neural nets change sample complexity.**
- **Sparse parity** (Barak et al. 2022): SGD needs $n^{O(k)}$ steps, matching statistical-query lower bounds. Progress hides in small Fourier gaps.
- **Leap complexity** (Abbe, Boix-Adserà & Misiakiewicz 2023, COLT). SGD time scales as $\tilde\Theta(d^{\max(\mathrm{Leap}-1,1)})$. Leap is the largest number of *new* coordinates any single monomial adds, in the best ordering.
- **Mapping to motifs [derived]:**
  - **Leap-1 staircases (cheap for SGD and for greedy search):** PWMs, and motifs extended one informative position at a time.
  - **Leap-2k:** a pure XOR-like interaction with **zero marginals**. SGD pays $d^{\Theta(k)}$, which is no better than exhaustive enumeration. **Exact gapped-pattern enumeration with interaction residuals (Comsequent's territory) is optimal here, not inferior.**
- **Minimum-degree bias** (Abbe et al. 2023). Trained networks interpolate unseen inputs with minimum-degree (low-order) functions. That is an implicit prior towards low-order epistasis, the same phenomenon as spectral bias.

### 7.7 Epistasis bases: Möbius, Walsh–Fourier and Shapley

**Reference (Möbius) basis:**

$$f(x)=\sum_{T}\sum_{a\in\Sigma^T}\beta_{T,a}\prod_{j\in T}\mathbb 1[x_j=a_j],\qquad \beta_T=\sum_{U\subseteq T}(-1)^{|T\setminus U|}g(U).$$

**Walsh–Fourier basis** (Brookes, Aghazadeh & Listgarten 2022, *PNAS*):

$$\hat f(T)=2^{-L}\sum_xf(x)(-1)^{\sum_{j\in T}x_j},\qquad \mathrm{Var}f=\sum_{T\ne\emptyset}\hat f(T)^2.$$

- With $q$ letters there are $(q-1)^r$ coefficients per order-$r$ interaction.
- Landscapes are approximately sparse in this basis. For tRNA ($L=8,q=4$), the top 150 of 65,536 coefficients explain more than 90% of the variance.
- LASSO recovers an $S$-sparse landscape from roughly $C\cdot S\log(\cdot)$ samples, with $C\approx2.62$ empirically.

**Background-averaged vs reference-based epistasis** (Poelwijk 2016 [unverified]). Reference-free analyses of protein landscapes find that low-order terms dominate (Park, Metzger & Thornton 2024 [unverified]).

**Möbius on the k-mer lattice: important for Comsequent [derived].** Order gapped patterns by specialisation: $Q\preceq P$ when $Q$ is $P$ with some fixed positions replaced by wildcards. Then

$$\iota(P)=\sum_{Q\preceq P}(-1)^{|P|-|Q|}\log c(Q)$$

is the **pointwise co-information**.
- At order 2 it is exactly PMI, $\log\frac{p(ab)}{p(a)p(b)}$.
- At order 3 and above it is the alternating-sign generalisation.
- It is **not** the log-linear or maximum-entropy interaction term, so it can be non-zero for patterns fully explained by their sub-patterns, and the sign can be misleading.

**The principled "excess beyond all sub-patterns" statistic** is the maximum-entropy residual:

$$\iota^*(P)=\log c(P)-\log\hat c_{\text{IPF}}\big(P\mid\{c(Q):Q\prec P\}\big),$$

where $\hat c_{\text{IPF}}$ is the maximum-entropy expectation fitted by iterative proportional fitting to every proper sub-pattern marginal.
- For contiguous k-mers, the maximum-entropy model given the $(k-1)$-mer marginals is the order-$(k-1)$ Markov chain. In that case $\iota^*$ is the familiar Markov-residual statistic.
- For gapped patterns, IPF generalises this correctly.

**Shapley link.** Harsanyi dividends $m(T)=\sum_{U\subseteq T}(-1)^{|T\setminus U|}v(U)$ are a Möbius transform, and the Shapley value is $\phi_i=\sum_{T\ni i}m(T)/|T|$. Shapley–Taylor and Faith-Shap truncate at order $k$. **Möbius on the pattern lattice and Shapley interaction attribution on a trained network are the same algebra applied to different set functions.** Results from the two can be reported in one shared language.

### 7.8 How much of TF specificity is non-additive?

| Study | Finding |
|---|---|
| Weirauch 2013 (DREAM5 PBM, 66 TFs) | Dinucleotide versions: 0.880 → 0.924 (FeatureREDUCE), 0.898 → 0.919 (BEEML); median per-TF gains 0.019 / 0.006; gains concentrated where PWMs fit poorly; **no gain in vivo** |
| Zhao & Stormo 2011 [unverified] | Mononucleotide energy models fit most PBM data |
| Siebert & Söding 2016 (BaMM, 446 ENCODE) | Order 1: +16% pAUC; order 5: +29%; order 5 with 8-bp flanks: +36%; beats PWM on 97% of sets and is never worse; Klf4 double-mutant affinity correlation rises from 0.26 to 0.64 |
| Mathelier & Wasserman 2013, TFFM [unverified] | Dinucleotide-dependent HMM motifs improve ChIP-seq prediction for many TFs |
| GRECO-BIT 2024 | Motif length, IC and GC only weakly related to performance; multi-PWM random forests consistently beat a single PWM |

**BaMM back-off** (a Bayesian hierarchy over context orders, closely analogous to CTW):

$$p_j^{(k)}(x_j\mid x_{j-k:j-1})=\frac{n_j(x_{j-k:j})+\alpha_k\,p_j^{(k-1)}(x_j\mid x_{j-k+1:j-1})}{n_{j-1}(x_{j-k:j-1})+\alpha_k},\qquad\alpha_k=20\cdot3^{k-1}.$$

**Synthesis.** Within-site epistasis is small, a few AUC points in vitro. In vivo gains come from multiple modes, flank and higher-order context, many-site avidity and co-binding grammar.

### 7.9 Inductive biases checklist

- **Translation equivariance:** convolution, or relative positional encodings in attention.
- **RC equivariance:** Mallet–Vert kernels, RCPS, or post-hoc conjoining; keep "−1" channels for strand-asymmetric signals.
- **Pooling:** global max or LSE gives shift invariance; the choice between them is a modelling claim (one site vs many).
- **Receptive field:** $r=1+\sum_\ell(k_\ell-1)d_\ell$. Exponential dilations give $r\sim2^{\text{depth}}$, the right bias for 10–1,000 bp spacing grammar.
- **Overparameterisation:** more first-layer filters (30 → 120) lowered the ground-truth match rate from 0.98 to 0.89 ("idle" filters). No genomic double-descent study was found.

---

## 8. Extracting motifs from trained models: algorithms

### 8.1 Attribution

**In-silico mutagenesis (ISM).** $\mathrm{ISM}(i,b)=f(X^{i\to b})-f(X)$, which takes $3L$ forward passes.
- Mean-centre each position across the four bases.
- ISM stays exact and on the simplex, but it is first-order only.
- Acceleration:
  - **Yuzu** (Schreiber & Singh 2022) uses compressed sensing. Its cost is roughly constant in $L$, with speedups growing with sequence length (>10× end to end).
  - **fastISM** caches activations outside the receptive field.

**Gradient methods.**
- Saliency: $\partial f/\partial X$.
- Gradient × input.
- **Integrated Gradients:**
  $$\mathrm{IG}_i=(x_i-x'_i)\int_0^1\partial_if\big(x'+\alpha(x-x')\big)\,d\alpha.$$
  It satisfies completeness, $\sum_i\mathrm{IG}_i=f(x)-f(x')$. It is computed with 20–100 steps from a shuffled baseline.

**DeepLIFT / DeepSHAP.**
- Contributions satisfy summation-to-delta, $\sum_iC_{\Delta x_i\Delta t}=\Delta t$, with multipliers $m=C/\Delta x$ combined by the chain rule.
- The **Rescale** rule replaces the local gradient with $\Delta y/\Delta x$.
- **Hypothetical contributions** are $\Delta x^h\cdot m$: what each unobserved base would contribute.
- DeepSHAP averages over references. The genomics convention is about 10–20 dinucleotide-shuffled references [count unverified], which removes composition from the attributions.

**Off-simplex gradient correction** (Majdandzic, Rajakumar & Koo 2023, *Genome Biol*). One-hot inputs live on a simplex, so the gradient component normal to the simplex is arbitrary noise. Correct it with
$$\tilde g_{i,b}=g_{i,b}-\tfrac14\sum_{b'}g_{i,b'}.$$
- Between 20% and 40% of positions have gradients more than 30° off the simplex, and 5–15% are more than 60° off.
- The correction consistently improves AUROC and AUPR against ground truth.
- **It costs nothing, so always apply it.**

**Robustness.**
- Attributions vary across random seeds.
- Ensembles of 10 models, or distilled students (DEGU, Zhou et al. 2025), reduce the variance.
- Accuracy does not imply interpretability: an accurate sigmoid CNN had poor attributions.

### 8.2 TF-MoDISco (Shrikumar et al., arXiv:1811.00416) step by step

1. **Seqlets.**
   - Score sliding-window sums of attributions against a null. The null is either window sums from shuffled sequences or a Laplace fit to the mode-centred distribution (percentiles 5–95).
   - Fit isotonic regression of real vs null scores (null weighted to a 50% prior) to get an FDR.
   - Clamp the passing fraction to [0.03, 0.20].
   - Take seqlets greedily from the top score, expand by a flank, and suppress anything with >50% overlap.
2. **Metaclusters.** Assign each seqlet a sign pattern across tasks ($\{-1,0,1\}^{n}$). Keep a pattern only if at least min(100, 1%·N) seqlets have it.
3. **Coarse similarity.** Embed the summed hypothetical contribution track as gapped k-mer scores. Compare with cosine similarity, taking the maximum over RC. Keep k-NN with k = 500.
4. **Fine similarity: continuous Jaccard.**
   - Define $x\cap y=\min(|x|,|y|)\,\mathrm{sign}(x)\,\mathrm{sign}(y)$ and $x\cup y=\max(|x|,|y|)$.
   - Then $J(v_1,v_2)=\sum_i(v_1^i\cap v_2^i)/\sum_i(v_1^i\cup v_2^i)\in[-1,1]$.
   - Compute it on L1-normalised, concatenated actual and hypothetical tracks, maximised over offsets with at least 70% overlap and over RC.
   - Drop seqlets whose coarse and fine similarities have Spearman correlation below 0.15.
5. **Density-adaptive affinity.**
   - Distance $d=\log(1/(0.5\max(x,0))-1)$.
   - $p_{j|i}\propto e^{-\beta_id_{ij}}$, with $\beta_i$ set to give perplexity 10 (t-SNE style).
   - Symmetrise.
6. **Cluster.** Use Louvain with 200-run consensus, or Leiden with 50 seeds.
7. **Aggregate.** Align greedily by best Jaccard offset and average the tracks.
8. **Trim.** Keep positions with support ≥ max(0.2·max, 30). Take the best 30-bp window and add flanks.
9. **Split spurious merges.** For clusters with more than 30 seqlets, split in two if the sub-motif Pearson correlation is below 0.8, recursively.
10. **Merge.** Merge clusters using cross-correlation plus density similarity, e.g. $(p>10^{-4},c>0.84)$ or $(p>10^{-6},c>0.90)$.
11. **Clean up.** Disband clusters with fewer than 30 seqlets and reassign members where J > 0.2.

**tfmodisco-lite** (Schreiber) is the same algorithm, faster. On 100k seqlets it ran in about 1.5 h and 8 GB, against 5 h and 35 GB for the original.

**CWM vs PWM.**
- The contribution weight matrix is $\mathrm{CWM}[p,b]=\frac1{|S|}\sum_{s\in S}C_s[p,b]$: what the model *uses*.
- The PWM is the log-frequency matrix of the aligned instances.
- The CWM down-weights positions that are conserved but uninformative.
- **CWM scanning correlates much better with ChIP-nexus signal than PWM scanning.** For Nanog, Spearman ρ is 0.36 vs 0.06 (BPNet).

**Hit calling with Fi-NeMo** (Kundaje lab). Reconstruct the contribution track sparsely as a non-negative combination of CWMs on both strands:
$$\min_{\beta\ge0}\ \|C-\textstyle\sum_k\beta_k*\mathrm{CWM}_k\|^2+\lambda\|\beta\|_1.$$
- Default λ = 0.7 (about 0.6 for ChIP). It is solved by proximal gradient on GPU.
- Overlapping motifs compete for the same signal, so redundant hits are suppressed.
- Hits are scored by shape (cosine similarity) as well as by coefficient.
- The objective here is paraphrased from the README.

**ARGMINN** (Tseng et al. 2024). An interpretable-by-design model with 8 filters of 10 bp, an overlap penalty plus L1, and a memory-query attention "syntax builder".
- In simulations it produced **0 extraneous motifs, against 7–10 for TF-MoDISco**.
- Instance precision was 0.953 vs 0.852 for FIMO with MoDISco motifs (SPI1).
- On E2F6, precision was 0.954 vs 0.777 and recall 0.816 vs 0.584.
- Predictive accuracy is on par with standard CNNs.

### 8.3 Filter-level and interpretable-by-design models

- **DeepBind (2015):** align the subsequences that activate each filter.
- **Basset:** about 45% of 300 filters match CIS-BP. Filter influence is measured by nulling a filter.
- **ExplaiNN** (Novakovsky 2023): N independent single-filter units (19 bp, BN → exp → max-pool → 2 FC layers) combined linearly.
  - Accuracy plateaus at about 100 units, close to DanQ.
  - For 85% of TFs, the best PWM came from one of the top 10 units.
  - It is more than 100× faster than STREME on large datasets, **but STREME was better on small datasets**.
  - It cannot represent motif interactions.

### 8.4 Global, epistatic and perturbation analyses

- **Global Importance Analysis (GIA)** (Koo 2021): $I(\phi)=\mathbb E_{x\sim p}[f(x^\phi)-f(x)]$, which embeds a hypothesised feature (motif, pair at spacing $d$, GC change) into sampled backgrounds.
- **Pair-insertion super-additivity:**
  $$f(AB)-f(A)-f(B)+f(\varnothing),$$
  computed over distance, orientation and strand. Applied in BPNet and DeepSTARR.
- **DFIM** (Greenside 2018): the feature interaction score is $\mathrm{FIS}(s\to t)=C_t(X^{s\text{ mut}})-C_t(X)$. It needs $O(F)$ backward passes rather than $O(F^2)$. Significance comes from a Gaussian null built on shuffles, then Benjamini–Hochberg.
  - GATA1–TAL1 interactions are strongest below 20 bp and detectable up to about 75 bp.
  - On Cbf1 BET-seq (over 1M variants) it found flank–flank epistasis such as T₋₁/A₊₁.
- **SQUID** (Seitz 2024, *Nat Mach Intell*): build an in-silico MAVE of 100k variants at a 10% mutation rate, then fit a MAVE-NN surrogate $y=g(\phi(x))+\text{noise}$.
  - $\phi$ is additive or pairwise; $g$ is a global-epistasis nonlinearity (a sum of sigmoids); the noise is a skew-t.
  - **Without $g$, saturation masquerades as uniform positive epistasis.**
  - Nanog–Sox2 interaction strength oscillates with spacing over 0–32 nt (helical).
  - Attribution variance is lower than for ISM or DeepSHAP.
  - The pairwise surrogate becomes intractable above about 100 nt.
- **SEAM** (Seitz 2025): cluster the attribution maps of a variant library to expose *mechanistic heterogeneity*.
  - A single substitution flips an AP-1 site from activator to repressor-like.
  - A double mutant creates a CAAT box that reverses transcription direction.
- **CREME** (Toneyan & Koo 2024): 5-kb tile perturbations of Enformer around 203 TSSs.
  - About 5 enhancers drive more than 80% of activity.
  - Enhancer sets split as 59 redundant, 28 cooperative and 65 additive.
  - **Necessity does not imply sufficiency.**

### 8.5 Transformers and genomic LMs

- **DNABERT attention motifs:**
  - Keep regions where attention is above the mean, more than 10× the minimum and at least 4 tokens long.
  - Test enrichment with a hypergeometric test (BH p < 0.005).
  - Merge by alignment in a 24-bp window with at least 3 instances, then build a PWM and run Tomtom.
  - Result: 1,595 of 1,999 motifs match JASPAR.
  - This is essentially k-mer enrichment seeded by attention. The "attention is not explanation" caveat applies.
- **DART-Eval:** gLM embeddings are near chance on motif identity.
- **Sparse autoencoders:**
  - Evo 2 SAEs give gene-structure features, but there is no TF-motif benchmark.
  - Decode-gLM SAEs on NT lean on 5-mer statistics.
  - **Current verdict: do not use gLM attention or embeddings as a motif discovery engine.**

---

## 9. What deep-model interpretation finds that PWMs miss

1. **Flanks change site strength.**
   - DeepSTARR GATAA sites: important instances have G at +1; unimportant ones have T at +1 and G at −1.
   - Swapping 2-bp flanks between strong and weak sites in 47 enhancers moved activity as predicted, and 5-bp swaps were slightly stronger.
   - Trl site importance depends on up to 5 bp of GA-repeat flank.
   - DFIM finds flank–flank epistasis (Cbf1).
2. **Instances of the same motif are not equivalent.**
   - Across 1,013 GATA instances, predicted vs measured mutation effects correlate at PCC 0.53.
   - 57% of multi-instance enhancers have a >2-fold difference between instances.
   - Syntax-aware linear models predict instance importance at PCC 0.50–0.66 in human HCT116 cells, beating PWM-only models.
3. **Syntax is soft: preferred ranges rather than strict grammar.**
   - BPNet: Nanog periodicity of about 10.5 bp (10.9 ± 0.6 in the Fourier analysis); protein-range interactions below 35 bp and nucleosome-range ones at 70–150 bp; directional cooperativity (Sox2 → Nanog, not the reverse); no co-occurrence enrichment beyond 150 bp.
   - Strict spacings are mostly ERV-derived (83% of the top 1% of pair distances).
   - DeepSTARR pair modes: synergy only below 25 bp (GATA/ETS); synergy only above 25 bp (GATA/GATA); close preference with a plateau beyond 75 bp (GATA/AP-1); monotone decay (housekeeping pairs).
   - ChromBPNet composite elements: FOS–TEAD at 6 bp, IRF1–SPI1, IRF1–AP1.
4. **Motif content does not determine activity.** Three designed enhancers with identical motif content (2 GATA, 2 AP-1) ranged in activity from 0.87 to 630.
5. **Additive motif counts explain about half.** ProCapNet: R² 0.25 vs 0.57.
6. **Low-affinity and homotypic sites tune dosage sensitivity.** Naqvi et al. 2025 (*Cell*; TWIST1/SOX9 dosage) found that distributed low-affinity and homotypic sites barely affect baseline accessibility but *sensitise* it to TF dosage. Central high-affinity heterotypic sites *buffer* against dosage changes. Both show purifying selection.
7. **New sub-motifs.** ProCapNet found cryptic Inr-like sub-motifs inside ETS, YY1, NRF1, AP1 and ATF1 motifs, and TCT motifs in 4 of 6 cell lines. 99.7% of PRO-cap peaks contain predictive motif instances, about 6 per peak.
8. **Bias motifs.** Tn5 and DNase enzyme-bias motifs appear unless a bias model is factored out (ChromBPNet). The ATAC vs DNase profile JSD fell from 0.58 to 0.26 after correction.

---

## 10. Benchmark compendium

Values marked [fig-only] are reported only in figures. Check the PDFs before quoting any of these numbers in a publication.

### 10.1 TF binding: PWM vs k-mer/gkm-SVM vs CNN

**DREAM5 PBM** (Weirauch et al. 2013, *Nat Biotechnol*). 66 mouse TFs, 26 algorithms, about 40k 60-mer probes per array. The score combines Pearson r and AUROC on bright probes.

| Algorithm | Model class | Final score |
|---|---|---|
| FeatureREDUCE | dinucleotide + k-mer | 0.997 |
| Team_D | k-mer | 0.984 |
| Team_E | PWM with a strong background model | 0.952 |
| Team_G | k-mer | 0.950 |
| FeatureREDUCE_dinuc | dinucleotide | 0.924 |
| BEEML-PBM_dinuc | dinucleotide | 0.919 |
| BEEML-PBM | PWM | 0.898 |
| FeatureREDUCE_PWM | PWM | 0.880 |
| MatrixREDUCE | PWM | 0.791 |
| Seed-and-Wobble | PWM | 0.647 |
| RankMotif++ | PWM | 0.582 |

- The median per-TF gap between the best k-mer model and the best PWM was 0.014.
- The best PWM beat the best k-mer model on 18 of 66 TFs.
- **On in vivo data (ChIP-seq/exo, 9 TFs), PWMs mostly beat the k-mer and dinucleotide models, and FeatureREDUCE was substantially worse.**
- The best PWMs had *low* information content.

**DeepBind** (Alipanahi 2015)

| Benchmark | Score |
|---|---|
| ChIP-seq, 137 TFs, mean AUC | DeepBind 0.90 (all peaks), 0.85 (top-500 setting); MEME-SUM 0.82; MEME best single motif 0.78 |
| HT-SELEX, 566 models | DeepBind 0.71; Jolma PWMs 0.63 |
| Trained on SELEX, tested on ChIP (35 TFs) | DeepBind 0.85; best PWM 0.82; MEME-SUM-5k 0.82 |

**Zeng, Edwards, Liu & Gifford 2016** (*Bioinformatics*; 690 ENCODE datasets, 101-bp inputs)

- A 1-layer CNN with global max-pooling matched DeepBind closely (R² of per-dataset AUCs = 0.886).
- Increasing the number of kernels helps up to about 128, then saturates.
- **Local pooling and extra conv layers made discovery worse.**
- On motif occupancy, CNNs beat gkm-SVM (AUC about 0.8).
- At 40k training sequences all architectures performed similarly. **At 10k and 2.5k, complex models were clearly worse**, with ANOVA p down to below 1e-16.

**Trabelsi, Chaabane & Ben-Hur 2019** (83 ChIP-seq and 31 CLIP datasets; median AUC)

| Setting | Results |
|---|---|
| All ChIP-seq | ECBLSTM 0.930; DeepBind 0.902 |
| More than 10k positives | ECBLSTM 0.993; DeepBind 0.967 |
| Fewer than 10k positives | DeepBind-E 0.879; DeepBind 0.864 |

On small datasets the architecture complexity gains only about 0.015.

**Variant effects** (Boytsov et al. 2022, *F1000Research*; SNP-SELEX, 129 TFs, 816,594 TF–SNP pairs)

| Comparison | Result |
|---|---|
| Best-chosen CIS-BP PWM vs ΔPWM (Yan et al.) | r = 0.531 vs 0.534 |
| Strong hits only (P < 1e-4) | r ≈ 0.83 |
| Top-bound TF per SNP: PWM vs deltaSVM | 0.711 vs \|r\| 0.777 |
| TFs with AUPRC > 0.75 | 72/129 (vs 24/129 for Yan et al.'s ΔPWM) |
| TFs where the PWM beats deltaSVM | 34/129 |

**ENCODE-DREAM reanalysis 2024** (12 TF–cell-type pairs; positive rate about 1%, so random AUPR ≈ 0.02)

| Model | AUPR |
|---|---|
| Best new model (TF + transformer) | 0.505 |
| J-Team (original top team) | 0.48 |
| FactorNet | 0.46 |
| autosome.ru | 0.46 |
| Unseen TFs (general model) | 0.36 |

- Removing ATAC features costs 0.12 AUPR; removing motif features costs 0.07.
- Per-TF AUPR ranges from 0.31 to 0.75.

### 10.2 Motif discovery tools

**STREME** (Bailey 2021). 40 ENCODE K562 datasets with SELEX reference motifs; 100-bp sequences; set sizes from 1,233 to 56,058. Accuracy is −log10 of the Tomtom p-value against the reference.

| Metric | STREME | HOMER | MEME | DREME |
|---|---|---|---|---|
| Datasets with accuracy ≥ 5 | 82.5% | ≥ 70% | ≥ 70% | lower |
| Datasets with accuracy ≥ 9 | 32.5% | 12.5% | about half of STREME's rate | — |
| Head-to-head wins against HOMER | 29/40 | — | — | — |
| Speed | 1× (about 1 min per motif for 10k × 100 bp; scales about n log n) | about 10× slower | about 10× slower at widths ≤ 30 | about 10× slower |

- **Small n (10–1,000 sequences): MEME is best overall.** STREME is slightly below MEME at n < 100 but beats the other tools.
- On dilution down to 1% purity, STREME ≥ HOMER > the rest.
- STREME's p-values were accurate within 2× on random data.

**BaMMmotif2** (Ge et al. 2021). 427 ENCODE ChIP-seq and 164 HT-SELEX datasets; the metric is AvRec, recall averaged over TP:FP ratios from 1:1 to 100:1.

| Comparison | Result |
|---|---|
| 5th-order BaMM vs InMoDe | +13.6% median; BaMM wins 415/427 |
| 5th-order BaMM vs ChIPMunk, cross-cell-line | +32.7% |
| Cross-cell-line drop | about −0.2 AvRec for every tool (PWM tools 0.5 → 0.3) |
| ±4-bp flanks, ChIP | +3% (PWM) vs +11.4% (5th order) |
| ±4-bp flanks, SELEX | +1.3% vs +8% |
| Unmasked vs masked refinement | +12% (ChIP), +5–8% (SELEX), but co-occurring motifs leak in; masking is about 10× faster |
| CTCF hg19 scan | BaMM: 1.8M sites at FDR < 0.5%; JASPAR PWM: 67k sites at about 15% FDR |

**GRECO-BIT / Codebook** (Vorontsov et al. 2024)

- **Scale:** 394 TFs; 4,237 experiments, of which 1,460 datasets were approved; 5 platforms (ChIP-seq, GHT-SELEX, HT-SELEX, SMiLE-seq, PBM); about 220k PWMs; 16.8M performance estimates.
- **Tools compared:** Dimont, ChIPMunk, ProBound, MEME, STREME, HOMER, Autoseed, ExplaiNN, RCade and gkmSVM.

| Platform | Best tool within platform | Best tool across platforms |
|---|---|---|
| ChIP-seq | Dimont | Dimont |
| GHT-SELEX | Dimont | ChIPMunk |
| HT-SELEX | ProBound | Dimont |
| SMiLE-seq | ChIPMunk | Dimont |
| PBM | Dimont | MEME |

- Dimont produced about half of all top motifs. HOMER, STREME, ExplaiNN and gkmSVM were not top performers.
- Dimont's median auROC exceeded 0.85 on ChIP and GHT-SELEX.
- **Moving from synthetic to genomic sequences costs every method 0.1–0.2 auROC.**
- **The first reported motif was the best in about 75% of cases (90–95% for Dimont and ChIPMunk).**
- Motif length, IC and GC content correlate only weakly with performance (ρ ≤ 0.2).
- **Dimont's lead depends on sum-occupancy scoring and largely vanishes under best-hit scoring.**
- Archipelago, a random forest (depth 6, 100 trees) over 2–4 PWMs, gave significant gains over the best single PWM [fig-only].

**IBIS challenge 2025** (40 TFs, 372 experiments, 5 platforms; leaderboard of 28 teams and 1,905 submissions; final of 19 teams on 30 TFs; tracks A2G, artificial → genomic, and G2A, genomic → artificial)

| Discipline | Winner | Method |
|---|---|---|
| PWM track | Bench Pressers | **ZMotif (CNN-based discovery) + STREME / MEME** |
| Any-architecture, A2G | Biology Impostor | CNN / CNN-RNN |
| Any-architecture, G2A | mj | dilated temporal conv net |

- Participant PWMs were on par with, or slightly below, the best GRECO-BIT single PWM (MEX1).
- The best any-architecture model beat the best PWM only "by a limited margin".
- Many deep-learning entries were worse than the best PWMs.
- **gLMs (DNABERT-2, GENA-LM, NT) were far below supervised CNNs and at or below the top PWM; none reached the final.**
- **Better TF-binding-site classification did not translate into better rSNP prediction.** LegNet-Max, added after the challenge, is better on binding sites but lags on rSNPs.
- All scores are [fig-only] or normalised.

**TF-MoDISco vs HOMER / MEME-ChIP / DiChIPMunk** (unpublished Kundaje-lab notebook, SPI1; preliminary)

- The score is −log10 of the Tomtom p-value against the most reproducible MoDISco motif.

| Tool | Score range across tasks |
|---|---|
| MEME-ChIP | 20.6–25.1 |
| HOMER | 12.2–16.9 |
| DiChIPMunk | 13.4–15.9 |

- Many HOMER and DiChIPMunk motifs had no MoDISco counterpart.
- **No peer-reviewed TF-MoDISco vs HOMER/MEME benchmark was found.**

### 10.3 Genomic LMs vs simple baselines

**DART-Eval** (Patel et al., NeurIPS 2024)

| Task | Best gLM | Simple or supervised baseline |
|---|---|---|
| cCRE vs shuffled background, fine-tuned (paired) | NT 0.976 | small CNN 0.932 (zero-shot Caduceus 0.971) |
| Cell-type-specific cCREs, probed (accuracy / K562 AUROC) | HyenaDNA 0.587 / 0.799 | ChromBPNet-like CNN **0.667 / 0.848** |
| Chromatin activity, Spearman (GM12878 / K562) | NT fine-tuned 0.515 / 0.583 | ChromBPNet **0.540 / 0.574** |
| caQTL effect, fine-tuned r / AUROC | Caduceus 0.471 / 0.650 | ChromBPNet **0.671 / 0.772** |
| dsQTL effect, fine-tuned r / AUROC | NT 0.557 / 0.697 | ChromBPNet **0.738 / 0.892** |
| Motif sensitivity from embeddings | 40–60% per motif | FIMO motif counts give the best embedding |

**BEND** (Marin et al., ICLR 2024). Selected rows; metrics are MCC for gene finding, AUPRC for enhancers and AUROC elsewhere.

| Model | Gene finding | Enhancer | Chromatin accessibility | Histone | CpG methylation | Variant effect (expression) |
|---|---|---|---|---|---|---|
| Expert (AUGUSTUS / Enformer / Basset / DeepSEA) | 0.80 | 0.07 | 0.85 | 0.74 | 0.93 | 0.70 |
| One-hot CNN | 0.00 | 0.03 | 0.75 | 0.76 | 0.84 | — |
| NT-MS | 0.68 | 0.06 | 0.79 | 0.78 | 0.92 | 0.54 |
| DNABERT-2 | 0.43 | 0.03 | 0.81 | 0.78 | 0.90 | 0.49 |
| HyenaDNA (large) | 0.35 | 0.03 | 0.84 | 0.76 | 0.91 | 0.51 |

**GUE** (DNABERT-2; MCC × 100). Transcription-factor rows only.

| Model | Params | GUE average | TF-H (5 tasks) | TF-M (5 tasks) |
|---|---|---|---|---|
| DNABERT (3-mer) | 86M | 61.62 | 68.0 / 70.9 / 60.5 / 53.0 / 69.8 | 42.3 / 79.1 / 69.9 / 55.4 / 42.0 |
| NT-2500M-multi | 2.5B | 66.93 | 66.6 / 70.3 / 58.7 / 51.7 / 69.3 | 63.3 / 83.8 / 71.5 / 69.4 / 47.1 |
| DNABERT-2 | 117M | 66.80 | 72.0 / 76.1 / 66.5 / 58.5 / 77.4 | 56.8 / 84.8 / 79.3 / 66.5 / 52.7 |

GUE includes no non-pretrained CNN baseline.

**Vishniakov et al. 2024 (pretrained vs randomly initialised)**

- A randomly initialised Caduceus (8M parameters) reached about 0.62 MCC on the NT benchmark, beating NT-500M, NTv2-50M and GENA-LM.
- On GUE the random initialisation scored **+0.114 MCC** higher.
- Pretraining gains are usually within 2–3%. The exception is very low label regimes: with 1% of labels, pretraining added about 0.25 MCC.
- Randomly initialised models with a tuned tokenizer and width beat the pretrained models in 5 of 7 cases.
- ClinVar log-likelihood-ratio AUROC was 0.345–0.536, i.e. chance.

**Lindsey et al. 2025** (44 tasks, mean MCC): Caduceus 0.674 > DNABERT-2 0.673 > HyenaDNA 0.643 > NT 0.634 > DNABERT 0.617 > CNN 0.594. Here the CNN baseline was lightly tuned. Contrast this with DART-Eval, where a tuned ChromBPNet wins.

### 10.4 Sequence-to-function

| Study | Model | Score | Baseline or ceiling |
|---|---|---|---|
| DeepSTARR (Drosophila S2; 242k training sequences of 249 bp) | CNN | PCC 0.68 (developmental) / 0.74 (housekeeping) | Replicate ceiling 0.73 / 0.76; gkm-SVM and lasso-on-motif-counts clearly lower [fig-only] |
| ProCapNet (PRO-cap) | BPNet-like | R² 0.57 | Linear model on motif counts: R² 0.25 |
| lentiMPRA (Agarwal 2025; HepG2 / K562 / WTC11) | EnformerMPRA | r ≈ 0.84 | Lasso on biochemical + motif features: r ≈ 0.72 |
| Random Promoter DREAM 2022 (yeast; 6.74M training sequences of 80 random bp; 71k test) | LegNet (EfficientNetV2-style CNN, about 2M parameters, soft classification over 18 bins) won | Penzar: r 0.9756 (own evaluation) | Prior state of the art, the Vaishnav transformer (r ≈ 0.96 [unverified]), beaten by 19 of 28 finalists |

DREAM 2022 detail:
- The final top 5 were LegNet (CNN), BHI (conv–biLSTM–conv), Unlock_DNA (transformer with masked-nucleotide auxiliary loss), Camformers (residual CNN), and NAD and WZTR tied.
- A human-pretrained DNABERT placed 16th.
- Spread between best and worst finalists in Pearson r²: 5% on random sequences, 17.6% on native sequences, 14.6% on SNVs, 6.5% on motif perturbation, and 17.7% on motif tiling.
- **Lesson:** with millions of labelled sequences, a well-engineered fully convolutional network wins, and pretraining does not help.

### 10.5 Bacterial promoters

**Promotech** (Chevez-Guardado & Peña-Castillo 2021, *Genome Biol*). Trained on 9 species with 27,766 promoters of 40 nt.

| Evaluation | Result |
|---|---|
| Genome-wide, mean over species, RF-HOT | **AUPRC 0.14** at AUROC 0.82 |
| Genome-wide, *M. smegmatis* | AUPRC 0.27 (random: 0.0003) |
| *E. coli* RegulonDB, RF-HOT | MCC 0.54, AUPRC 0.845 |
| *E. coli*-specific tools on RegulonDB | MCC 0.498–0.763 |
| Recall of bTSSfinder and BPROM | 59% and 49% |

**Bacterial σ70 promoter tool benchmark 2020** (865 strong-evidence *E. coli* σ70 promoters, 81 bp each, vs 1,000 composition-matched negatives)

| Tool | Sensitivity |
|---|---|
| CNNProm | 94.8% |
| iPro70-FMWin | 94.5% (fewest false positives in the top group) |
| 70ProPred | 89.7% |
| iPromoter-2L | 83.8% |
| MULTiPly, bTSSfinder, BPROM, IBBP | MCC ≈ 0 |

**Promoter Calculator** (LaFleur, Hossain & Salis 2022, *Nat Commun*) [not readable; all details unverified]. A statistical-thermodynamic model of σ70 promoters with free-energy terms for the UP element, −35 box, spacer length, −10 box, extended −10 and discriminator. Its key property is extrapolation to designed sequences outside the training distribution.

**Evo / Evo 2 on prokaryotes.** Genomic LMs carry the most signal in compact prokaryotic genomes (Evo: essentiality AUROC 0.86). **No neutral head-to-head exists between a σ70 PWM or biophysical model and a gLM on genome-wide bacterial TSS calling.** That is a publishable gap.

**Note on metrics for bacterial promoters.** Balanced-set accuracies of 0.85–0.95 coexist with genome-wide AUPRC around 0.14. At a realistic prior of roughly one TSS per few kb, AUROC means little; precision at low recall is what separates tools.

---

## 11. Failure modes and pitfalls

| Failure mode | Evidence | Mitigation |
|---|---|---|
| Off-simplex gradient noise | 20–40% of positions more than 30° off the simplex | Mean-centre across bases |
| Assay-bias "motifs" (Tn5, DNase) | ChromBPNet | Factor in a bias model; use shuffled or decoy controls |
| Repeats masquerading as motifs | BPNet: more than 80% of long composites in ERV LTRs; ProCapNet: 19 "esoteric" motifs, each at least 60% transposable element | Annotate RepeatMasker overlap; report repeat-free statistics |
| Strict spacing that is really repeat-derived | 83% of the top 1% of pair distances are ERV | Spacing nulls that exclude repeats |
| Redundant or extraneous motifs | TF-MoDISco gives 7–10 extras on simulations | Sparse competitive hit calling (Fi-NeMo); overlap penalties; merge rules |
| Composition / GC driving signal | SEAM backgrounds track GC; 85% of the gLM probe signal comes from 5-mer statistics | Dinucleotide-shuffled or Markov backgrounds; k-mer-preserving shuffle tests |
| Seed instability | DEGU | Ensembles or distillation; report motif reproducibility across seeds or folds |
| Hard low-complexity TFs (SP1, STAT1, ARID3A; C2H2 zinc fingers in PBM) | Koo & Eddy; DREAM5 | Multi-mode models; zinc-finger-specific handling |
| Saturation mistaken for epistasis | SQUID | Global-epistasis output nonlinearity |
| Balanced-benchmark inflation | Promotech (AUROC 0.82 vs AUPRC 0.14 genome-wide) | Evaluate at realistic priors |
| In-distribution evaluation | GRECO-BIT (−0.1 to −0.2 auROC synthetic → genomic); BaMM (−0.2 cross-cell-line) | Always include a cross-platform or cross-condition split |
| TFBS accuracy not equal to variant-effect accuracy | IBIS | Score the two separately |
| Attention mistaken for attribution | DNABERT pipeline, DART-Eval | Require enrichment statistics |
| Max-pooling / best-hit blind to avidity | §7.2 theory; GRECO-BIT scoring result | LSE or sum-occupancy scoring |

---

## 12. Design roadmap for a frontier motif detector

*Context for the coding agent.* HASTUR is near STREME in accuracy, makes far fewer false calls and is faster. It is strongest on short-core motifs and still catching up on longer ones. Comsequent does exact GPU counting of sparse/gapped k-mers (k ≤ 5–7, span w ≤ 13+). It uses Möbius inversion to drop patterns explained by lower-order sub-patterns. The recommendations below build on those strengths.

### 12.1 Core thesis

Large models win through a handful of identifiable mechanisms (§2, §7, §9), and each one has an explicit, exact and interpretable counterpart:

| What the deep model does | Explicit counterpart to implement |
|---|---|
| Layer-1 conv + exp | PSAM/PWM energy with Boltzmann weighting |
| Sum or LSE pooling | Occupancy score $\log\sum_x e^{E(S_x)}$ over offsets and strands |
| Several filters for one TF | Mixture of 2–4 modes per motif family |
| Deeper receptive field on flanks | Higher-order (BaMM-style, order 2–5) context in and around the core, with back-off |
| Attention or dilated layers over motif channels | Pairwise spacing/orientation kernel $\omega_{a:b}(\delta,\text{strand})$ (ProBound term) |
| Learned background via LM pretraining | Variable-order (CTW / interpolated Markov) background, optionally an LM null |
| Bias factorisation | Explicit composition/assay-bias nuisance model |

### 12.2 Pipeline sketch

1. **Background model.**
   - Fit an interpolated variable-order Markov model to the control set or to shuffled positives. Use CTW, or BaMM back-off with $\alpha_k=20\cdot3^{k-1}$, up to order 5–6.
   - Use it for both expected counts and site scoring. It is the main defence against GC, low-complexity and repeat "motifs".
   - Run a test where motif calls must survive a k-mer-preserving shuffle of the positives (k = 2–5).

2. **Exact candidate enumeration (Comsequent's engine).**
   - Count all gapped patterns (k fixed positions, span ≤ w) in positives and background, on both strands, canonicalised under reverse complement.
   - Compute enrichment with a proper test: Fisher or binomial on sequence presence, plus a count-based Poisson test on occurrences.
   - **Swap the interaction filter.**
     - Möbius on log counts gives *pointwise co-information*, which is not "excess beyond sub-patterns".
     - Instead compute the **max-ent residual** $\iota^*(P)=\log c(P)-\log\hat c_{\text{IPF}}(P\mid\text{all proper sub-pattern marginals})$.
     - For contiguous patterns this reduces to an order-(k−1) Markov expectation, which is cheap. For gapped patterns run a few IPF sweeps over the sub-pattern lattice; that is GPU-friendly because patterns are small.
     - Keep Möbius co-information as a *diagnostic*, since its sign and order structure are informative. Use $\iota^*$ for selection.
   - **Significance and information floor.**
     - Report E-values against the size of the search space actually explored.
     - Prefer STREME-style held-out p-values: discover on one half, test on the other.
     - Compute the planted-motif floor $4^lP^t$ (§7.6) for the dataset's $(t,n)$, and suppress or flag calls below it.
     - This keeps HASTUR's low-false-call advantage principled.

3. **Seed → motif model (fixing the long-motif weakness).**
   - Turn significant patterns into seeds, cluster them by RC-aware alignment, and initialise PSAMs.
   - **Long motifs are often composites or dimers.** Model them as *two short cores plus a spacing/orientation kernel*: $E=E_a(S_{x})+E_b(S_{x+\delta})+\log\omega_{a:b}(\delta,\text{orient})$, as in ProBound's pairwise term.
     - This plays to the strength on short cores.
     - It explains long motifs through their sub-structure, the same logic as the Möbius/max-ent lattice but at motif level.
     - It captures soft spacing (helical ~10.5 bp periodicity; ranges below 25 bp and above 25 bp) that a monolithic long PWM blurs.
   - For genuinely long monolithic motifs (CTCF, zinc fingers), extend seeds greedily one informative column at a time. That is the leap-1 staircase, which is efficient (§7.6). Stop when the held-out likelihood gain is no longer significant.

4. **Refinement objective.**
   - Use a discriminative or likelihood fit with **occupancy (LSE/sum) pooling** over all offsets and strands. Fit a temperature $\tau$ so the data choose between best-hit ($\tau\to0$) and avidity behaviour.
   - If counts are available (SELEX rounds, ChIP read depth, CUT&RUN), use a Poisson/multinomial count likelihood rather than binary labels, as ProBound does.
   - Add **BaMM-style higher-order terms** (order 1–3 within the core, then test up to 5) plus a **±4–8 bp flank model**. This is the best-documented gain: +11–36% across studies.
   - Provide **masked and unmasked** modes. Unmasked gains up to about 12% but absorbs co-occurring motifs.

5. **Multi-mode and ensemble.**
   - Allow 2–4 modes per TF family: a PSAM mixture with shared occupancy.
   - Also provide an "Archipelago" random forest (depth 6, 100 trees) over the motif scores as a cheap ensemble ceiling.
   - Order the output so that **the first reported motif is the best**: rank by held-out performance, not training enrichment. GRECO-BIT treats this as a quality metric (75% overall; 90–95% for the leaders).

6. **Optional neural proposal engine** (the IBIS-winning recipe).
   - Train a *small* CNN on positives vs background: exponential first-layer activation, first-layer pool size ≥ 10, RC-equivariant or augment-and-conjoin, about 5–10 seeds.
   - Attribute with DeepSHAP or ISM, using dinucleotide-shuffled references and simplex correction.
   - Extract seqlets as in TF-MoDISco, and feed the resulting CWMs as *additional seeds* into steps 3–5.
   - This catches weak and composite signals that enumeration ranks low, while the final deliverable stays a calibrated motif.
   - Keep it optional and CPU-light: the CNN is small, and training sets up to roughly 100k sequences of 200 bp are fine.

7. **Instance calling.**
   - Use sparse, non-negative, competitive decomposition (Fi-NeMo style) so overlapping or redundant motifs compete for signal.
   - Report each instance's occupancy contribution, not only whether it matched.
   - Report sub-threshold low-affinity hits separately. They matter for dosage sensitivity.

8. **Syntax report** (where tools are weakest and deep models are best).
   - For each motif pair, give a co-occurrence vs spacing histogram against a spacing null built from position shuffles with repeats excluded.
   - Run a Fourier test for 10–11 bp periodicity.
   - Report orientation preference.
   - Classify each pair into the DeepSTARR modes: synergy only below 25 bp; synergy only above 25 bp; close preference with a plateau; monotone decay.
   - If a small CNN is trained, add in-silico pair-insertion super-additivity $f(AB)-f(A)-f(B)+f(\varnothing)$.

9. **Hygiene.**
   - Annotate the overlap of each motif's instances with RepeatMasker and transposable elements.
   - Flag motifs explained by composition alone.
   - Handle reverse complement throughout.
   - For bacterial promoters, model the σ factor explicitly: −35 box + spacer (15–19 bp, the ω kernel) + −10 box + extended −10 + discriminator. This is exactly the "two cores + spacing kernel" template above and matches the biophysics of the Promoter Calculator.

### 12.3 Where *not* to spend effort

- **Genomic LM embeddings or attention as discovery engines.** The evidence says they are near chance on motif identity and below k-mers on most tasks (§6, §10.3).
- If a pretrained component is used at all, use **LM likelihoods as a background or conservation prior**: weight each site by $\log q_{\text{LM}}(x_i\mid\text{ctx})-\log p_B(x_i\mid\text{ctx})$.
  - The candidates are SpeciesLM-style or GPN-style models, which do best on compact genomes.
  - For prokaryotes, Evo-family likelihoods are the most plausible useful prior, but this is untested against PWMs.
- **Deep within-core epistasis.** Gains are small (§7.8). Prioritise occupancy, multiple modes, flanks and spacing.
- **Large CNNs on small n.** Below about 10k positives, complex models lose (Zeng 2016; Trabelsi 2019). Below about 100 sequences, MEME beats STREME and HOMER. Make the small-n path a first-class mode.

### 12.4 Evaluation protocol (what "beating the frontier" should mean)

**Baselines:**
- STREME, MEME, HOMER
- BaMMmotif2
- Dimont, ChIPMunk / DiChIPMunk
- ProBound (SELEX)
- TF-MoDISco on a BPNet/ChromBPNet model
- ExplaiNN

**Datasets:**
1. **Synthetic planted motifs** (JASPAR PWMs in Markov backgrounds). Sweep IC, width, number of sequences, purity (100% down to 1%), dimers with set spacings, and many-weak-site regimes. Show that performance tracks the theoretical floor $4^lP^t$.
2. **GRECO-BIT / Codebook and IBIS datasets.** These give multi-platform data and published comparisons. Use the cross-platform splits.
3. **ENCODE ChIP-seq.** Use the STREME set of 40 K562 datasets with SELEX references, and the 427-set BaMM collection.
4. **Bacterial:** RegulonDB σ70 (genome-wide, not balanced), plus other σ factors and species, for HASTUR/Comsequent's promoter interest.

**Metrics:**
- Tomtom accuracy against the reference: share of datasets with −log10 p ≥ 5 and ≥ 9 (STREME's yardstick: 82.5% and 32.5%).
- AvRec, and AUPRC at realistic priors.
- auROC under **both** best-hit and sum-occupancy scoring.
- First-motif-is-best rate.
- Number of extraneous motifs on synthetic data, which shows the false-call advantage.
- Calibration of E-values on shuffled or random data (type I error ≤ 5% at 0.05).
- Cross-cell-line and cross-platform drop.
- Runtime and memory against n.
- Small-n curve (n = 10 to 1,000).
- Separate variant-effect evaluation (SNP-SELEX: Boytsov's best-PWM r = 0.531 to beat; deltaSVM 0.777 on the top-bound subset).

**Ablations to report:**
- Mono vs higher-order core.
- With vs without flanks.
- Best-hit vs LSE vs sum pooling.
- One mode vs 2–4 modes.
- Möbius vs max-ent residual filtering.
- With vs without CNN seeds.
- With vs without spacing kernel.

### 12.5 Concrete quick wins, ordered by expected gain per unit effort

1. **Sum/LSE occupancy scoring.** A small change, with documented gains of +0.04 AUC and a large effect on tool rankings in GRECO-BIT.
2. **Order 2–5 Markov context with back-off, plus ±4–8 bp flanks.** +11–36% in BaMM studies.
3. **Variable-order background**, to cut false calls further.
4. **2–4 modes per TF**, plus Archipelago as a cheap ensemble.
5. **Two-core + spacing-kernel model for long motifs.** This targets HASTUR's current weakness directly.
6. **Max-ent (IPF) residuals in place of raw Möbius co-information** in Comsequent.
7. **Held-out E-values plus a planted-motif information floor.**
8. **Optional small-CNN seed proposals** (exp activation, pool size ≥ 10, simplex-corrected DeepSHAP) → seqlets → seeds.
9. **Syntax report** with periodicity and spacing modes.

---

## 13. Open gaps in the literature (opportunities)

1. No published head-to-head of gLM perplexity against order-k Markov or CTW baselines on the same held-out DNA, so we don't know how much of gLM likelihood goes beyond k-mer statistics.
2. No neutral genome-wide benchmark comparing σ70 PWMs, biophysical models (Promoter Calculator) and Evo/Evo 2 likelihoods for bacterial TSS/promoter calling, scored by AUPRC at realistic priors.
3. No peer-reviewed quantitative comparison of TF-MoDISco motif recovery against HOMER, MEME or STREME on the same peaks with ground truth.
4. No rigorous study of how max-entropy interaction residuals over gapped-pattern lattices compare with PMI or Möbius co-information for motif discovery. A methods paper on this would be novel.
5. Evo 2 SAE features for TF motifs have not been benchmarked against JASPAR.
6. No study of double descent or overparameterisation in genomic CNNs for motif recovery.
7. GRECO-BIT's per-tool numbers, the Archipelago gains and the IBIS scores are reported only in figures or as normalised values. Pull them from the supplementary tables or data repositories before quoting them.

---

## 14. Bibliography (primary sources used)

**Models and tokenisation**
- Ji et al. 2021, DNABERT, *Bioinformatics*: https://www.biorxiv.org/content/10.1101/2020.09.17.301879.full.pdf
- Zhou et al. 2024, DNABERT-2 / GUE, ICLR: https://arxiv.org/abs/2306.15006
- Dalla-Torre et al. 2024, Nucleotide Transformer, *Nat Methods*: https://www.biorxiv.org/content/10.1101/2023.01.11.523679.full.pdf
- NTv3 (InstaDeep 2025): https://www.biorxiv.org/content/10.64898/2025.12.22.695963.full.pdf
- Nguyen et al. 2023, HyenaDNA, NeurIPS: https://arxiv.org/abs/2306.15794
- Schiff et al. 2024, Caduceus, ICML: https://arxiv.org/abs/2403.03234
- Nguyen et al. 2024, Evo, *Science*: https://arcinstitute.org/manuscripts/Evo.pdf
- Brixi et al. 2025/26, Evo 2: https://arcinstitute.org/manuscripts/Evo2.pdf ; https://www.goodfire.com/research/interpreting-evo-2
- Sanabria et al. 2024, GROVER, *Nat Mach Intell*: https://www.biorxiv.org/content/10.1101/2023.07.19.549677v2
- Fishman et al. 2025, GENA-LM, *NAR*: https://academic.oup.com/nar/article/53/2/gkae1310/7954523
- Cornman et al. 2025, gLM2 / OMG, ICLR: https://huggingface.co/tattabio/gLM2_650M
- Benegas et al. 2023, GPN, *PNAS*: https://www.biorxiv.org/content/10.1101/2022.08.22.504706.full.pdf
- Benegas et al. 2025, GPN-MSA, *Nat Biotechnol*
- Karollus et al. 2024, SpeciesLM, *Genome Biol*: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-024-03221-x
- AIDO.DNA 2024: https://www.biorxiv.org/content/10.1101/2024.12.01.625444.full.pdf
- Wu et al. 2025, GENERator: https://arxiv.org/abs/2502.07272
- Li et al. 2024, VQDNA, ICML: https://arxiv.org/abs/2405.10812
- Qiao et al. 2024, MxDNA, NeurIPS: https://arxiv.org/abs/2412.13716
- dnaHNet 2026: https://arxiv.org/pdf/2602.10603
- Lindsey et al. 2025, *Bioinformatics*: https://pmc.ncbi.nlm.nih.gov/articles/PMC12453675
- Gong et al. 2025: https://arxiv.org/abs/2507.15087
- Avsec et al. 2021, Enformer, *Nat Methods*; Linder et al. 2025, Borzoi, *Nat Genet*: https://www.biorxiv.org/content/10.1101/2023.08.30.555582.full.pdf
- Avsec et al. 2026, AlphaGenome, *Nature*: https://deepmind.google/discover/blog/alphagenome-ai-for-better-understanding-the-genome/

**Representations and benchmarks**
- Tang, Somia, Yu & Koo 2025, *Genome Biol*: https://link.springer.com/article/10.1186/s13059-025-03674-8
- Civale, Bagdanov & Magi 2026, NeurIPS: https://neurips.cc/virtual/2026/poster/139816
- Datta et al. 2026: https://www.alphaxiv.org/abs/2608.05329
- Patel et al. 2024, DART-Eval, NeurIPS D&B: https://arxiv.org/abs/2412.05430
- Marin et al. 2024, BEND, ICLR: https://arxiv.org/abs/2311.12570
- Vishniakov et al. 2024: https://www.biorxiv.org/content/10.1101/2024.12.18.628606.full.pdf
- Weirauch et al. 2013, DREAM5, *Nat Biotechnol*: https://pmc.ncbi.nlm.nih.gov/articles/PMC3687085
- Alipanahi et al. 2015, DeepBind, *Nat Biotechnol*
- Zeng et al. 2016, *Bioinformatics*: https://academic.oup.com/bioinformatics/article/32/12/i121/2240609
- Trabelsi, Chaabane & Ben-Hur 2019: https://arxiv.org/abs/1901.10526
- Boytsov et al. 2022, *F1000Research*: https://f1000research.com/articles/11-33/v2
- ENCODE-DREAM reanalysis 2024: https://www.biorxiv.org/content/10.1101/2024.11.08.622635.full.pdf
- Bailey 2021, STREME, *Bioinformatics*: https://www.biorxiv.org/content/10.1101/2020.11.23.394619.full.pdf
- Ge et al. 2021, BaMMmotif2: https://www.biorxiv.org/content/10.1101/2020.07.12.197053.full.pdf
- Siebert & Söding 2016, BaMM, *NAR*: https://epub.ub.uni-muenchen.de/37677/1/10.1093_nar_gkw521.pdf
- Vorontsov et al. 2024, GRECO-BIT / Codebook: https://www.biorxiv.org/content/10.1101/2024.11.11.619379.full.pdf
- IBIS challenge 2025: https://www.biorxiv.org/content/10.1101/2025.11.16.688692.full.pdf
- de Almeida et al. 2022, DeepSTARR, *Nat Genet*: https://www.biorxiv.org/content/10.1101/2021.10.05.463203v1.full
- Rafi et al. 2024, Random Promoter DREAM, *Nat Biotechnol*: https://www.biorxiv.org/content/10.1101/2023.04.26.538471v2
- Penzar et al. 2023, LegNet: https://www.biorxiv.org/content/10.1101/2022.12.22.521582.full.pdf
- Agarwal et al. 2025, lentiMPRA, *Nature*: https://www.biorxiv.org/content/10.1101/2023.03.05.531189.full.pdf
- Chevez-Guardado & Peña-Castillo 2021, Promotech, *Genome Biol*: https://link.springer.com/article/10.1186/s13059-021-02514-9
- Bacterial promoter tool benchmark 2020: https://www.biorxiv.org/content/10.1101/2020.05.05.079335.full.pdf
- LaFleur, Hossain & Salis 2022, Promoter Calculator, *Nat Commun*: https://pmc.ncbi.nlm.nih.gov/articles/PMC9440211 (not read)

**Interpretation**
- Shrikumar et al., TF-MoDISco: https://arxiv.org/abs/1811.00416 ; tfmodisco-lite: https://github.com/jmschrei/tfmodisco-lite
- Fi-NeMo: https://pypi.org/project/finemo/
- Avsec et al. 2021, BPNet, *Nat Genet*: https://www.biorxiv.org/content/10.1101/737981v1.full-text
- Pampari et al. 2024, ChromBPNet: https://www.biorxiv.org/content/10.1101/2024.12.25.630221.full
- Cochran et al. 2024, ProCapNet: https://www.biorxiv.org/content/10.1101/2024.05.28.596138.full
- Koo & Eddy 2019, *PLoS Comput Biol*: https://www.biorxiv.org/content/10.1101/362756v2.full
- Koo & Ploenzke 2021, *Nat Mach Intell*: https://www.biorxiv.org/content/10.1101/2020.06.14.150706v2.full-text
- Majdandzic, Rajakumar & Koo 2023, *Genome Biol*: https://www.biorxiv.org/content/10.1101/2022.04.29.490102v2.full
- Novakovsky et al. 2023, ExplaiNN, *Genome Biol*: https://www.biorxiv.org/content/10.1101/2022.05.20.492818v2.full
- Tseng et al. 2024, ARGMINN: https://arxiv.org/abs/2410.06211
- Koo et al. 2021, GIA, *PLoS Comput Biol*: https://www.biorxiv.org/content/10.1101/2020.09.08.288068.full.pdf
- Greenside et al. 2018, DFIM, *Bioinformatics*: https://www.biorxiv.org/content/10.1101/302711v1.full
- Seitz et al. 2024, SQUID, *Nat Mach Intell*: https://www.biorxiv.org/content/10.1101/2023.11.14.567120.full
- Seitz et al. 2025, SEAM: https://www.biorxiv.org/content/10.1101/2025.10.07.681052.full
- Toneyan & Koo 2024, CREME, *Nat Genet*: https://repository.cshl.edu/id/eprint/40963/
- Zhou et al. 2025, DEGU: https://www.biorxiv.org/content/10.1101/2024.11.13.623485.full
- Schreiber & Singh 2022, Yuzu: https://pubmed.ncbi.nlm.nih.gov/35678521/
- Decode-gLM 2025: https://www.biorxiv.org/content/10.1101/2025.10.31.685860v4.full-text
- Naqvi et al. 2025, *Cell*: https://pubmed.ncbi.nlm.nih.gov/40020686/

**Theory**
- Rube et al. 2022, ProBound, *Nat Biotechnol*: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9546773/
- Zhao, Granas & Stormo 2009, BEEML: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2777355/
- Ghandi et al. 2014, gkm-SVM, *PLoS Comput Biol*: https://pmc.ncbi.nlm.nih.gov/articles/PMC4102394
- Mallet & Vert 2021, RC-equivariant networks, NeurIPS: https://proceedings.nips.cc/paper_files/paper/2021/file/706608cfdbcc1886bb7eea5513f90133-Paper.pdf
- Zhou, Shrikumar & Kundaje 2022, MLCB: https://proceedings.mlr.press/v165/zhou22a.html
- Edelman et al. 2024, statistical induction heads, NeurIPS: https://arxiv.org/abs/2402.11004
- Akyürek et al. 2024, in-context language learning, ICML: https://arxiv.org/abs/2401.12973
- Zhou, Tian & Diggavi 2024, transformers and VOMC/CTW: https://arxiv.org/abs/2410.05493
- Nichani, Damian & Lee 2024, ICML: https://arxiv.org/abs/2402.14735
- Bietti et al. 2023, NeurIPS: https://arxiv.org/abs/2306.00802
- Hahn 2020, *TACL*: https://aclanthology.org/2020.tacl-1.11
- Dao & Gu 2024, Mamba-2 / SSD: https://arxiv.org/abs/2405.21060
- Barak et al. 2022, sparse parity, NeurIPS: https://arxiv.org/abs/2207.08799
- Abbe, Boix-Adserà & Misiakiewicz 2023, leap complexity, COLT: https://arxiv.org/abs/2302.11055
- Abbe et al. 2023, min-degree bias: https://arxiv.org/abs/2301.13105
- Brookes, Aghazadeh & Listgarten 2022, *PNAS*: https://www.biorxiv.org/content/10.1101/2021.05.24.445506.full.pdf
- Pevzner & Sze 2000, ISMB; Buhler & Tompa 2002, *J Comput Biol*; Schneider et al. 1986, *J Mol Biol*; Leslie et al. 2002/2004 (spectrum/mismatch kernels) [cited from background knowledge]
