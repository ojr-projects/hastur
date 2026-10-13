# HASTUR ROADMAP v2: real-data evidence plan

*Started 2026-10-12. This file is the forward plan. `ROADMAP.md` keeps the dated development log (§5) up to the
FR1 confirmation and is the record of how we got here. Lessons that carry over are in `docs/DEV_INSIGHTS.md`.*

---

> **2026-10-12: E-series stopped; focus is PROMETHEUS only** (user decision, after HASTUR-first quick looks).
> All quick looks ran HASTUR first; comparators were added only where the comparison looked winnable. The plain
> conclusion: HASTUR is calibrated everywhere tested, but has no sensitivity advantage on small sets.
>
> | experiment | data | HASTUR result | comparator result | verdict |
> |---|---|---|---|---|
> | E1 | *E. coli* iModulons, 32 evaluable modules | 8/32 native, 0.044 false calls/run | STREME 30/32 native but 9.6 false calls/run; level with HASTUR at matched rates | calibrated, no frontier win |
> | E1, known-motif arm | same | own known-motif test 16/32 | SEA 29/32 at 0.059 false calls/run | SEA owns known motifs |
> | E1b | *P. putida* | 0.032 false calls/run; 12/48 modules called, including the FNR box | not run | calibration only |
> | E6 | 100 null screens | per-run calibration predicts screen counts (4.7 per screen) | not run | p_adj is anti-conservative in the deep tail, so BH over runs fails (46% of null screens) |
> | E5 | *P. falciparum*, 13% GC nulls | 0.04 false calls/run | STREME calibrated at N ≥ 50; fails only at N = 25 (S-mode) | composition does not separate the tools; small N does |
> | E2 | 23 TF knockouts | own motif in 1/14 scoreable contrasts | not run | no |
>
> Records: `splits/protocol_E1*.md`, `splits/protocol_E6.md`, `bench/{e2_make,promoter_nulls_make,imodulon_other_make,known_motif_test}.py`.
> **Remaining direction:** `docs/PROMETHEUS_PLAN.md`. **Development paused 2026-10-12** after MA1/MA2 and the deep-model
> literature review: no large expected gain was identified (see PROMETHEUS_PLAN §6).

## 1. Where HASTUR stands

**Method (frozen F2, `splits/freeze_v18.md`).**
- **Detection:** the gapped k-of-w lattice and the L-mer table, both double-centred against each set's own Markov
  chain.
- **Ranking:** v8's statistics rank the seeds.
- **Significance:** exact tests on the seeds, calibrated by **permutation max-T** (20 label permutations, Gumbel tail,
  no genomic control on the tested statistics).
- **Representation:** frame or EM, chosen by held-out evidence, then a site refit.
- **Implementations:**
  - `proto/hlat.py`: the Python reference;
  - `rust/hastur_core`: production, byte-identical outputs on every validation unit.

**Evidence** (STREME = permissive STREME scored by its E-value; all at matched precision unless noted;
`results/`):

| claim | evidence | strength |
|---|---|---|
| Calibrated on exchangeable real sequence | 0.00–0.15 false calls per motif-free run on fresh human and bacterial nulls (STREME 0.00, old v8 0.5–1.8) | registered confirmation |
| Calibrated under composition shift | simulated CpG shift: 0.05 per run vs STREME 3.45 | registered confirmation |
| Calibrated on small sets | random E. coli promoter sets (10–206): 0.08 vs STREME 7.65; fresh small-N panel 0.05 | dev + verified-fresh panel |
| Better frontier with composition-shifted nulls | +0.24 [+0.20, +0.28] sensitivity at 0.1 false calls per run (ChIP +0.24, SELEX +0.25); +0.27 repeat-masked | registered FR1 confirmation |
| Level with STREME on exchangeable data | +0.01 overall; SELEX +0.10, unmasked ChIP −0.06, masked ChIP −0.02 | post-hoc decomposition |
| Long motifs | 0.31 vs 0.09 at matched precision (> 10 informative columns) | confirmation, secondary |
| Fast | 0.8–7 s per full run on 12 cores; 15–40× STREME per core, 100–330× wall-clock | measured |

**Known weak spots:**
- GC-selected human windows: both tools are miscalibrated, even repeat-masked.
- Very small sets (about 25 sequences) are calibrated but conservative.
- The ChIP motif statement (P5) is slightly behind STREME unmasked.
- STREME is the only comparator so far.

**Positioning:** *a calibrated, very fast motif finder for high-throughput and small-set discovery, matching STREME's
accuracy.* The case for it is trustworthy output at screen scale: many small or composition-shifted contrasts, where
false calls from other tools pile up. It is not "finds more motifs in one clean ChIP experiment".

---

## 2. Rules carried forward (from DEV_INSIGHTS §3.7 and §8)

1. **Register every experiment before running it** (`splits/protocol_*.md`, sha256 in `SHA256SUMS`): question, data,
   truth, comparators, readouts, decision rule.
2. **Oracle-test every truth source before registering it.** If the truth itself does not score as a hit, and
   shuffled truths do not score as misses, there is no verdict.
3. **Nulls are real, exchangeable and fresh.**
   - Real sequence at every deployment size, with tiled (non-overlapping) windows.
   - Tuple seeds, plus a measured freshness check against earlier panels.
   - A composition-shifted null wherever the application has one.
4. **Report the frontier decomposed by null family,** together with each tool's native operating point. No
   headline number without its exchangeable-only counterpart.
5. **Comparators are fixed per experiment, run on the same inputs and machine,** at both native and matched
   thresholds. Runtime is reported.
6. **HASTUR runs come from the Rust binary.** The Python reference re-checks identity on a sample of units.
7. **No real-data CRISPR analyses** (standing user instruction).

---

## 3. Experiment plan

Each experiment lists which measured advantage it plays to. **A** = calibration under composition shift;
**S** = small-set calibration; **X** = speed and scale; **L** = long motifs; **E** = strong-selection data (SELEX-like).

### E1. E. coli iModulon screen (headline; A, S, X)
- **Question:** run on every co-regulated module, does HASTUR recover each module's regulator motif without the false
  calls that make other tools' output need curation?
- **Data:** iModulonDB, E. coli PRECISE-1K (about 200 iModulons with curated regulator annotations). Upstream regions
  of each module's operons, against all other promoters (the D1 promoter construction).
- **Nulls:** random gene sets matched to each module's size and operon structure; at least 20 per module.
- **Truth:** each module's annotated regulator, with its site matrix from RegulonDB. The readout is site recovery
  (oracle-validated in D1). Modules with no annotated sequence-specific regulator are scored for false calls.
- **Comparators:**
  - STREME (default and permissive);
  - HOMER de novo (fasta mode, matched background);
  - MEME (ZOOPS, small sets);
  - XSTREME;
  - for the known-motif arm (E1b): SEA and HOMER known.
- **Readouts:** regulator recovery at native thresholds and at matched false calls per null run; screen-level FDR;
  share of "no motif" calls on modules without a cis regulator; total runtime.
- **Extension E1b:** the same modules in a second organism from iModulonDB (B. subtilis, or P. aeruginosa if curated
  enough).

### E2. Direct vs indirect targets in TF perturbation data (S, A)
- **Question:** given the genes that respond to deleting or overexpressing a TF, does a calibrated motif call in their
  upstream regions separate direct from indirect regulation?
- **Data:** E. coli TF knockout and overexpression transcriptomes (PRECISE-1K perturbation samples, GEO knockout
  series). Differentially expressed genes per perturbation.
- **Truth:** direct binding from RegulonDB high-throughput datasets (ChIP-exo, gSELEX, DAP-seq), oracle-tested.
- **Comparators:** STREME, HOMER, MEME; SEA with the perturbed TF's known motif as a positive control.
- **Readouts:**
  - is the perturbed TF's motif found (sensitivity);
  - share of perturbations with indirect-only responses correctly called "no motif";
  - false calls on matched random sets.

### E3. Pangenome-scale phylogenetic footprinting (X, S, L)
- **Question:** run on the upstream regions of every orthologous gene group across many genomes, does HASTUR recover
  known regulon motifs at scale, with a usable global FDR?
- **Data:** a clade with many genomes (e.g. Enterobacterales or Bacillales; 100–500 RefSeq genomes). Orthogroups from
  an orthology tool, and upstream regions per genome. Includes the 21 local genomes as a pilot.
- **Truth:** RegPrecise comparative regulons (curated motifs and regulons across many bacteria), oracle-tested.
- **Comparators:** STREME and HOMER on a stratified subsample (full runs are infeasible for them; that is part of the
  point); MEME on the smallest orthogroups.
- **Readouts:** regulon motif recovery; false calls on shuffled-orthogroup nulls; runtime at full scale (HASTUR) vs
  extrapolated (comparators).

### E4. Gene-set annotation in other organisms (S, A)
- **Question:** which gene sets carry a cis-regulatory signature, in yeast and plant co-expression clusters and GO
  sets?
- **Data:**
  - yeast co-expression modules and GO sets (promoters from SGD);
  - Arabidopsis co-expression clusters (promoters from TAIR/Araport).
- **Truth:** yeast regulator assignments from YEASTRACT, with CIS-BP or JASPAR motifs; plant motifs from DAP-seq
  (O'Malley et al.), oracle-tested.
- **Comparators:** STREME, HOMER, MEME; AME or SEA for the known-motif arm.
- **Readouts:** as E1, per organism.

### E5. Extreme-composition genomes (A, S)
- **Question:** does calibration under composition shift translate into usable motif discovery in genomes far from
  average composition, where the other tools' false calls should be worst?
- **Data:**
  - Plasmodium falciparum (about 80% AT): stage-specific co-expression clusters; ApiAP2 TFs.
  - A high-GC actinomycete (Streptomyces, about 72% GC): co-expression modules or regulons.
- **Truth:**
  - Plasmodium: ApiAP2 binding motifs (protein-binding-microarray motifs, CIS-BP) and published ChIP-seq where
    available;
  - Streptomyces: curated regulons (RegPrecise / literature).
- **Comparators:** STREME, HOMER, MEME.
- **Readouts:** false calls on random gene sets of matched size (the calibration test that matters here); recovery of
  known ApiAP2 motifs.

### E6. Screen-level FDR (D6; X, S)
- **Question:** because each run is calibrated, can a global FDR over thousands of runs be estimated and controlled?
- **Data:** the E1, E3 and E4 screens, plus fully null screens (thousands of random sets).
- **Readouts:** estimated vs realised FDR at screen scale (realised from the null screens); the same exercise with
  STREME's E-values and HOMER's p-values, showing whether their per-run statistics support a screen-level FDR.

### E7. Many-contrast differential accessibility (A, X)
- **Question:** in many cell-type or condition contrasts, where every background is composition-shifted, does HASTUR
  give calibrated differential motifs where others pile up false calls?
- **Data:** single-cell ATAC peaks per cluster against the rest (e.g. 10x PBMC multiome or a public atlas). Each
  contrast = cluster-specific peaks vs the other clusters' peaks. Plus tiled-window nulls from the same clusters.
- **Truth:** cluster identities with known master regulators (e.g. SPI1 in myeloid, PAX5 in B cells, TCF7 in T
  cells), oracle-tested. Per-contrast false calls from label-swapped contrasts.
- **Comparators:** HOMER (the field default here, with GC-matched background), STREME, chromVAR (deviation scores, as
  the non-discovery baseline).
- **Readouts:** master-regulator recovery per contrast; false calls on label-swapped contrasts; total runtime across
  all contrasts.

### E8. Low-input binding assays (S, X)
- **Question:** at the peak counts of CUT&RUN / CUT&Tag or low-cell ChIP (100–1,000 peaks), does calibration hold with
  usable sensitivity?
- **Data:** public CUT&RUN / CUT&Tag sets for TFs with known motifs (ENCODE / 4DN / GEO), and ChIP subsampled to the
  same counts. Use matched repeat-masked backgrounds.
- **Truth:** the TF's known matrix (oracle-tested).
- **Comparators:** STREME, HOMER, MEME.
- **Readouts:** frontier and native operating point by peak count; replicate reproducibility (calibration should make
  calls reproduce).

### E9. Calibrated "nothing found" (D2; S, A)
- **Question:** does a calibrated empty result track failed or non-specific experiments?
- **Data:** ENCODE TF ChIP-seq for sequence-specific factors, with and without audit flags (low enrichment, poor
  replicate concordance).
- **Readouts:** the rate of empty results in flagged vs unflagged experiments, per tool. The confounder (factors
  without a direct motif) is controlled by restricting to factors with strong known motifs.

### E10. Long and complex motifs (L)
- **Question:** is the long-motif advantage real on curated long sites?
- **Data:**
  - KRAB zinc-finger ChIP-exo sets (long motifs; repeat-derived; report masked and unmasked);
  - bacterial 16–22 bp operator sites (from E1 and E3);
  - CTCF and REST.
- **Comparators:** STREME, MEME (designed for long motifs), HOMER with -len up to 20.
- **Readouts:** accuracy by motif length (P5 and PWM similarity at matched precision), and the representation gap
  (column accuracy).

### E11. Auditing deep-model motifs (D7)
- **Question:** which motifs reported by attribution methods (TF-MoDISco on BPNet / ChromBPNet models) carry
  calibrated enrichment in the same peaks?
- **Data:** public BPNet / ChromBPNet models with TF-MoDISco reports, and their peaks.
- **Comparators:** TF-MoDISco's own motif set (the object audited), plus STREME and HOMER on the same peaks.
- **Readouts:** the share of deep-model motifs with calibrated enrichment; motifs found de novo but missed by
  attribution, and the reverse.

### E12. A calibrated motif atlas (D4; X)
- **Question:** can a calibrated motif-presence matrix across every ENCODE TF ChIP set be built cheaply, and does it
  recover cofactor modules?
- **Data:** all ENCODE TF ChIP-seq with conservative IDR peaks.
- **Comparators:** STREME on a subsample (full runs are expensive); HOMER known-motif enrichment as the reference
  atlas.
- **Readouts:** cofactor-module recovery against known complexes; total compute.

---

## 4. Engineering work these experiments need

1. **Batch mode** (E1–E7): many foreground sets against one shared background in a single invocation. Background
   counts and Markov fits are computed once. This is the main speed lever for screens.
2. **Known-motif test mode** (E1b, E2, E4, small-N power): test a motif library with the same permutation max-T
   calibration. The family is hundreds of motifs instead of about 10⁶ cells, so the thresholds and the small-N power
   are far better. This answers the N = 25 limit.
3. **Background builders:** promoter, operon and orthogroup upstream sets; GC- and repeat-matched genomic
   backgrounds; label-swapped contrast nulls.
4. **Repeat masking** (soft-masked input → N) as a documented option; masking is an input choice.
5. **Screen-level FDR utility** over per-run adjusted p-values (E6).
6. **Rust CLI and packaging:** documented flags (frozen F2 as the default), MEME-format output, a container image.
   Identity tests against the Python reference in CI.
7. **Optional speed steps that change numerics** (each validated on accuracy, not identity): smaller backgrounds for
   the frame/EM choice, precomputed background window sums, a vectorised exp, tiled counting. Together about 2–4× at
   5k peaks.

## 5. Open method questions (dev, not confirmation)

- **Composition beyond order-1 Markov** (GC-selected windows, CpG islands): matched backgrounds vs a higher-order or
  mixture background model, judged on the GC-selected nulls.
- **Representation gap on ChIP (P5):** the LN5 round-2 candidate L12 (+0.014 PWM similarity); discriminative polish.
- **Small-N power beyond known-motif mode.**

## 6. Order of work and compute

| step | what | depends on | compute (box) |
|---|---|---|---|
| 1 | batch mode + known-motif mode + background builders (engineering) | — | none |
| 2 | **E1 iModulon screen** (register, oracle, run, read) | 1 | HASTUR minutes; comparators about 1–2 box-days |
| 3 | E6 screen-level FDR on E1, plus null screens | 2 | small |
| 4 | E2 perturbation screen | 1 | about 1 box-day |
| 5 | E5 extreme-composition genomes | 1 | about 1 box-day |
| 6 | E3 pangenome footprinting (pilot on the 21 local genomes, then a clade) | 1 | HASTUR hours; comparators subsampled |
| 7 | E7 single-cell ATAC contrasts, E8 low-input assays | 3 | 1–2 box-days each |
| 8 | E4, E9–E12 | as capacity allows | |

**Decision points:**
- **After E1:** if HASTUR's advantage at the native threshold on small bacterial sets does not hold against HOMER
  and MEME as well, re-scope the positioning before E2–E5.
- **After E6:** if screen-level FDR is not well estimated, the "trust at scale" claim is limited to per-run
  calibration.
