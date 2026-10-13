# Registered protocol E1: E. coli iModulon screen (registered 2026-10-12, after the oracle test, before any scored tool run)

**Question (ROADMAP_v2 E1).** Run on every co-regulated module from an independent source (expression-derived
iModulons, not RegulonDB regulons), does HASTUR F2 recover each module's regulator motif, and how many false calls do
it and the comparators make on the way? This extends D1, which defined regulons from RegulonDB itself.

**Disclosure.**
- D1 (RegulonDB regulons, same promoter population) was a readout during LN7 development. E1 modules share regulators
  with D1, but not region sets: the best Jaccard of any evaluable module with any D1 regulon is ≤ 0.46.
- E1 nulls are fresh: there are 1,880, the largest Jaccard with a same-size D1 null is 0.11, and none is identical.
- A pipeline and timing test ran every tool on three modules (SoxS, Crp-1, Fatty_Acid) in a scratch directory. Those
  outputs are discarded; every unit is rerun from scratch.

## Data (sha256 in `~/data/imodulon/raw/SHA256SUMS`, recorded below)
- **PRECISE-1K** (github.com/SBRG/precise1k, commit 4829b83ace51ee7980de6950077f1cf8f3746f88, fetched 2026-10-12):
  `imodulon_table.csv`, `M.csv`, `precise1k.json.gz` (thresholds), `gene_info.csv`.
- **RegulonDB v12 GFF3:** `PromoterSet.gff3`, `RISet.gff3` (as D1) and `TUSet.gff3`.
- **Genome:** NC_000913.3, as in D1.

## Units (`bench/imodulon_make.py`)
- **Promoter regions:** the D1 construction, 1,942 regions.
- **Module membership:** genes with |M| above the module's PRECISE-1K threshold. Each gene maps to the regions that
  hold the TSS of every RegulonDB transcription unit containing it (on its strand). Median gene coverage is 0.80.
- **Positive unit:** fg = the module's regions, bg = every other region. All 94 modules with ≥ 5 regions are included,
  out of 201.
- **Nulls:** 20 per module (1,880 runs), each a random set of regions of the same count, with bg = every other region.
  Seeds come from sha256 of `E1|<module>`.
- **Truths:** each annotated regulator token (split on + and /) that is a RegulonDB TF with ≥ 8 sites at its modal
  length. The truth matrix is built as in D1.

## Oracle test (`bench/imodulon_eval.py oracle`; `e1_oracle.txt`)
- **Readout:** D1's site recovery (oracle-validated in D1), over the module's fg regions that hold a site of that TF.
  Pairs with ≥ 4 such regions are tested: 49 pairs.
- **Chance threshold τ:** the larger of the 95th percentile of 20 column-shuffled truths and the 95th percentile of
  every other E1 truth TF.
- **Result:** truth median 0.91, against chance 95th percentiles of 0.50 (other TFs) and 0.55 (shuffled).
- **Evaluable pair:** truth recovery ≥ τ + 0.15. There are 37 evaluable pairs, in **32 evaluable modules**:
  Arginine, CpxR, Cra, Crp-1, Crp-2, Cysteine-1, DNA_Damage, EvgA, FDH-N_NRA, Fatty_Acid, Fnr-1, Fnr-2, Fnr-3,
  Fur-1, Fur-2, GadXW, Galactose, GlcNAc, Gluconate, Glutamine, Magnesium, Nitrate_Nitrite, NtrC-1, NtrC-2, NtrC-3,
  Nucleoside, OxyR, Phosphate-1, Purine, Pyrimidine, Thr_Ser-1, Tyr_Trp_Phe.
- **Hit:** a call on an evaluable module whose recovery reaches τ on any evaluable pair of that module.
- **Mismatched-truth chance:** the same rule against a fixed other evaluable module.

## Tools (all on one machine, the local workstation; `bench/run_tools.py`)
- **`hlat_f2r`:** frozen F2 from the Rust binary (`rust/hastur_core`, byte-identical to `proto/hlat.py` F2). Its
  native threshold is `sig` (p_adj ≤ 0.05).
- **`streme`:** default settings, `--nmotifs 10`. Native threshold E ≤ 0.05, or S ≤ 0.05 when STREME holds out no
  test set (D1 amendment 1).
- **`streme_all`:** `--thresh 1 --nmotifs 25`, for the frontier.
- **`homer`:** `findMotifs.pl` in fasta mode, `-len 8,10,12`, frontier only.
- **`meme`:** ZOOPS, `-nmotifs 5 -minw 6 -maxw 25 -revcomp`, with an order-1 background fitted on bg. Native
  threshold E ≤ 0.05; it also enters the frontier.
- **`xstreme`:** default, with no motif database (`combined.meme`). Native threshold E ≤ 0.05.
- **Scale:** scores are on the FR1 scale (−log10 of the tool's significance).
- **Parity:** every tool runs on every unit (94 + 1,880), with the same inputs, the same machine and the same
  `MIN_FREE_GB` memory gate.

## Readouts (`bench/imodulon_eval.py tools`, then `bench/imodulon_report.py`)
1. **Native operating points** (hlat_f2r, streme, meme, xstreme):
   - sensitivity on the 32 evaluable modules;
   - false calls per null run, and the share of runs with any;
   - native calls on all 94 modules, with those on the 62 non-evaluable ones shown separately;
   - estimated screen FDR = false calls per null run × 94 / native calls on the 94 modules;
   - mismatched-truth chance.
2. **Frontier** (hlat_f2r, streme_all, homer, meme): sensitivity at 0.05, 0.1 and 0.5 false calls per null run, with
   bootstrap CIs over modules. Paired differences HASTUR − each comparator.
3. **"Nothing found" on modules without a cis truth:** native calls per module on the Genomic, Technical, Unknown and
   Biological categories. These are descriptive, since some of these modules may still carry real motifs.
4. **Cost:** total and median runtime per tool.
5. **Strata (descriptive):**
   - sensitivity by module size (5–9, 10–19, ≥ 20 regions);
   - by site length (≤ 16 vs ≥ 17 bp);
   - on the regulators that were not among D1's 21 evaluable TFs (OxyR, EvgA, GadW, GalR, GalS, GntR, TyrR).

## What may be claimed
- **"Calibrated screen":** only if HASTUR's native false calls per null run are ≤ 0.1.
- **"Lowest screen FDR at the native point":** only if HASTUR's estimated screen FDR is below that of every native
  comparator (STREME, MEME, XSTREME), with a calibrated screen as above.
- **"Better frontier on iModulons":** only against a named comparator, where HASTUR's sensitivity at 0.1 false calls
  per run exceeds that comparator's by ≥ 0.05 with the bootstrap CI excluding 0.
- **ROADMAP_v2 decision point:** if HASTUR's frontier at 0.1 is not at least level (CI including or above 0) with
  HOMER and MEME as well as STREME, the positioning is re-scoped before E2–E5.
- Otherwise results are reported as they stand. n = 32, so the CIs are wide, and that is stated with every result.
