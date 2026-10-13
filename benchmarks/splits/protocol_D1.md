# Registered protocol D1: E. coli regulon screen (registered 2026-10-11, after the oracle test, before any tool run)

**Question.** Run as a screen (every regulon, plus many motif-free random sets of the same promoter population), how
many real regulon motifs does each tool find, and how many false calls does it make on the way? This is the headline
demonstration for "calibrated, fast discovery for high-throughput use" (ROADMAP §3b, D1). No bacterial positive data
has been used in HASTUR's development; the configuration is frozen v8 (tool key `hlat_a3`).

## Data (sha256)
- RegulonDB v12 GFF3 releases, fetched 2026-10-11 from `regulondb.ccg.unam.mx/media/raw/gff3/`:
  `RISet.gff3` 2152a3e38f74bcf3529bc26116df7a1283b5e1716390431d832e6a77a90b158b,
  `PromoterSet.gff3` 50344119ecad07a24b724b9ffeab70e18324cff77413849b6ac930b8c68ca026.
- Genome NC_000913.3 (GCF_000005845.2) 53bb6a51b6e92139ced1e38f74b7938781027c52200922ff03718c2237d23bb4.

## Units (`bench/regulon_make.py`)
- Promoter regions: [TSS − 250, TSS + 50] around every RegulonDB TSS, overlapping windows merged: 1,942 regions,
  median 305 bp, 718 kb.
- Regulon of TF t: the regions holding at least half of one of t's annotated sites (RISet, riType tf-*). Positive
  unit: fg = the regulon, bg = every other region.
- Included: TFs with ≥ 10 regions and ≥ 8 sites at their modal site length: 34 TFs.
- Nulls: 10 per TF (340 runs): fg = random regions, the same count as the regulon; bg = every other region.
- Truth matrix: t's sites at their modal length, strand-oriented, +0.25 pseudocount.

## Oracle test (`bench/regulon_eval.py oracle`, results `regulon_oracle.txt`)
The truth matrix, 20 column-shuffled truths per TF and every other TF's truth were scored with every readout.
- **IUPAC readouts (contains, P5) fail the test.** Weak truths have near-empty consensuses (Lrp `W`), so contains
  passes trivially, and shuffled matrices hit 9% of the time.
- **PWM similarity fails the test.** The shuffled-truth 95th percentile is 1.0, because truths with few informative
  columns match anything, and five truths have no informative column at all.
- **Site recovery passes the test:** the fraction of regulon regions whose best-scoring window (both strands)
  overlaps an annotated site of t by at least half the shorter of the two. Truth median 0.66 vs shuffled 95th
  percentile 0.40.
- **Primary readout: site recovery.** Per-TF chance threshold τ_t = max(95th percentile of the shuffled truths,
  95th percentile of the other TFs' truths). A call is a **hit** if its site recovery ≥ τ_t.
- **Evaluable TFs:** those whose own truth reaches τ_t + 0.15. These are 21 TFs: ArgR, CRP, CpxR, Cra, CysB, FNR,
  FadR, Fis, Fur, LexA, MarA, NagC, NsrR, NtrC, PdhR, PhoB, PhoP, PurR, RcsB-BglJ, Rob, ZraR.
- **Not evaluable:** the other 13 TFs, whose curated sites do not localise from their own matrix. These are ArcA,
  CsgD, GadX, GlaR, H-NS, IHF, Lrp, Nac, NarL, NarP, OmpR, OxyR, SoxS. They run in the screen, and their calls are
  counted but not scored.

## Tools (all on one machine, the box)
- `hlat_a3`: frozen v8.
- `streme`: default settings, the native operating point, E ≤ 0.05.
- `streme_all`: `--thresh 1 --nmotifs 25`, for the frontier.
- `homer`: `-len 8,10,12`, frontier only, since HOMER has no native threshold.

Scores are on the FR1 scale (`bench/frontier.py`).

## Readouts (`bench/regulon_eval.py tools`, then `bench/regulon_report.py`)
1. **Native operating points.** HASTUR uses its Yellow Sign calls (z ≥ zcrit); STREME uses E ≤ 0.05 from its
   default run. Report:
   - sensitivity: the share of the 21 evaluable regulons with at least one hit;
   - false calls per null run, over 340 runs;
   - calls on the 13 non-evaluable regulons;
   - the estimated screen FDR: false calls per null run × 34, divided by the native calls on the 34 regulons;
   - mismatched-truth chance: the hit rate when a unit's calls are scored against a fixed other TF's sites.
2. **Frontier**, for HASTUR, permissive STREME and HOMER: sensitivity on the evaluable regulons against false calls
   per null run, with the threshold swept. Report sensitivity at 0.05, 0.1 and 0.5 false calls per run, with
   bootstrap CIs over regulons. n = 21, so the CIs are wide, and that is stated with the result.
3. **Cost:** total and median runtime per tool over the 374 runs.
4. **Secondary:** P5 and PWM similarity on the evaluable regulons (oracle-weak, so descriptive only), and
   sensitivity by site length (≤ 16 bp vs ≥ 17 bp).

## What may be claimed
- **"Calibrated screen"**: only if HASTUR's native false calls per null run are ≤ 0.1.
- **"Better frontier on bacterial regulons"**: only if HASTUR's sensitivity at 0.1 false calls per run exceeds
  STREME's by ≥ 0.05, with the CI excluding 0.
- **"Matches STREME at the native point"**: only if HASTUR's native sensitivity is within 0.10 of STREME's while
  its false calls per run are lower.
- Otherwise the result is reported as it stands. A sensitivity gap on long (≥ 17 bp) sites is the expected failure
  mode (THEORY §11b) and is reported, not averaged away.
