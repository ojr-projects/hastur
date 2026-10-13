# Matched genome backgrounds (Mode A)

This tool builds a background whose tiles match the foreground's composition. Use it when the foreground set is
composition-selected (for example GC-rich windows). In that case a random genomic background makes any motif finder,
HASTUR included, call composition as motif.

```bash
python matched_index.py index    hg38.fa.gz rmsk.txt.gz gencode.gtf.gz idx/   # once; ~17 GB index
python matched_index.py annotate idx/ fg.bed fg_cov.tsv                         # optional: inspect covariates
python match_sets.py idx/ fg.bed out/ --k 10 [--tss] [--exclude peaks.bed.gz] [--masked]
hastur out/fg.fa out/bg.fa -o out/hastur
```

**How tiles are matched:**
- **Buckets:** GC (0.05 bins), CpG o/e (4 bins), soft-masked repeat share (3 bins, with SINE vs other at the centre),
  and optionally TSS distance.
- **Sparse buckets** are coarsened in a fixed order.
- **Exclusions:** tiles that overlap the foreground or `--exclude` are never drawn.

`strata.txt` is also written. It was the input to an experimental stratified permutation (`--strata`), which MA1 showed
added nothing over the matched background alone, so it was not carried into the engine.

**Evidence and cost:** see the main README and `benchmarks/splits/protocol_MA1.md`.

**Requirements:** Python 3 with numpy.
