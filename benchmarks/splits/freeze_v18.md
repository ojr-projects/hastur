# Freeze v18: HASTUR configuration for the FR1 confirmation (2026-10-11, before any confirmation output is read)

**Tool key:** `hlat_f2` (`bench/run_tools.py`).

**Flags:** `--sites --gtrim 0.01 --plainz --exact-all --lmer-ball --perm 20 --rank-v8 --no-gc-test`

**What each part does:**
- **Detection and ranking:** v8's lattice statistics (Wald gapped z, soft L-mer z) rank the seeds, so the motif list is
  v8's.
- **Significance:**
  - **Tested statistics:** the exact conditional binomial test on gapped cells and the exact radius-1 Hamming-ball
    test on L-mer cells, with no genomic control.
  - **Thresholds:** they are calibrated by permutation max-T. There are 20 label permutations of the pooled
    sequences, and each family's maximum is recorded. A Gumbel distribution is fitted by moments, giving a
    family-wise α of 0.025 per family.
  - **Score per call:** a family-wise adjusted p, with significance at p_adj ≤ 0.05.
- **Representation:** v8's site refit with LN5's L2 (plain fg/bg site threshold).

**Prototype:** `proto/hlat.py` VERSION v18, sha256 181475d992ce4e4a1659e6ba475c306582d09faa678b191461607ae9be286193 (the box copy is identical). Git commit 2737ab1.
`bench/run_tools.py` sha256 f8e114c9db76cc6ebc8e986a515855016a47822276f3ec66dc707d8b05d01acd.

**Dev evidence** (LN7 cycle e; dev panels, seen; `results/ln7e_report.txt`, `results/frontier_ln7e.txt`):
- **Gates:**
  - small-N nulls 0.038 false calls per run;
  - exchangeable real-genome nulls 0.03;
  - simulated nulls ≤ 0.05 (CpG shift 0.05, at the gate's bound);
  - ChIP top-10 list identical to v8.
- **Sensitivity at the native threshold** equals v8's at N = 250–5,000.
- **FR1 dev decision cell:** F2 − STREME +0.190 [+0.142, +0.237].

**Combination check (FR1 amendment 2):** L2 leaves the calls unchanged, and its best-similarity gain is the same on
top of P5 as alone (+0.007 / +0.010).
