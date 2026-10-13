# Registered protocol LN6 round 3 (registered 2026-10-11, before any run on the decision panels). The last round:
# if it does not pass, development pauses (user decision, ROADMAP §5).

**Round 2** (`results/ln6r2_report.txt`): R2 (exact tests, per-sequence presence counts, exact radius-1 ball test for
L-mers) passed both real-data null gates:
- 0.05 false calls per run on exchangeable real-genome nulls (v8 0.71, STREME 0.00);
- 0.019 on small-N sets (v8 1.76).

It failed the simulated CpG-shift null (4.45 per run vs v8 0.05), and located at N = 5,000 (−0.105, n = 19).

**Diagnosis.** The false calls are GC-rich words. Presence expectations from 1 − exp(−λ) assume Poisson
occurrences, while CpG-rich chains clump them, differently in fg and bg, so the double-centring no longer cancels.

**Fix.** Expected presence counts are computed per set from synthetic replicates of that set, drawn from its own
Markov chain with the same lengths (≥ 2,000 sequences, at most 50 replicates). The empirical/analytic ratio is shrunk
to 1 where synthetic counts are small: E = E_an (S + 5) / (M E_an + 5). In a sanity run, this removed every false
call on three CpG-shift nulls and kept CTCF at N = 250.

## Arm (prototype v14; base frozen v8 `--sites --gtrim 0.01`)
| arm | tool key | flags added |
|---|---|---|
| R4 | `hlat_z1` | `--exact-all --seqcount --lmer-ball --presence-emp` |

## Decision data, all fresh (no run of any decision panel has been looked at)
- **Small-N null panel** (`bench/upnull_make.py <root> 2`): same three genomes and design as rounds 1 and 2, new
  sets. 210 runs.
- **Simulated nulls** (`bench/null_make.py`, seed offset 920000): iid, Markov, CpG shift, 20 each.
- **Real-genome nulls** (`bench/realnull_make.py`, `SEED_OFFSET=7000`): human_exch N = 500 and 2,000, and
  {bceno, cjejuni, ecoli}_coding, 20 each (100 runs). The intergenic and GC-selected kinds are built but not gated,
  for the reasons in round 2's amendment.
- **ChIP dev panels as in rounds 1 and 2** (N = 25 to 5,000, matrix readout).

## Keep rule (R4 against frozen v8 `hlat_a3`). All must hold:
1. **Small-N nulls:** ≤ 0.1 false calls per run pooled, and ≤ 0.2 in every size bin.
2. **Real-genome nulls:** ≤ 0.1 false calls per run pooled.
3. **Simulated nulls:** false calls per run not above v8's in any kind (iid, Markov, CpG shift), within one call in
   20 runs (+0.05).
4. **ChIP, N ≥ 250:** P5 top-10 and located top-10 within −0.03 of v8 at every N. N = 25–100 reported, not gated.

If R4 passes, it is combined with LN5's L2 and frozen (FR1 amendment 2), then the FR1 frontier is rerun for it.
If it fails, development pauses.
