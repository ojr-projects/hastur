# Registered protocol LN3: one configuration for ChIP and SELEX (2026-10-10, before any LN3 run)

**Context.** The LN2-frozen configuration (v8, `--sites --gtrim 0.01`) was developed on ChIP alone. On SELEX dev it
regresses against plain HASTUR (P5 top-10 0.64 vs 0.72; located 0.22 vs 0.25), and HASTUR reports few motifs on
SELEX (median 2 vs STREME's 10). LN2 Part B (ChIP test) proceeds as registered and is unaffected.

## Part A: development (dev units only: ChIP dev low-N + N = 5,000, SELEX dev 213 libraries)

- **C1 (representation by evidence):** the site matrix becomes a third candidate in the held-out representation
  choice (frame / EM / sites), as THEORY §6.2 prescribes, instead of a blanket refit.
- **C2 (more distinct calls):** 20 seeds per family instead of 12 (de-duplication and fragment absorption unchanged).
- C1 + C2 combined.

**Acceptance** (paired, against frozen v8 on ChIP and against plain v6 on SELEX): SELEX located and P5 top-10 must
improve on v8 and not fall below plain v6 by more than 0.03; ChIP located / P5 top-10 must not drop by more than
0.03 against v8 at any N; motif-free false `sig` calls per run must not rise. At most 3 rounds; stop at the first
round whose best SELEX change is < 0.02.

## Part B: confirmation (run once, configuration frozen and hashed before the run)

**Primary:** the locked HT-SELEX **test** split (`splits/selex_split.tsv`, never run by this project), built with
`data/prep_selex.py test` (20k final-cycle reads vs 100k cycle-0 reads; JASPAR truths). Tools: frozen HASTUR, STREME,
HOMER (LN1 settings). Criteria: as LN2 Part B (paired located and P5 top-10 ≥ −0.05 with lower bound ≥ −0.10 vs
STREME; PWM similarity not lower by > 0.02; runtime ≥ 3× lower), with the chance rates reported.
**Secondary (not independent):** the ChIP test factors re-run with the LN3 configuration, reported beside LN2 Part B.
