# D1 amendment 1 (2026-10-11, before any D1 tool run on the box; found in a 3-unit pipeline test)

STREME holds out no test set when the input is small (e.g. PurR, 18 promoters). It then reports `S=` (a
training-set p-value) instead of `E=`. Its `--thresh` applies to whichever it reports. The FR1 parser read only
`E=`, so it would have scored such runs as having no calls.

**Change:** `bench/frontier.py` and `bench/null_eval.py` read `E=` or `S=`.
- STREME's native operating point is E ≤ 0.05 or S ≤ 0.05, whichever STREME reports, which is its own threshold rule.
- Null units have the same sizes as their regulons, so the false-call side is treated the same way.
- No earlier panel (low-N, ChIP 5k, SELEX, nulls) has any `S=` output, so earlier results are unaffected.

Also reported: the number of D1 units in S mode, per tool run.
