# FR1 amendment 2 (2026-10-11, before any confirmation output is read)

D1 exposed a calibration failure at small N, which is now registered as LN6. The configuration confirmed in FR1 is
therefore:
- **Frozen v8 plus the kept arms:** whichever LN5 arm (matrix accuracy) and LN6 arm (small-N calibration) are kept
  under their registered rules, combined.
- **Frozen together:** recorded in `splits/freeze_v12.md` (flags and git hash), after one dev sanity run of the
  combination. That run checks that its dev numbers match each arm's own within 0.01, so the effects do not
  interact. If they do interact, the combination is not used, and the more important arm alone (LN6, for
  calibration) is confirmed.
- **v8 is also reported** on the confirmation data. If no arm is kept, v8 alone is confirmed (amendment 1, item 5).

The confirmation runs for the combined configuration are separate tool runs on the same confirmation inputs and fresh
nulls. Comparator runs are not repeated.
