# E1 amendment 2 (2026-10-12): comparator runs stopped; E1 closed as an indicative result

The E1 comparators were stopped by the user's decision. The front-loaded comparison already answered the question:
HASTUR is calibrated on small bacterial sets but has no frontier advantage there.

**State at the stop**
- **Complete on all 1,974 units:** `hlat_f2r`, `streme`, `meme` and `xstreme`. `sea` is also complete; it ran on
  the workstation as a known-motif reference, so it falls outside E1's one-machine rule.
- **Partial:** `homer`, 740 units. Not scored.
- **Not run:** `streme_all`.

**Indicative readouts (HASTUR and STREME; the registered report was not produced)**
- **HASTUR F2:**
  - native: 8/32 evaluable modules at 0.044 false calls per null run;
  - at matched rates of 0.05 / 0.1 / 0.5 / 1.0 false calls per run: 8 / 13 / 17 / 19.
- **STREME:**
  - native: 30/32 at 9.6 false calls per run;
  - at the same matched rates (default run, E-sorted): 10 / 11 / 18 / 19.
- **SEA (37-TF RegulonDB bank, circular with the truth):** 29/32 at 0.059 false calls per run.
- **HASTUR known-motif quick test:** 16/32 at 0.055.

**Conclusion:** HASTUR's distinctive property on small bacterial sets is calibration at its native threshold. Its
sensitivity at matched false-call rates is level with STREME's, and SEA owns the known-motif case.
