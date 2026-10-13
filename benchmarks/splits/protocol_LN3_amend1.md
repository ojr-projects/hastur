# LN3 amendment 1 (2026-10-10, before any LN3 Part A result was produced or seen)

1. **Tie-break.** If more than one C arm meets the acceptance rule, the arm frozen for Part B is the one with the highest
   SELEX dev P5 top-10; ties broken by SELEX located top-10, then by lower median runtime.
2. **Venue.** Part A and Part B run on the workstation (the EC2 box was terminated). Tool parallelism is capped by a
   memory gate (MIN_FREE_GB = 5); every tool run that is killed or fails writes a `failed` marker and is rerun, never
   scored (an OOM kill during LN2 Part B would otherwise have scored 24 STREME runs as misses).
3. **Timing.** Runtimes measured under contention are reported but the speed criterion is re-measured on an idle
   machine with one job at a time on a fixed subset (10 units per N), before any claim.
