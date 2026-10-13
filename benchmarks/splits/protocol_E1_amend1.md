# E1 amendment 1 (2026-10-12, before any E1 output was scored or read)

**Change: machine.** All tools run on the box (64 cores) instead of the workstation. The box was idle while the
workstation run needed about 15 h. The parity rule is unchanged: every tool runs on every unit on one machine.

**Discarded.** The partial workstation run (HASTUR on every unit, STREME on part) is discarded unread. No output of
it was passed to `bench/imodulon_eval.py` or inspected. Units are copied to the box without tool outputs, and every
tool is rerun there.

**Runtime.** Runtimes are measured on the box. HASTUR runs with `--threads 2`, as on the workstation.
