# DEV_INSIGHTS: general learnings from Comsequent and HASTUR, for future projects

> **This copy builds on Comsequent's.** Sections 0–7 are Comsequent's learnings (2026-09-28 → 2026-10-10), with HASTUR
> additions marked **[HASTUR]** inside them. Section 8 is the HASTUR development cycle (2026-10-10 → ), whose headline
> lesson is a calibration claim that held only on simulated nulls.

> **General learnings, not project status.** This file distils what the Comsequent development cycle (2026-09-28 → 2026-10-10)
> taught about developing and benchmarking bioinformatics tools quickly. It is written as the **starting point for the next
> project** in the same setup: one human directing, Claude Code agents doing the work, a WSL2 workstation (Ryzen 7900X,
> RTX 4070, ~15 GB RAM visible to WSL, disk image on C:) plus an EC2 spot box.
> For Comsequent's own status, see `docs/ROADMAP.md` and `docs/REAL_DATA_PLAN.md`. For the narrative, see `docs/COLLABORATION_STORY.md`.

The goal is **speed to a defensible endpoint**: reaching the true conclusion with the fewest wasted cycles.

---

## 0. The lesson in one paragraph

Comsequent's final, correct conclusion was a **domain-of-applicability rule**: the method wins when the specificity is exact
and discrete *and* the contrast is near-binary and matched by construction (methylation), and it loses on PWM-like
signals (ChIP, SELEX, RBNS, promoters). That conclusion arrived on **day 13** of the collaboration, after about 8 active
days, 2 full simulation grids (1,050 and 2,720 datasets), 20 ordering variants and 6 evaluator versions. Almost all of it could
have been predicted or measured by **day 3–4**. Three things would have achieved that:

1. **A mechanistic "why would this win or lose" analysis on day 0.** The conjunctive-exact vs additive-soft argument in
   `UNIFIED_MOTIF_MODEL.md` §1 predicts the PWM loss. It was written on day 10.
2. **A tiny real-data pilot of every candidate application on days 1–2**, scored by one registered rule against the
   standard tools.
3. **A cheap pre-check that predicts the outcome** before any discovery run. The signal-to-drift composition ratio was
   built on day 13 and separates every win from every loss so far.

Simulation should then have been used to **engineer inside the winning domain**, not to discover the domain.

---

## 1. Where the time went

| Period | What happened | What it bought | What was avoidable |
|---|---|---|---|
| 2025, ~10 weeks solo | Promoter scoring, validated against shuffled negatives | The core method | The negatives *were the method's own null*: a dinucleotide Markov baseline matched the tool on day 1 of review |
| 09-28 (1 day) | Stock-take: confound + 6 real bugs | Correct code; honest baseline | Nothing. This was the highest-value day. Do it first, always |
| 10-02 → 10-08 (~4 active days) | Engine rewrite, Phase 0 simulation grids, 6 evaluator versions, gapped panel, U-series | Calibrated engine; simulated wins on gapped/bipartite | Evaluator bugs re-scored grids 4+ times; the gapped comparator panel came only after an outside critique; full-grid dev cycles took hours until the 95-dataset dev panel arrived on day 11 |
| 10-02 → 10-06 | Raw methylation data (146 GB) sat on disk; GPU idle across a holiday | — | The lead application's basecalling (~15 h unattended) could have run then. "Finish Phase 0 first" was applied to all compute, not just to dev attention |
| 10-08 → 10-09 | u8–u20 ordering variants | u8/u9 real gains, then **trade-offs within ±0.02** | ~6 variants after the plateau. A stopping rule would have cut them |
| 10-09 | One real-data defect fix (round cap): 0.61 → 0.83 | Biggest gain of the project | It was found on real data. It would have been found days earlier with real data in the loop from the start |
| 10-09 → 10-10 | RBNS, ChIP, HT-SELEX, CAP-SELEX, PAM, promoters | The domain map | Most outcomes were predictable from motif class + drift pre-check; several runs used untested settings (HOMER L ≤ 20 on 50k reads: ~4 days → restarted; STREME on 200k reads stalled the box) |
| Throughout | Disk-full crashes ×2, OOM storms, self-killed shells, oversized commit, lost box, idle local machine | Guard scripts | Nearly all avoidable with the day-0 infrastructure in §5 |

**Late outside critiques that cost rework**, each preventable by a day-0 hostile-review pass (§3.10):
- "STREME/HOMER aren't built for gapped motifs; use dyad tools."
- "They won't switch at parity."
- "Your real data became your dev set; declare held-out data."
- "Check that comparators reproduce their published performance."

---

## 2. The faster path, in hindsight (template for the next tool)

| Day | Do | Exit criterion |
|---|---|---|
| 0 | Code review for bugs; **trivial baselines** (composition, Markov LLR, the obvious field tool) on **realistic negatives**; mechanistic win/lose theory; candidate applications listed with motif/signal class; **dev/test splits fixed**; prediction registry started; hostile-review pass on the plan; infrastructure from §5 installed | Written predictions per application; the tool beats the trivial baseline on *something*, or you stop |
| 1–2 | **Real-data smoke panel:** 2–5 items per application, matched N, every tool at defaults, one registered scorer, a pre-check per item; GPU and box queues fed with the long data-prep jobs for the most promising application | A domain map: where it wins, ties, loses, and whether the pre-check predicted it |
| 2–4 | Develop **only in the winning domain**, on dev items. Use a ≤ 15-minute simulation dev panel for mechanisms; keep real dev data in the loop | Gains stop exceeding replicate noise (§3.8) |
| 4–5 | Freeze, tag, run held-out sets once; bounds (losing domains) run once as reported limits | Held-out within the pre-registered tolerance |

Full simulation grids are a **reporting artefact**: run them once, on the frozen binary, for the paper. They are not a development loop.

---

## 3. Science-process rules (speed-oriented)

### 3.1 Baselines and negatives before anything else
- First experiment: the simplest baselines (composition, k-th order Markov LLR, the field-standard tool) on the intended
  benchmark. If the negative set shares the method's null model (e.g. shuffles for a shuffle-null method), the benchmark is invalid.
- Negatives must be **realistic** (genomic background, matched composition), not synthetic.
- **[HASTUR]** This applies to *nulls* as much as to negatives. A synthetic null drawn from the method's own background
  model satisfies the test's assumptions by construction, so it cannot expose the test's failures (§8).

### 3.2 Real data in the loop from day 1
- In Comsequent, simulations plateaued into reshuffling replicates, and the decisive defect was found on real data.
  Keep 2–5 real dev items per application running alongside simulation work from the start.
- Start the long data-prep pipelines (downloads, basecalling, alignment) for the lead application **immediately and
  unattended**. They gate everything downstream and need no dev attention.

### 3.3 Write the mechanism first
- Before benchmarking, write one page: what structure in the data does the method exploit that the standard tools cannot,
  and what structure does it *discard*? Turn this into per-application predictions (win / parity / loss).
- Comsequent: "exact conjunctive matches vs soft additive PWMs" predicts both the gapped-motif win and the PWM loss.

### 3.4 Build the outcome pre-check early
- Find a cheap statistic, computed *before* any discovery run, that predicts whether the method will win on a dataset.
  For contrast methods: **signal vs drift**. Divergence with the reference motif masked is the drift; the total minus that is
  the signal. Signal/drift > ~1 tracked every win; < 0.4 every loss (`benchmarks/tools/drift_check.py`).
- The pre-check is also what you tell users: a domain rule they can apply before running the tool.

### 3.5 Validate the evaluator before the first grid
The evaluator went through v1 → v6. Every version invalidated earlier numbers and forced re-scores. Bugs found:
comparator PWM calls scored as false, lenient credit for background sub-patterns, a capped comparator, IUPAC truths
unscoreable, **every tool scoring 0 on PWM truths**, spurious zeros from unrenamed side outputs, wrong truth motifs (2 of 4 ChIP),
and a match rule revised after seeing results. Before any grid:
- **Oracle test:** feed the truth itself (and its IUPAC members, its PWM form, a near-miss) as calls. Scores must be
  1 / 1 / 1 / 0.
- **Comparator positive control:** reproduce each comparator's published or obvious performance on an easy case. If a
  well-established tool scores 0, the evaluator is wrong until proven otherwise. Comsequent's owner asked this question,
  and it caught the bug that would have inflated every early win.
- **One scoring rule, defined from the reference alone**, used for every phase (Comsequent ended up with
  `registered_scoring.py`: cover every reference position with IC ≥ 1 bit). Have it from day 0.
- Truth motifs from a curated database, checked for orientation and form before scoring.

### 3.6 Validate the simulator
- Check that planted-site counts and flank composition match the intent. A context-aware planting mode planted almost nothing;
  the default mode created a Markov-context artefact. Both were caught late.
- A quick "planted vs detected by a trivial counter" table per class catches most of this.

### 3.7 Calibrate nulls before measuring power
- λ ≈ 1 on motif-free splits for every statistic, before reading any recall figure. Comsequent found an SE overstated
  ~9× at l = 7, a Wald z collapsing under separation, anti-conservative background statistics, and a pooled epistasis
  "inroad" that turned out to be anti-conservative. Each would have been seen on the first null run.

- **[HASTUR] Calibrate on real motif-free sequence, at every size the tool will be used at.**
  - **Use exchangeable real splits:** two disjoint random sets of genomic windows, or a random subset of a fixed
    population (promoters, peaks) against the rest. Include the smallest deployment N (10–60 sequences for bacterial
    regulons), not just "typical" N.
  - **Tile the windows; never sample them with overlap.** Overlapping windows duplicate sequence and inflate every
    statistic (one HASTUR panel had 376 of 500 foreground windows sharing a 20-mer).
  - **Give the null panel its own oracle test:** a deliberately uncorrected statistic (e.g. a plain Wald z) must fail
    on it. A null that nothing can fail proves nothing.
- **[HASTUR] Occurrence counts are clumpy on real sequence** (low complexity, repeats, codon structure), so Poisson
  tests on occurrences are anti-conservative in the tail. Median-based genomic control (λ) rescales the bulk but cannot
  fix a tail. Count sequences (presence), as STREME and HOMER do, or test on held-out data.
- **[HASTUR] Wald z on log ratios with observed-count variance is badly anti-conservative at small expected counts.**
  One observed vs 0.09 expected gives z ≈ 4.2, where the exact tail gives 1.4. Use exact conditional tests wherever
  expected counts are small. Note that λ, estimated only on cells with large expectations, never sees these cells.
- **[HASTUR] Max-over-many designs need extreme-tail accuracy.** Taking the best of ~10⁶ cells with Bonferroni needs
  the per-cell null tail to be right at p ≈ 1e-8, which real data rarely gives. A held-out test of a few preselected
  candidates needs p ≈ 0.05. Weigh this before choosing all-data detection: it is a power-vs-robustness trade, not a
  free gain.

### 3.8 Size the dev loop, and set a stopping rule
- **Dev cycle ≤ 15 minutes.** Comsequent's 95-dataset dev panel turned hours into 15 minutes, but only arrived on day 11.
  Build it with the first grid.
- **Stopping rule, written in advance:** compute the replicate noise of the headline metric (e.g. ±0.02 on 1,050
  datasets). When two consecutive variants move the net metric by less than that, or only trade wins between classes,
  stop and freeze. u10–u20 were mostly trade-offs.

### 3.9 Splits, registry and claims at day 0
- Fix dev/test splits (strains, TFs, RBPs, libraries) with a seed and a committed hash **before looking at any data**.
  Comsequent declared them on day 12, by which time every strain had been inspected, so a new 112 GB held-out set was needed.
- Keep a **prediction registry**: motif class, contrast type, pre-check, predicted outcome and readout, registered before
  each run. Mark retrospective entries honestly.
- **No outcome statements before the registered scorer has run on the whole set.** A "looks like a win" from one PAM
  sample was wrong.
- Any post hoc choice (command, match rule) turns that phase into development data. Label it.

### 3.10 Hostile-review pass at day 0 and at every milestone
Ask an agent, explicitly, "what would a reviewer use to dismiss this?" before the plan runs:
- the right comparators for each motif class, including specialist tools (dyad/gapped tools here), at their defaults;
- matched N and matched inputs across tools;
- held-out data;
- adoption: why would anyone switch at parity?
Every one of these arrived late in Comsequent as a relayed outside comment, and each forced rework.

### 3.11 Time every tool before launching a batch
- Run one unit at the planned N and settings, measure wall time and peak RSS (`/usr/bin/time -v`), and multiply out.
  Abort the plan if the total is unreasonable. HOMER at motif length 20 on 50k reads was ~4 days for one phase;
  STREME on 200k reads stalled the box for hours. Both were restarted, and 31 finished libraries were discarded.
- Hard per-run timeouts that kill the **whole process group** (tools outlived their timeouts as orphans).
- Matched N across all tools from the first run.

---

## 4. The plan document itself

- **One master status file** (here `docs/ROADMAP.md`) that supersedes every other doc's "status" section. Fixed-input
  briefs stay read-only in `docs/briefs/`.
- Dated update logs, decision pages (`DECISION.md`, `SIM_DEV_FREEZE.md`) and a **"State at <date> — read this first
  after a restart or compaction"** block. These survived three context compactions and a VS Code crash.
- Keep the status snapshot current or delete it. The roadmap's "Status snapshot (2026-10-06)" was stale within a day.
- Verify large doc edits: a silent non-matching replace once dropped whole findings sections (`bd7f5db`).
- Too many parallel planning docs (benchmark plan, addendum, differential brief, theory addendum, schedule,
  differential plan, unified model, real-data plan v1–v3) make "what is the current plan?" expensive. Prefer one plan
  with a dated changelog.

---

## 5. Compute and infrastructure (same setup next time)

**Install all of this on day 0, before the first long job.** Comsequent's worst lost time came from the infrastructure
crises, and every guard that fixed them would have been cheap to set up front. Templates to copy are listed in §7.

### 5.1 Local WSL disk: budget against the image's high-water mark
- WSL's `/` is a virtual disk (`ext4.vhdx`) on C:. `df -h /` inside WSL reports free space inside the image (~400 GB):
  **that figure is fiction.** C: free space is not the right check either: deleting inside WSL never returns space to C:,
  and C: fell 88 → 79 GB on a day when space had been freed inside WSL.
- **The rule (adopted 2026-10-10):**
  1. Know the image's physical size, its **high-water mark** (576 GB here). Measure it from Windows (WinDirStat or
     Explorer on `ext4.vhdx`), or `stat` the vhdx via `/mnt/c/...` from inside WSL.
  2. Usage inside WSL comes from `df` on `/` (140 GB here).
  3. **Usable budget = high-water mark − used − margin** (~436 GB here). Writes inside that budget reuse space the
     image already holds. Writes beyond it grow the image and eat C:.
  4. Cycle data through (download → process → delete). Never accumulate.
- **Treat any growth of the vhdx file as an alarm.** Have the disk guard check the vhdx file size itself, not only C: free
  (ext4 does not promise to reuse previously touched blocks before untouched ones).
- Keep a guard that stops local writers below a C: floor (80 GB here) as a backstop.
- Bulk data (multi-GB downloads, large work dirs) goes to the box by default. One 36 GB local download filled C: and
  crashed WSL, killing every local job.
- **[HASTUR] Fixed on 2026-10-12.**
  - **Before:** the image was 619 GB apparent and 537 GB allocated, while WSL used 170 GB.
  - **Why diskpart failed:** the image was already sparse, and diskpart refuses sparse disks.
  - **What fixed it:** a one-off `sudo fstrim -v /` inside WSL. The root fs is mounted with `discard`, so a sparse
    vhdx hands trimmed blocks back to Windows.
  - **After:** the image is 179 GB allocated, and C: went from 113 to 487 GB free.
  - **Keep sparse on.** Re-run `fstrim` after large deletions. Check the allocated size, not the apparent one: `du`
    on the vhdx via `/mnt/c`, not `ls`.
- Worth checking once: WSL's sparse-VHD option (`wsl --manage <distro> --set-sparse true` on recent WSL versions) and
  periodic compaction let freed space return to C:. Verify against current WSL docs before relying on it.

### 5.2 Local memory
- WSL saw ~15 GB. Every heavy launch reserves from a **shared memory ledger** (`memory_slot()` in
  `benchmarks/shared_eval.py`; 60 % of RAM), and a **watchdog** (`memwatch.sh`) kills the newest heavy child when
  MemAvailable < 1.2 GB. Before these existed, OOM storms killed WSL repeatedly.
- **Reservations from measured peak RSS, never from guesses.** An XXmotif guess of 6 GB (measured 3.6) throttled a grid
  ~10×; k = 7 runs measured 8.1 GB against a 7 GB slot and starved.
- Check whether `.wslconfig` `memory=` can be raised (WSL defaults to a fraction of host RAM). More visible RAM is the
  cheapest local throughput gain.

### 5.3 Keeping the local machine busy (the persistent problem)
The box was kept busy by a cron-driven queue plus a half-hourly utilisation watcher. **The local machine never got the
equivalent**, and it idled repeatedly ("nothing running locally"). Causes, from the logs and docs:

1. **Scheduling depended on the agent's attention.** Local jobs started only when the agent was awake and remembered.
   The local wake-up watcher (`watch_once.sh`) fired every 110 minutes or on a queue event. The box had cron every 10 minutes.
2. **"Nothing restarts automatically."** This was the right reaction to the first OOM storm (everything restarted at once),
   but it turned every failure, refusal or watchdog kill into idle time until someone noticed. Serial queues
   **skip** failed items and move on, or stop.
3. **Over-conservative guards with no fallback work.** Memory reservations were inflated, launchers refused below
   4 GB free, and the disk guard stopped local SELEX and made local "write-free". With nothing guard-compatible queued,
   the CPU idled.
4. **Shallow backlog.** Queues held the next 1–3 items, not a day's worth. Once data-heavy work moved to the box, local
   had little left to chew on.
5. **GIL-bound Python.** 4 threads of evaluation ran at ~1.3 CPUs until it was moved to process pools.
6. **The GPU was not treated as a separate always-fed resource** until basecalling got its own queue (it sat idle
   10-02 → 10-06).

**Recipe for next time:**
- Run a **local clone of `box_queue.sh`**: a jobs file (`name | workdir | after | command`), done markers, idempotent
  launcher, started every 5–10 minutes by cron inside WSL (enable systemd in `/etc/wsl.conf` so cron runs). The setsid loop
  must keep the WSL VM alive.
- Give the launcher a **utilisation target**: if load average < ~0.8 × threads and the memory/disk guards pass, start the
  next eligible job. Have separate classes for CPU-heavy, memory-heavy, GPU and network jobs.
- Keep a **deep, pre-approved backlog** (≥ 12 h) ordered by value, with low-priority filler at `nice 19`: extra
  replicates, sweeps, re-scores on the latest evaluator, timing runs, held-out prep. Filler is what absorbs idle cycles.
- **Retry through the ledger:** a watchdog-killed or refused job goes back on the queue and restarts one at a time when
  resources allow, with a retry cap. That is safe because units are atomic and resumable.
- Plan local-friendly work explicitly: CPU-bound, small-write jobs (re-scores, simulations that stream and delete, dev
  panels) stay local; data-heavy jobs go to the box; the GPU queue always holds the next basecall.
- **Before the human leaves (holiday, overnight), fill every queue** with long, safe jobs, and confirm load on all three
  resources (local CPU, local GPU, box).
- One status command shows both machines: local `uptime` + `nvidia-smi` + queue depth, and the same on the box.
  The agent's periodic check uses it.

**The agent's code-writing time is not a scheduling constraint (HASTUR, 2026-10-12).** The agent justified keeping
the paid box on a workstation experiment with "the next PROMETHEUS compute needs code written before it needs the
box, so the box would otherwise sit mostly idle". That had it backwards. For these tools, writing a component (a
builder, a Rust flag, an evaluator) takes the agent minutes. Evaluating it (registered rounds, null panels, every
comparator) takes hours to days.
- **Rule:** treat code time as ≈ 0 when deciding what runs where. Never leave a machine idle "while the code is
  written", and never justify an allocation by it.
- **Do instead:** write the next job's code while the current jobs run, so that every machine always has its next
  job ready. Move work between machines freely when that keeps both busy.
- **Idle capacity is the cost**, especially on the paid box. The utilisation watcher
  (`infra/local/util_watch.sh`, every 15 min on both machines) flags it on every check.

### 5.4 EC2 spot box: setup, monitoring, interruption, recovery
- **Idempotent setup script** (`box_setup.sh`: Go, micromamba, one env per tool, symlinks, engine build). A fresh box
  is productive in minutes, which matters because boxes do get lost (one was terminated by accident; spot can reclaim).
- **Boot-time, cron-driven job queue** (`box_queue.sh`): `@reboot` and every 10 min; jobs file with dependencies;
  a job runs unless it is running (`BOXJOB=<name>` in the environment) or done (marker on exit 0).
- **Spot interruption handler** (`spot_notice.sh`, IMDSv2, 5 s poll): on notice it writes a flag (no new jobs start),
  flushes, and snapshots all small result tables into one tarball.
- **Pull-back** (`box_pull.sh`): a local job pulls the result tables every 30 min and checks the notice flag every 5 min.
  A lost box then costs at most 30 min of results. This pull-back was added hours before the accidental termination,
  and it is why that cost almost nothing.
- **Utilisation watcher every 30 min** (fixed the box's idle stalls). Run the same for local (§5.3).
- Network: limit parallel downloads per host (EBI refused 20 workers; 3 with retries worked), and verify completeness
  afterwards. A job "finished" with 281 failed downloads.
- Keep ≥ 50 GiB free on the box volume; compress intermediates; delete scratch.
- Host/IP change on every stop/start: read them from environment variables, never from the repo.

### 5.5 Job design rules (all long jobs)
- One **atomic part file per unit**: write to temp, then rename. Zero-byte crash leftovers must count as "not done";
  they once passed as finished.
- Expensive tool outputs are cached per unit, so re-scoring never re-runs tools.
- **Private work directories per job/tag.** Concurrent jobs deleted each other's files, and an empty-assembly crash
  lost 19 datasets.
- Launch detached (`setsid nohup … &`). A VS Code crash killed every non-detached job.
- Timeouts kill the process group; check for orphaned tools after stops.
- Compress or drop regenerable intermediates as soon as a unit is scored.
- Evaluator and engine versions recorded in every result row, so stale rows can be found and re-scored.

### 5.6 Agent shell hazards
- **Never `pkill -f` / `kill $(pgrep -f PATTERN)` when PATTERN appears in the same command line.** The Bash tool runs
  `bash -c "<whole command>"`, so it matches itself; this killed the agent's shell or SSH session several times. List PIDs
  in one call (`ps -eo pid,cmd | grep "[p]attern"`), then kill exact PIDs in the next.
- `git stash push` on a clean tree saves nothing, and the following `stash pop` applies an older stash. Check
  `git stash list` first.

### 5.7 Git hygiene
- Never `git commit -a`. Add paths explicitly. Regenerable result tables went over GitHub's 100 MB limit and blocked a push.
- Gitignore regenerable large outputs from day 0. Before each push, check blob sizes:
  `git rev-list --objects origin/<branch>..HEAD | git cat-file --batch-check='%(objecttype) %(objectsize) %(rest)' | awk '$1=="blob" && $2>50000000'`
- Decide on day 0 who pushes (here the agent, from day 12).

### 5.8 Secrets
- Keep keys and box briefs **outside the repo tree** (e.g. `~/.ssh/`). In Comsequent the `.pem` and box brief sat in `docs/`,
  excluded only via the local `.git/info/exclude`. A fresh clone, another worktree, or a broad `git add` would not be protected.

---

## 6. Working with agents

- **Rigour as a standing instruction.** "We need to find out rigorously, even if the method isn't great" licensed plain
  negative reporting. That made the final conclusions credible. Say it on day 0.
- **The human's highest-leverage moves were short questions:** "are the comparators really that bad?",
  "is the bridge uncrossable?", "competitors winning is fair game", "are we finding trade-offs rather than gains?".
  Schedule these as explicit checkpoints, not luck: a comparator sanity check after the first grid, and a plateau check
  after every 3 variants.
- **Opus to design, a cheaper model to implement** (schedule with acceptance tests per task, "implement as written,
  escalate deviations") worked well.
- **Standing rules go in agent memory**, one rule per file with the reason. Examples: disk method, no `pkill -f`,
  no `commit -a`, keep machines busy, downloads to the box. They carry across sessions and compactions.
- **Autonomy needs fences first:** grant "continue without my say-so" only after the resource guards (§5) exist. The
  first autonomous day here filled the disk.
- **Parallel agents:** give each one a written brief (e.g. `FINDINGS_BRIEF.md`) rather than shared chat context, and
  keep a single owner for the master status file.
- **Outside critiques are the best input.** Relay them early and turn each into a concrete task or protocol rule.

---

## 7. Starter kit: files worth copying from Comsequent

| Need | File (in this repo) |
|---|---|
| Resumable job launcher with disk/memory refusal | `benchmarks/tools/resume_all.sh` |
| Serial / named queues | `benchmarks/tools/queue.sh` |
| Memory ledger | `memory_slot()` in `benchmarks/shared_eval.py` |
| OOM watchdog | `benchmarks/tools/memwatch.sh` |
| Agent wake-up on job events | `benchmarks/tools/watch_once.sh` (shorten the interval) |
| Box setup | `benchmarks/tools/box_setup.sh` |
| Box cron queue (also the template for a local one) | `benchmarks/tools/box_queue.sh` |
| Spot interruption handler | `benchmarks/tools/spot_notice.sh` |
| Result pull-back | `benchmarks/tools/box_pull.sh`, `sync_box.sh` |
| Registered, tool-agnostic scorer | `benchmarks/tools/registered_scoring.py` |
| Outcome pre-checks | `benchmarks/tools/composition_check.py`, `drift_check.py` |
| Shared evaluation + simulation templates | `benchmarks/shared_eval.py`, `shared_sim.py`, `tools/inject_shared.py` |
| Decision / freeze page formats | `benchmarks/phase0_sim/DECISION.md`, `docs/SIM_DEV_FREEZE.md` |
| Prediction registry + domain rule format | `docs/REAL_DATA_PLAN.md` §0 |

### Day-0 checklist (copy into the new project's ROADMAP)
- [ ] Code review for bugs; trivial baselines on realistic negatives
- [ ] One-page mechanism: what the method exploits and what it discards → per-application predictions
- [ ] Dev/test splits fixed and hashed; prediction registry started
- [ ] Single registered scorer; oracle test passes; each comparator reproduces an easy known result
- [ ] Hostile-review pass: comparator panel per class, matched N, held-out, adoption argument
- [ ] Tool timing + peak RSS at planned N for every tool; batch totals estimated
- [ ] Null calibration (λ ≈ 1) before any power figure
- [ ] **[HASTUR]** Real-sequence exchangeable null panel (tiled windows; human + AT-/GC-rich bacterial; every deployment
      size incl. the smallest), with its own oracle test, before any calibration claim
- [ ] **[HASTUR]** Oracle test for **every** truth source before it is registered (not just the first)
- [ ] Dev panel ≤ 15 min; stopping rule with replicate-noise threshold written down
- [ ] WSL disk high-water mark measured; vhdx-growth guard and C: floor guard running
- [ ] Memory ledger + watchdog; reservations from measured RSS
- [ ] Local cron queue with utilisation target, deep backlog and retry-through-ledger; GPU queue fed
- [ ] Box: setup script, cron queue, spot handler, 30-min pull-back, 30-min utilisation watcher
- [ ] `.gitignore` for regenerable outputs; push policy; secrets outside the repo tree
- [ ] Master status file with a "read this first after compaction" block
- [ ] Real-data smoke panel of every candidate application scheduled for days 1–2

---

## 8. HASTUR (2026-10-10 → ): learnings added to Comsequent's

### 8.1 The headline: a calibration claim that held only on simulated nulls
- **The claim.** From day 1, HASTUR was "calibrated": λ ≈ 1 and **0.0 false calls per run** on simulated motif-free
  pairs (iid, Markov, CpG shift), against STREME's 3.8. By day 2 that claim carried the positioning ("calibrated, fast
  discovery you can run 10,000 times").
- **What tested it.** All the calibration panels drew from order-≤1 Markov chains, which is exactly the background
  model HASTUR double-centres against. Their counts were Poisson-like by construction. They tested the arithmetic,
  not the assumptions.
- **What broke it, in order:**
  1. **D1 (E. coli regulon screen):** 1.5 false calls per run on random promoter sets of 10–206 regions.
  2. **The small-N diagnosis:** a Wald z at small expected counts.
  3. **The FR1 preliminary on real-genome nulls:** **0.45–1.05 false calls per run at N = 500–2,000 on exchangeable
     human and bacterial sequence, where STREME made 0.00.** The matched-threshold frontier against STREME came out
     at −0.31.
- **What would have caught it on day 1.** The real-genome null panel takes about an hour to build. It was added on day
  2 only because an outside review asked for "realistic nulls", and it was read on day 2 afternoon.
- **The rule:** *no calibration claim from nulls generated by the method's own model.* Real exchangeable splits at
  every deployment size, before the word "calibrated" appears in any plan or summary.
- **Outcome, part 2: resumed, and fixed by changing the design (LN7, same day).**
  - **The change:** permutation max-T calibration of the maximum statistic replaced per-cell analytic calibration.
    Ranking and testing were separated, and genomic control was dropped on the tested statistics.
  - **Confirmed result:** in the registered FR1 confirmation, the frozen configuration beat STREME by +0.24 in
    sensitivity at matched false-call rates, on fresh nulls and unseen positives.
  - **The decomposition matters:** post hoc, the advantage comes from GC-shifted nulls. On exchangeable real nulls
    the two tools are at parity.
  - **Rule:** decompose any headline frontier by null family before stating what it means.
  - **Process:** five quick dev cycles, about 30 minutes each, HASTUR-only on the box, each preceded by a written
    hypothesis.
- **Outcome, part 1: development paused.** Three registered repair rounds (LN6) fixed one regime at a time.
  - **Exact tests** fixed small N.
  - **Per-sequence presence counts** fixed repeats.
  - **An exact ball test** fixed the bacterial L-mer table.
  - **Synthetic-replicate expectations** fixed CpG composition shift, but added estimation noise that the test
    ignored.
  - **No configuration passed all regimes at once.** The best (round 2) reached 0.05 false calls per run on real
    nulls but lost composition robustness. Round 3 restored composition robustness but regressed to 0.35 on fresh
    real nulls.
  - **The deeper lesson (§3.7):** an all-data, max-over-a-million-cells design needs every regime's tail right at
    p ≈ 1e-8. A held-out test of a few candidates needs about 0.05. Choose the significance design from the
    deployment data's messiness, not from the cleanest null.
- **Once the significance design is right, the evaluation shows it.** Every LN7 step was decided by which
  hypothesis survived a cheap dev cycle.
  - The SELEX ceiling was first misdiagnosed as permutation contamination; the test of that refuted it in one
    cycle.
  - The real cause (genomic control on the tested statistic) was then found by reading the outputs, not by more
    compute.
- **Patching calibration one regime at a time is a warning sign.** When each fix moves the failure to a different
  null family, the problem is the significance design, not the statistic. Stop and change the design.

### 8.2 Other HASTUR learnings
- **Dev-to-test optimism, measured.**
  - **What happened:** on dev, HASTUR was at parity with STREME. On unseen factors (LN2 Part B, post hoc) it was
    behind by 0.09–0.19 on P5.
  - **Why:** dev truths were short canonical cores, which reward crisp short calls. Test truths were full matrices.
  - **Rule:** keep dev and test truths in the same form, and expect the dev margin to shrink.
- **Oracle-test every truth source, every time.**
  - **What happened:** a new truth source registered without the oracle test (test-factor consensus strings) made a
    whole confirmation uninformative. Its own matrix passed `contains` only 6% of the time.
  - **D1:** for bacterial regulons, the oracle test showed that IUPAC and PWM-similarity readouts fail outright. Only
    site recovery separated true from shuffled matrices, so it became the registered readout.
- **Diagnose before designing.**
  - **What happened:** the planned fix for the matrix gap was width search, because "motifs too short" was the
    story.
  - **What the anatomy showed:** a ten-minute anatomy (coverage, softening, phase blur, site count) refuted it. The
    gap was column accuracy, fixed in part by one more refit round.
  - **Rule:** measure the failure's anatomy before choosing the fix.
- **Run the cheapest decisive readout first.** The FR1 frontier against STREME needed only HASTUR and STREME outputs.
  Running it before HOMER finished answered the question about 5 hours early and saved HOMER's slowest arm (SELEX,
  up to 70 min per unit).
- **Small-N and large-N are different regimes.** Small sets (10–60 sequences) are where bacterial regulons live, and
  they were absent from every early panel. A tool's deployment sizes belong in the null panel from day 0.
- **Comparator output formats change with input size.** STREME reports `S=` (a training p-value) instead of `E=` when
  the input is too small for a hold-out set. A parser that reads only `E=` scores those runs as having no calls. Feed
  each comparator the smallest planned input and check its output parses.
- **Evaluators must discover tool outputs, not list them.** Hard-coded tool lists silently skipped new arms (LN5),
  and the report showed only the baseline. Discover output dirs, and fail loudly on an unknown tool.
- **Never let one job's cleanup touch another job's in-progress output.** A concurrent job's audit deleted an
  in-progress run (no completion marker yet) belonging to a different job. Audits remove only finished-but-incomplete
  outputs; anything without a completion marker is left alone.
- **`pkill -f` with a pattern that appears in your own command line kills your own shell** (twice in one day). Kill by
  PID or process group, found with a bracketed pattern (`[w]atch_box`) in a separate command.
- **"Fresh" panels must be checked for freshness, not trusted.** Additive seeds (`base + offset + rep`, or
  `+ 1000·rep`) collide across offsets. Consecutive "fresh" small-N panels shared 9 of 10 replicates, and a
  "confirmation" panel was 90% a dev panel.
  - **How it was caught:** identical false calls at shifted replicate indices.
  - **Rule:** seed with tuples (`default_rng([base, offset, rep, ...])`). Before calling a panel fresh, measure its
    sequence overlap with every earlier panel.
- **Simulated-null artefacts are real artefacts.** A uniform first base in the Markov null generator (a
  non-stationary start) produced the only residual false calls on simulated nulls. Stationary burn-in removed them.
  Validate the null generator as §3.6 asks of the simulator.
- **Public repos.** AWS identifiers (volume ids, availability zone, role names) were pushed once and scrubbed. Keep
  them in a gitignored `*.local.env`, and grep every commit for id patterns.
- **Persistent EBS for spot boxes.** A persistent EBS volume, attached by user data that installs its own CLI,
  identifies the device by volume id and refuses to format a non-blank disk, makes spot loss cost only in-flight
  work (`infra/box/BOX.md`).
