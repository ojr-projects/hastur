# Final reflections: why Comsequent did not make the grade

**Status:** written 2026-10-10, at the point where simulation development was frozen and the real-data programme had
returned its verdict on every phase. Phase 1's held-out confirmation was still running (5 of 25 strains), so the
one genuinely open question is flagged as such throughout.

**Who this is for:** us, and anyone who picks this repository up later. It assumes no knowledge of the sessions that
produced the numbers. Every figure is sourced to a script or a section of `REAL_DATA_PLAN.md` so it can be checked
or overturned.

---

## 1. The short version

Comsequent enumerates sparse k-mer patterns over a Möbius lattice, scores them with an exact log-linear
interaction statistic (θ) and an enrichment statistic (ENR), and reports a single merged ranked list. On simulated
data it beat every comparator we could run, often by a wide margin. On real data it beat none of them.

| ground | Comsequent | best comparator | verdict |
|---|---|---|---|
| simulated core grid (1,040 datasets) | **0.412 / 0.475** @1/@10 | 0.31 / 0.32 (STREME 0.30, HOMER 0.14) | **large win** |
| simulated dev panel (ctx2) | **0.59 / 0.68** | 0.37 / 0.39 | **large win** |
| simulated gapped, vs 8 gapped-capable tools | **+0.40 / +0.73 / +1.00** | — | **large win** |
| ChIP-seq (24 TFs) | 0.500 net | STREME 0.663 (HOMER 0.246) | loss to STREME, **beats HOMER 2:1** |
| RNA Bind-n-Seq (53 RBPs) | 0.372 net | HOMER 0.551 | loss |
| HT-SELEX (255 libraries) | 0.133 net | HOMER 0.291 | loss |
| bacterial promoters (16 species) | 0.188 net | STREME 0.238 | near-parity on −10, nothing on −35 |
| CAP-SELEX half-sites (135 libraries) | 0.170 / 0.067 | HOMER 0.504 / 0.444 | loss |
| CRISPR PAM depletion (12 samples) | — | — | **no valid readout** |
| CAP-SELEX composite sites | — | — | **no valid readout** |
| nanopore methylation (10 dev strains) | 0.707 | Nanomotif 0.714 | **parity** (see §4) |
| nanopore methylation (held-out, 5 of 25) | 0.857 lenient | Nanomotif 0.829 | **provisional**, see §4 |

The tool is not broken. It finds motifs, it is exactly calibrated where its assumptions hold, and it is fast on some
problems. It is simply not better than the established tools at the thing those tools are for; the one class where we
predicted it would dominate turned out to be its worst real-data class; and its assumptions do not hold on any real
dataset we ran.

What follows is why — eight reasons, roughly in order of how much each one cost us.

---

## 2. The central defect: it locates motifs but will not commit to them

This is the one finding that holds everywhere, and it is the thing to fix if anyone resumes this work.

Two scoring tiers, applied identically to every tool (`benchmarks/tools/registered_scoring.py`):

- **strict**: the call must cover every informative reference position with a letter set at least as specific as the
  reference's (`match_ic`).
- **located**: the call must agree with the reference wherever both specify, in either direction of specificity, with
  no conflicts (`match_contains`). A call that carries the truth but states it too vaguely passes here and fails
  strict.

| phase | ours, strict | ours, located | **precision gap** | best comparator, strict | its gap |
|---|---|---|---|---|---|
| RBNS | 0.151 | **0.453** | **0.302** | HOMER 0.604 | 0.000 |
| HT-SELEX | 0.133 | 0.227 | 0.094 | HOMER 0.357 | 0.047 |
| ChIP | 0.500 | 0.583 | 0.083 | STREME 0.667 | 0.000 |
| promoters (−10) | 0.188 | 0.188 | 0.000 | STREME 0.250 | 0.000 |

**The asymmetry is one-sided.** STREME and HOMER have precision gaps of 0.000–0.047: when they find a motif they
state it at full specificity. Ours run to 0.302. On RBNS we *locate* 0.453 of the truths and state only 0.151 of them
exactly — two thirds of the apparent gap to HOMER is vagueness, not blindness.

A continuous metric says the same (E4, `pwm_similarity`: call and reference as probability matrices, mean
1 − total-variation distance over the reference's columns at the best offset, with letter-shuffled truths as the
null):

| phase | null | STREME | HOMER | best Comsequent | our share of the best net |
|---|---|---|---|---|---|
| ChIP | 0.501 | 0.938 | 0.844 | **0.858** | **82%** (and above HOMER) |
| RBNS | 0.615 | 0.872 | 0.886 | **0.827** | **78%** |

On exact match we reach 25% of HOMER's RBNS net and 75% of STREME's ChIP net. On representation, 78% and 82%. The
lattice is locating these motifs in roughly the region of sequence space the PWM tools do; it fails at collapsing
that region into one committed string.

**The mechanism is measurable.** The share of top-10 calls carrying under 6 bits of information is 0.19–0.58 for
Comsequent arms against 0.00–0.02 for STREME and HOMER. And the arms trade specificity against coverage rather than
achieving both: on RBNS the Markov-background arm (`cs_mkv`) has a precision gap of **0.000** at 0.377 strict, while
the all-in-one arm has **0.302** at 0.151 strict but locates 0.453. One states crisply and finds less; the other
finds more and states it vaguely. **The all-in-one ranking prefers the vague one.**

So the defect lives in the readout — the ranking rule and the letter-collapse step — not in the enumeration or the
statistics. That is the most tractable negative finding here, and it is the only one that points at a specific fix.

---

## 3. The simulator was ours, and it certified our own assumptions

Every win was on data we generated. The simulator encoded a model of what a motif *is* — a discrete IUPAC pattern
planted at a rate, against a background — which is precisely the model the lattice is built to invert. It was not a
test of the method so much as a consistency check on it, and we treated it as a test for far too long.

The sharpest single illustration is bipartite (gapped) motifs, where simulation promised the most and real data
delivered the least:

| | Comsequent | Nanomotif / best comparator | margin |
|---|---|---|---|
| **simulated** bipartite (ctx2 panel) | 1.00 @1 | 0.47 / 0.80 | **+0.20 to +0.53** |
| **simulated** bipartite (vs Nanomotif) | 0.70 | 0.56 | **+0.14** |
| **simulated** gapped vs 8 gapped-capable tools | — | — | **+0.40 / +0.73 / +1.00** |
| **real** bipartite methylation, dev 40× (45 motifs) | 0.533 | **0.800** | **−0.267** |
| **real** bipartite methylation, dev 40×, exact rule | 0.378 | **0.800** | **−0.422** |
| **real** bipartite methylation, held-out (10 motifs) | 0.800 | **1.000** | **−0.200** |
| **real** gapped HT-SELEX truths (31 libraries) | 0.065 net | HOMER 0.113 | −0.048 |

A swing of roughly 0.41 on the same motif class between our simulator and real data. Our best real-data methylation
strata are the *contiguous* ones (non-palindrome +0.091, palindrome +0.154 lenient on held-out), the exact opposite
of the simulated profile.

Why the simulation was wrong about gaps: a planted `GCNNNANNNNGC` has an exactly constant gap, exactly uniform
filler and a clean planting rate. Real bipartite methylation motifs have variable spacing, biased filler and partial
methylation, and Nanomotif models their structure explicitly. We built a lattice that enumerates the idealisation
and a simulator that produces it, then measured one against the other.

**The lesson is not "simulation is useless".** It is that a simulator written by the same people, from the same model
of the problem, as the method under test cannot adjudicate that method. It can debug it. Ours was treated as
evidence, and it carried the project for months on a margin that did not exist.

---

## 4. The one real-data win was a property of which configuration we quoted

Phase 1 (nanopore methylation) was the designated win: discrete enzymatic specificity, a contrast matched by
construction (methylated vs unmethylated positions in the same genome), and the best signal-to-drift ratio of any
phase (3.5, against HT-SELEX's 0.39). The reported result was 0.79 against Nanomotif's 0.71.

It reproduces exactly — but only in one configuration. Scored on the **same seven strains**, same truth, same rule,
same denominators, with Nanomotif's figure identical in both rows:

| configuration | Comsequent | Nanomotif | margin |
|---|---|---|---|
| `r10m4` (the configuration quoted) | **0.799** | 0.714 | **+0.085** |
| `rwcf` (cross-fitted walk) | **0.707** | 0.714 | **−0.007** |

Cross-fitting — which exists precisely to stop the walk fitting its own selection — costs 0.092 and erases the
margin. The three extra strains in the `rwcf` run are not the cause; restricting to the shared seven isolates the
configuration as the only variable.

We had registered the *protocol* (strains, truth, comparator, rule) but not the *configuration*. With several arms
and several walk variants live at once, the reported number drifted to whichever combination scored best, without
anyone choosing it. That is selection without intent, and the registry did not catch it because it was not watching
that degree of freedom.

**The held-out run is incomplete and is the one open question.** At 5 of 25 strains, under four rules applied
identically to both tools: the four pre-amendment strains have Nanomotif ahead by 0.069–0.138 under *every* rule,
while strain 5 (*Sphingomonas echinoides*, 6 motifs) went 6/6 to us against its 2/6 and swings the pooled figure by
0.13 on its own. Pooled, 0.857 vs 0.829 lenient and 0.571 vs 0.714 exact. At n = 6 strain 5 is noise. **The remaining
20 strains decide whether Phase 1 is parity or a loss; they are very unlikely to make it a win.**

---

## 5. Two phases produced numbers that measured nothing

We ran both to completion and reported both before checking whether their readouts had any power. Neither did.

**CRISPR PAM depletion (P7).** The counts reference — which *is* the truth by construction — matches letter-shuffled
truths **0.942** of the time under the uniform rule and **0.992** under the phase's own rule. A PAM reference like
`DGVG` is 4 nt of which `D` and `V` carry 0.415 bits each, leaving **two informative positions**, and two positions
cannot distinguish a motif from a permutation of itself. HOMER's apparent 0.75 and our 0.50 were both almost
entirely chance. This was reported as "parity-to-loss" for a while. It is neither: it is unmeasurable.

**CAP-SELEX composite sites (P4b).** One composite hit in roughly **21,000 calls** across every arm and both runs —
HOMER, on a single library. The rule requires two non-overlapping JASPAR half-sites inside one pattern, needing
~15–25 bp; HOMER's median call specifies 12 positions and still never contains both. A rule no tool can reach has no
power to separate them.

These are mirror images — one readout accepts everything, the other accepts nothing — and both are disqualifying.
The fix is cheap and we should have had it from the start: **score shuffled truths alongside real ones, and refuse to
state an outcome for any phase where the chance rate exceeds a threshold or where the oracle itself fails.** That is
now `registered_scoring.py`'s E1 baseline, implemented after the fact.

Related: once each phase's own registered rule was given a chance baseline, **three of five could not separate
tools** — RBNS (STREME chance 0.626), promoters (HOMER 0.581) and PAM (0.992). Only ChIP and HT-SELEX were sound
under both their own rule and the uniform one. We had been comparing tools through rules we had never validated.

---

## 6. The comparisons were asymmetric in our favour, repeatedly

None of these was deliberate. All of them pointed the same way, which is the part worth noticing.

- **Slot asymmetry.** Comsequent's reported recall was the union of its arms' top-10 lists — up to 30 slots — against
  Nanomotif's single list (median 4–7 calls). Collapsing the arms into one de-duplicated list of 10 removes this.
- **Rule asymmetry.** Our calls were matched with a lenient occurrence-based rule while the comparator was matched on
  exact rank. Giving Nanomotif the same lenient rule lifts it 0.793 → 0.931 (+0.138) against our 0.793 → 0.828
  (+0.035). It gains nearly four times as much, because our advantage under lenient matching *was* the asymmetry.
- **Denominator asymmetry.** Rows where a tool produced no calls were dropped from that tool's denominator rather
  than counted as misses. Nanomotif returns nothing for 15 of 207 dev rows, so this silently added ~7% to its
  measured recall — in that instance against us, which is how it was found.
- **Offset shopping.** RBNS's registered rule accepted 4 of 5 positions with one conflict tolerated at any offset on
  either strand. Under it our headline was 0.792; under the uniform rule, 0.151. The ordering did not change (we lose
  either way), but a rule that generous mostly measures itself: its chance rate is 0.35–0.63.

The general shape: a method with many arms, many parameters and a flexible output format will, in the hands of people
who want it to work, be compared in whatever way shows it best, unless the comparison is fixed in advance and applied
mechanically to everyone. Ours was not, until late.

---

## 7. Two harness faults manufactured results, and the readout never flagged either

- **PAM run through the wrong command.** The first PAM pass used `discover` (offset-agnostic) on data where position
  relative to the protospacer is the entire signal, scoring 0.08. The `anchored` command on identical data scored
  0.50. A 6× swing from a command-line choice.
- **The factorial arm could not express what it was scored on.** `43_capselex_factorial.py` runs only
  `discover -k 5 -w 13` and scores the raw output, while the plain arm runs `discover` *and* `extend -kmax 16 -wmax
  30 -join 30`. Maximum call length: **13 bp** for the factorial arm against 25 (plain), 20 (HOMER), 16 (STREME). A
  composite site needs ~15–25 bp. Its `composite_top10 = 0.000` across 135 libraries carries no information about
  the factorial design at all, and the pre-registered prediction it was built to test is **untested**, not failed.

Both were found by reading the code, not by anything in the output. A 0.000 and a 0.08 look exactly like negative
results. **The pipeline had no way to tell a method failure from a harness failure**, which means every null it
produced needed manual verification to be worth anything — and nulls are most of what a project like this generates.

---

## 8. The deeper reason: the niche the method fits is already better served

Strip away the process failures and a substantive problem remains.

The lattice's premise is that biological specificity is **discrete and exact**: a site either matches `GATC` or it
does not, and the interesting structure is the interaction between positions, recoverable by Möbius inversion over
the pattern lattice. Where that premise holds — restriction and methyltransferase recognition, enzymatic
specificity — the premise is right. Where it does not — transcription factor binding, which is graded and energetic
and genuinely well described by a position weight matrix — exhaustive exact enumeration buys nothing a PWM method
does not already get, and costs the commitment to a single string that §2 shows we cannot make.

That alone would leave a real niche. The problem is that **the niche is narrow and already occupied**. In the one
place our premise holds and the contrast is matched by construction, Nanomotif — which models bipartite methylation
motif structure explicitly, and uses the modification-probability signal directly rather than through a
discretisation — matches or beats us, and beats us worst on exactly the structured motifs we were built for (§3).
A general lattice over discrete patterns loses to a specialised model of the actual generating process, which is
not surprising in retrospect and would have been worth predicting in advance.

The honest statement of the domain of applicability is therefore much smaller than the one we developed: not
"discrete specificity with a matched contrast", which we thought was validated, but something closer to "discrete
specificity with a matched contrast, *and no specialised tool for that assay*". That is a thin claim on which to
rest a method paper.

A secondary blow: the "enough motif structure" clause of our domain rule is **falsified** (§3), and the
signal-to-drift ratio that was supposed to predict wins (> ~1) was fitted post hoc to five phases and has never been
confirmed on a locked phase.

The prediction registry (§0a of `REAL_DATA_PLAN.md`) divides cleanly, and the split is informative: **every
prediction of a loss held, and every prediction of a win or parity failed.**

| prediction | phases | outcome |
|---|---|---|
| **loss** | P5 ChIP (pre-registered), P6 promoters | **correct** |
| **parity** | P3 RBNS, P4a HT-SELEX, P4b CAP-SELEX half-sites (pre-registered) | all **losses** |
| **win** | P1 methylation, P1b REBASE | parity, configuration-dependent (§4) |
| **win or parity** | P7 PAM | **no valid readout** (§5) |

So the domain model was not useless — it reliably identified where we would lose. It was systematically optimistic
wherever it expected us to do well, which is the more dangerous failure, and the only *pre-registered* prediction it
got right was a prediction against ourselves.

---

## 9. The significance machinery dissolves exactly where real data lives

This belongs in the list of reasons, and it was found while checking a claim in the next section, which is why it directly
precedes it.

The same null calibration run tells two very different stories depending on how the background was built
(`14_allinone_null.py`; motif-free pairs, so *every* significant call is a false positive):

| null design | significance rate | patterns passing BH at 0.05 |
|---|---|---|
| i.i.d. split of one sequence set | **0.000** | 0 |
| codon-preserving split | 0.007 | 0.3 |
| Markov-2 split | 0.013 | 0.3 |
| GC-content shift | **1.000** | **20,000 (all of them)** |
| different species | **1.000** | **20,000** |
| order mismatch | **1.000** | **20,000** |
| coding vs Markov background | **1.000** | 12,557 |

When foreground and background come from the same distribution, the statistics are exact and the gate holds. **Under
any composition mismatch it saturates: every call is declared significant and tens of thousands of patterns pass
multiple-testing correction.**

Real foreground/background pairs are *never* drawn from the same distribution. ChIP peaks differ from genomic
background in GC and repeat content; SELEX cycles differ by PCR and selection drift; promoter sets differ from coding
sequence in every compositional respect. The measured drift in our own phases runs from KL₃ 0.019 (ChIP) to 0.211
(HT-SELEX). So on every real dataset we ran, the `sig` column was somewhere between uninformative and meaningless,
and the ranking it was meant to guard was being driven partly by composition rather than by the motif.

This is the mechanism behind the HT-SELEX result (drift KL₃ 0.104, five times methylation's, our worst phase) and it
explains why methylation was the best case: an anchored contrast within one genome is the one design where
foreground and background genuinely *are* exchangeable apart from the modification. The method's exactness is
conditional on an assumption that only a matched anchored contrast satisfies — which is a much sharper statement of
the domain limit than §8's, and it was available from a simulation we had already run.

The remedy is not cosmetic. It means either a background model rich enough to absorb composition (we tried a Markov
background arm: `cs_mkv`, which is indeed our best-calibrated and best-stated arm on RBNS, with a precision gap of
0.000 — a hint worth following), or a factorial/difference-in-difference design that cancels shared drift (built, but
never validly tested, §7). Both point the same way: **the exact lattice statistics need a drift-robust contrast to
mean anything, and supplying one is a harder problem than the enumeration we spent the project optimising.**

## 10. What does survive

Stated conservatively, because everything above is a reason to distrust a positive claim here.

1. **ChIP-seq: we beat HOMER about 2:1.** Net 0.500 against 0.246, chance-corrected, with the same ordering under
   both the uniform rule and the phase's own rule, and 0.858 against 0.844 on the representation metric. STREME
   beats us (0.663). "Competitive on ChIP, ahead of one of the two standard tools" is defensible and robust.
2. **Speed on short unanchored problems.** RBNS: 22 s against HOMER's 33 s and STREME's 109 s. Not general —
   methylation runs 871 s per strain against Nanomotif's 166 s.
3. **Calibration under an exchangeable background.** On motif-free pairs drawn from the *same* distribution the
   `sig` gate is essentially exact: significance rate 0.000 on an i.i.d. split, 0.013 on a Markov-2 split, 0.007 on
   a codon split (`14_allinone_null.py`, 1,077 calls over 7 null designs × 3 replicates). This is a real property
   but a narrow one — see §9 for why it does not transfer.
4. **The evaluation infrastructure.** The chance baseline, the rule-sensitivity analysis, the precision tier, the
   representation metric and the readout-validity criterion caught a retracted phase, an invalid composite readout,
   three non-separating rules, a configuration-specific win and two harness faults — in our own work, before review.
   This is reusable by anyone benchmarking motif finders, and it is arguably the most transferable thing the project
   produced.

---

## 11. What we would do differently

1. **Validate the readout before running the method.** Score shuffled truths first. If the oracle cannot beat its own
   permutations, there is no experiment. This would have killed PAM and the CAP-SELEX composite readout before either
   consumed a box-week.
2. **Never let your own simulator adjudicate your own method.** Use it to debug. Require at least one real,
   independently-generated dataset before claiming a capability — and specifically before claiming a *structural*
   capability like gapped motifs, which is where an idealised generator flatters you most.
3. **Register the configuration, not just the protocol.** Binary, flags, arm set, walk variant — pinned by hash,
   before the run. With several live arms, "the method" is not a well-defined object and the reported number will
   drift to the best one without anyone deciding.
4. **Fix the comparison and apply it mechanically.** Same slots, same rule, same denominators, missing output counts
   as a miss, both tools' own call lists. Write it as one function all phases call. Every asymmetry we found pointed
   our way.
5. **Hold out before developing, not after.** Phase 1 became a development set through a series of individually
   reasonable tuning decisions. The held-out set was declared only after an outside critic pointed this out, and the
   first results from it do not reproduce the development margin.
6. **Report the located/stated split from the start.** Binary recall conflates "did not find it" with "found it and
   described it loosely". These have different fixes, and we spent months of development on ranking rules without
   knowing which one we were failing at.
7. **Build harness self-checks.** Assert that each arm can in principle express the target (length, alphabet,
   positional capability) before scoring it. Both manufactured results in §7 would have tripped a two-line assertion.

---

## 12. What is still genuinely unknown

In fairness to the method, and so nobody over-reads this document:

- **Phase 1's held-out verdict.** 5 of 25 strains. Parity or a loss, most likely; the data will say.
- **The factorial interaction arm.** Never tested (§7). Testing it needs the arm given a walk *and* a readout that
  scores half-sites and spacing separately rather than demanding both in one pattern. θ-interaction between two
  TFs' contributions is a question no comparator answers, so this is the one place a novel claim could still live.
- **Whether the representation defect is fixable.** §2 localises it to the ranking and letter-collapse steps and
  shows a 0.08–0.30 recoverable margin if calls could be stated as specifically as they are located. Nobody has
  tried. Unlike the u8–u20 sweep that plateaued on trade-offs, this target has a named cause, a measured size and a
  comparator that demonstrably does not share it.
- **Whether a drift-robust contrast rescues the statistics.** §9 narrows this to a concrete question: the
  Markov-background arm is both our best-stated and best-calibrated arm on RBNS. Nobody has tested whether a
  properly-run factorial design, or a richer background model, restores the `sig` gate on real data.

If anyone resumes this: start at §2, not at the lattice.
