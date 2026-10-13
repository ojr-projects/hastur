# LN6 round 2, amendment 1 (2026-10-11, before any round-2 result is read)

The preliminary FR1 readout (`results/frontier_prelim.txt`) showed that frozen v8 is miscalibrated on the real-genome
nulls even at N = 500–2,000, where STREME makes 0.00 false calls per run:
- 0.45–0.75 false calls per run on bacterial coding windows;
- 1.05 on human random windows at N = 500, and 0.80 at N = 2,000.

So round 2 is also gated on these nulls.

**Added gate:** on the exchangeable real-genome nulls of the FR1 dev panel (`runs/realnull`), native false calls per
run must be ≤ 0.1 pooled:
- `human_exch_n500` and `human_exch_n2000`;
- `{bceno, cjejuni, ecoli}_coding`.
That is 100 runs in all.

The intergenic panels are excluded because of a construction flaw: 100-bp windows were sampled with overlap from
small intergenic spaces, and 376 of 500 C. jejuni foreground windows share a 20-mer with another foreground window.
They are rebuilt by non-overlapping tiling before any further use. GC-selected human sets are a composition
difference, not an exchangeable null, so they are reported but not gated.
