# Secondary analysis: repeat-masked FR1 confirmation (2026-10-12; not a registered confirmation)

**Why.** F2's false calls on the GC-selected human nulls are Alu elements (ROADMAP §5). The registered FR1 advantage
(+0.241) comes from those nulls. Repeat masking is the standard remedy.

**What.** The same confirmation windows, with soft-masked hg38 bases set to N:
- test2 ChIP, low-N panel and N = 5,000;
- the human real-genome nulls.

Each window is built with decisions made on the unmasked sequence, and the identity of the windows is verified.
SELEX and the bacterial nulls are not genomic-repeat data, so they are unchanged. Tools: frozen F2 (`hlat_f2`) and
STREME (`streme_all`). Readouts: the registered frontier cell, and the exchangeable-only cell.

**Caveat.** KRAB zinc-finger factors (many in test2) bind transposable elements, so masking removes some true sites.
Masked and unmasked results are reported side by side. Masking is an input choice, not a property of either tool.
