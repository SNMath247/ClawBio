# Examples — aso-off-target-screening

Regenerate the reference outputs any time with:

```bash
python skills/aso-off-target-screening/aso_off_target_screening.py --demo --output /tmp/aso_demo
# or
python clawbio.py run aso-screen --demo --output /tmp/aso_demo
```

## What the demo demonstrates

The synthetic target (`../demo_target.json`) is a muscle (`Muscle_Skeletal`)
cryptic-exon pre-mRNA carrying a patient SNV (T>A@60). 261 candidate 15–20mer
ASOs are tiled across the target window and pushed through the full funnel.

| Funnel stage | Demo evidence |
|---|---|
| Tissue-weighted penalty | `ESSGENE` exonic off-target is the top liability (penalty ≈ 10.6, TPM 55.3, essential ×3.4) |
| Tissue zeroing | `SILENTGENE` (TPM 0 in muscle) contributes exactly 0 — `log2(0+1)=0` |
| Pre-mRNA / intron search | `INTRONGENE` **intronic** off-target overlaps a PTBP1 CLIP peak → collision liability |
| Intended-vs-unintended CLIP | on-target overlaps a `HNRNPA1` mechanistic silencer → asset, not a liability |
| Sequence-intrinsic tox | G-quadruplex / CpG candidates are excluded even at zero off-target penalty |
| Multi-objective arbitration | a Pareto front (not a single scalar) of clean, mechanism-engaging 15-mers |

Expected top-line summary: `261 candidates, ~84 excluded, ~132 intended-mechanism,
2 Pareto-optimal`. The report leads with a **"triage, not clearance"** caveat
because a large fraction of PS-ASO toxicity is hybridization-independent and
therefore invisible to sequence screening.
