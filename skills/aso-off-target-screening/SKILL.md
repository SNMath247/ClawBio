---
name: aso-off-target-screening
description: >-
  Deterministic multi-stage off-target and toxicity screening funnel for
  antisense oligonucleotide (ASO) candidates — tiles 15-20mers across an
  in-silico mutant transcript, then applies a tissue-weighted off-target
  penalty, CLIP-collision oracle, and sequence-intrinsic tox motifs, surfacing
  a Pareto front (never a single safety scalar).
license: MIT
metadata:
  version: "0.1.0"
  author: ClawBio
  domain: oligonucleotide-therapeutics
  tags:
    - antisense-oligonucleotide
    - aso
    - off-target
    - toxicity
    - splice-switching
    - gtex
    - clip
    - depmap
  inputs:
    - name: input_file
      type: file
      format:
        - json
      description: >-
        Target spec — baseline pre-mRNA (introns included), optional patient
        SNV, target window, tissue, and modality
      required: true
  outputs:
    - name: report
      type: file
      format:
        - md
      description: Ranked screening report with the Pareto front and per-axis liabilities
    - name: result
      type: file
      format:
        - json
      description: Machine-readable per-candidate scores and off-target hits
  dependencies:
    python: ">=3.10"
    packages:
  demo_data:
    - path: demo_target.json
      description: Synthetic cryptic-exon target + bundled off-target/GTEx/CLIP/DepMap tracks
  endpoints:
    cli: python skills/aso-off-target-screening/aso_off_target_screening.py --input {input_file} --output {output_dir}
  openclaw:
    requires:
      bins:
        - python3
    always: false
    emoji: "🧬"
    homepage: https://github.com/ClawBio/ClawBio
    os:
      - darwin
      - linux
    install:
    trigger_keywords:
      - ASO off-target screening
      - antisense oligonucleotide off-target
      - splice-switching ASO design
      - ASO toxicity screen
      - gapmer off-target
      - weighted off-target score
---

# 🧬 ASO Off-Target & Toxicity Screening

You are **ASO Off-Target & Toxicity Screening**, a specialised ClawBio agent that triages antisense-oligonucleotide (ASO) candidates through a deterministic, computationally-cheapest-first funnel and reports a multi-objective Pareto front — never a single "safety" score.

## Trigger

**Fire this skill when the user says any of:**
- "screen ASO candidates for off-targets"
- "antisense oligonucleotide off-target / toxicity screening"
- "design / triage a splice-switching ASO for a cryptic exon"
- "which ASO / gapmer has the lowest off-target liability"
- "tissue-weighted off-target score" or "weighted off-target penalty"
- "tile 15-20mers across a transcript and rank them"
- "check my ASO for CpG / G-quadruplex / hepatotoxic / CLIP-collision liabilities"

**Do NOT fire when:**
- The user wants siRNA seed-region design only with no ASO/gapmer context → generic RNA tools.
- The user wants CRISPR guide off-targets → this is a nucleic-acid *hybridization* screen, not a Cas nuclease screen.
- The user wants variant pathogenicity / ACMG calls → `clinical-variant-reporter` / `cnv-acmg-classifier`.
- The user wants clinical trial matching or drug pharmacology → other skills.

**Design notes:** Loud triggers on "ASO", "antisense oligonucleotide", "gapmer", "splice-switching", "off-target", and "weighted off-target score". This skill is the only ClawBio skill that scores oligo hybridization liabilities against tissue expression, CLIP peaks, and gene essentiality together.

## Why This Exists

- **Without it**: Running deep-learning models on every 15-20mer across a transcript is expensive and unstructured; naive BLAST off-target lists ignore tissue expression, gene essentiality, intronic pre-mRNA collisions, and sequence-intrinsic immunogenicity, so teams over- or under-flag candidates.
- **With it**: A single deterministic funnel tiles candidates, ranks off-targets by a tissue- and essentiality-weighted penalty, flags RBP CLIP collisions on the **pre-mRNA (introns included)**, screens for toxic/immunogenic motifs, and surfaces a defensible Pareto front for human review — in seconds, fully locally.
- **Why ClawBio**: Every score traces to a documented formula and to bundled table schemas that map 1:1 onto real production oracles (GTEx v10 API, ENCODE/POSTAR3 eCLIP, DepMap Chronos, gnomAD LOEUF). Nothing is hallucinated.

## Core Capabilities

1. **In-silico mutant transcript**: merge a baseline pre-mRNA with a patient SNV.
2. **Candidate tiling**: enumerate 15-20mer antisense candidates across the target window.
3. **Pre-mRNA off-target alignment**: modality-aware near-match search over a transcript atlas that **includes introns**.
4. **Tissue-weighted off-target penalty**: `sum_i P(bind_i) · log2(TPM_tissue,i + 1) · E(essentiality_i)`.
5. **CLIP-collision oracle**: genomic-interval intersect vs RBP peaks, separating **intended mechanism** (asset) from **unintended collision** (liability).
6. **Sequence-intrinsic tox**: CpG/TLR9, G-quadruplex, hepatotoxic and GU-rich motif screen on the ASO itself, with chemistry-mitigation flags.
7. **Multi-objective arbitration**: a Pareto front over the axes; no single scalar safety score.

## Scope

**One skill, one task.** This skill triages a supplied ASO target and its candidate panel against **bundled or user-supplied reference tracks**. It does not: run deep-learning phenotype predictors, call variants, fetch live GTEx/DepMap/ENCODE data, model phosphorothioate backbone chemistry, or clear a candidate for the clinic. It is an in-silico triage layer that feeds wet-lab, not a replacement for it.

## Input Formats

| Format | Extension | Required Fields | Example |
|--------|-----------|-----------------|---------|
| Target spec | `.json` | `tissue`, `pre_mrna`, `target_window` | `demo_target.json` |

Optional spec fields: `variant` (SNV `{pos,ref,alt}` to build the mutant), `modality` (`gapmer` \| `steric-blocker`), `aso_length_min`/`max` (default 15/20), `chrom`, `strand`, `pre_mrna_genomic_start`, `max_mismatch`. Reference tracks default to `data/` and follow these schemas: `gtex_median_tpm.csv` (gene × tissue), `essentiality.csv` (gene,chronos,loeuf,common_essential), `clip_peaks.csv` (chrom,start,end,strand,rbp,cell_line,idr,mechanistic_silencer), `offtarget_atlas.json` (transcripts with genomic_start + pre-mRNA sequence), `tox_motifs.json`.

## Workflow

1. **Validate** (prescriptive): confirm the spec has `tissue`, `pre_mrna` (ACGT only), and `target_window`; fail cleanly on missing/invalid fields.
2. **Build mutant** (prescriptive): apply the SNV to the baseline pre-mRNA; verify the `ref` base matches before substituting.
3. **Tile** (prescriptive): enumerate every 15-20mer antisense candidate across the window; each ASO is the reverse complement of its sense target site.
4. **Align off-targets** (prescriptive): scan the atlas (introns included) for near-matches at the modality-specific mismatch ceiling (gapmer ≤3, steric-blocker ≤1).
5. **Weight** (prescriptive): score each hit `P(bind)·log2(TPM+1)·E`; sum per candidate.
6. **CLIP intersect** (prescriptive): intersect on-target and off-target footprints with CLIP peaks; classify intended-mechanism vs unintended collision.
7. **Tox scan** (prescriptive): count CpG/G-quad/hepatotoxic/GU-rich motifs on the ASO; attach chemistry-mitigation flags.
8. **Arbitrate** (flexible narrative): compute the Pareto front, then explain the trade-offs and next wet-lab steps.

## CLI Reference

```bash
# Standard usage
python skills/aso-off-target-screening/aso_off_target_screening.py \
  --input <target_spec.json> --output <report_dir>

# Demo mode (synthetic data, no user files needed)
python skills/aso-off-target-screening/aso_off_target_screening.py --demo --output /tmp/aso_demo

# Overrides
python skills/aso-off-target-screening/aso_off_target_screening.py \
  --input spec.json --tissue Brain_Cortex --modality steric-blocker --max-mismatch 1 --output <dir>

# Via ClawBio runner
python clawbio.py run aso-screen --demo
python clawbio.py run aso-screen --input spec.json --output <dir>
```

## Demo

```bash
python clawbio.py run aso-screen --demo
```

Expected output: a report over **261** tiled 15-20mers against a synthetic muscle (`Muscle_Skeletal`) cryptic-exon target. It flags an ESSGENE exonic off-target as the top weighted liability, **zeroes** a silent-in-muscle off-target (`log2(0+1)=0`), catches an intronic PTBP1 CLIP collision, excludes G-quadruplex/CpG candidates, and surfaces a Pareto front of clean 15-mers that engage the intended mechanistic silencer.

## Algorithm / Methodology

An agent can apply this without the script:

1. **Mutant transcript** = baseline pre-mRNA with the SNV applied (single-base, `ref` must match).
2. **Tiling**: for each length L∈[15,20], for each offset in the window, `ASO = reverse_complement(mutant[offset:offset+L])`; the ASO's genomic footprint is colinear with the pre-mRNA.
3. **P(bind)** = `identity³ · (0.5 + 0.5·GC)` where `identity = (L − mismatches)/L`. Monotonically decreasing in mismatches, increasing in GC (a thermodynamic proxy for ΔΔG°37 ranking).
4. **E(essentiality)** = `1 + max(0,−Chronos) + max(0, 0.6−LOEUF)/0.6`. A neutral gene returns exactly 1.0; essential/constrained genes multiply the penalty.
5. **Off-target penalty** = `Σ_i P(bind_i)·log2(TPM_tissue,i + 1)·E_i` over off-target hits (the intended on-target is scored separately).
6. **CLIP**: half-open genomic-interval intersect, strand-matched; IDR peaks weight 1.0, relaxed 0.5. An on-target overlap with a `mechanistic_silencer` peak is an **asset**, not a liability.
7. **Pareto front**: non-dominated candidates over benefits `{efficacy, intended_mechanism}` and liabilities `{offtarget_penalty, clip_collision, essentiality_hit, seq_tox}`.

**Key thresholds / parameters**:
- Gapmer off-target mismatch ceiling: ≤3 (source: Kamola et al., NAR 2015 — validated exonic + intronic off-targets at 2–3 mismatches).
- "Not expressed" TPM: <1 (field convention); exactly 0 → penalty contribution 0.
- LOEUF loss-of-function-intolerant regime: <0.6 (gnomAD v4).
- DepMap common-essential proxy: Chronos ≤ −0.5.
- Hepatotoxic trinucleotides: TGC / TCC (source: Burdick et al., NAR 2014). CpG→TLR9 (Bauer 2001). G-quadruplex: ≥4 consecutive G.

## Example Queries

- "Screen these 15-20mer ASOs against skeletal muscle for off-target and tox liabilities."
- "My patient has a deep-intronic SNV creating a cryptic exon — tile and rank splice-switching ASOs."
- "Compute the weighted off-target penalty for a gapmer and show me the Pareto front."
- "Does any candidate collide with a PTBP1 CLIP peak on the pre-mRNA?"

## Example Output

```markdown
# ASO Off-Target & Toxicity Screening Report

**Target**: SYN-CRYPTIC-01 (SYNTGT) — tissue Muscle_Skeletal, modality gapmer
**Candidates tiled**: 261 | **Pareto-optimal**: 2 | **Excluded**: 84

**Penalty formula**: `sum_i P(bind_i) * log2(TPM_tissue,i + 1) * E(essentiality_i)`

## Pareto-optimal candidates (non-dominated)

| Rank | ASO ID | ASO (5'->3') | Efficacy | Off-tgt penalty | CLIP collision | Ess. hit | Seq tox | Intended mech. | Verdict |
|---:|---|---|---:|---:|---:|---:|---:|:---:|---|
| 1 | ASO_40_15 | `CTAGAGCCCTGTTGG` | 0.9 | 0.0 | 0.0 | 0.0 | 0.0 | yes | prioritise-for-wet-lab |

## Highest weighted off-target liability
`ASO_45_18` — penalty 10.58, driven by ESSGENE (exon, 2 mismatch, TPM 55.3).
```

## Output Structure

```
output_directory/
├── report.md                     # Ranked report + Pareto front + interpretation
├── result.json                   # Per-candidate scores, off-target hits, provenance
├── tables/
│   ├── candidates.csv            # Every tiled ASO with per-axis scores + verdict
│   ├── offtarget_hits.csv        # Per-hit alignment: gene, region, mismatches, P(bind), TPM, contribution
│   └── pareto_front.csv          # Ranked non-dominated candidates
└── reproducibility/
    └── commands.sh               # Exact command to reproduce
```

<!-- Output contract: the demo run (aso_off_target_screening.py --demo) writes every
     file listed above. No figures are produced (stdlib-only, deterministic). -->

## Dependencies

**Required**: Python 3.10+ standard library only. No third-party packages, no network.

**Optional (production swap-in, not required by this skill)**: GTEx Portal API v2, ENCODE/POSTAR3 eCLIP BED tracks, DepMap Chronos CSV, gnomAD LOEUF — each maps onto the bundled table schemas above.

## Gotchas

- **You will want to search only the spliced mRNA. Do not.** RNase H1 acts on nuclear pre-mRNA; the atlas and CLIP intersect **must include introns**. Many real off-target toxicity events are intronic cryptic-splice disruptions on unrelated genes (Kamola 2015). The demo deliberately seeds an intronic PTBP1 collision.
- **You will want to treat any CLIP overlap as a liability. Do not.** An on-target overlap with a mechanistic silencer is the *mechanism* (nusinersen works BY blocking hnRNP A1 at ISS-N1). Only **off-target** footprint overlaps are collisions; on-target mechanistic overlaps are surfaced as an asset benefit axis.
- **You will want to down-weight or drop off-targets by bulk TPM. Be careful.** Bulk GTEx masks rare cell types; a gene essential in a rare neuron/satellite cell can show low bulk TPM. Never hard-filter on bulk TPM — it only *weights* the penalty here, and exactly-zero TPM (untranscribed) is the only value that fully zeroes a hit.
- **You will want to emit one "safety score". Do not.** A large fraction of PS-ASO toxicity is hybridization-**independent** (backbone protein binding, complement, thrombocytopenia) and is invisible to sequence screening. This skill reports a Pareto front and an explicit "triage, not clearance" caveat; it never claims a candidate is safe.
- **Cancer-cell essentiality ≠ neuron/muscle essentiality.** DepMap Chronos is a proliferating-cancer signal; pair it with gnomAD LOEUF (bundled) and, in production, neuron/muscle CRISPRi screens.

## Safety

- **Local-first**: No network calls, no data upload; all processing is deterministic and local.
- **Disclaimer**: Every report includes the ClawBio medical disclaimer — *"ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions."*
- **Triage, not clearance**: The report states plainly that in-silico screening cannot clear a candidate on safety and that wet-lab tolerability cascades are mandatory.
- **Audit trail**: The exact command is written to `reproducibility/commands.sh`; every threshold traces to a cited source.
- **No hallucinated science**: Parameters trace to the formula and cited databases; the agent must not invent gene-drug or gene-essentiality associations.

## Agent Boundary

The agent (LLM) dispatches, explains trade-offs, and recommends next steps. The skill (Python) tiles, aligns, scores, and writes outputs. The agent must NOT override the mismatch ceilings, penalty weights, or verdict thresholds, nor claim a candidate is clinically safe.

## Integration with Bio Orchestrator

**Trigger conditions**: the orchestrator routes here on "ASO / antisense oligonucleotide off-target", "splice-switching ASO", "gapmer", "weighted off-target score", or a target-spec JSON with `pre_mrna` + `target_window`.

## Chaining Partners

- `variant-annotation` / `gi-splice`: upstream — identify the cryptic-exon-creating variant and splice consequence that defines the target window.
- `omics-target-evidence-mapper` / `target-validation-scorer`: cross-reference off-target genes for essentiality/druggability context.
- `struct-predictor`: downstream — model an ASO:RNA or protein interaction for a prioritised candidate.

> Output is structured JSON + CSV with stable headers, so it chains cleanly into downstream evidence aggregation.

## Maintenance

- **Review cadence**: Re-evaluate quarterly, or whenever GTEx (v10→), ENCODE eCLIP, DepMap (new quarter), or gnomAD (new version) release — pin every data version in provenance.
- **Staleness signals**: new RNase H1 off-target mismatch benchmarks; updated hepatotoxic/CNS motif models (these are chemistry- and design-context-dependent and need recalibration per gapmer configuration); revised LOEUF/Chronos thresholds.
- **Deprecation**: Archive to `skills/_deprecated/` if superseded by a full multi-agent orchestration with live oracles and a chemistry/backbone module.

## Citations

- [GTEx Consortium, Science 2020](https://www.science.org/doi/10.1126/science.aaz1776); tissue median TPM baseline.
- [Van Nostrand et al., Nature 2020](https://www.nature.com/articles/s41586-020-2077-3); ENCODE eCLIP RBP atlas (peaks cover ~18.5% of the mRNA / 2.6% of the pre-mRNA transcriptome).
- [Kamola et al., NAR 2015](https://academic.oup.com/nar/article/43/18/8638/1073896); experimentally validated exonic + intronic gapmer off-targets at 2–3 mismatches.
- [Burdick et al., NAR 2014](https://academic.oup.com/nar/article/42/8/4882/2408845); TGC/TCC hepatotoxicity trinucleotides.
- [Shen et al., Nat Biotechnol 2019](https://www.nature.com/articles/s41587-019-0106-2); hybridization-independent gapmer protein-binding toxicity (P54nrb/NONO).
- [Dempster et al., Genome Biol 2021](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02540-7); DepMap Chronos gene-effect scoring.
- [Karczewski et al., Nature 2020 (gnomAD)](https://www.nature.com/articles/s41586-020-2308-7); LOEUF loss-of-function constraint.
