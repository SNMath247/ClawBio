---
name: hpo-differential
description: >-
  Rank OMIM diseases, Orphanet diseases and disease-associated genes against a list of
  HPO phenotype terms, fully offline, using information-content-weighted Resnik semantic
  similarity plus a discriminating-feature analysis that does not depend on the score.
license: MIT
metadata:
  version: "0.1.0"
  author: ClawBio contributors
  domain: clinical-genetics
  tags:
    - hpo
    - phenotype
    - rare-disease
    - differential-diagnosis
    - semantic-similarity
    - resnik
    - offline
  trigger_keywords:
    - HPO terms
    - phenotype-driven diagnosis
    - rare disease differential
    - HPO to gene
    - phenotype similarity
    - candidate gene ranking
    - undiagnosed case
    - HP:0000000
  inputs:
    - name: terms
      type: string
      format: [text]
      description: Comma/whitespace/newline separated HPO IDs; trailing labels tolerated
      required: false
    - name: terms_file
      type: file
      format: [txt]
      description: File of HPO IDs, one per line, labels tolerated
      required: false
  outputs:
    - name: report
      type: file
      format: [md]
      description: Human-readable ranked differential with discriminating-feature analysis
    - name: result
      type: file
      format: [json]
      description: Machine-readable rankings, term IC table and provenance
    - name: tables
      type: directory
      format: [csv]
      description: Ranked OMIM, Orphanet and gene tables
  dependencies:
    python: ">=3.10"
    packages:
      - pyhpo>=4,<5
  demo_data:
    - path: skills/hpo-differential/examples/demo_terms.txt
      description: Synthetic 6-term ciliopathy-flavoured phenotype set (never real patient data)
  endpoints:
    cli: python skills/hpo-differential/hpo_differential.py --demo --output /tmp/hpo_diff_demo
  openclaw:
    requires:
      bins: []
    always: false
    emoji: "🧬"
    homepage: https://clawbio.ai
    os: [linux, darwin]
    install:
      - kind: pip
        package: "clawbio[hpo]"
---

# 🧬 HPO Differential

## Trigger

**Fire this skill when the user says any of:**

- "here are the HPO terms for a case, what could it be"
- "rank candidate genes / diseases from this phenotype list"
- "HP:0001250, HP:0001263, ..." (a bare list of HPO IDs)
- "phenotype-driven prioritisation", "which disease explains these features"
- "undiagnosed case", "what should we sequence for"
- "which feature here is most discriminating"

**Do NOT fire when:**

- The user has a VCF and wants variant classification → `clinical-variant-reporter`
- The user has CNV/SV calls → `cnv-acmg-classifier`
- The user wants RNA expression outliers → `rare-disease-rnaseq`
- The user wants free text converted into HPO terms — this skill consumes HPO IDs, it does not extract them
- The user wants a clinical diagnosis. This skill never diagnoses.

## Why This Exists

ClawBio had no phenotype-driven prioritisation anywhere: no HPO term handling, no
Exomiser or LIRICAL integration, no phenotype-similarity scoring. Gene panels existed
only as flat symbol lists. That left the most common starting point for an undiagnosed
case — a clinician's HPO list, with no sequencing data yet — entirely unserved.

It also runs offline on purpose. The public phenotype-ranking services (JAX HPO,
Monarch, PubCaseFinder, PhenoBrain, Orphanet) are frequently unreachable from
governed or air-gapped analysis nodes, and sending phenotype data to a third-party
API is exactly what a PHI boundary is meant to prevent. The HPO release ships inside
`pyhpo`, so the whole computation stays on the machine.

## Core Capabilities

1. Validate every supplied HPO ID against the ontology; abort rather than guess.
2. Rank all OMIM diseases, Orphanet diseases and disease-associated genes by
   phenotypic similarity to the query set.
3. Report, per candidate, which query terms it explains and which it does not.
4. Compute the information content of each query term, so the user can see which
   features are doing the discriminating work.
5. Run a discriminating-feature analysis: list every disease carrying the rarest
   query terms together. This result does not depend on the scoring formula.

## Scope

One skill, one task: **HPO term list in, ranked candidate list out.** It does not
extract HPO terms from text, does not read VCFs, does not classify variants, and
does not diagnose.

## Input Formats

| Format | Extension | Required fields | Example |
|---|---|---|---|
| Inline term list | — | HPO IDs | `--terms "HP:0001263,HP:0000252"` |
| Term file | `.txt` | One HPO ID per line, trailing label optional | `--terms-file case.txt` |
| Demo | — | none | `--demo` |

A line such as `HP:0000252 Microcephaly` is accepted — the label is ignored and the
ontology's own label is used in the report.

## Workflow

1. **Parse and validate the term list.** (Prescriptive — no latitude.) Accept only
   `HP:` + exactly 7 digits. On `HPO:` or a wrong digit count, raise. Deduplicate,
   preserving first-occurrence order.
2. **Resolve every ID against the ontology.** (Prescriptive.) An ID that does not
   resolve aborts the run with exit code 3. Never drop it, never substitute.
3. **Load the annotation corpus** — OMIM, Orphanet and gene term sets from the HPO
   release bundled with pyhpo.
4. **Compute information content** per term from annotation frequency, with ancestor
   closure, separately for each corpus.
5. **Score every entity** with symmetric best-match-average Resnik similarity and
   IC-weighted recall; rank by their product.
6. **Run the discriminating-feature analysis** on the three highest-IC query terms
   individually and on the top pair.
7. **Emit the output contract** (below). (Flexible — the prose framing of the report
   may adapt to the case; the artifact set may not.)

## CLI Reference

```bash
# Rank from an inline term list
python skills/hpo-differential/hpo_differential.py \
  --terms "HP:0001263,HP:0000252,HP:0002170" --output /tmp/hpo_diff

# Rank from a file (labels tolerated on each line)
python skills/hpo-differential/hpo_differential.py \
  --terms-file case_terms.txt --output /tmp/hpo_diff --top 25

# Synthetic demo
python skills/hpo-differential/hpo_differential.py --demo --output /tmp/hpo_diff_demo

# Options
#   --top N               rows per ranking (default 25)
#   --min-annotations N   skip entities with fewer annotations (default 3)
#   --force               overwrite a non-empty output directory
```

Requires the `hpo` extra: `pip install 'clawbio[hpo]'` (or `uv run --extra hpo ...`).

## Demo

```bash
python skills/hpo-differential/hpo_differential.py --demo --output /tmp/hpo_diff_demo
```

Runs on `examples/demo_terms.txt` — a synthetic six-term set. Never real patient data.

## Algorithm / Methodology

```
IC(t)     = -log( |entities annotated with t or a descendant of t| / N )
Resnik(a,b) = max IC over the common ancestors of a and b   (the MICA)
funSimAvg = ½ · ( mean_q max_d Resnik(q,d) + mean_d max_q Resnik(q,d) )
recall    = Σ IC(q) over directly-matched q  +  ½ · Σ IC(q) over ancestor-matched q
            ────────────────────────────────────────────────────────────────────
                                   Σ IC(q) over all q
composite = funSimAvg × recall
```

**Match semantics — note the direction, it is the easiest thing to get wrong:**

- **Direct match**: the entity is annotated with the query term *or a more specific
  descendant of it*. The entity has that feature. Counts fully.
- **Ancestor match**: the entity is annotated only with a *more general ancestor* of
  the query term. It may have that feature. Counts half.

Key thresholds and parameters:

| Parameter | Default | Source |
|---|---|---|
| Minimum annotations per entity | 3 | Avoids stub entries scoring spuriously high on one term |
| Ancestor-match weight | 0.5 | Partial credit convention; a more general annotation is weaker evidence than a specific one |
| Rows per ranking | 25 | Presentation only |
| Ontology release | bundled with `pyhpo` | Reported verbatim in `result.json` as `hpo_release` |

## Example Queries

- "Rank these HPO terms: HP:0002170, HP:0000573, HP:0001263"
- "Which gene best explains this phenotype list?"
- "Which of these features is most discriminating?"
- "What diseases carry both of the rare features in this list?"

## Example Output

```
# HPO Differential — phenotype-driven candidate ranking

- HPO release: **2025-01-16**
- Corpus scored: 8359 OMIM · 4281 Orphanet · 5132 genes

## 1. Query term verification

| HPO ID     | Ontology label             | IC   | # OMIM diseases |
|------------|----------------------------|------|-----------------|
| HP:0000573 | Retinal hemorrhage         | 6.33 | 14              |
| HP:0002170 | Intracranial hemorrhage    | 4.79 | 65              |
| HP:0001263 | Global developmental delay | 1.43 | 1882            |

## Ranked OMIM diseases

| # | ID          | Name                            | Gene(s) | composite | direct/total |
|---|-------------|---------------------------------|---------|-----------|--------------|
| 1 | OMIM:620371 | Neurodevelopmental disorder ... | ESAM    | 1.843     | 12/15        |
| 2 | OMIM:610759 | Cornelia de Lange syndrome 3    | SMC3    | 1.373     | 10/15        |

## Discriminating-feature analysis

**Retinal hemorrhage + Intracranial hemorrhage** — 4 disease(s):
`OMIM:175780`, `OMIM:615368`, `OMIM:620371`, `OMIM:177850`
```

## Output Structure

```
<output_dir>/
├── report.md                          # ranked differential + discriminating analysis
├── result.json                        # rankings, term IC table, corpus sizes, release
├── tables/
│   ├── ranked_diseases.csv            # OMIM
│   ├── ranked_genes.csv               # genes
│   └── ranked_orphanet.csv            # Orphanet
└── reproducibility/
    ├── commands.sh
    └── environment.yml
```

<!--
This tree is the skill's output contract. TestOutputContract in tests/ runs the skill
in --demo mode and asserts every file listed here is produced. If an artifact is only
produced conditionally, add "(optional)" to its comment to exempt it, or remove it
from this tree.
-->

## Dependencies

**Required:** `pyhpo>=4,<5` (bundles the HPO ontology, `phenotype.hpoa`,
`genes_to_phenotype.txt` and `phenotype_to_genes.txt`, ~60 MB).

**Optional:** none. No network access is used at any point.

## Gotchas

1. **The model will want to write an HPO ID from memory when one is missing. Do not.**
   A wrong HPO ID is worse than no HPO ID, because nothing downstream can distinguish
   a typo from a genuinely rare term, and the wrong ID will silently reshape the
   ranking. This skill raises on a malformed ID and exits 3 on an unresolvable one.
   `HPO:` is the specific typo that recurs — the prefix is `HP:`.

2. **The model will want to read the composite score as a probability. It is not.**
   HPO annotation depth is deeply uneven. A disorder described in 2023 from 13
   patients carries an annotation set tightly concentrated on exactly the features
   that defined it, which flatters a best-match-average metric. A syndrome described
   in 1933 carries hundreds of annotations that dilute the same metric. The numeric
   margin between rank 1 and rank 2 routinely overstates the clinical margin. Read
   the discriminating-feature analysis alongside the ranking, not after it.

3. **The model will want to treat an unexplained term as evidence against a candidate.
   Do not.** A feature that a disease's HPO annotation does not mention may simply
   never have been recorded in the handful of published cases. `unexplained_terms` in
   the output means "not explained by the annotation", never "contradicted". Only an
   explicitly assessed-and-absent feature is contradictory evidence, and HPO's
   `negative_hpo` is the only place that is recorded.

4. **The model will want to run this and stop. Do not.** A ranked list is a starting
   point for examination and testing, not an answer. The highest-value output is
   usually the discriminating-feature analysis, because it tells the clinician which
   single observation would most change the ranking.

5. **The model will want to skip `--min-annotations`. Do not lower it to 0.** Entities
   with one or two annotations score spuriously high on best-match-average — a disease
   annotated only with "seizure" gets a perfect symmetric match against a seizure query.

## Safety

ClawBio is a research and educational tool. It is not a medical device and does not
provide clinical diagnoses. Consult a healthcare professional before making any
medical decisions. This exact disclaimer is emitted at the top and bottom of every
generated `report.md` and is asserted by the test suite.

Phenotype data never leaves the machine: there are no network calls in this skill.
The skill refuses to overwrite a non-empty output directory without `--force`.

## Agent Boundary

The agent dispatches and explains. The skill executes. The agent chooses the term
list, interprets the ranking in clinical context, and decides what to do next; the
skill computes the ranking and never editorialises about which candidate is correct.

## Integration with Bio Orchestrator

**Trigger conditions:** a phenotype-only case; a list of HPO IDs; a request for
candidate genes with no sequencing data in hand.

**Chaining partners:**

- → `clinical-variant-reporter`: pass the top-ranked gene symbols as `--genes` once a
  VCF exists.
- → `cnv-acmg-classifier`: the same gene list, for the CNV/SV callset.
- → `rare-disease-rnaseq`: build a custom `--panel` CSV from the ranked genes.
- → `clinical-trial-finder`, `pubmed-summariser`: per-candidate evidence gathering.
- ← upstream: any step that produces an HPO term list.

## Maintenance

**Review cadence:** on each `pyhpo` release (it tracks the HPO monthly releases).

**Staleness signals:** `result.json.hpo_release` more than ~6 months behind the
current HPO release; newly described disorders missing from the corpus; a candidate
the user expects failing to appear at all.

**Deprecation criteria:** retire in favour of a full Exomiser or LIRICAL integration
if one lands in ClawBio, since those add variant-level evidence this skill cannot.
Until then this is the only phenotype-driven prioritisation in the repo.

## Citations

- Köhler S, *et al.* The Human Phenotype Ontology in 2021. *Nucleic Acids Res.*
  2021;49(D1):D1207–D1217. doi:10.1093/nar/gkaa1043
- Resnik P. Using information content to evaluate semantic similarity in a taxonomy.
  *IJCAI* 1995;448–453.
- Schlicker A, *et al.* A new measure for functional similarity of gene products based
  on Gene Ontology. *BMC Bioinformatics* 2006;7:302. doi:10.1186/1471-2105-7-302
  (funSimAvg / best-match-average combination)
- pyhpo — https://pypi.org/project/pyhpo/
