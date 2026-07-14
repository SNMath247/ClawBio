#!/usr/bin/env python3
"""ASO Off-Target & Toxicity Screening — deterministic multi-stage funnel.

Implements the computable, local core of the antisense-oligonucleotide (ASO)
off-target screening blueprint:

  1. In-silico mutant transcript  = baseline pre-mRNA + patient variant
  2. Candidate tiling             = 15-20mer ASOs across the target window
  3. Transcriptome alignment      = near-match search over a pre-mRNA atlas
                                    (INTRONS INCLUDED) with a modality-aware
                                    mismatch ceiling
  4. Tissue-weighted penalty      = sum_i P(bind_i) * log2(TPM_tissue,i + 1)
                                    * E(essentiality_i)
  5. CLIP collision oracle        = genomic-interval intersect vs an RBP peak
                                    atlas, separating INTENDED mechanism (asset)
                                    from UNINTENDED collision (liability)
  6. Sequence-intrinsic tox       = CpG/TLR9, G-quadruplex, hepatotoxic and
                                    GU-rich motif screen on the ASO itself
  7. Multi-objective arbitration  = a Pareto front over the axes; NO single
                                    scalar "safety score" is emitted

This is a deterministic, stdlib-only triage helper. It ships with small
synthetic reference tracks; it does not clear a candidate on safety. Production
oracles (GTEx v10 API, ENCODE/POSTAR3 eCLIP, DepMap Chronos, gnomAD LOEUF) are
documented in SKILL.md and plug in via the same table schemas.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path


SKILL_DIR = Path(__file__).resolve().parent
DATA_DIR = SKILL_DIR / "data"
DISCLAIMER = (
    "ClawBio is a research and educational tool. It is not a medical device "
    "and does not provide clinical diagnoses. Consult a healthcare professional "
    "before making any medical decisions."
)

COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G"}
VALID_BASES = set(COMPLEMENT)
# Modality-specific mismatch tolerance for the off-target search.
#   gapmers cleave via RNase H1 and tolerate 1-3 mismatches (Kamola 2015);
#   steric blockers need near-complete complementarity to act.
DEFAULT_MAX_MISMATCH = {"gapmer": 3, "steric-blocker": 1}
REQUIRED_SPEC_FIELDS = {"tissue", "pre_mrna", "target_window"}

# Multi-objective arbitration axes. Benefits are maximised (negated into costs),
# liabilities are minimised. Intended mechanistic-silencer engagement is treated
# as a benefit (asset), per the nusinersen/ISS-N1 precedent.
BENEFIT_AXES = ("efficacy", "intended_mechanism")
LIABILITY_AXES = ("offtarget_penalty", "clip_collision", "essentiality_hit", "seq_tox")


# --------------------------------------------------------------------------- #
# Sequence primitives
# --------------------------------------------------------------------------- #

def reverse_complement(seq: str) -> str:
    seq = seq.upper()
    return "".join(COMPLEMENT[b] for b in reversed(seq))


def gc_fraction(seq: str) -> float:
    if not seq:
        return 0.0
    seq = seq.upper()
    return sum(1 for b in seq if b in ("G", "C")) / len(seq)


def _validate_dna(seq: str, label: str) -> str:
    seq = seq.upper().strip()
    bad = set(seq) - VALID_BASES
    if bad:
        raise ValueError(f"{label} contains non-ACGT characters: {sorted(bad)}")
    if not seq:
        raise ValueError(f"{label} is empty")
    return seq


def build_mutant_transcript(pre_mrna: str, variant: dict | None) -> str:
    """Merge the baseline pre-mRNA with a patient SNV to get the mutant transcript."""
    pre_mrna = _validate_dna(pre_mrna, "pre_mrna")
    if not variant:
        return pre_mrna
    pos = int(variant["pos"])
    ref = str(variant["ref"]).upper()
    alt = str(variant["alt"]).upper()
    if pos < 0 or pos >= len(pre_mrna):
        raise ValueError(f"variant pos {pos} out of range for pre_mrna length {len(pre_mrna)}")
    if len(ref) != 1 or len(alt) != 1 or ref not in VALID_BASES or alt not in VALID_BASES:
        raise ValueError("only single-nucleotide variants (SNVs) are supported")
    if pre_mrna[pos] != ref:
        raise ValueError(f"variant ref '{ref}' does not match pre_mrna base '{pre_mrna[pos]}' at pos {pos}")
    return pre_mrna[:pos] + alt + pre_mrna[pos + 1:]


# --------------------------------------------------------------------------- #
# Stage 2 — candidate tiling
# --------------------------------------------------------------------------- #

def tile_candidates(mutant: str, window: dict, length_min: int, length_max: int,
                    genomic_start: int) -> list[dict]:
    """Tile ASO candidates across the target window on the mutant pre-mRNA.

    Each ASO is antisense to the sense target site (its reverse complement).
    Its genomic footprint is colinear with the pre-mRNA (introns included).
    """
    start = int(window["start"])
    end = min(int(window["end"]), len(mutant))
    if start < 0 or start >= end:
        raise ValueError("target_window is empty or out of range")
    length_min = max(1, int(length_min))
    length_max = max(length_min, int(length_max))
    candidates: list[dict] = []
    for length in range(length_min, length_max + 1):
        for offset in range(start, end - length + 1):
            sense = mutant[offset: offset + length]
            aso = reverse_complement(sense)
            candidates.append({
                "aso_id": f"ASO_{offset}_{length}",
                "start": offset,
                "length": length,
                "target_sense": sense,
                "aso": aso,
                "gc": round(gc_fraction(aso), 3),
                "genomic_start": genomic_start + offset,
                "genomic_end": genomic_start + offset + length,
            })
    if not candidates:
        raise ValueError("no ASO candidates could be tiled from the target window")
    return candidates


# --------------------------------------------------------------------------- #
# Stage 3 — transcriptome / pre-mRNA off-target alignment
# --------------------------------------------------------------------------- #

def p_bind(length: int, mismatches: int, gc: float) -> float:
    """Heuristic binding-likelihood score in [0, 1] (NOT a thermodynamic ΔG).

    Monotonically DECREASING in mismatch count and INCREASING in GC content.
    identity**3 gives a steep-but-nonzero falloff so 1-3 mismatches still bind
    at reduced likelihood. This is a coarse ranking proxy only: it is blind to
    mismatch position/identity and has no nearest-neighbour term, so it is not
    a ΔΔG°37 calculation. A real thermodynamic model (OligoWalk/RIsearch2) is
    the production swap-in.
    """
    if length <= 0:
        return 0.0
    identity = max(0.0, (length - mismatches) / length)
    return round((identity ** 3) * (0.5 + 0.5 * max(0.0, min(1.0, gc))), 6)


def align_offtargets(aso_target: str, atlas: dict, max_mismatch: int) -> list[dict]:
    """Find near-matches of ``aso_target`` (= reverse complement of the ASO)
    across every atlas transcript, scanning all offsets (introns included)."""
    aso_target = aso_target.upper()
    L = len(aso_target)
    gc = gc_fraction(aso_target)
    hits: list[dict] = []
    for tx in atlas.get("transcripts", []):
        seq = str(tx["sequence"]).upper()
        gstart = int(tx["genomic_start"])
        for offset in range(0, len(seq) - L + 1):
            window = seq[offset: offset + L]
            mism = sum(1 for a, b in zip(aso_target, window) if a != b)
            if mism <= max_mismatch:
                hits.append({
                    "transcript_id": tx["transcript_id"],
                    "gene": tx["gene"],
                    "chrom": tx.get("chrom", ""),
                    "strand": tx.get("strand", "+"),
                    "region": tx.get("region", "unknown"),
                    "offset": offset,
                    "mismatches": mism,
                    "genomic_start": gstart + offset,
                    "genomic_end": gstart + offset + L,
                    "p_bind": p_bind(L, mism, gc),
                })
    hits.sort(key=lambda h: (h["mismatches"], h["gene"], h["offset"]))
    return hits


# --------------------------------------------------------------------------- #
# Stage 4 — tissue-weighted off-target penalty
# --------------------------------------------------------------------------- #

def essentiality_multiplier(chronos: float, loeuf: float) -> float:
    """E(essentiality): highly essential / constrained genes multiply the penalty.

    Inverse DepMap Chronos (more negative -> more essential) plus a gnomAD LOEUF
    constraint boost (LOEUF < 0.6 is the loss-of-function-intolerant regime).
    A neutral gene (chronos 0, loeuf >= 0.6) returns exactly 1.0.
    """
    chronos = float(chronos)
    loeuf = float(loeuf)
    essential_term = max(0.0, -chronos)
    loeuf_boost = max(0.0, (0.6 - loeuf)) / 0.6
    return 1.0 + essential_term + loeuf_boost


def penalty_contribution(hit: dict, tpm: dict, tissue: str, essentiality: dict) -> float:
    """P(bind_i) * log2(TPM_tissue,i + 1) * E(essentiality_i)."""
    gene = hit["gene"]
    tissue_tpm = float(tpm.get(gene, {}).get(tissue, 0.0))
    chronos, loeuf = essentiality.get(gene, (0.0, 1.0))
    e = essentiality_multiplier(chronos, loeuf)
    return float(hit["p_bind"]) * math.log2(tissue_tpm + 1.0) * e


# --------------------------------------------------------------------------- #
# Stage 5 — CLIP collision oracle
# --------------------------------------------------------------------------- #

def intervals_overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> bool:
    """Half-open [start, end) overlap test (touching endpoints do not overlap)."""
    return a_start < b_end and b_start < a_end


def clip_hits(chrom: str, start: int, end: int, strand: str, peaks: list[dict]) -> list[dict]:
    """Return CLIP peaks overlapping a genomic footprint (strand-matched)."""
    out = []
    for peak in peaks:
        if peak["chrom"] != chrom:
            continue
        if peak.get("strand", strand) != strand:
            continue
        if intervals_overlap(start, end, int(peak["start"]), int(peak["end"])):
            out.append(peak)
    return out


def _clip_weight(peak: dict) -> float:
    # IDR-reproducible peaks weigh above relaxed peaks.
    return 1.0 if peak.get("idr") else 0.5


# --------------------------------------------------------------------------- #
# Stage 6 — sequence-intrinsic toxicity / immunogenicity motifs
# --------------------------------------------------------------------------- #

def _count_occurrences(seq: str, pattern: str) -> int:
    count = start = 0
    while True:
        idx = seq.find(pattern, start)
        if idx == -1:
            return count
        count += 1
        start = idx + 1   # overlapping occurrences count (e.g. GGGGG -> two GGGG)


def scan_tox_motifs(aso: str, rules: list[dict]) -> dict:
    aso = aso.upper()
    flags: list[dict] = []
    penalty = 0.0
    severity = "none"
    chemistry: list[str] = []
    for rule in rules:
        total = sum(_count_occurrences(aso, p.upper()) for p in rule["patterns"])
        if total:
            weight = float(rule.get("weight", 1.0))
            penalty += total * weight
            flags.append({
                "id": rule["id"],
                "count": total,
                "severity": rule.get("severity", "flag"),
                "mechanism": rule.get("mechanism", ""),
            })
            if rule.get("mitigation"):
                chemistry.append(f"{rule['id']}: {rule['mitigation']}")
            if rule.get("severity") == "exclude":
                severity = "exclude"
            elif severity != "exclude":
                severity = "flag"
    return {
        "flags": flags,
        "tox_penalty": round(penalty, 3),
        "severity": severity,
        "chemistry_flags": chemistry,
    }


# --------------------------------------------------------------------------- #
# Stage 7 — multi-objective Pareto arbitration
# --------------------------------------------------------------------------- #

def _cost_vector(cand: dict) -> tuple:
    # All axes expressed as costs to MINIMISE (benefits negated).
    return (
        -float(cand["efficacy"]),
        -1.0 if cand.get("intended_mechanism") else 0.0,
        float(cand["offtarget_penalty"]),
        float(cand["clip_collision"]),
        float(cand["essentiality_hit"]),
        float(cand["seq_tox"]),
    )


def _dominates(a: dict, b: dict) -> bool:
    """True if a dominates b: no worse on every axis and strictly better on one."""
    av, bv = _cost_vector(a), _cost_vector(b)
    if any(x > y for x, y in zip(av, bv)):
        return False
    return any(x < y for x, y in zip(av, bv))


def pareto_front(candidates: list[dict]) -> list[dict]:
    front = []
    for cand in candidates:
        if not any(other is not cand and _dominates(other, cand) for other in candidates):
            front.append(cand)
    return front


# --------------------------------------------------------------------------- #
# On-target efficacy (transparent proxy)
# --------------------------------------------------------------------------- #

def on_target_efficacy(aso: str, gc: float, tox_severity: str) -> float:
    """Balanced-GC, low-liability ASOs score higher. A convenience proxy only."""
    gc_score = max(0.0, 1.0 - 2.0 * abs(gc - 0.55))     # peaks at ~55% GC
    penalty = 0.3 if tox_severity == "exclude" else (0.1 if tox_severity == "flag" else 0.0)
    return round(max(0.0, min(1.0, gc_score - penalty)), 3)


# --------------------------------------------------------------------------- #
# Reference loaders
# --------------------------------------------------------------------------- #

def load_target_spec(path: Path) -> dict:
    spec = json.loads(Path(path).read_text(encoding="utf-8"))
    missing = REQUIRED_SPEC_FIELDS - set(spec)
    if missing:
        raise ValueError(f"target spec is missing required fields: {sorted(missing)}")
    return spec


def load_atlas(path: Path) -> dict:
    atlas = json.loads(Path(path).read_text(encoding="utf-8"))
    if "transcripts" not in atlas:
        raise ValueError("off-target atlas must contain a 'transcripts' list")
    return atlas


def load_tpm(path: Path) -> dict:
    with Path(path).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        tissues = [c for c in (reader.fieldnames or []) if c != "gene"]
        out: dict = {}
        for row in reader:
            out[row["gene"]] = {t: float(row[t]) for t in tissues}
    return out


def load_essentiality(path: Path) -> dict:
    with Path(path).open(newline="", encoding="utf-8") as handle:
        out = {}
        for row in csv.DictReader(handle):
            out[row["gene"]] = (float(row["chronos"]), float(row["loeuf"]))
    return out


def load_clip(path: Path) -> list[dict]:
    with Path(path).open(newline="", encoding="utf-8") as handle:
        peaks = []
        for row in csv.DictReader(handle):
            peaks.append({
                "chrom": row["chrom"],
                "start": int(row["start"]),
                "end": int(row["end"]),
                "strand": row.get("strand", "+"),
                "rbp": row.get("rbp", ""),
                "cell_line": row.get("cell_line", ""),
                "idr": str(row.get("idr", "")).lower() == "true",
                "mechanistic_silencer": str(row.get("mechanistic_silencer", "")).lower() == "true",
            })
    return peaks


def load_motifs(path: Path) -> list[dict]:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    if isinstance(data, list):
        return data
    return data.get("rules", [])


# --------------------------------------------------------------------------- #
# Orchestration
# --------------------------------------------------------------------------- #

def _verdict(cand: dict) -> str:
    """Triage tier (NOT a safety clearance)."""
    if cand["tox_severity"] == "exclude":
        return "exclude"
    if cand["clip_collision_essential"]:
        return "exclude"
    if cand["tox_flags"] or cand["clip_collision"] > 0 or cand["essentiality_hit"] > 0:
        return "chemistry-redesign"
    if cand["offtarget_penalty"] > 0:
        return "chemistry-redesign"
    return "prioritise-for-wet-lab"


def screen(spec: dict, atlas: dict, tpm: dict, essentiality: dict,
           clip_peaks: list[dict], motif_rules: list[dict],
           tissue_override: str | None = None) -> dict:
    tissue = tissue_override or spec["tissue"]
    modality = spec.get("modality", "gapmer")
    max_mismatch = int(spec.get("max_mismatch", DEFAULT_MAX_MISMATCH.get(modality, 3)))
    lmin = int(spec.get("aso_length_min", 15))
    lmax = int(spec.get("aso_length_max", 20))
    chrom = spec.get("chrom", "chr")
    strand = spec.get("strand", "+")
    genomic_start = int(spec.get("pre_mrna_genomic_start", 0))

    mutant = build_mutant_transcript(spec["pre_mrna"], spec.get("variant"))
    candidates = tile_candidates(mutant, spec["target_window"], lmin, lmax, genomic_start)

    scored: list[dict] = []
    for cand in candidates:
        # Stage 3: off-target hits (aso_target = revcomp(aso) = the sense the ASO recognises).
        hits = align_offtargets(cand["target_sense"], atlas, max_mismatch)
        # Stage 4: tissue-weighted penalty per hit.
        penalty = 0.0
        essentiality_hit = 0.0
        for h in hits:
            contrib = penalty_contribution(h, tpm, tissue, essentiality)
            h["penalty_contribution"] = round(contrib, 4)
            h["tissue_tpm"] = round(float(tpm.get(h["gene"], {}).get(tissue, 0.0)), 3)
            h["common_essential"] = h["gene"] in essentiality and essentiality[h["gene"]][0] <= -0.5
            penalty += contrib
            if contrib > 0:
                essentiality_hit = max(essentiality_hit, essentiality_multiplier(*essentiality.get(h["gene"], (0.0, 1.0))))
        hits.sort(key=lambda x: -x["penalty_contribution"])

        # Stage 5: CLIP collision — intended mechanism vs unintended collision.
        intended = clip_hits(chrom, cand["genomic_start"], cand["genomic_end"], strand, clip_peaks)
        intended_mech = [p for p in intended if p["mechanistic_silencer"]]
        collision = 0.0
        collision_essential = False
        collision_peaks = []
        for h in hits:
            for peak in clip_hits(h["chrom"], h["genomic_start"], h["genomic_end"], h["strand"], clip_peaks):
                collision += _clip_weight(peak)
                collision_peaks.append({"gene": h["gene"], "region": h["region"], "rbp": peak["rbp"],
                                        "cell_line": peak["cell_line"], "idr": peak["idr"]})
                if h.get("common_essential"):
                    collision_essential = True

        # Stage 6: sequence-intrinsic toxicity / immunogenicity.
        tox = scan_tox_motifs(cand["aso"], motif_rules)

        efficacy = on_target_efficacy(cand["aso"], cand["gc"], tox["severity"])
        record = {
            **cand,
            "efficacy": efficacy,
            "offtarget_penalty": round(penalty, 4),
            "offtarget_hits": hits,
            "clip_collision": round(collision, 3),
            "clip_collision_peaks": collision_peaks,
            "clip_collision_essential": collision_essential,
            "intended_mechanism": bool(intended_mech),
            "intended_mechanism_rbps": sorted({p["rbp"] for p in intended_mech}),
            "essentiality_hit": round(essentiality_hit, 3),
            "seq_tox": tox["tox_penalty"],
            "tox_severity": tox["severity"],
            "tox_flags": tox["flags"],
            "chemistry_flags": tox["chemistry_flags"],
        }
        record["verdict"] = _verdict(record)
        scored.append(record)

    front = pareto_front(scored)
    front_ids = [c["aso_id"] for c in front]
    # Convenience ordering WITHIN the front (transparent; not a clearance score).
    front_sorted = sorted(
        front,
        key=lambda c: (c["verdict"] != "prioritise-for-wet-lab",
                       c["offtarget_penalty"] + c["clip_collision"] + c["seq_tox"],
                       -c["efficacy"]),
    )

    n_excluded = sum(1 for c in scored if c["verdict"] == "exclude")
    return {
        "skill": "aso-off-target-screening",
        "target": {
            "target_id": spec.get("target_id", ""),
            "gene": spec.get("gene", ""),
            "tissue": tissue,
            "modality": modality,
            "variant": spec.get("variant"),
            "target_window": spec["target_window"],
            "mutant_length": len(mutant),
            "n_candidates": len(candidates),
        },
        "parameters": {
            "max_mismatch": max_mismatch,
            "aso_length_range": [lmin, lmax],
            "penalty_formula": "sum_i P(bind_i) * log2(TPM_tissue,i + 1) * E(essentiality_i)",
            "pareto_axes": {"benefits": list(BENEFIT_AXES), "liabilities": list(LIABILITY_AXES)},
        },
        "summary": {
            "n_candidates": len(candidates),
            "n_pareto": len(front_ids),
            "n_excluded": n_excluded,
            "n_intended_mechanism": sum(1 for c in scored if c["intended_mechanism"]),
            "n_with_offtarget": sum(1 for c in scored if c["offtarget_penalty"] > 0),
            "n_clip_collision": sum(1 for c in scored if c["clip_collision"] > 0),
            "best_candidate": front_sorted[0]["aso_id"] if front_sorted else None,
        },
        "candidates": scored,
        "pareto_front": front_ids,
        "pareto_front_ranked": [c["aso_id"] for c in front_sorted],
        "disclaimer": DISCLAIMER,
    }


def _default_refs() -> dict:
    return {
        "atlas": load_atlas(DATA_DIR / "offtarget_atlas.json"),
        "tpm": load_tpm(DATA_DIR / "gtex_median_tpm.csv"),
        "essentiality": load_essentiality(DATA_DIR / "essentiality.csv"),
        "clip_peaks": load_clip(DATA_DIR / "clip_peaks.csv"),
        "motif_rules": load_motifs(DATA_DIR / "tox_motifs.json"),
    }


def screen_from_spec(spec_path: Path, tissue_override: str | None = None,
                     modality_override: str | None = None,
                     max_mismatch: int | None = None) -> dict:
    spec = load_target_spec(Path(spec_path))
    if modality_override:
        spec["modality"] = modality_override
    if max_mismatch is not None:
        spec["max_mismatch"] = max_mismatch
    refs = _default_refs()
    return screen(spec, refs["atlas"], refs["tpm"], refs["essentiality"],
                  refs["clip_peaks"], refs["motif_rules"], tissue_override=tissue_override)


# --------------------------------------------------------------------------- #
# Output pack
# --------------------------------------------------------------------------- #

def _write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({f: row.get(f, "") for f in fieldnames})


def write_outputs(result: dict, input_path: Path, output_dir: Path, command: list[str], demo: bool) -> None:
    if output_dir.exists() and any(output_dir.iterdir()):
        print(f"WARNING: output directory already exists and files may be overwritten: {output_dir}", file=sys.stderr)
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "tables").mkdir(exist_ok=True)
    (output_dir / "reproducibility").mkdir(exist_ok=True)

    cand_rows = [{
        "aso_id": c["aso_id"], "start": c["start"], "length": c["length"], "aso": c["aso"],
        "gc": c["gc"], "efficacy": c["efficacy"], "offtarget_penalty": c["offtarget_penalty"],
        "clip_collision": c["clip_collision"], "essentiality_hit": c["essentiality_hit"],
        "seq_tox": c["seq_tox"], "intended_mechanism": c["intended_mechanism"],
        "tox_flags": ";".join(f["id"] for f in c["tox_flags"]),
        "in_pareto_front": c["aso_id"] in set(result["pareto_front"]), "verdict": c["verdict"],
    } for c in result["candidates"]]
    _write_csv(output_dir / "tables" / "candidates.csv", cand_rows, list(cand_rows[0]))

    hit_rows = []
    for c in result["candidates"]:
        for h in c["offtarget_hits"]:
            hit_rows.append({
                "aso_id": c["aso_id"], "gene": h["gene"], "transcript_id": h["transcript_id"],
                "region": h["region"], "mismatches": h["mismatches"], "p_bind": h["p_bind"],
                "tissue_tpm": h["tissue_tpm"], "penalty_contribution": h["penalty_contribution"],
                "common_essential": h["common_essential"],
                "genomic": f"{h['chrom']}:{h['genomic_start']}-{h['genomic_end']}",
            })
    _write_csv(output_dir / "tables" / "offtarget_hits.csv", hit_rows,
               ["aso_id", "gene", "transcript_id", "region", "mismatches", "p_bind",
                "tissue_tpm", "penalty_contribution", "common_essential", "genomic"])

    front = {c["aso_id"]: c for c in result["candidates"]}
    front_rows = [{
        "rank": i + 1, "aso_id": aid, "aso": front[aid]["aso"], "efficacy": front[aid]["efficacy"],
        "offtarget_penalty": front[aid]["offtarget_penalty"], "clip_collision": front[aid]["clip_collision"],
        "essentiality_hit": front[aid]["essentiality_hit"], "seq_tox": front[aid]["seq_tox"],
        "intended_mechanism": front[aid]["intended_mechanism"], "verdict": front[aid]["verdict"],
    } for i, aid in enumerate(result["pareto_front_ranked"])]
    _write_csv(output_dir / "tables" / "pareto_front.csv", front_rows,
               ["rank", "aso_id", "aso", "efficacy", "offtarget_penalty", "clip_collision",
                "essentiality_hit", "seq_tox", "intended_mechanism", "verdict"])

    (output_dir / "result.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    (output_dir / "report.md").write_text(_render_report(result, input_path, demo), encoding="utf-8")
    (output_dir / "reproducibility" / "commands.sh").write_text(
        "#!/usr/bin/env bash\n" + " ".join(command) + "\n", encoding="utf-8")


def _render_report(result: dict, input_path: Path, demo: bool) -> str:
    t = result["target"]
    s = result["summary"]
    lines = [
        "# ASO Off-Target & Toxicity Screening Report",
        "",
        f"**Input**: `{input_path}`",
        f"**Mode**: {'Synthetic demo data' if demo else 'User-provided local data'}",
        f"**Target**: {t['target_id']} ({t['gene']}) — tissue **{t['tissue']}**, modality **{t['modality']}**",
        f"**Mutant transcript**: {t['mutant_length']} nt"
        + (f", variant {t['variant']['ref']}>{t['variant']['alt']}@{t['variant']['pos']}" if t.get("variant") else ""),
        f"**Candidates tiled**: {s['n_candidates']} | "
        f"**Pareto-optimal**: {s['n_pareto']} | **Excluded**: {s['n_excluded']}",
        "",
        "> **This is triage, not clearance.** A large fraction of PS-ASO toxicity is "
        "hybridization-INDEPENDENT (backbone protein binding, complement, thrombocytopenia) "
        "and is NOT addressable by sequence screening. No single scalar clears a candidate on "
        "safety; the Pareto front below preserves the trade-offs for human review and wet-lab.",
        "",
        f"**Penalty formula**: `{result['parameters']['penalty_formula']}`",
        "",
        "## Pareto-optimal candidates (non-dominated)",
        "",
        "| Rank | ASO ID | ASO (5'->3') | Efficacy | Off-tgt penalty | CLIP collision | Ess. hit | Seq tox | Intended mech. | Verdict |",
        "|---:|---|---|---:|---:|---:|---:|---:|:---:|---|",
    ]
    ranked = {c["aso_id"]: c for c in result["candidates"]}
    display_front = result["pareto_front_ranked"][:15]
    for i, aid in enumerate(display_front, 1):
        c = ranked[aid]
        lines.append(
            f"| {i} | {aid} | `{c['aso']}` | {c['efficacy']} | {c['offtarget_penalty']} | "
            f"{c['clip_collision']} | {c['essentiality_hit']} | {c['seq_tox']} | "
            f"{'yes' if c['intended_mechanism'] else '-'} | {c['verdict']} |"
        )
    if len(result["pareto_front_ranked"]) > len(display_front):
        lines.append(f"| ... | _+{len(result['pareto_front_ranked']) - len(display_front)} "
                     f"more non-dominated_ | | | | | | | | |")

    # Highlight the worst weighted off-target offender.
    worst = max(result["candidates"], key=lambda c: c["offtarget_penalty"])
    lines += ["", "## Highest weighted off-target liability", ""]
    if worst["offtarget_penalty"] > 0 and worst["offtarget_hits"]:
        h = worst["offtarget_hits"][0]
        lines.append(
            f"`{worst['aso_id']}` ({worst['aso']}) — penalty **{worst['offtarget_penalty']}**, "
            f"driven by **{h['gene']}** ({h['region']}, {h['mismatches']} mismatch, "
            f"TPM {h['tissue_tpm']}, contribution {h['penalty_contribution']})."
        )
    else:
        lines.append("No candidate carried a non-zero weighted off-target penalty.")

    excluded = [c for c in result["candidates"] if c["verdict"] == "exclude"]
    if excluded:
        lines += ["", "## Excluded on sequence-intrinsic / essential-collision grounds", ""]
        for c in excluded[:10]:
            reasons = ";".join(f["id"] for f in c["tox_flags"]) or ("essential-gene CLIP collision" if c["clip_collision_essential"] else "")
            lines.append(f"- `{c['aso_id']}` ({c['aso']}): {reasons}")

    lines += [
        "",
        "## Interpretation",
        "",
        "1. **Tissue weighting** zeroes off-targets that are untranscribed in the target tissue "
        "(`log2(0+1)=0`), while essential, constrained genes multiply the penalty.",
        "2. **Pre-mRNA search** includes introns, so intronic collisions with RBP CLIP peaks are "
        "caught — a common blind spot for spliced-transcript-only pipelines.",
        "3. **Intended-vs-unintended CLIP**: an on-target overlap with a mechanistic silencer is an "
        "asset (e.g. nusinersen blocks hnRNP A1 at ISS-N1), not a liability.",
        "4. **Next step**: escalate Pareto-optimal candidates carrying any essentiality hit or "
        "hepatotoxic/CNS motif to wet-lab tolerability cascades. In-silico triage cannot clear safety.",
        "",
        DISCLAIMER,
        "",
    ]
    return "\n".join(lines)


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #

def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="ASO Off-Target & Toxicity Screening")
    parser.add_argument("--input", type=Path, help="Target spec JSON")
    parser.add_argument("--output", type=Path, default=Path("aso_off_target_out"))
    parser.add_argument("--demo", action="store_true", help="Run with bundled synthetic data")
    parser.add_argument("--tissue", type=str, default=None, help="Override the target tissue")
    parser.add_argument("--modality", type=str, default=None, choices=["gapmer", "steric-blocker"])
    parser.add_argument("--max-mismatch", type=int, default=None, help="Off-target mismatch ceiling")
    args = parser.parse_args(argv)

    input_path = SKILL_DIR / "demo_target.json" if args.demo else args.input
    if input_path is None:
        parser.error("--input is required unless --demo is used")
    try:
        result = screen_from_spec(input_path, tissue_override=args.tissue,
                                  modality_override=args.modality, max_mismatch=args.max_mismatch)
    except (ValueError, KeyError, json.JSONDecodeError, OSError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    write_outputs(result, input_path, args.output, [sys.executable, __file__, *sys.argv[1:]], args.demo)
    print(f"aso-off-target-screening wrote {args.output / 'report.md'} "
          f"({result['summary']['n_candidates']} candidates, "
          f"{result['summary']['n_pareto']} Pareto-optimal)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
