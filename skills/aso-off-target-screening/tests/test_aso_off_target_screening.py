import importlib.util
import json
import math
import subprocess
import sys
from pathlib import Path


SKILL_DIR = Path(__file__).resolve().parents[1]
MODULE_PATH = SKILL_DIR / "aso_off_target_screening.py"
DEMO_SPEC = SKILL_DIR / "demo_target.json"
DISCLAIMER_FRAGMENT = "not a medical device"


def load_module():
    spec = importlib.util.spec_from_file_location("aso_off_target_screening", MODULE_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# --------------------------------------------------------------------------- #
# Pure-function unit tests
# --------------------------------------------------------------------------- #

def test_reverse_complement():
    m = load_module()
    assert m.reverse_complement("ATGC") == "GCAT"
    assert m.reverse_complement("AAAA") == "TTTT"
    assert m.reverse_complement(m.reverse_complement("ACGTACGT")) == "ACGTACGT"


def test_gc_fraction():
    m = load_module()
    assert m.gc_fraction("GGCC") == 1.0
    assert m.gc_fraction("ATAT") == 0.0
    assert abs(m.gc_fraction("ATGC") - 0.5) < 1e-9


def test_build_mutant_transcript_applies_snv():
    m = load_module()
    pre = "AAAAGAAAA"
    mut = m.build_mutant_transcript(pre, {"pos": 4, "ref": "G", "alt": "A"})
    assert mut == "AAAAAAAAA"
    assert len(mut) == len(pre)


def test_build_mutant_transcript_rejects_ref_mismatch():
    m = load_module()
    try:
        m.build_mutant_transcript("AAAAGAAAA", {"pos": 4, "ref": "C", "alt": "A"})
    except ValueError as exc:
        assert "ref" in str(exc).lower()
    else:
        raise AssertionError("ref mismatch should raise")


def test_build_mutant_transcript_none_is_identity():
    m = load_module()
    assert m.build_mutant_transcript("ACGTACGT", None) == "ACGTACGT"


def test_p_bind_monotonic_decreasing_in_mismatches():
    m = load_module()
    p0 = m.p_bind(18, 0, 0.5)
    p1 = m.p_bind(18, 1, 0.5)
    p3 = m.p_bind(18, 3, 0.5)
    assert 0.0 <= p3 < p1 < p0 <= 1.0


def test_p_bind_increases_with_gc():
    m = load_module()
    assert m.p_bind(18, 1, 0.7) > m.p_bind(18, 1, 0.3)


def test_essentiality_multiplier_ranks_essential_high():
    m = load_module()
    ess = m.essentiality_multiplier(-1.62, 0.10)   # strongly essential + constrained
    mid = m.essentiality_multiplier(-0.15, 0.75)
    neutral = m.essentiality_multiplier(0.0, 1.0)
    assert ess > mid > neutral
    assert abs(neutral - 1.0) < 1e-9


def test_tiling_candidate_count():
    m = load_module()
    mutant = "N" * 120
    cands = m.tile_candidates(mutant.replace("N", "A"), {"start": 40, "end": 100}, 18, 18, 1000)
    # window is 60 nt wide, 18-mers -> 43 candidates
    assert len(cands) == 43
    first = cands[0]
    assert first["length"] == 18
    assert first["genomic_end"] - first["genomic_start"] == 18


def test_align_offtargets_finds_planted_site_with_exact_mismatches():
    m = load_module()
    query = "ACGTACGTACGTACGTAC"           # 18 nt aso_target
    planted = list(query)
    planted[5] = "A" if planted[5] != "A" else "T"   # 1 mismatch
    transcript = "TTTTT" + "".join(planted) + "GGGGG"
    atlas = {"transcripts": [{
        "transcript_id": "T1", "gene": "G1", "chrom": "c", "strand": "+",
        "genomic_start": 100, "region": "exon", "sequence": transcript,
    }]}
    hits = m.align_offtargets(query, atlas, max_mismatch=3)
    assert any(h["gene"] == "G1" and h["mismatches"] == 1 and h["offset"] == 5 for h in hits)


def test_align_offtargets_respects_mismatch_ceiling():
    m = load_module()
    query = "ACGTACGTACGTACGTAC"
    transcript = "TTTTTGGGGGCCCCCAAAAA" * 2     # no near-match
    atlas = {"transcripts": [{
        "transcript_id": "T1", "gene": "G1", "chrom": "c", "strand": "+",
        "genomic_start": 100, "region": "exon", "sequence": transcript,
    }]}
    assert m.align_offtargets(query, atlas, max_mismatch=2) == []


def test_weighted_penalty_zeroes_when_untranscribed():
    m = load_module()
    # Two hits, identical p_bind and essentiality, differ only in tissue TPM.
    hit_expr = {"gene": "EXPR", "p_bind": 0.5, "mismatches": 1}
    hit_silent = {"gene": "SILENT", "p_bind": 0.5, "mismatches": 1}
    tpm = {"EXPR": {"Muscle_Skeletal": 50.0}, "SILENT": {"Muscle_Skeletal": 0.0}}
    ess = {"EXPR": (0.0, 1.0), "SILENT": (0.0, 1.0)}
    c_expr = m.penalty_contribution(hit_expr, tpm, "Muscle_Skeletal", ess)
    c_silent = m.penalty_contribution(hit_silent, tpm, "Muscle_Skeletal", ess)
    assert c_silent == 0.0
    assert c_expr > c_silent


def test_weighted_penalty_essentiality_multiplier_increases_penalty():
    m = load_module()
    hit = {"gene": "G", "p_bind": 0.5, "mismatches": 1}
    tpm = {"G": {"Muscle_Skeletal": 50.0}}
    ess_high = {"G": (-1.6, 0.1)}
    ess_low = {"G": (0.0, 1.0)}
    assert (m.penalty_contribution(hit, tpm, "Muscle_Skeletal", ess_high)
            > m.penalty_contribution(hit, tpm, "Muscle_Skeletal", ess_low))


def test_penalty_matches_formula():
    m = load_module()
    hit = {"gene": "G", "p_bind": 0.4, "mismatches": 1}
    tpm = {"G": {"Muscle_Skeletal": 30.0}}
    ess = {"G": (-0.5, 0.4)}
    expected = 0.4 * math.log2(30.0 + 1.0) * m.essentiality_multiplier(-0.5, 0.4)
    assert abs(m.penalty_contribution(hit, tpm, "Muscle_Skeletal", ess) - expected) < 1e-9


def test_interval_overlap():
    m = load_module()
    assert m.intervals_overlap(10, 20, 15, 25)
    assert m.intervals_overlap(15, 25, 10, 20)
    assert not m.intervals_overlap(10, 20, 20, 30)   # half-open, touching is not overlap
    assert not m.intervals_overlap(10, 20, 25, 30)


def test_scan_tox_motifs_flags_gquad_and_cpg():
    m = load_module()
    rules = m.load_motifs(SKILL_DIR / "data" / "tox_motifs.json")
    res_gquad = m.scan_tox_motifs("AAGGGGAATT", rules)
    ids = {f["id"] for f in res_gquad["flags"]}
    assert "g_quadruplex" in ids
    assert res_gquad["severity"] == "exclude"
    res_cpg = m.scan_tox_motifs("AAACGAAA", rules)
    assert "cpg_tlr9" in {f["id"] for f in res_cpg["flags"]}
    clean = m.scan_tox_motifs("AATTAATTAA", rules)
    assert clean["flags"] == []
    assert clean["tox_penalty"] == 0.0


def test_pareto_front_excludes_dominated():
    m = load_module()
    cands = [
        {"aso_id": "good", "efficacy": 0.9, "offtarget_penalty": 0.0, "clip_collision": 0.0,
         "essentiality_hit": 0.0, "seq_tox": 0.0},
        {"aso_id": "dominated", "efficacy": 0.5, "offtarget_penalty": 5.0, "clip_collision": 1.0,
         "essentiality_hit": 2.0, "seq_tox": 1.0},
        {"aso_id": "tradeoff", "efficacy": 0.95, "offtarget_penalty": 8.0, "clip_collision": 0.0,
         "essentiality_hit": 0.0, "seq_tox": 0.0},
    ]
    front = {c["aso_id"] for c in m.pareto_front(cands)}
    assert "good" in front
    assert "tradeoff" in front       # higher efficacy, non-dominated
    assert "dominated" not in front


# --------------------------------------------------------------------------- #
# Demo / integration tests
# --------------------------------------------------------------------------- #

def _expected_candidate_count():
    spec = json.loads(DEMO_SPEC.read_text())
    w = spec["target_window"]
    width = w["end"] - w["start"]
    return sum(width - L + 1 for L in range(spec["aso_length_min"], spec["aso_length_max"] + 1))


def test_screen_demo_end_to_end():
    m = load_module()
    result = m.screen_from_spec(DEMO_SPEC)
    assert result["skill"] == "aso-off-target-screening"
    assert result["target"]["tissue"] == "Muscle_Skeletal"
    assert result["summary"]["n_candidates"] == _expected_candidate_count()

    cands = {c["aso_id"]: c for c in result["candidates"]}

    # (1) Weighted off-target penalty: the worst offender's top hit is the
    #     broadly-expressed essential gene.
    worst = max(result["candidates"], key=lambda c: c["offtarget_penalty"])
    assert worst["offtarget_penalty"] > 0
    assert worst["offtarget_hits"][0]["gene"] == "ESSGENE"

    # (2) Tissue weighting: a silent-in-muscle off-target never contributes.
    for c in result["candidates"]:
        for h in c["offtarget_hits"]:
            if h["gene"] == "SILENTGENE":
                assert h["penalty_contribution"] == 0.0

    # (3) CLIP collision on pre-mRNA (intron) is caught as an unintended liability.
    intronic = [c for c in result["candidates"]
                if any(h["gene"] == "INTRONGENE" and h["region"] == "intron"
                       for h in c["offtarget_hits"])]
    assert intronic
    assert any(c["clip_collision"] > 0 for c in intronic)

    # (4) Intended-mechanism nuance: some on-target sites overlap a mechanistic
    #     silencer (asset), which is NOT counted as a collision liability.
    assert any(c["intended_mechanism"] for c in result["candidates"])

    # (5) Sequence-intrinsic tox can exclude a candidate even at zero penalty.
    assert any(c["verdict"] == "exclude" for c in result["candidates"])

    # (6) Pareto front preserves genuine trade-offs (>1 non-dominated candidate),
    #     is a proper subset, and no member is dominated.
    front_ids = set(result["pareto_front"])
    assert 1 < len(front_ids) < len(result["candidates"])
    front = [cands[i] for i in front_ids]
    for other in result["candidates"]:
        for f in front:
            assert not m._dominates(other, f) or other["aso_id"] == f["aso_id"]


def test_demo_cli_writes_output_contract(tmp_path):
    out = tmp_path / "aso_out"
    completed = subprocess.run(
        [sys.executable, str(MODULE_PATH), "--demo", "--output", str(out)],
        text=True, capture_output=True, check=True,
    )
    assert "aso-off-target-screening" in completed.stdout.lower() or "ASO" in completed.stdout
    report = (out / "report.md").read_text(encoding="utf-8")
    assert DISCLAIMER_FRAGMENT in report
    assert "Pareto" in report
    result = json.loads((out / "result.json").read_text(encoding="utf-8"))
    assert result["skill"] == "aso-off-target-screening"
    assert result["summary"]["n_candidates"] == _expected_candidate_count()
    assert (out / "tables" / "candidates.csv").exists()
    assert (out / "tables" / "offtarget_hits.csv").exists()
    assert (out / "tables" / "pareto_front.csv").exists()
    assert (out / "reproducibility" / "commands.sh").exists()


def test_cli_rejects_malformed_input(tmp_path):
    bad = tmp_path / "bad.json"
    bad.write_text('{"target_id": "x"}', encoding="utf-8")   # missing required fields
    completed = subprocess.run(
        [sys.executable, str(MODULE_PATH), "--input", str(bad), "--output", str(tmp_path / "o")],
        text=True, capture_output=True, check=False,
    )
    assert completed.returncode == 2
    assert "ERROR:" in completed.stderr
    assert "Traceback" not in completed.stderr


def test_cli_rejects_missing_input_file(tmp_path):
    completed = subprocess.run(
        [sys.executable, str(MODULE_PATH), "--input", str(tmp_path / "nope.json"),
         "--output", str(tmp_path / "o")],
        text=True, capture_output=True, check=False,
    )
    assert completed.returncode == 2
    assert "ERROR:" in completed.stderr
    assert "Traceback" not in completed.stderr


def test_load_motifs_accepts_list_form(tmp_path):
    m = load_module()
    obj_form = m.load_motifs(SKILL_DIR / "data" / "tox_motifs.json")
    assert isinstance(obj_form, list) and obj_form
    list_file = tmp_path / "motifs_list.json"
    list_file.write_text(json.dumps([{"id": "cpg_tlr9", "patterns": ["CG"], "weight": 1.0}]),
                         encoding="utf-8")
    list_form = m.load_motifs(list_file)
    assert isinstance(list_form, list) and list_form[0]["id"] == "cpg_tlr9"


def test_cli_tissue_override_changes_penalty(tmp_path):
    # SILENTGENE is silent in muscle but highly expressed in liver, so switching
    # the tissue to Liver should surface a non-zero SILENTGENE contribution.
    m = load_module()
    liver = m.screen_from_spec(DEMO_SPEC, tissue_override="Liver")
    contribs = [h["penalty_contribution"]
                for c in liver["candidates"] for h in c["offtarget_hits"]
                if h["gene"] == "SILENTGENE"]
    assert any(x > 0 for x in contribs)
