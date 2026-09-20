"""Tests for the hpo-differential skill.

Red/green TDD contract (ClawBio CLAUDE.md, Development Rules): these tests were
written before the implementation and define the expected behaviour.

Scientific framework encoded here:
  * HPO IDs are exactly "HP:" + 7 digits. An ID that does not resolve against the
    ontology is NEVER silently carried - a wrong ID is worse than none.
  * Information content IC(t) = -log(|entities annotated with t or a descendant| / N).
  * Disease/gene similarity = symmetric best-match-average (funSimAvg) Resnik,
    where Resnik(a,b) = IC of the most informative common ancestor.
  * A "direct" phenotype hit means the entity is annotated with the query term or
    a MORE SPECIFIC descendant of it; an "ancestor" hit (entity annotated with a
    more general term) counts half in IC-weighted recall.

The ontology-dependent maths is tested against a small synthetic DAG so the tests
run without the 60 MB pyhpo data bundle. Tests that need the real ontology are
marked `integration` and skipped when pyhpo is unavailable.
"""

import json
import math
import subprocess
import sys
from pathlib import Path

import pytest

SKILL_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(SKILL_DIR))

from hpo_differential import (  # noqa: E402
    DISCLAIMER,
    MiniOntology,
    best_match_average,
    build_ic,
    demo_query,
    discriminating_entities,
    ic_weighted_recall,
    parse_terms,
    rank_entities,
    resnik,
)


# --------------------------------------------------------------------------
# A small synthetic DAG standing in for HPO:
#
#            root
#           /    \
#      abnorm_eye  abnorm_brain
#        /    \          |
#   retina   lens     seizure
#      |                  |
#  ret_hem             focal_seizure
# --------------------------------------------------------------------------
EDGES = {
    "HP:0000001": [],
    "HP:0000002": ["HP:0000001"],  # abnormal eye
    "HP:0000003": ["HP:0000001"],  # abnormal brain
    "HP:0000004": ["HP:0000002"],  # retinal abnormality
    "HP:0000005": ["HP:0000002"],  # lens abnormality
    "HP:0000006": ["HP:0000003"],  # seizure
    "HP:0000007": ["HP:0000004"],  # retinal hemorrhage
    "HP:0000008": ["HP:0000006"],  # focal-onset seizure
}


@pytest.fixture
def onto():
    return MiniOntology(EDGES)


@pytest.fixture
def annotations():
    """Five synthetic 'diseases'."""
    return {
        "D1": ["HP:0000007", "HP:0000008"],  # ret hem + focal seizure
        "D2": ["HP:0000007"],                # ret hem only
        "D3": ["HP:0000005"],                # lens only
        "D4": ["HP:0000008"],                # focal seizure only
        "D5": ["HP:0000004", "HP:0000006"],  # general retina + general seizure
    }


# --------------------------------------------------------------------------
# Term parsing and validation
# --------------------------------------------------------------------------
def test_parse_terms_accepts_canonical_ids():
    assert parse_terms("HP:0001263,HP:0000252") == ["HP:0001263", "HP:0000252"]


def test_parse_terms_handles_newlines_labels_and_whitespace():
    raw = "HP:0001263 Global developmental delay\n  HP:0000252 Microcephaly \n"
    assert parse_terms(raw) == ["HP:0001263", "HP:0000252"]


def test_parse_terms_deduplicates_preserving_first_occurrence_order():
    assert parse_terms("HP:0000252,HP:0001263,HP:0000252") == [
        "HP:0000252",
        "HP:0001263",
    ]


def test_parse_terms_rejects_hpo_prefix_typo():
    """'HPO:' is the documented failure mode - it must raise, never be coerced."""
    with pytest.raises(ValueError, match="HP:"):
        parse_terms("HPO:0001263")


def test_parse_terms_rejects_wrong_digit_count():
    with pytest.raises(ValueError):
        parse_terms("HP:12345")


def test_parse_terms_rejects_empty_input():
    with pytest.raises(ValueError):
        parse_terms("   ")


# --------------------------------------------------------------------------
# Ontology closure
# --------------------------------------------------------------------------
def test_ancestors_include_self_and_transitive_parents(onto):
    assert onto.ancestors("HP:0000007") == {
        "HP:0000007", "HP:0000004", "HP:0000002", "HP:0000001",
    }


def test_descendants_include_self_and_transitive_children(onto):
    assert onto.descendants("HP:0000002") == {
        "HP:0000002", "HP:0000004", "HP:0000005", "HP:0000007",
    }


def test_leaf_descendants_is_just_itself(onto):
    assert onto.descendants("HP:0000007") == {"HP:0000007"}


# --------------------------------------------------------------------------
# Information content
# --------------------------------------------------------------------------
def test_ic_is_zero_for_the_root(onto, annotations):
    ic, n = build_ic(annotations.values(), onto)
    assert n == 5
    assert ic["HP:0000001"] == pytest.approx(0.0)


def test_ic_counts_propagate_up_the_dag(onto, annotations):
    """A disease annotated with a leaf also counts towards every ancestor."""
    ic, n = build_ic(annotations.values(), onto)
    # retinal hemorrhage annotated in D1, D2 only -> 2/5
    assert ic["HP:0000007"] == pytest.approx(-math.log(2 / 5))
    # retinal abnormality inherited from D1, D2 plus direct in D5 -> 3/5
    assert ic["HP:0000004"] == pytest.approx(-math.log(3 / 5))


def test_rarer_term_has_strictly_higher_ic(onto, annotations):
    ic, _ = build_ic(annotations.values(), onto)
    assert ic["HP:0000007"] > ic["HP:0000004"] > ic["HP:0000002"]


# --------------------------------------------------------------------------
# Resnik similarity
# --------------------------------------------------------------------------
def test_resnik_of_a_term_with_itself_is_its_own_ic(onto, annotations):
    ic, _ = build_ic(annotations.values(), onto)
    assert resnik("HP:0000007", "HP:0000007", ic, onto) == pytest.approx(ic["HP:0000007"])


def test_resnik_uses_most_informative_common_ancestor(onto, annotations):
    ic, _ = build_ic(annotations.values(), onto)
    # ret_hem vs lens -> MICA is abnormal eye
    assert resnik("HP:0000007", "HP:0000005", ic, onto) == pytest.approx(ic["HP:0000002"])


def test_resnik_across_disjoint_branches_is_root_level(onto, annotations):
    ic, _ = build_ic(annotations.values(), onto)
    assert resnik("HP:0000007", "HP:0000008", ic, onto) == pytest.approx(ic["HP:0000001"])


def test_resnik_is_symmetric(onto, annotations):
    ic, _ = build_ic(annotations.values(), onto)
    a = resnik("HP:0000007", "HP:0000005", ic, onto)
    b = resnik("HP:0000005", "HP:0000007", ic, onto)
    assert a == pytest.approx(b)


# --------------------------------------------------------------------------
# Best-match average
# --------------------------------------------------------------------------
def test_bma_of_identical_sets_equals_mean_ic(onto, annotations):
    ic, _ = build_ic(annotations.values(), onto)
    q = ["HP:0000007", "HP:0000008"]
    score, _ = best_match_average(q, q, ic, onto)
    expected = (ic["HP:0000007"] + ic["HP:0000008"]) / 2
    assert score == pytest.approx(expected)


def test_bma_is_symmetric(onto, annotations):
    ic, _ = build_ic(annotations.values(), onto)
    q, d = ["HP:0000007"], ["HP:0000005", "HP:0000008"]
    assert best_match_average(q, d, ic, onto)[0] == pytest.approx(
        best_match_average(d, q, ic, onto)[0]
    )


def test_bma_returns_zero_for_empty_sets(onto, annotations):
    ic, _ = build_ic(annotations.values(), onto)
    assert best_match_average([], ["HP:0000007"], ic, onto)[0] == 0.0


def test_bma_reports_the_matching_partner_per_query_term(onto, annotations):
    ic, _ = build_ic(annotations.values(), onto)
    _, matches = best_match_average(["HP:0000007"], ["HP:0000004"], ic, onto)
    assert matches[0][0] == "HP:0000007"
    assert matches[0][1] == "HP:0000004"


# --------------------------------------------------------------------------
# IC-weighted recall: direct vs ancestor hits
# --------------------------------------------------------------------------
def test_recall_is_one_when_every_query_term_is_directly_annotated(onto, annotations):
    ic, _ = build_ic(annotations.values(), onto)
    q = ["HP:0000007", "HP:0000008"]
    rec, direct, anc = ic_weighted_recall(q, annotations["D1"], ic, onto)
    assert rec == pytest.approx(1.0)
    assert set(direct) == set(q)
    assert anc == []


def test_descendant_annotation_counts_as_a_direct_hit(onto, annotations):
    """Query 'retinal abnormality'; disease annotated with the more specific
    'retinal hemorrhage' - that IS the feature, so it is a direct hit."""
    ic, _ = build_ic(annotations.values(), onto)
    rec, direct, anc = ic_weighted_recall(["HP:0000004"], ["HP:0000007"], ic, onto)
    assert direct == ["HP:0000004"]
    assert rec == pytest.approx(1.0)


def test_ancestor_only_annotation_counts_half(onto, annotations):
    """Query 'retinal hemorrhage'; disease annotated only with the more general
    'retinal abnormality' - partial credit, not full."""
    ic, _ = build_ic(annotations.values(), onto)
    rec, direct, anc = ic_weighted_recall(["HP:0000007"], ["HP:0000004"], ic, onto)
    assert direct == []
    assert anc == ["HP:0000007"]
    assert rec == pytest.approx(0.5)


def test_recall_is_zero_when_nothing_matches(onto, annotations):
    ic, _ = build_ic(annotations.values(), onto)
    rec, direct, anc = ic_weighted_recall(["HP:0000007"], ["HP:0000006"], ic, onto)
    assert rec == 0.0 and direct == [] and anc == []


def test_recall_weights_by_ic_not_by_term_count(onto, annotations):
    """Matching the one rare term beats matching the one common term."""
    ic, _ = build_ic(annotations.values(), onto)
    q = ["HP:0000007", "HP:0000006"]  # rare, common
    rare_only, _, _ = ic_weighted_recall(q, ["HP:0000007"], ic, onto)
    common_only, _, _ = ic_weighted_recall(q, ["HP:0000006"], ic, onto)
    assert rare_only > common_only


# --------------------------------------------------------------------------
# Ranking
# --------------------------------------------------------------------------
def test_ranking_puts_the_fully_matching_entity_first(onto, annotations):
    q = ["HP:0000007", "HP:0000008"]
    results = rank_entities(q, annotations, onto, min_annotations=1)
    assert results[0]["id"] == "D1"
    assert results[0]["n_direct"] == 2


def test_ranking_is_sorted_by_descending_composite(onto, annotations):
    q = ["HP:0000007", "HP:0000008"]
    results = rank_entities(q, annotations, onto, min_annotations=1)
    scores = [r["composite"] for r in results]
    assert scores == sorted(scores, reverse=True)


def test_ranking_respects_min_annotations_filter(onto, annotations):
    q = ["HP:0000007", "HP:0000008"]
    results = rank_entities(q, annotations, onto, min_annotations=2)
    assert {r["id"] for r in results} == {"D1", "D5"}


def test_ranking_composite_is_similarity_times_recall(onto, annotations):
    q = ["HP:0000007", "HP:0000008"]
    for r in rank_entities(q, annotations, onto, min_annotations=1):
        assert r["composite"] == pytest.approx(r["funsimavg"] * r["ic_recall"], abs=1e-6)


def test_ranking_never_returns_an_unresolvable_term(onto, annotations):
    with pytest.raises(KeyError):
        rank_entities(["HP:9999999"], annotations, onto, min_annotations=1)


# --------------------------------------------------------------------------
# Discriminating-feature analysis
# --------------------------------------------------------------------------
def test_discriminating_entities_requires_all_terms_present(onto, annotations):
    hits = discriminating_entities(["HP:0000007", "HP:0000008"], annotations, onto)
    assert hits == ["D1"]


def test_discriminating_entities_matches_via_descendants(onto, annotations):
    hits = discriminating_entities(["HP:0000004"], annotations, onto)
    assert set(hits) == {"D1", "D2", "D5"}


def test_discriminating_entities_returns_empty_when_no_entity_has_all(onto, annotations):
    assert discriminating_entities(["HP:0000005", "HP:0000008"], annotations, onto) == []


# --------------------------------------------------------------------------
# Demo data
# --------------------------------------------------------------------------
def test_demo_query_is_valid_and_non_trivial():
    terms = demo_query()
    assert len(terms) >= 5
    assert all(t.startswith("HP:") and len(t) == 10 for t in terms)
    assert parse_terms(",".join(terms)) == terms


def test_disclaimer_is_the_clawbio_standard_text():
    assert "not a medical device" in DISCLAIMER
    assert "Consult a healthcare professional" in DISCLAIMER


# --------------------------------------------------------------------------
# Output contract (see SKILL.md ## Output Structure)
# --------------------------------------------------------------------------
def _interpreter_with_pyhpo() -> str | None:
    """Find a Python that can import pyhpo.

    `uv run pytest` executes pytest from an isolated tool environment that does
    not carry the project's optional deps, so sys.executable is usually NOT the
    project interpreter. Fall back to the project venv before giving up.
    """
    candidates = [sys.executable]
    for rel in (".venv/bin/python", ".venv/bin/python3"):
        p = SKILL_DIR.parent.parent / rel
        if p.exists():
            candidates.append(str(p))
    for exe in candidates:
        try:
            probe = subprocess.run(
                [exe, "-c", "import pyhpo"], capture_output=True, timeout=120
            )
        except (OSError, subprocess.SubprocessError):
            continue
        if probe.returncode == 0:
            return exe
    return None


@pytest.mark.integration
def test_output_contract_demo_produces_every_declared_artifact(tmp_path):
    exe = _interpreter_with_pyhpo()
    if exe is None:
        pytest.skip("pyhpo not installed (pip install 'clawbio[hpo]')")
    out = tmp_path / "demo"
    proc = subprocess.run(
        [exe, str(SKILL_DIR / "hpo_differential.py"),
         "--demo", "--output", str(out), "--top", "10"],
        capture_output=True, text=True, timeout=1800,
    )
    assert proc.returncode == 0, proc.stderr
    for rel in ("report.md", "result.json",
                "tables/ranked_diseases.csv", "tables/ranked_genes.csv",
                "reproducibility/commands.sh", "reproducibility/environment.yml"):
        assert (out / rel).exists(), f"missing declared artifact: {rel}"
    assert DISCLAIMER in (out / "report.md").read_text()
    payload = json.loads((out / "result.json").read_text())
    assert payload["hpo_release"]
    assert all(t["verified"] for t in payload["query_verification"])
