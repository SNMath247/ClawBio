#!/usr/bin/env python3
"""hpo-differential - offline phenotype-driven rare disease and gene ranking.

Takes a list of HPO term IDs and ranks every OMIM disease, Orphanet disease and
disease-associated gene in the Human Phenotype Ontology annotation corpus by
phenotypic similarity to that list. Runs entirely offline against the HPO
release bundled with pyhpo - no API calls, no patient data leaves the machine.

Scoring
-------
IC(t)     = -log(|entities annotated with t or a descendant| / N)
Resnik    = IC of the most informative common ancestor of two terms
funSimAvg = symmetric best-match-average Resnik over (query, target) term sets
recall    = IC-weighted fraction of query terms the target explains, where a
            descendant annotation counts fully and an ancestor annotation counts half
composite = funSimAvg * recall        <- the ranking key

See SKILL.md for methodology, gotchas and the output contract.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import subprocess
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Set, Tuple

DISCLAIMER = (
    "ClawBio is a research and educational tool. It is not a medical device and "
    "does not provide clinical diagnoses. Consult a healthcare professional before "
    "making any medical decisions."
)

HPO_ID_RE = re.compile(r"^HP:\d{7}$")
_TOKEN_RE = re.compile(r"\b([A-Za-z]+):(\d+)\b")

# Synthetic demo phenotype set: a textbook ciliopathy-flavoured vignette.
# Synthetic only - never real patient data.
DEMO_TERMS = [
    "HP:0000510",  # Rod-cone dystrophy
    "HP:0001263",  # Global developmental delay
    "HP:0001829",  # Postaxial polydactyly (demo)
    "HP:0000135",  # Hypogonadism
    "HP:0001513",  # Obesity
    "HP:0000105",  # Enlarged kidney
]


# ---------------------------------------------------------------------------
# Term parsing
# ---------------------------------------------------------------------------
def parse_terms(raw: str) -> List[str]:
    """Extract and validate HPO IDs from free-form text.

    Accepts comma-, whitespace- or newline-separated IDs, with or without
    trailing labels. Deduplicates preserving first-occurrence order.

    Raises ValueError on a malformed ID rather than dropping it: a wrong HPO ID
    is worse than no HPO ID, because downstream vocabulary validation cannot
    tell a typo from a real rare term.
    """
    if not raw or not raw.strip():
        raise ValueError("no HPO terms supplied")

    out: List[str] = []
    seen: Set[str] = set()
    for prefix, digits in _TOKEN_RE.findall(raw):
        candidate = f"{prefix}:{digits}"
        if prefix.upper() == "HPO":
            raise ValueError(
                f"{candidate!r}: HPO IDs use the prefix 'HP:' followed by exactly "
                "7 digits (e.g. HP:0001250), never 'HPO:'"
            )
        if prefix.upper() != "HP":
            continue  # not an HPO token - e.g. a PMID or OMIM id in a pasted line
        if not HPO_ID_RE.match(candidate):
            raise ValueError(
                f"{candidate!r}: HPO IDs are 'HP:' + exactly 7 digits "
                f"(got {len(digits)} digits)"
            )
        if candidate not in seen:
            seen.add(candidate)
            out.append(candidate)

    if not out:
        raise ValueError(f"no valid HP: IDs found in input: {raw[:80]!r}")
    return out


# ---------------------------------------------------------------------------
# Ontology adapters
# ---------------------------------------------------------------------------
class MiniOntology:
    """In-memory DAG adapter. Used by the test suite and as the reference
    implementation of the closure semantics the pyhpo adapter must match."""

    def __init__(self, parents: Dict[str, Sequence[str]]):
        self._parents = {k: list(v) for k, v in parents.items()}
        self._children: Dict[str, List[str]] = defaultdict(list)
        for child, plist in self._parents.items():
            for p in plist:
                self._children[p].append(child)
        self._anc: Dict[str, Set[str]] = {}
        self._desc: Dict[str, Set[str]] = {}

    def __contains__(self, tid: str) -> bool:
        return tid in self._parents

    def ancestors(self, tid: str) -> Set[str]:
        if tid not in self._parents:
            raise KeyError(tid)
        if tid in self._anc:
            return self._anc[tid]
        out, stack = {tid}, list(self._parents[tid])
        while stack:
            cur = stack.pop()
            if cur not in out:
                out.add(cur)
                stack.extend(self._parents.get(cur, []))
        self._anc[tid] = out
        return out

    def descendants(self, tid: str) -> Set[str]:
        if tid not in self._parents:
            raise KeyError(tid)
        if tid in self._desc:
            return self._desc[tid]
        out, stack = {tid}, list(self._children.get(tid, []))
        while stack:
            cur = stack.pop()
            if cur not in out:
                out.add(cur)
                stack.extend(self._children.get(cur, []))
        self._desc[tid] = out
        return out

    def label(self, tid: str) -> str:
        return tid


class PyHPOOntology:
    """Adapter over the HPO release bundled with pyhpo."""

    def __init__(self):
        try:
            from pyhpo import Ontology  # noqa: WPS433
        except ImportError as exc:  # pragma: no cover - environment dependent
            raise SystemExit(
                "hpo-differential requires the 'pyhpo' package, which bundles the "
                "HPO ontology and annotations.\n  pip install pyhpo"
            ) from exc
        self._ont = Ontology()
        self._anc: Dict[str, Set[str]] = {}
        self._desc: Dict[str, Set[str]] = {}

    @property
    def raw(self):
        return self._ont

    def _term(self, tid: str):
        try:
            return self._ont.get_hpo_object(tid)
        except Exception as exc:
            raise KeyError(tid) from exc

    def __contains__(self, tid: str) -> bool:
        try:
            self._term(tid)
            return True
        except KeyError:
            return False

    def ancestors(self, tid: str) -> Set[str]:
        if tid in self._anc:
            return self._anc[tid]
        t = self._term(tid)
        out = {t.id} | {p.id for p in t.all_parents}
        self._anc[tid] = out
        return out

    def descendants(self, tid: str) -> Set[str]:
        if tid in self._desc:
            return self._desc[tid]
        t = self._term(tid)
        out, stack = {t.id}, list(t.children)
        while stack:
            cur = stack.pop()
            if cur.id not in out:
                out.add(cur.id)
                stack.extend(cur.children)
        self._desc[tid] = out
        return out

    def label(self, tid: str) -> str:
        return self._term(tid).name

    def release(self) -> str:
        import pyhpo  # noqa: WPS433

        obo = Path(pyhpo.__file__).parent / "data" / "hp.obo"
        try:
            with obo.open() as fh:
                for line in fh:
                    if line.startswith("data-version:"):
                        return line.split("hp/releases/")[-1].strip()
        except OSError:  # pragma: no cover
            pass
        return "unknown"


# ---------------------------------------------------------------------------
# Scoring primitives
# ---------------------------------------------------------------------------
def build_ic(annotation_sets: Iterable[Sequence[str]], onto) -> Tuple[Dict[str, float], int]:
    """Information content per term from annotation frequency, with closure."""
    counts: Dict[str, int] = defaultdict(int)
    n = 0
    for terms in annotation_sets:
        n += 1
        closed: Set[str] = set()
        for tid in terms:
            closed |= onto.ancestors(tid)
        for tid in closed:
            counts[tid] += 1
    if n == 0:
        return {}, 0
    return {tid: -math.log(c / n) for tid, c in counts.items()}, n


def resnik(a: str, b: str, ic: Dict[str, float], onto) -> float:
    """IC of the most informative common ancestor of two terms."""
    common = onto.ancestors(a) & onto.ancestors(b)
    if not common:
        return 0.0
    return max(ic.get(t, 0.0) for t in common)


def best_match_average(
    query: Sequence[str],
    target: Sequence[str],
    ic: Dict[str, float],
    onto,
    cache: Dict[Tuple[str, str], float] | None = None,
) -> Tuple[float, List[Tuple[str, str | None, float]]]:
    """Symmetric best-match-average (funSimAvg) Resnik similarity."""
    if not query or not target:
        return 0.0, []
    if cache is None:
        cache = {}

    def sim(q: str, d: str) -> float:
        key = (q, d)
        if key not in cache:
            cache[key] = resnik(q, d, ic, onto)
        return cache[key]

    best_q, matches = [], []
    for q in query:
        best, arg = 0.0, None
        for d in target:
            v = sim(q, d)
            if v > best:
                best, arg = v, d
        best_q.append(best)
        matches.append((q, arg, round(best, 4)))

    best_d = [max((sim(q, d) for q in query), default=0.0) for d in target]
    score = 0.5 * (sum(best_q) / len(best_q) + sum(best_d) / len(best_d))
    return score, matches


def ic_weighted_recall(
    query: Sequence[str],
    target: Sequence[str],
    ic: Dict[str, float],
    onto,
    q_desc: Dict[str, Set[str]] | None = None,
) -> Tuple[float, List[str], List[str]]:
    """Fraction of the query's information content the target explains.

    A target annotated with the query term or a MORE SPECIFIC descendant counts
    fully (it has that feature). A target annotated only with a MORE GENERAL
    ancestor counts half (it may have that feature).
    """
    if not query:
        return 0.0, [], []
    if q_desc is None:
        q_desc = {q: onto.descendants(q) for q in query}

    tset = set(target)
    # direct: the target carries the query term or a MORE SPECIFIC descendant.
    direct = [q for q in query if tset & q_desc[q]]
    # ancestor: the target carries only a MORE GENERAL form of the query term,
    # i.e. one of the query term's proper ancestors is in the target's own set.
    # (Note the direction: it is the QUERY term's ancestors we intersect with the
    # target annotations, not the target's ancestors with the query.)
    ancestor = [
        q for q in query
        if q not in direct and (onto.ancestors(q) - {q}) & tset
    ]

    total = sum(ic.get(q, 0.0) for q in query)
    if total <= 0:
        return 0.0, direct, ancestor
    got = sum(ic.get(q, 0.0) for q in direct) + 0.5 * sum(ic.get(q, 0.0) for q in ancestor)
    return got / total, direct, ancestor


# ---------------------------------------------------------------------------
# Ranking
# ---------------------------------------------------------------------------
def rank_entities(
    query: Sequence[str],
    annotations: Dict[str, Sequence[str]],
    onto,
    min_annotations: int = 3,
    names: Dict[str, str] | None = None,
    top: int | None = None,
) -> List[dict]:
    """Rank annotated entities (diseases or genes) against the query term set."""
    for q in query:
        onto.ancestors(q)  # raises KeyError on an unresolvable ID - never silently skip

    usable = {k: list(v) for k, v in annotations.items() if len(v) >= min_annotations}
    ic, _ = build_ic(usable.values(), onto)
    q_desc = {q: onto.descendants(q) for q in query}
    cache: Dict[Tuple[str, str], float] = {}

    results = []
    for key, terms in usable.items():
        sim, matches = best_match_average(query, terms, ic, onto, cache)
        rec, direct, ancestor = ic_weighted_recall(query, terms, ic, onto, q_desc)
        results.append(
            {
                "id": key,
                "name": (names or {}).get(key, key),
                "n_annotations": len(terms),
                "funsimavg": round(sim, 4),
                "ic_recall": round(rec, 4),
                "composite": round(sim * rec, 4),
                "n_direct": len(direct),
                "n_ancestor": len(ancestor),
                "direct_terms": direct,
                "ancestor_terms": ancestor,
                "unexplained_terms": [q for q in query if q not in direct and q not in ancestor],
                "best_matches": [m for m in matches if m[2] > 0][:20],
            }
        )
    results.sort(key=lambda r: (-r["composite"], -r["funsimavg"], r["id"]))
    return results[:top] if top else results


def discriminating_entities(
    terms: Sequence[str],
    annotations: Dict[str, Sequence[str]],
    onto,
) -> List[str]:
    """Entities annotated with EVERY given term (or a descendant of each).

    This is the annotation-independent half of the method: when a term
    combination is rare enough, the carrier list is short enough to read, and
    that list does not depend on any scoring choice.
    """
    dsets = [onto.descendants(t) for t in terms]
    return sorted(
        key for key, ann in annotations.items()
        if all(set(ann) & ds for ds in dsets)
    )


def demo_query() -> List[str]:
    return list(DEMO_TERMS)


# ---------------------------------------------------------------------------
# Corpus loading
# ---------------------------------------------------------------------------
def load_corpus(onto: PyHPOOntology):
    """Build {entity_id: [hpo ids]} maps for OMIM, Orphanet and genes."""
    raw = onto.raw
    omim, omim_names = {}, {}
    for d in raw.omim_diseases:
        omim[f"OMIM:{d.id}"] = [t.id for t in d.hpo_set()]
        omim_names[f"OMIM:{d.id}"] = d.name
    orpha, orpha_names = {}, {}
    for d in raw.orpha_diseases:
        orpha[f"ORPHA:{d.id}"] = [t.id for t in d.hpo_set()]
        orpha_names[f"ORPHA:{d.id}"] = d.name
    genes, gene_names = {}, {}
    for g in raw.genes:
        genes[g.name] = [t.id for t in g.hpo_set()]
        gene_names[g.name] = g.name
    return (omim, omim_names), (orpha, orpha_names), (genes, gene_names)


def disease_gene_map() -> Dict[str, List[str]]:
    """disease_id -> gene symbols, from the bundled phenotype_to_genes.txt."""
    import pyhpo  # noqa: WPS433

    path = Path(pyhpo.__file__).parent / "data" / "phenotype_to_genes.txt"
    out: Dict[str, Set[str]] = defaultdict(set)
    try:
        with path.open() as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                if row.get("disease_id"):
                    out[row["disease_id"]].add(row["gene_symbol"])
    except OSError:  # pragma: no cover
        return {}
    return {k: sorted(v) for k, v in out.items()}


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------
def write_outputs(out_dir: Path, payload: dict, argv: List[str]) -> None:
    (out_dir / "tables").mkdir(parents=True, exist_ok=True)
    (out_dir / "reproducibility").mkdir(parents=True, exist_ok=True)

    with (out_dir / "result.json").open("w") as fh:
        json.dump(payload, fh, indent=2)

    for kind, fname in (("omim", "ranked_diseases.csv"), ("gene", "ranked_genes.csv"),
                        ("orpha", "ranked_orphanet.csv")):
        rows = payload["rankings"].get(kind, [])
        with (out_dir / "tables" / fname).open("w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["rank", "id", "name", "genes", "composite", "funsimavg",
                        "ic_recall", "n_direct", "n_query", "unexplained_terms"])
            for i, r in enumerate(rows, 1):
                w.writerow([i, r["id"], r["name"], ";".join(r.get("genes", [])),
                            r["composite"], r["funsimavg"], r["ic_recall"],
                            r["n_direct"], len(payload["query_ids"]),
                            ";".join(r["unexplained_terms"])])

    with (out_dir / "reproducibility" / "commands.sh").open("w") as fh:
        fh.write("#!/usr/bin/env bash\nset -euo pipefail\n")
        fh.write(f"# generated {payload['generated_utc']}\n")
        fh.write("python " + " ".join(argv) + "\n")

    with (out_dir / "reproducibility" / "environment.yml").open("w") as fh:
        fh.write("name: hpo-differential\nchannels: [conda-forge]\ndependencies:\n")
        fh.write(f"  - python={sys.version_info.major}.{sys.version_info.minor}\n")
        fh.write("  - pip\n  - pip:\n")
        fh.write(f"      - pyhpo=={payload.get('pyhpo_version', 'unknown')}\n")
        fh.write(f"# HPO release: {payload['hpo_release']}\n")

    (out_dir / "report.md").write_text(render_report(payload))


def render_report(p: dict) -> str:
    lab = {t["hpo_id"]: t["ontology_label"] for t in p["query_verification"]}
    L = []
    L.append("# HPO Differential — phenotype-driven candidate ranking\n")
    L.append(f"> ⚠️ **VERIFY BEFORE CLINICAL USE.** {DISCLAIMER}\n")
    L.append(f"- HPO release: **{p['hpo_release']}**")
    L.append(f"- Generated (UTC): {p['generated_utc']}")
    L.append(f"- Corpus scored: {p['corpus_sizes']['omim']} OMIM · "
             f"{p['corpus_sizes']['orpha']} Orphanet · {p['corpus_sizes']['gene']} genes\n")

    L.append("## 1. Query term verification\n")
    L.append("| HPO ID | Ontology label | IC | # OMIM diseases |")
    L.append("|---|---|---|---|")
    for t in p["query_verification"]:
        ic = p["term_ic"].get(t["hpo_id"], {})
        L.append(f"| {t['hpo_id']} | {t['ontology_label']} | "
                 f"{ic.get('ic', 0):.2f} | {ic.get('n_omim', 0)} |")
    L.append("\nEvery ID above resolved against the ontology. "
             "An ID that does not resolve aborts the run rather than being dropped.\n")

    for kind, title in (("omim", "OMIM diseases"), ("orpha", "Orphanet diseases"),
                        ("gene", "Genes")):
        rows = p["rankings"].get(kind, [])
        if not rows:
            continue
        L.append(f"## Ranked {title}\n")
        L.append("| # | ID | Name | Gene(s) | composite | direct/total | Unexplained |")
        L.append("|---|---|---|---|---|---|---|")
        for i, r in enumerate(rows, 1):
            unexp = ", ".join(lab.get(x, x) for x in r["unexplained_terms"]) or "—"
            L.append(f"| {i} | {r['id']} | {r['name'][:56]} | "
                     f"{', '.join(r.get('genes', [])) or '—'} | {r['composite']:.3f} | "
                     f"{r['n_direct']}/{len(p['query_ids'])} | {unexp} |")
        L.append("")

    if p.get("discriminating"):
        L.append("## Discriminating-feature analysis\n")
        L.append("The rarest query terms, and every disease in the corpus carrying "
                 "them together. This result does not depend on the scoring formula.\n")
        for combo, hits in p["discriminating"].items():
            names = ", ".join(lab.get(x, x) for x in combo.split("+"))
            L.append(f"**{names}** — {len(hits)} disease(s): "
                     + (", ".join(f"`{h}`" for h in hits[:25]) or "none"))
            L.append("")

    L.append("## Method\n")
    L.append("- `IC(t) = -log(|entities annotated with t or a descendant| / N)`")
    L.append("- `Resnik(a,b)` = IC of the most informative common ancestor")
    L.append("- `funSimAvg` = symmetric best-match-average Resnik over the term sets")
    L.append("- `recall` = IC-weighted fraction of query terms explained "
             "(descendant match counts fully, ancestor match counts half)")
    L.append("- `composite = funSimAvg × recall` — the ranking key\n")
    L.append("**Annotation-depth caveat.** Recently described disorders have "
             "annotation sets tightly concentrated on the features that defined "
             "them, which flatters a best-match-average metric; long-established "
             "syndromes carry broad annotation sets that dilute it. Read the "
             "discriminating-feature analysis alongside the ranking, not after it.\n")
    L.append(f"> ⚠️ **VERIFY BEFORE CLINICAL USE.** {DISCLAIMER}")
    return "\n".join(L)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main(argv: List[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--terms", help="comma/space separated HPO IDs (labels tolerated)")
    src.add_argument("--terms-file", help="file containing HPO IDs, one per line")
    src.add_argument("--demo", action="store_true", help="run on the synthetic demo term set")
    ap.add_argument("--output", required=True, help="report directory")
    ap.add_argument("--top", type=int, default=25, help="rows per ranking (default 25)")
    ap.add_argument("--min-annotations", type=int, default=3,
                    help="skip entities with fewer annotations (default 3)")
    ap.add_argument("--force", action="store_true", help="overwrite a non-empty output directory")
    args = ap.parse_args(argv)

    if args.demo:
        query = demo_query()
    elif args.terms_file:
        query = parse_terms(Path(args.terms_file).read_text())
    else:
        query = parse_terms(args.terms)

    out_dir = Path(args.output)
    if out_dir.exists() and any(out_dir.iterdir()) and not args.force:
        print(f"refusing to overwrite non-empty {out_dir} (pass --force)", file=sys.stderr)
        return 2
    out_dir.mkdir(parents=True, exist_ok=True)

    onto = PyHPOOntology()
    verification = []
    for q in query:
        if q not in onto:
            print(f"HPO ID {q} does not resolve against release "
                  f"{onto.release()} — aborting rather than guessing", file=sys.stderr)
            return 3
        verification.append({"hpo_id": q, "ontology_label": onto.label(q), "verified": True})

    (omim, omim_names), (orpha, orpha_names), (genes, gene_names) = load_corpus(onto)
    d2g = disease_gene_map()

    rankings = {}
    for kind, (ann, names) in (("omim", (omim, omim_names)),
                               ("orpha", (orpha, orpha_names)),
                               ("gene", (genes, gene_names))):
        rows = rank_entities(query, ann, onto, args.min_annotations, names, top=args.top)
        for r in rows:
            r["genes"] = d2g.get(r["id"], []) if kind != "gene" else [r["id"]]
        rankings[kind] = rows

    ic_omim, n_omim = build_ic(
        [v for v in omim.values() if len(v) >= args.min_annotations], onto
    )
    counts: Dict[str, int] = defaultdict(int)
    for terms in omim.values():
        if len(terms) < args.min_annotations:
            continue
        closed: Set[str] = set()
        for t in terms:
            closed |= onto.ancestors(t)
        for t in closed:
            counts[t] += 1
    term_ic = {q: {"label": onto.label(q), "ic": round(ic_omim.get(q, 0.0), 3),
                   "n_omim": counts.get(q, 0)} for q in query}

    rarest = sorted(query, key=lambda q: -term_ic[q]["ic"])[:3]
    discriminating = {}
    if len(rarest) >= 2:
        pair = rarest[:2]
        discriminating["+".join(pair)] = discriminating_entities(pair, omim, onto)
    for q in rarest:
        discriminating[q] = discriminating_entities([q], omim, onto)

    try:
        import pyhpo  # noqa: WPS433
        pyhpo_version = getattr(pyhpo, "__version__", "unknown")
    except Exception:  # pragma: no cover
        pyhpo_version = "unknown"

    payload = {
        "skill": "hpo-differential",
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "hpo_release": onto.release(),
        "pyhpo_version": pyhpo_version,
        "query_ids": query,
        "query_verification": verification,
        "term_ic": term_ic,
        "corpus_sizes": {"omim": len(omim), "orpha": len(orpha), "gene": len(genes)},
        "rankings": rankings,
        "discriminating": discriminating,
        "disclaimer": DISCLAIMER,
    }

    argv_used = [str(Path(__file__).name)] + (argv or sys.argv[1:])
    write_outputs(out_dir, payload, argv_used)
    print(f"wrote {out_dir}/report.md  (top hit: "
          f"{rankings['omim'][0]['id']} {rankings['omim'][0]['name']})"
          if rankings.get("omim") else f"wrote {out_dir}/report.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
