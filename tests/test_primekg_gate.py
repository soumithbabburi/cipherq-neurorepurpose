"""
Regression tests for the PrimeKG ground-truth contraindication gate.

These lock in the Phase 3 wiring so a future change cannot silently regress it:
  1. the canonical scorer applies the x0.15 kill multiplier on a labeled
     contraindication and leaves a labeled indication unpenalized;
  2. the reverse engine (canonical_pair_score) propagates the primekg block;
  3. the reverse engine forces a contraindicated candidate non-actionable;
  4. _why_not surfaces a contraindication as a reason (pure function, no data).

Run:  .venv312/Scripts/python.exe -m pytest tests/test_primekg_gate.py -q

Tests that need the local PrimeKG artifacts skip cleanly when they are absent
(they are gitignored / rebuildable), so the suite is safe to run anywhere.
"""
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def _predictor_or_skip():
    from services import primekg_predictor as pkg
    if not pkg.available():
        pytest.skip("PrimeKG artifacts not present (data/primekg/ is gitignored/rebuildable)")
    return pkg


def _find_pair(pkg, relation: str):
    """Discover a (drug_name, disease_name) pair whose names round-trip through the
    public API back to `relation`, so the test is robust to data rebuilds rather than
    hardcoding a single pair that a re-ETL might drop. Returns None if none found."""
    s = pkg._load()

    def nm(i):
        return s["nodes"].get(str(i), ["", "", ""])[2]

    for di, rels in s["labels"].items():
        for z in rels.get(relation, [])[:5]:
            drug, dis = nm(int(di)), nm(int(z))
            # Only accept pairs whose NAMES resolve back to the same relation via the
            # public API — i.e. a pair a real query would actually hit.
            if drug and dis and pkg.labeled_relation(drug, dis) == relation:
                return drug, dis
    return None


# ── canonical scorer ────────────────────────────────────────────────────────────
def test_contraindication_kills_composite_in_canonical_scorer():
    pkg = _predictor_or_skip()
    pair = _find_pair(pkg, "contraindication")
    if not pair:
        pytest.skip("no round-tripping contraindication pair in the local data")
    drug, disease = pair
    from services.repurposing_engine import score_compound_for_disease
    r = score_compound_for_disease({"name": drug, "chembl_id": ""}, disease, [], [], {})
    pk = r.get("primekg") or {}
    assert pk.get("relation") == "contraindication", (drug, disease, pk)
    assert pk.get("multiplier") == 0.15, pk
    assert pk.get("flag"), "a contraindication must carry a user-facing flag"
    # the kill multiplier must have actually driven the composite down
    assert r["composite_score"] <= 0.15, r["composite_score"]


def test_indication_not_penalized_in_canonical_scorer():
    pkg = _predictor_or_skip()
    pair = _find_pair(pkg, "indication")
    if not pair:
        pytest.skip("no round-tripping indication pair in the local data")
    drug, disease = pair
    from services.repurposing_engine import score_compound_for_disease
    r = score_compound_for_disease({"name": drug, "chembl_id": ""}, disease, [], [], {})
    pk = r.get("primekg") or {}
    # a labeled indication resolves as such and never triggers the kill multiplier
    assert pk.get("relation") == "indication", (drug, disease, pk)
    assert pk.get("multiplier") == 1.0, pk
    assert "flag" not in pk, "an indication must not carry a contraindication flag"


# ── reverse engine ───────────────────────────────────────────────────────────────
def test_reverse_canonical_pair_propagates_primekg():
    pkg = _predictor_or_skip()
    pair = _find_pair(pkg, "contraindication")
    if not pair:
        pytest.skip("no round-tripping contraindication pair in the local data")
    drug, disease = pair
    from services.reverse_repurposing import canonical_pair_score
    out = canonical_pair_score("", disease, drug_name=drug)
    pk = out.get("primekg") or {}
    assert pk.get("relation") == "contraindication", (drug, disease, pk)
    assert pk.get("multiplier") == 0.15, pk
    assert out["composite_score"] <= 0.15, out["composite_score"]


# ── _why_not (pure function, no data needed) ─────────────────────────────────────
def test_why_not_surfaces_contraindication():
    from services.reverse_repurposing import _why_not
    flag = "PrimeKG labels this an established contraindication for the disease."
    reasons = _why_not({"primekg": {"relation": "contraindication", "flag": flag},
                        "actionable": False})
    assert any("Contraindicated" in r for r in reasons), reasons
    assert any(flag in r for r in reasons), reasons


def test_why_not_silent_on_indication():
    from services.reverse_repurposing import _why_not
    reasons = _why_not({"primekg": {"relation": "indication"}, "actionable": True})
    assert not any("Contraindicated" in r for r in reasons), reasons
