"""
Unit tests for services.offtarget — the off-target / polypharmacology inference.

Pure (no DB): they lock in the promiscuity guard so a future change cannot
silently let a weak predicted off-target be treated like a curated primary target
or let a single-assay outlier through.

Run:  python -m pytest tests/test_offtarget.py -q
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from services import offtarget as ot
from services.repurposing_engine import _score_target_overlap, _score_target_overlap_weighted


# Imatinib-like measured profile: (gene, max_pchembl, n_measurements).
IMATINIB = [
    ("ABL1", 9.00, 75), ("KIT", 7.90, 45), ("PDGFRA", 8.70, 12),
    ("ERBB2", 10.22, 1),   # single-assay outlier — must be rejected
    ("DDR1", 9.15, 8), ("CSF1R", 8.00, 7), ("DDR2", 7.82, 5),
    ("NQO2", 7.41, 5), ("SLC47A1", 7.40, 3),
]
PRIMARY = {"ABL1", "KIT", "PDGFRA"}


def test_single_assay_outlier_rejected():
    preds = ot.weight_candidates(IMATINIB, PRIMARY)
    assert "ERBB2" not in {p["gene"] for p in preds}


def test_primary_targets_excluded():
    genes = {p["gene"] for p in ot.weight_candidates(IMATINIB, PRIMARY)}
    assert genes.isdisjoint(PRIMARY)


def test_offtarget_weight_below_primary():
    # No predicted off-target may reach a curated primary's weight (1.0).
    for p in ot.weight_candidates(IMATINIB, PRIMARY):
        assert 0.0 < p["weight"] <= ot.MAX_WEIGHT < 1.0


def test_count_capped():
    assert len(ot.weight_candidates(IMATINIB, PRIMARY, max_offtargets=3)) <= 3


def test_weaker_offtarget_gets_lower_weight():
    preds = {p["gene"]: p["weight"] for p in ot.weight_candidates(IMATINIB, PRIMARY)}
    # DDR1 (9.15, near the ABL1 anchor) must outweigh the weaker DDR2 (7.82).
    assert preds["DDR1"] > preds["DDR2"]


def test_provenance_present():
    for p in ot.weight_candidates(IMATINIB, PRIMARY):
        assert p["provenance"] and str(p["max_pchembl"]) in p["provenance"]


def test_scorer_backward_compatible():
    # drug_weights=None must reproduce the original unweighted score exactly.
    dg, dis = ["ABL1", "KIT"], ["ABL1", "TP53", "MYC"]
    assert _score_target_overlap(dg, dis) == _score_target_overlap(dg, dis, drug_weights=None)


def test_weak_offtarget_cannot_outrank_primary_hit():
    # A drug whose only disease hit is a WEAK off-target must score strictly below
    # a drug whose hit is a genuine primary target — the core promiscuity guard.
    disease = ["FOO", "BAR", "TP53"]
    primary_hit = _score_target_overlap(["TP53"], disease, drug_weights={"TP53": 1.0})
    offtarget_hit = _score_target_overlap(
        ["AAA", "BBB", "TP53"], disease,
        drug_weights={"AAA": 1.0, "BBB": 1.0, "TP53": 0.15})
    assert offtarget_hit < primary_hit


def test_weighted_scorer_backward_compatible():
    dw_dis = {"ABL1": 0.8, "TP53": 0.5}
    a = _score_target_overlap_weighted(["ABL1", "KIT"], dw_dis)
    b = _score_target_overlap_weighted(["ABL1", "KIT"], dw_dis, drug_weights=None)
    assert a == b
