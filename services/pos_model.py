"""
Probability-of-Success (PoS) / phase-progression model
═══════════════════════════════════════════════════════════════════════════════
Turns "we display the current clinical phase" into "we *predict* the likelihood
of a molecule progressing through each remaining development phase to approval"
— the claim the platform makes for Asset Scoring.

This is a phase-transition SURVIVAL model, not a black box. Every number is
traceable to the published clinical-development literature, then conditioned on
two things the platform actually knows about a *repurposing* candidate:

  1. It is an already-approved molecule  → Phase-1 safety/PK is largely cleared,
     which is precisely why repurposing de-risks development (the flyer's
     "70% vs 90% failure" is this effect). We raise the Phase 1→2 transition
     accordingly and document it.
  2. Our evidence/relevance score for the NEW indication is an estimate of
     efficacy odds there — so it modulates the efficacy-determining Phase 2→3
     transition. Weak mechanistic/clinical evidence lowers it; strong evidence
     raises it (bounded).

An optional fitted scikit-learn calibrator (data/pos_model.pkl, produced by
services/train_pos_model.py from real ChEMBL outcomes) is blended in when
present. With no model file and no sklearn, the analytic model stands alone and
still returns real, citable numbers — so the scoring path never hard-fails.

Sources
  • Thomas D.W. et al., "Clinical Development Success Rates 2006-2015",
    BIO / Biomedtracker / Amplion, 2016.   (overall + therapeutic-area LoA)
  • Wong C.H., Siah K.W., Lo A.W., "Estimation of clinical trial success rates
    and related parameters", Biostatistics 20(2):273-286, 2019.
  • Hay M. et al., "Clinical development success rates for investigational
    drugs", Nat. Biotechnol. 32:40-51, 2014.
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)

_ROOT = Path(__file__).parent.parent
_MODEL_FILE = _ROOT / "data" / "pos_model.pkl"

# ── Industry phase-transition success rates (BIO 2016, all indications) ────────
# Probability of advancing from each stage to the next.
P_1_2   = 0.632     # Phase 1 → Phase 2
P_2_3   = 0.307     # Phase 2 → Phase 3   (efficacy gate — the usual graveyard)
P_3_NDA = 0.581     # Phase 3 → regulatory filing
P_NDA_A = 0.853     # filing → approval
# Sanity: 0.632*0.307*0.581*0.853 = 0.0961  → matches BIO's 9.6% Phase-1 LoA.
_OVERALL_LOA = P_1_2 * P_2_3 * P_3_NDA * P_NDA_A

# For an already-approved molecule entering a NEW indication, dedicated Phase-1
# safety/PK is largely cleared. Published repurposing analyses put de-risking at
# ~2-3x the de-novo LoA; we encode it conservatively at the safety transition.
P_1_2_REPURPOSE = 0.90

# Therapeutic-area Phase-1 Likelihood-of-Approval (BIO 2016). Used as a relative
# multiplier vs the all-indication overall; missing areas fall back to overall.
AREA_LOA: Dict[str, float] = {
    "hematology":           0.261,
    "infectious disease":   0.191,
    "ophthalmology":        0.171,
    "metabolic":            0.153,
    "endocrine":            0.153,
    "autoimmune":           0.151,
    "immunology":           0.151,
    "respiratory":          0.131,
    "gastroenterology":     0.122,
    "allergy":              0.122,
    "cardiovascular":       0.087,
    "neurology":            0.084,
    "psychiatry":           0.083,
    "urology":              0.082,
    "oncology":             0.051,
}

_STAGES = ["Phase 1", "Phase 2", "Phase 3", "Filing", "Approved"]


def _clamp(x: float, lo: float, hi: float) -> float:
    return max(lo, min(hi, x))


def _area_multiplier(area: Optional[str]) -> float:
    if not area:
        return 1.0
    key = area.strip().lower()
    loa = None
    for k, v in AREA_LOA.items():
        if k in key or key in k:
            loa = v
            break
    if loa is None:
        return 1.0
    # Relative to the all-indication baseline, gently bounded.
    return _clamp(loa / _OVERALL_LOA, 0.35, 2.4)


# ── Optional fitted calibrator (lazy, fail-soft) ───────────────────────────────
_ml_cache: Optional[Dict] = None
_ml_loaded = False


def _load_ml() -> Optional[Dict]:
    """Load the documented sklearn calibrator if it exists. Never raises."""
    global _ml_cache, _ml_loaded
    if _ml_loaded:
        return _ml_cache
    _ml_loaded = True
    try:
        if _MODEL_FILE.exists():
            import joblib
            bundle = joblib.load(_MODEL_FILE)
            ok = (isinstance(bundle, dict) and bundle.get("model") is not None
                  and bundle.get("features"))
            # Only blend a calibrator the training step judged non-circular and
            # genuinely informative (see train_pos_model.py). A model that scores
            # near-perfect on this label is leaking (annotation richness is a
            # CONSEQUENCE of approval, not a predictor) and must NOT touch a
            # customer-facing score — the analytic model stands alone instead.
            trustworthy = bool(bundle.get("meta", {}).get("trustworthy")) if ok else False
            if ok and trustworthy:
                _ml_cache = bundle
                logger.info("PoS calibrator loaded (%s features, cv_auc=%s)",
                            len(bundle["features"]), bundle["meta"].get("cv_auc"))
            elif ok:
                logger.info("PoS calibrator present but flagged non-trustworthy "
                            "(cv_auc=%s) — using analytic model only.",
                            bundle.get("meta", {}).get("cv_auc"))
    except Exception as e:                       # missing sklearn/joblib, schema drift…
        logger.debug("PoS calibrator unavailable: %s", e)
        _ml_cache = None
    return _ml_cache


def _ml_probability(features: Dict[str, float]) -> Optional[float]:
    bundle = _load_ml()
    if not bundle:
        return None
    try:
        x = [[float(features.get(f, 0.0)) for f in bundle["features"]]]
        return float(bundle["model"].predict_proba(x)[0][1])
    except Exception as e:
        logger.debug("PoS calibrator predict failed: %s", e)
        return None


# ── Public API ─────────────────────────────────────────────────────────────────
def predict_progression(
    current_phase: float = 0.0,
    evidence_score: Optional[float] = None,
    therapeutic_area: Optional[str] = None,
    is_repurposing: bool = True,
    features: Optional[Dict[str, float]] = None,
) -> Dict:
    """
    Predict the likelihood of clearing each remaining development phase.

    current_phase   : platform max_phase (0 preclinical … 4 approved in original use)
    evidence_score  : 0-1 relevance/evidence for the NEW indication (efficacy proxy)
    therapeutic_area: e.g. "neurology" — selects published area LoA when known
    is_repurposing  : if True, treat Phase-1 safety as de-risked (known molecule)
    features        : optional feature dict for the fitted calibrator blend

    Returns per-transition probabilities, the cumulative Likelihood of Approval
    from the current phase, the dominant residual risk, and full provenance.
    """
    phase = current_phase or 0.0
    ev = None if evidence_score is None else _clamp(float(evidence_score), 0.0, 1.0)

    # A repurposing candidate is entering a NEW indication: its efficacy there is
    # unproven no matter how far it advanced in its ORIGINAL use, so the
    # new-indication program starts at Phase 1 (safety de-risked via
    # P_1_2_REPURPOSE below). The phase shortcut therefore applies ONLY to a
    # de-novo program, where current_phase is the molecule's progress in THIS
    # indication — otherwise an approved (Phase-4 original-use) molecule would be
    # scored as if it had already cleared Phase 1-3 of the new indication.
    if is_repurposing:
        start_idx = 0                  # new indication starts fresh at Phase 1
    elif phase >= 3:
        start_idx = 2                  # in/through Phase 3
    elif phase == 2:
        start_idx = 1                  # in Phase 2
    else:
        start_idx = 0                  # preclinical / Phase 1

    p12 = P_1_2_REPURPOSE if is_repurposing else P_1_2
    # Evidence modulates the efficacy gate (Phase 2 → 3): the single transition our
    # mechanistic + clinical-trial evidence actually informs for a new indication.
    if ev is None:
        p23 = P_2_3
    else:
        p23 = _clamp(P_2_3 * (0.55 + 0.90 * ev), 0.05, 0.85)

    ladder = [
        ("Phase 1", "Phase 2", p12),
        ("Phase 2", "Phase 3", p23),
        ("Phase 3", "Filing",  P_3_NDA),
        ("Filing",  "Approved", P_NDA_A),
    ]
    remaining = ladder[start_idx:]

    cumulative = 1.0
    transitions: List[Dict] = []
    for frm, to, p in remaining:
        cumulative *= p
        transitions.append({"from": frm, "to": to, "p": round(p, 4)})

    cumulative *= _area_multiplier(therapeutic_area)
    analytic_loa = _clamp(cumulative, 0.0, 0.95)

    # Blend with fitted calibrator when both the model and features are available.
    model_kind = "phase-transition survival (analytic)"
    loa = analytic_loa
    ml_p = _ml_probability(features) if features else None
    if ml_p is not None:
        loa = _clamp(0.5 * analytic_loa + 0.5 * ml_p, 0.0, 0.95)
        model_kind = "phase-transition survival + fitted calibrator"

    dominant = min(transitions, key=lambda t: t["p"]) if transitions else None
    band = "High" if loa >= 0.30 else "Moderate" if loa >= 0.12 else "Low"

    return {
        "loa": round(loa, 4),
        "loa_pct": round(loa * 100, 1),
        "from_phase": _STAGES[start_idx],
        "transitions": transitions,
        "dominant_risk": (f'{dominant["from"]} to {dominant["to"]}' if dominant else None),
        "band": band,
        "is_repurposing": is_repurposing,
        "evidence_score": ev,
        "therapeutic_area": therapeutic_area,
        "model": model_kind,
        "baseline_de_novo_loa": round(_OVERALL_LOA, 4),
        "citations": [
            "BIO/Biomedtracker 2016 (Thomas et al.)",
            "Wong, Siah & Lo, Biostatistics 2019",
            "Hay et al., Nat. Biotechnol. 2014",
        ],
    }


def progression_curve(current_phase: float = 0.0, **kw) -> List[Dict]:
    """Cumulative survival across stages — convenient for a small chart later."""
    res = predict_progression(current_phase, **kw)
    pts, cum = [{"stage": res["from_phase"], "survival": 1.0}], 1.0
    for t in res["transitions"]:
        cum *= t["p"]
        pts.append({"stage": t["to"], "survival": round(cum, 4)})
    return pts


if __name__ == "__main__":
    import json
    print("De-novo baseline Phase-1 LoA:", round(_OVERALL_LOA, 4), "\n")
    cases = [
        ("Approved drug, strong evidence, neurology",
         dict(current_phase=4, evidence_score=0.82, therapeutic_area="neurology")),
        ("Approved drug, weak evidence, oncology",
         dict(current_phase=4, evidence_score=0.20, therapeutic_area="oncology")),
        ("In Phase 2, moderate evidence, metabolic",
         dict(current_phase=2, evidence_score=0.55, therapeutic_area="metabolic",
              is_repurposing=False)),
        ("De-novo Phase 1 (no repurposing), avg evidence",
         dict(current_phase=1, evidence_score=0.5, is_repurposing=False)),
    ]
    for label, kw in cases:
        r = predict_progression(**kw)
        print(f"• {label}")
        print(f"    LoA from {r['from_phase']}: {r['loa_pct']}%  ({r['band']}), "
              f"dominant risk {r['dominant_risk']}  [{r['model']}]")
    print("\nProgression curve (Phase 2, ev=0.55, metabolic):")
    print(json.dumps(progression_curve(2, evidence_score=0.55, therapeutic_area="metabolic"), indent=2))
