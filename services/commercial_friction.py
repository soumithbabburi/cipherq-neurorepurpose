"""
Commercial & Off-Label Friction (Layer 2).
════════════════════════════════════════════════════════════════════════════════
Repurposing platforms are great at biology and blind to development ECONOMICS. Two
frictions a BD / search-and-evaluate team gates on before funding a 505(b)(2):

 1. OFF-LABEL CANNIBALIZATION — if the candidate population is already reachable
    OFF-LABEL under the drug's EXISTING approval (same mechanism, same organ), nobody
    funds a $30M+ dedicated trial; clinicians already prescribe it. The classic case:
    ABPA scored "High Priority" for mepolizumab while doctors already treat it off-label
    under the severe-asthma approval. Lower the score.

 2. MARKET VIABILITY — a very small, NON-orphan market can't support development
    economics; tag "low commercial viability". (A small ORPHAN market is fine — premium
    orphan pricing offsets the small population.)

This is a BOUNDED ranking adjustment + honest flags, NOT a kill: label expansion still
carries reimbursement/coding value, and niche orphan markets can be very profitable.
"""
from __future__ import annotations

import re
import logging
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)

_STOP = {"disease", "disorder", "syndrome", "with", "the", "of", "and", "chronic", "acute"}


def _tok(s: str) -> set:
    return {w for w in re.split(r"[^a-z0-9]+", (s or "").lower()) if len(w) > 3 and w not in _STOP}


def _approved_organs(approved_indications: List[str]) -> set:
    from services.safety_filter import disease_organ_profile
    organs = set()
    for nm in approved_indications or []:
        organs |= set(disease_organ_profile(nm).keys())
    return organs


def _norm01(v) -> Optional[float]:
    try:
        v = float(v)
    except (TypeError, ValueError):
        return None
    return v / 100.0 if v > 1.0 else v


def assess(candidate_disease: str, approved_indications: List[str],
           disease_value: Optional[Dict] = None) -> Dict:
    """Commercial Friction Index for a repurposing candidate. Returns a bounded ranking
    multiplier + off-label / market flags. Fail-soft."""
    from services.safety_filter import disease_organ_profile
    flags: List[str] = []
    off_label = 0.0
    market = 0.0
    market_note = ""

    # 1) Off-label cannibalization — candidate organ overlaps an APPROVED-indication organ.
    cand_organs = set(disease_organ_profile(candidate_disease).keys())
    appr_organs = _approved_organs(approved_indications)
    shared = cand_organs & appr_organs
    if shared:
        ct = _tok(candidate_disease)
        close = ct and any(len(ct & _tok(nm)) / max(1, len(ct)) >= 0.5
                           for nm in (approved_indications or []))
        off_label = 0.65 if close else 0.45
        organ = ", ".join(sorted(shared)).replace("_", "/")
        flags.append(
            f"Off-label reachable — same {organ} organ + mechanism as an already-approved "
            f"indication, so clinicians can prescribe off-label today. This weakens the "
            f"505(b)(2) incentive (label expansion, not a new market).")

    # 2) Market viability — small NON-orphan market is low commercial viability.
    if disease_value:
        prev = (disease_value.get("prevalence_class") or "").lower()
        orphan = disease_value.get("is_orphan")
        mfit = _norm01((disease_value.get("pillars") or {}).get("market_fit"))
        small = (prev in ("ultra-rare", "very rare", "rare")) or (mfit is not None and mfit < 0.30)
        if small and orphan is False:
            market = 0.5
            market_note = ("Small, non-orphan market — low commercial viability without "
                           "orphan pricing to offset the population size.")
            flags.append(market_note)
        elif small and orphan:
            market_note = "Small market but orphan-eligible — premium pricing can make it viable."

    friction = min(0.75, 0.7 * off_label + 0.5 * market)
    multiplier = round(1.0 - friction * 0.5, 3)          # ≥ ~0.625: a down-weight, never a kill
    return {"friction": round(friction, 3), "multiplier": multiplier,
            "off_label_risk": round(off_label, 2), "market_note": market_note,
            "flags": flags, "penalized": friction > 0.05}
