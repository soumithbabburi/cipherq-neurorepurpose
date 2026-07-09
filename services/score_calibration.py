"""
Score calibration (P1)  —  make a raw composite score interpretable.
═══════════════════════════════════════════════════════════════════════════════
A raw 0-1 composite is meaningless without a baseline: a genuinely strong
repurposing lead for a NOVEL indication legitimately scores low in absolute terms
(little direct clinical evidence), so "0.28" reads as weak when it may be a
top-few-percent signal. This module calibrates a score against a NULL distribution
of random drug-disease pairs (data/score_null.json, built offline by
services.build_score_null) and reports:

  • percentile — where the score sits vs the random-pair background (0..1)
  • enrichment — fold over the null median (how far above background)
  • tier       — Strong / Promising / Moderate / Weak, by percentile

so the platform can say "top 3% of all drug-disease pairs — pursue" instead of
showing a weak-looking raw number. FAIL-SOFT: with no (or too small) a null, it
falls back to the documented heuristic bands so nothing breaks.
"""
from __future__ import annotations

import bisect
import json
import logging
from pathlib import Path
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)

_NULL_FILE = Path(__file__).parent.parent / "data" / "score_null.json"
_MIN_NULL = 30          # need at least this many null scores to calibrate

_null: Optional[List[float]] = None
_null_median: float = 0.0
_loaded = False


def _load() -> Optional[List[float]]:
    global _null, _null_median, _loaded
    if _loaded:
        return _null
    _loaded = True
    try:
        if _NULL_FILE.exists():
            data = json.loads(_NULL_FILE.read_text())
            ss = sorted(data.get("sorted") or data.get("scores") or [])
            if len(ss) >= _MIN_NULL:
                _null = ss
                _null_median = ss[len(ss) // 2] or 0.0001
                logger.info("Score calibration: null n=%d median=%.3f", len(ss), _null_median)
    except Exception as e:
        logger.debug("null load failed: %s", e)
    return _null


# Tier from the RAW composite bands. This is the primary tier basis: after the
# 2026-07 scoring rework the composite is well-separated in ABSOLUTE terms (random
# pairs ~0.03, a pathway/novel-target hypothesis ~0.15-0.30, a mechanistically or
# clinically grounded lead ~0.4-0.6, an approved/trial-backed pair ~0.6-0.85), so the
# raw value is directly interpretable. Percentile-vs-null is NOT used for the tier:
# the fix widened the signal/background gap so far that the null collapsed near zero,
# which SATURATES any percentile mapping (0.19 and 0.83 both read "top 98%"). The
# percentile/enrichment are still reported as background context, just not the tier.
def _tier_from_score(score: float) -> str:
    if score >= 0.60:
        return "Strong"
    if score >= 0.40:
        return "Promising"
    if score >= 0.18:
        return "Moderate"
    return "Weak"


def calibrate(score: Optional[float]) -> Dict:
    """Map a raw composite score to {percentile, enrichment, tier, basis}.
    Calibrated against the null when available, else heuristic (basis says which)."""
    if score is None:
        return {"percentile": None, "enrichment": None, "tier": None, "basis": "none"}
    s = float(score)
    tier = _tier_from_score(s)          # tier from the raw composite bands (primary)
    null = _load()
    if null:
        # Percentile + enrichment vs the random-pair background — reported as context
        # (how far above noise), but NOT the tier basis (see _tier_from_score).
        pct = bisect.bisect_right(null, s) / len(null)
        enrichment = round(s / _null_median, 2) if _null_median else None
        return {"percentile": round(pct, 4), "enrichment": enrichment,
                "tier": tier, "basis": "calibrated"}
    return {"percentile": None, "enrichment": None,
            "tier": tier, "basis": "heuristic"}
