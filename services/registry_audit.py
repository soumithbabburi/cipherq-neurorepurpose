"""
Registry audit (CTPA Rule 2)  —  kill the historical ghosts.
═══════════════════════════════════════════════════════════════════════════════
Temporal-regulatory blindness: the platform credits a drug's highest clinical
phase as positive evidence even when that program DEFINITIVELY FAILED and was
abandoned years ago. Tarenflurbil reached Phase 3 for Alzheimer's, missed its
endpoint, and was permanently discontinued ~2010 — yet a naive "Phase 3" tag
ranked it a top lead. A drug you would recommend a pharma company RE-RUN a failed
human trial on is a false positive of the worst kind.

This audit queries ClinicalTrials.gov for the drug IN the indication and flags a
GHOST when the program is dead: a late phase was reached, there is no active
trial, the newest trial is years old, the drug is not approved for it — or a
late-phase trial was terminated for efficacy/safety. A ghost's score is zeroed.

FAIL-SOFT: no trials / network error → not a ghost (multiplier 1.0). We only zero
on positive evidence of a dead program, never on absence of data.
"""
from __future__ import annotations

import datetime
import logging
import re
from typing import Dict

from services import http_client

logger = logging.getLogger(__name__)

CT_BASE = "https://clinicaltrials.gov/api/v2/studies"
_CURRENT_YEAR = datetime.date.today().year
_STALE_YEARS = 7                       # newest trial older than this = dormant
_ACTIVE = {"RECRUITING", "ACTIVE_NOT_RECRUITING", "ENROLLING_BY_INVITATION",
           "NOT_YET_RECRUITING"}
_STOPPED = {"TERMINATED", "WITHDRAWN", "SUSPENDED"}
_PHASE_NUM = {"PHASE4": 4, "PHASE3": 3, "PHASE2/PHASE3": 3, "PHASE2": 2,
              "PHASE1/PHASE2": 2, "PHASE1": 1, "EARLY_PHASE1": 1}
_EFF_SAF = ("efficac", "futil", "did not meet", "endpoint", "lack of response",
            "no benefit", "ineffective", "safety", "adverse", "toxic")


def _phase_num(phases) -> int:
    return max((_PHASE_NUM.get((p or "").upper(), 0) for p in (phases or [])), default=0)


def _year(struct) -> int:
    d = (struct or {}).get("date", "") if isinstance(struct, dict) else ""
    m = re.match(r"(\d{4})", d or "")
    return int(m.group(1)) if m else 0


def audit(drug: str, disease: str, *, approved_for_it: bool = False) -> Dict:
    """Flag a discontinued 'ghost' program for a drug-disease pair. Returns
    {ghost, multiplier, max_phase, latest_year, active, reason}. Fail-soft."""
    out = {"ghost": False, "multiplier": 1.0, "max_phase": 0,
           "latest_year": 0, "active": False, "reason": ""}
    if not drug or not disease or approved_for_it:
        return out
    try:
        r = http_client.get(CT_BASE, params={
            "query.intr": drug, "query.cond": disease,
            "pageSize": 30, "format": "json"}, timeout=15)
        studies = r.json().get("studies", []) if (r and r.ok) else []
    except Exception as e:
        logger.debug(f"registry audit fetch failed for {drug}/{disease}: {e}")
        return out
    if not studies:
        return out

    max_phase, latest_year, active, failed_late = 0, 0, False, False
    for s in studies:
        ps = s.get("protocolSection", {})
        st = ps.get("statusModule", {})
        status = (st.get("overallStatus") or "").upper()
        ph = _phase_num(ps.get("designModule", {}).get("phases", []))
        max_phase = max(max_phase, ph)
        if status in _ACTIVE:
            active = True
        yr = max(_year(st.get("startDateStruct")),
                 _year(st.get("primaryCompletionDateStruct")),
                 _year(st.get("completionDateStruct")))
        latest_year = max(latest_year, yr)
        why = (st.get("whyStopped") or "").lower()
        if status in _STOPPED and ph >= 3 and any(k in why for k in _EFF_SAF):
            failed_late = True

    out.update(max_phase=max_phase, latest_year=latest_year, active=active)
    dormant = (latest_year and latest_year < _CURRENT_YEAR - _STALE_YEARS)
    # Ghost: a late-phase program that is dead (dormant, no active trials, not
    # approved) OR a late-phase trial explicitly failed for efficacy/safety.
    if failed_late or (max_phase >= 2 and not active and dormant):
        out["ghost"] = True
        out["multiplier"] = 0.05
        out["reason"] = (
            f"Phase {max_phase} program discontinued (last trial {latest_year}, "
            "no active trials, not approved) - recommending it would repeat a "
            "failed/abandoned trial." if not failed_late else
            f"A Phase {max_phase} trial in this indication was terminated for "
            "efficacy/safety - the program failed in humans.")
    return out
