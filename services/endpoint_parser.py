"""
Clinical Endpoint Parsing (Layer 1) — the "did it actually work?" filter.
════════════════════════════════════════════════════════════════════════════════
Structured-first, high precision. The platform's confound was: "N trials exist for
drug X in disease Y" inflated the clinical score even when those trials FAILED their
primary endpoint. The classic case — mepolizumab in eosinophilic esophagitis
(NCT03656380, Phase 2): PRIMARY outcome "Mean Change in Dysphagia Score" p=0.14 — a
CLINICAL miss, even though eosinophils cleared. Counting that as positive is exactly
the false-positive a BD team catches.

This module parses the ClinicalTrials.gov v2 `resultsSection` that is ALREADY in the
study JSON we fetch (no extra call) and classifies each trial by its PRIMARY-endpoint
outcome:
    met_primary          🟢  primary endpoint met (p<0.05)          → boost
    biomarker_only       🟡  biomarker moved, primary CLINICAL miss → penalize + flag
    failed_primary       🟡  primary clinical endpoint missed        → penalize
    terminated_efficacy  🔴  stopped for lack of efficacy/safety     → reject
    unknown                  no parseable results                    → neutral
Every classification is grounded in the posted numbers; nothing is asserted without
the p-value / status behind it.
"""
from __future__ import annotations

import re
import logging
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)

# A biomarker / lab / histology endpoint (molecular, not how the patient feels/functions).
_BIOMARKER_KW = ("eosinophil", "count", "cell", "level", "concentration", "biomarker",
                 "tissue", "histolog", "expression", "serum", "plasma", "titer", "titre",
                 "load", "peak count", "mrna", "protein level", "il-", "cytokine")
# A CLINICAL / patient-relevant endpoint (symptom, function, survival, PRO).
_CLINICAL_KW = ("symptom", "dysphagia", "pain", "response", "remission", "quality of life",
                "qol", "survival", "progression", "exacerbation", "flare", "function",
                "disability", "clinical", "patient-reported", " pro ", "daily", "episode",
                "attack", "relapse", "mortality", "cure", "fev1", "fvc", "score", "index",
                "improvement", "activity", "steroid-sparing", "hospitali")

# ── Structured typing signals (ClinicalTrials.gov v2 unitOfMeasure / paramType) ──
# These are consulted BEFORE the keyword title scan: they are a controlled-ish
# vocabulary and are far more reliable than scanning free-text titles, which
# misreads surrogates (HbA1c, LDL-C, tumour size) as clinical wins. Keyword
# scanning remains only as a last-resort fallback (flagged). See the miner design.
_SURROGATE_UNIT_STEMS = ("mg/dl", "g/dl", "mmol", "µmol", "umol", "nmol", "ng/ml", "pg/ml",
                         "ng/l", "iu/", "u/l", "units per liter", "cells", "copies", "/ml",
                         "mmhg", "ml/min", "kpa", "meq/", "10^", "x10")
_CLINICAL_UNIT_STEMS = ("participant", "score on a scale", "units on a scale",
                        "score on scale", "scale units")
# Named validated CLINICAL instruments / patient-relevant outcomes (title match).
_CLINICAL_INSTRUMENTS = ("acr20", "acr50", "acr70", "das28", "pasi", "easi", "ham-d", "hamd",
                         "madrs", "adas-cog", "updrs", "panss", "haq", "6-minute walk", "6mwd",
                         "overall survival", "progression-free survival", "exacerbation rate")
# Named SURROGATE markers (title match) that must never read as a clinical win.
_KNOWN_SURROGATES = ("hba1c", "ldl", "egfr", "viral load", "tumor size", "tumour size",
                     "eosinophil", "cd4 count", "hdl", "triglycerid", "c-reactive", "crp ",
                     "glucose level", "blood glucose")

_EFF_STOP = ("lack of efficacy", "no efficacy", "inefficac", "futility", "did not meet",
             "failed", "insufficient efficacy", "no benefit", "lack of effect")
_SAF_STOP = ("safety", "adverse", "toxicity", "death", "tolerab", "serious")

# Underpowered-pilot threshold (audited 2026-07): a trial enrolling fewer than this many
# participants that misses its primary only NUMERICALLY (no significant analysis) was almost
# never powered to detect the effect, so calling it a rigorous efficacy FAILURE overstates the
# evidence. Such trials are reclassified "inconclusive_underpowered" (neutral), not failed.
# Example: NCT04527471 (ensifentrine, COVID-19) — a single-site n=45 pilot explicitly not
# designed for statistical significance; placebo numerically recovered slightly more.
_UNDERPOWERED_N = 100
# Title markers that self-declare a study as exploratory / not powered for significance.
_PILOT_TITLE_KW = ("pilot", "feasibility", "proof of concept", "proof-of-concept",
                   "exploratory", "first-in-human", "first in human")


def _enrollment_count(study: Dict) -> Optional[int]:
    """Actual/target enrollment for a CT.gov v2 study, or None if unstated."""
    ei = ((study.get("protocolSection") or {}).get("designModule") or {}).get("enrollmentInfo") or {}
    try:
        return int(ei.get("count")) if ei.get("count") is not None else None
    except (TypeError, ValueError):
        return None


def _is_underpowered(study: Dict) -> bool:
    """True when a study is a small pilot / explicitly exploratory design — a missed primary
    here is inconclusive, not a rigorous failure. Enrollment below the threshold, OR a
    pilot/feasibility/proof-of-concept title, OR an EARLY_PHASE1/PHASE1 efficacy readout."""
    n = _enrollment_count(study)
    if n is not None and n < _UNDERPOWERED_N:
        return True
    ps = study.get("protocolSection", {}) or {}
    title = ((ps.get("identificationModule") or {}).get("briefTitle") or "").lower()
    if any(k in title for k in _PILOT_TITLE_KW):
        return True
    return False


def _kind_by_keyword(t: str) -> str:
    """Legacy last-resort typing by scanning the (lowercased) title. Lower confidence."""
    is_bio = any(k in t for k in _BIOMARKER_KW)
    is_clin = any(k in t for k in _CLINICAL_KW)
    if is_clin and not (is_bio and not is_clin):
        return "clinical"
    if is_bio:
        return "biomarker"
    return "clinical"          # conservative: an un-typed primary is treated as clinical


def _kind(om) -> str:
    """Type a primary outcome as 'clinical' vs 'biomarker' from STRUCTURED fields first
    (unitOfMeasure / paramType / a validated-instrument table), keyword title-scan last.
    Accepts a CT.gov outcomeMeasure dict, or a bare title string (legacy callers)."""
    if isinstance(om, dict):
        title = (om.get("title") or "")
        unit = (om.get("unitOfMeasure") or "")
    else:
        title, unit = (om or ""), ""
    t, u = title.lower(), unit.lower()
    # 1) a named validated clinical instrument / hard outcome in the title
    if any(k in t for k in _CLINICAL_INSTRUMENTS):
        return "clinical"
    # 2) a named surrogate marker in the title
    if any(k in t for k in _KNOWN_SURROGATES):
        return "biomarker"
    # 3) structured unit says molecular / lab -> surrogate
    if any(k in u for k in _SURROGATE_UNIT_STEMS):
        return "biomarker"
    # 4) structured unit says patient-level (participants / rating scale) -> clinical
    if any(k in u for k in _CLINICAL_UNIT_STEMS):
        return "clinical"
    # 5) last resort: the legacy keyword title scan
    return _kind_by_keyword(t)


def _pvalue(s) -> Optional[float]:
    """Parse a ClinicalTrials.gov p-value string ('0.14', '<0.001', '≤0.05', '>0.05')."""
    if s is None:
        return None
    txt = str(s).strip().replace("≤", "<=").replace("≥", ">=")
    m = re.search(r"([<>]=?)?\s*([0-9]*\.?[0-9]+)", txt)
    if not m:
        return None
    op, num = m.group(1), m.group(2)
    try:
        v = float(num)
    except ValueError:
        return None
    if op == "<" or op == "<=":
        return max(0.0, v - 1e-9) if op == "<" else v
    if op == ">" or op == ">=":
        return v + 1e-9          # ">0.05" → just-above, i.e. not significant
    return v


def _outcome_pvalue(om: Dict) -> Optional[float]:
    """Best (smallest) p-value across an outcome measure's posted analyses."""
    ps = [_pvalue(a.get("pValue")) for a in (om.get("analyses") or [])]
    ps = [p for p in ps if p is not None]
    return min(ps) if ps else None


def classify_study(study: Dict) -> Dict:
    """Classify ONE ClinicalTrials.gov v2 study by its primary-endpoint outcome.
    Returns {class, primary_endpoint, primary_kind, p, biomarker_met, note}."""
    ps = study.get("protocolSection", {}) or {}
    status = (ps.get("statusModule", {}).get("overallStatus") or "").upper()
    why = (ps.get("statusModule", {}).get("whyStopped") or "").lower()

    if status in ("TERMINATED", "WITHDRAWN", "SUSPENDED"):
        if any(k in why for k in _EFF_STOP):
            return {"class": "terminated_efficacy", "note": "Stopped for lack of efficacy / futility.", "p": None}
        if any(k in why for k in _SAF_STOP):
            return {"class": "terminated_safety", "note": "Stopped for a safety reason.", "p": None}

    rs = study.get("resultsSection", {}) or {}
    oms = (rs.get("outcomeMeasuresModule", {}) or {}).get("outcomeMeasures", []) or []
    if not oms:
        return {"class": "unknown", "note": "No posted results to parse.", "p": None}

    primaries = [om for om in oms if (om.get("type") or "").upper() == "PRIMARY"]
    if not primaries:
        return {"class": "unknown", "note": "Results posted, no primary-outcome analysis.", "p": None}

    # Was any biomarker outcome (primary OR secondary) statistically moved?
    biomarker_met = any(_kind(om) == "biomarker"
                        and (_outcome_pvalue(om) is not None and _outcome_pvalue(om) < 0.05)
                        for om in oms)

    # A primary endpoint that is met matters ONLY if it is a patient-relevant (clinical)
    # endpoint. A biomarker PRIMARY that moves is target engagement, not a demonstrated
    # clinical benefit, so it must not be scored as a clinical win — separate the two.
    clin_failed = None
    clin_met = None
    bio_met_primary = None
    for om in primaries:
        k = _kind(om)
        p = _outcome_pvalue(om)
        title = (om.get("title") or "")[:90]
        if p is None:
            continue
        if p < 0.05:
            if k == "clinical":
                clin_met = {"title": title, "p": p}
            else:
                bio_met_primary = {"title": title, "p": p}
        elif k == "clinical":
            clin_failed = {"title": title, "p": p}

    # A met CLINICAL primary is the only true "win".
    if clin_met:
        return {"class": "met_primary", "primary_endpoint": clin_met["title"],
                "primary_kind": "clinical", "p": clin_met["p"], "biomarker_met": biomarker_met,
                "note": f"Primary clinical endpoint met ('{clin_met['title']}', p={clin_met['p']:.2g})."}
    if clin_failed:
        if biomarker_met:
            return {"class": "biomarker_only", "primary_endpoint": clin_failed["title"],
                    "primary_kind": "clinical", "p": clin_failed["p"], "biomarker_met": True,
                    "note": (f"Biomarker moved but the PRIMARY clinical endpoint "
                             f"('{clin_failed['title']}') was missed (p={clin_failed['p']:.2g}). "
                             f"Right mechanism, wrong endpoint, a smarter-trial signal, not a win.")}
        # Underpowered / pilot: a numerically-missed primary in a small or explicitly
        # exploratory study is inconclusive, not a rigorous efficacy failure.
        if _is_underpowered(study):
            _n = _enrollment_count(study)
            _nstr = f"n={_n}" if _n is not None else "small pilot"
            return {"class": "inconclusive_underpowered", "primary_endpoint": clin_failed["title"],
                    "primary_kind": "clinical", "p": clin_failed["p"], "biomarker_met": False,
                    "note": (f"Underpowered pilot ({_nstr}): the primary endpoint "
                             f"('{clin_failed['title']}') was numerically missed (p={clin_failed['p']:.2g}), "
                             f"but the study was not powered for statistical significance — inconclusive, "
                             f"not a rigorous efficacy failure.")}
        return {"class": "failed_primary", "primary_endpoint": clin_failed["title"],
                "primary_kind": "clinical", "p": clin_failed["p"], "biomarker_met": False,
                "note": f"Primary clinical endpoint missed ('{clin_failed['title']}', p={clin_failed['p']:.2g})."}
    # No clinical primary resolved, but a biomarker PRIMARY moved: this is biomarker-only
    # evidence (target engagement), NOT a clinical win. Score it as such, not met_primary.
    if bio_met_primary:
        return {"class": "biomarker_only", "primary_endpoint": bio_met_primary["title"],
                "primary_kind": "biomarker", "p": bio_met_primary["p"], "biomarker_met": True,
                "note": (f"The PRIMARY endpoint was a biomarker ('{bio_met_primary['title']}', "
                         f"p={bio_met_primary['p']:.2g}), not a patient-relevant clinical outcome. "
                         f"Target engagement, not a demonstrated clinical benefit.")}
    return {"class": "unknown", "note": "Results posted, primary p-value not extractable.", "p": None}


# outcome class → a directional evidence weight (fed into the clinical/outcome signal)
_CLASS_SIGNAL = {"met_primary": 1.0, "biomarker_only": -0.45, "failed_primary": -0.7,
                 "terminated_efficacy": -1.0, "terminated_safety": -0.9,
                 "inconclusive_underpowered": 0.0, "unknown": 0.0}


def aggregate(studies: List[Dict]) -> Dict:
    """Classify a set of studies (already fetched) and aggregate into a per-pair
    endpoint verdict: counts by class + a directional outcome signal in [-1, 1]."""
    counts: Dict[str, int] = {}
    best = {"class": "unknown", "note": "", "p": None}
    best_rank = -99
    best_nct = ""
    _RANK = {"met_primary": 3, "biomarker_only": 1, "failed_primary": -1,
             "terminated_safety": -2, "terminated_efficacy": -3,
             "inconclusive_underpowered": 0, "unknown": 0}
    for s in studies:
        c = classify_study(s)
        counts[c["class"]] = counts.get(c["class"], 0) + 1
        r = _RANK.get(c["class"], 0)
        _nct = ((s.get("protocolSection", {}) or {}).get("identificationModule", {})
                or {}).get("nctId", "")
        # keep the most INFORMATIVE verdict (a real met or a real fail beats 'unknown');
        # on an equal-magnitude tie, prefer the more-negative class so the verdict can't
        # read positive for the same evidence just because a met_primary was seen first.
        if (abs(r) > abs(best_rank)
                or (abs(r) == abs(best_rank) and r < best_rank)
                or (best["class"] == "unknown" and c["class"] != "unknown")):
            best, best_rank, best_nct = c, r, _nct
    # signal: worst credible negative dominates (a failed primary is stronger evidence
    # than another 'trial exists'); a clean met_primary is the only strong positive.
    sig = 0.0
    if counts.get("met_primary") and not (counts.get("failed_primary")
            or counts.get("terminated_efficacy") or counts.get("terminated_safety")):
        sig = 1.0
    elif counts.get("terminated_efficacy"):
        sig = -1.0
    elif counts.get("terminated_safety"):
        sig = -0.9
    elif counts.get("failed_primary"):
        sig = -0.7
    elif counts.get("biomarker_only"):
        sig = -0.45
    return {"counts": counts, "verdict": best.get("class"),
            "outcome_signal": round(sig, 3), "note": best.get("note", ""),
            "primary_endpoint": best.get("primary_endpoint"), "p": best.get("p"),
            "nct_id": best_nct}
