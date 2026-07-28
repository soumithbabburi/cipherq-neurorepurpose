"""
Coverage-aware Verification stage
════════════════════════════════════════════════════════════════════════════════
Every negative gate in the canonical scorer FAILS OPEN: it starts at a neutral pass
and only moves off it when a source is both queried AND returns a hit. So a pair
whose negative sources were never actually consulted (an obscure drug for a rare
disease with no trials, no FAERS mass, no PrimeKG edge, no genetic weights) reads
identically to a pair whose negative sources were all checked and came back clean.
Absence of evidence is being rewarded as if it were verified-clean evidence.

This module attaches a separate `evidence_balance` object that:
  1. records which negative sources were actually QUERIED vs SILENT for the pair,
  2. summarises support vs contradiction from gate outputs that ALREADY ran, and
  3. yields a verdict that lets the caller demote the TIER / lane of a low-coverage
     candidate.

Hard rules honoured here:
  • It is PURE and DETERMINISTIC. No network I/O, no re-fetch, no re-derivation of any
    penalty. It only RE-READS the return dict of `score_compound_for_disease` and a
    small `coverage_probe` of fetch-time markers the caller supplies.
  • It NEVER writes `composite` / `final_score` and introduces NO seventh multiplier.
    `contradiction_score` merely summarises the six hard gates that already fired.
  • `insufficient_coverage` is an UNKNOWN, not a negative: it must route a genuinely
    sparse candidate to an "unverified" lane, never sink it below junk. A source
    withheld by the own therapy / approval guards is SUPPRESSED (excluded from the
    applicable set), not SILENT, so an approved indication always reads verified.
"""
from __future__ import annotations

import logging
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)

try:
    from services.provenance import stamp as _prov_stamp
except Exception:                                   # pragma: no cover - fail soft
    _prov_stamp = None

# Fixed constants (design section B/G). These are DELIBERATELY fixed, not learned:
# calibration is out of scope, so the coverage threshold is a documented constant, not
# an empirically tuned cutoff. 0.34 = fewer than roughly one of three negative sources
# actually checked.
COVERAGE_THRESHOLD = 0.34
CONTRADICTION_MIN = 0.30
SUPPORT_MIN = 0.30

# Best-matching provenance registry key per source (fail soft: an unknown key still
# produces a usable row from provenance.stamp).
_SOURCE_PROV_KEY = {
    "clinicaltrials": "clinicaltrials",
    "faers": "openfda_faers",
    "primekg": "primekg",
    "target_coverage": "open_targets",
}

_SOURCE_LABEL = {
    "clinicaltrials": "clinical trials",
    "faers": "FAERS adverse events",
    "primekg": "PrimeKG ground truth",
    "target_coverage": "target coverage",
}


def _f(x, default: float = 0.0) -> float:
    try:
        return float(x)
    except (TypeError, ValueError):
        return default


# ── Support summary (display / tie-break only; NEVER multiplies the composite) ────
def _support(result: Dict, probe: Dict) -> float:
    """Positive-evidence mass in [0, 1], read straight from signals already in the
    result. Additive bumps on the mechanistic base, capped at 1. Used only to
    distinguish verified_supported / verified_clean / mixed - it gates no score."""
    base = min(1.0, max(0.0, _f(result.get("mechanistic_plausibility"))))
    s = base
    if (result.get("primekg") or {}).get("relation") == "indication":
        s += 0.5
    if probe.get("endpoint_verdict") == "met_primary":
        s += 0.4
    if _f((result.get("directional_evidence") or {}).get("signal")) > 0:
        s += 0.2
    if (result.get("actionability") or {}).get("clinical_present"):
        s += 0.2
    return round(min(1.0, s), 4)


# ── Contradiction summary (RE-READS already-applied gates; adds NO new penalty) ───
def _contradiction(result: Dict):
    """Max severity of any HARD negative gate that actually FIRED (moved off its
    neutral default), read from the result. This introduces no new penalty - the
    composite was already multiplied by these gates upstream."""
    fired = []   # (human_note, severity)
    tf = result.get("trial_failure") or {}
    if tf.get("penalized"):
        fired.append(("a prior trial failed for efficacy or safety here",
                      1.0 - _f(tf.get("multiplier"), 1.0)))
    d = result.get("direction") or {}
    if d.get("net") == "harmful":
        fired.append(("the mechanism direction would worsen the disease",
                      1.0 - _f(d.get("factor"), 1.0)))
    ap = result.get("appropriateness") or {}
    if ap.get("appropriate") is False:
        fired.append(("the disease is an adverse event the drug causes, or a direction mismatch",
                      1.0 - _f(ap.get("factor"), 1.0)))
    rg = result.get("registry") or {}
    if rg.get("ghost"):
        fired.append(("a late phase program here was abandoned",
                      1.0 - _f(rg.get("multiplier"), 1.0)))
    pk = result.get("primekg") or {}
    if pk.get("relation") == "contraindication":
        fired.append(("PrimeKG labels this an established contraindication", 0.85))
    cp = result.get("ctpa") or {}
    if cp.get("phantom"):
        fired.append(("no functional cohesion beyond a single target match",
                      1.0 - _f(cp.get("multiplier"), 1.0)))
    if not fired:
        return 0.0, []
    score = max(sev for _, sev in fired)
    return round(max(0.0, min(1.0, score)), 4), [n for n, _ in fired]


# ── Per-source coverage ledger (queried vs silent vs suppressed) ──────────────────
def _coverage_ledger(result: Dict, probe: Dict):
    """Build the per-source coverage dict + applicable / checked counts.

    checked   = the source was actually consulted for THIS pair (fetch-time marker or
                self-evident data present).
    suppressed= the source was deliberately withheld by an own therapy / approval guard
                -> NOT applicable (so an approved indication is never punished).
    silent    = applicable, but never run -> lowers coverage_fraction.
    """
    coverage: Dict[str, Dict] = {}
    n_applicable = 0
    n_checked = 0
    silent: List[str] = []

    # ---- ClinicalTrials.gov -------------------------------------------------------
    tc = int(_f(probe.get("trial_count"), 0))
    ctq = bool(probe.get("ctgov_queried"))
    ct_applicable = (tc > 0) or ctq            # only applicable when actually consulted
    ct_checked = (tc > 0) or ctq
    ct_verdict = probe.get("endpoint_verdict") or ("unknown" if ctq else None)
    coverage["clinicaltrials"] = {"checked": bool(ct_checked), "data": bool(tc > 0),
                                  "verdict": ct_verdict}
    if ct_applicable:
        n_applicable += 1
        if ct_checked:
            n_checked += 1
        else:
            silent.append(_SOURCE_LABEL["clinicaltrials"])

    # ---- FAERS / adverse events ---------------------------------------------------
    faers_suppressed = bool(probe.get("faers_suppressed"))
    faers_queried = bool(probe.get("faers_queried"))
    faers_total = probe.get("faers_total")
    if faers_total is None:
        faers_total = (result.get("appropriateness") or {}).get("faers_total")
    faers_checked = faers_queried and (faers_total is not None)
    coverage["faers"] = {"checked": bool(faers_checked),
                         "data": bool((faers_total or 0) > 0),
                         "serious_n": (int(faers_total) if faers_total is not None else None),
                         "suppressed": faers_suppressed}
    if not faers_suppressed:
        n_applicable += 1
        if faers_checked:
            n_checked += 1
        else:
            silent.append(_SOURCE_LABEL["faers"])

    # ---- PrimeKG ground truth -----------------------------------------------------
    pk = result.get("primekg") or {}
    pk_relation = pk.get("relation")
    pk_suppressed = bool(probe.get("primekg_suppressed"))
    pk_checked = pk_relation is not None       # relation None == no edge == not covered
    coverage["primekg"] = {"checked": bool(pk_checked), "data": bool(pk_checked),
                           "relation": pk_relation, "suppressed": pk_suppressed}
    if not pk_suppressed:
        n_applicable += 1
        if pk_checked:
            n_checked += 1
        else:
            silent.append(_SOURCE_LABEL["primekg"])

    # ---- target coverage (Open Targets genetics) ---------------------------------
    cov = result.get("coverage") or {}
    scope = result.get("mechanism_scope") or {}
    target_mediated = bool(scope.get("target_mediated"))
    cov_suppressed = ("suppressed" in cov) or bool(probe.get("coverage_suppressed"))
    cov_val = cov.get("coverage")
    cov_checked = cov_val is not None
    tc_applicable = target_mediated and not cov_suppressed
    coverage["target_coverage"] = {"checked": bool(cov_checked), "data": bool(cov_checked),
                                   "fraction": cov_val, "suppressed": cov_suppressed,
                                   "applicable": tc_applicable}
    if tc_applicable:
        n_applicable += 1
        if cov_checked:
            n_checked += 1
        else:
            silent.append(_SOURCE_LABEL["target_coverage"])

    coverage_fraction = 1.0 if n_applicable == 0 else round(n_checked / n_applicable, 4)
    return coverage, n_applicable, n_checked, coverage_fraction, silent


def _verdict(support: float, contradiction: float, coverage_fraction: float) -> str:
    """First-match-wins verdict order (design section B). insufficient_coverage is an
    UNKNOWN, never a positive; when nothing is applicable coverage_fraction is 1.0 so a
    fully-suppressed approved pair can never read insufficient_coverage."""
    if contradiction >= CONTRADICTION_MIN:
        return "contradicted"
    if coverage_fraction < COVERAGE_THRESHOLD:
        return "insufficient_coverage"
    if support >= SUPPORT_MIN and contradiction > 0:
        return "mixed"
    if support >= SUPPORT_MIN:
        return "verified_supported"
    return "verified_clean"


def _reason(verdict: str, fired: List[str], silent: List[str]) -> str:
    if verdict == "contradicted":
        return "Contradicted: " + ("; ".join(fired) if fired else "a hard negative gate fired.")
    if verdict == "insufficient_coverage":
        if silent:
            return ("Negative sources not yet checked for this pair: "
                    + ", ".join(silent) + ". Hypothesis only.")
        return "Too few negative sources were checked to verify this pair."
    if verdict == "mixed":
        return "Both supporting and contradicting evidence are present."
    if verdict == "verified_supported":
        return "Negative sources checked and clean, with positive supporting evidence."
    return "Negative sources checked and came back clean (mechanism only support)."


def _source_trace(coverage: Dict) -> List[Dict]:
    rows: List[Dict] = []
    for key in ("clinicaltrials", "faers", "primekg", "target_coverage"):
        c = coverage.get(key) or {}
        checked = bool(c.get("checked"))
        if c.get("suppressed"):
            note = "not applicable (withheld by the own therapy / approval guard)"
        elif checked:
            note = "checked for this pair"
        else:
            note = "not consulted for this pair"
        prov_key = _SOURCE_PROV_KEY.get(key, key)
        row = None
        if _prov_stamp:
            try:
                row = _prov_stamp(prov_key, source_key=key, checked=checked, note=note)
            except Exception:
                row = None
        if row is None:
            row = {"source_key": key, "checked": checked, "note": note}
        rows.append(row)
    return rows


def evidence_balance(result: Dict, coverage_probe: Optional[Dict] = None) -> Dict:
    """Build the coverage-aware evidence_balance object (design section B) from the
    already-computed `result` and a small `coverage_probe` of fetch-time markers.

    Deterministic and side-effect free. Fail soft: on any error it returns a neutral
    verified_clean object so a fault here can never demote a candidate."""
    probe = coverage_probe or {}
    try:
        support = _support(result, probe)
        contradiction, fired = _contradiction(result)
        coverage, n_app, n_chk, cov_frac, silent = _coverage_ledger(result, probe)
        verdict = _verdict(support, contradiction, cov_frac)
        return {
            "support_score": support,
            "contradiction_score": contradiction,
            "coverage": coverage,
            "coverage_fraction": cov_frac,
            "n_checked": n_chk,
            "n_applicable": n_app,
            "verdict": verdict,
            "verdict_reason": _reason(verdict, fired, silent),
            "source_trace": _source_trace(coverage),
        }
    except Exception as e:                              # pragma: no cover - fail soft
        logger.debug(f"evidence_balance failed: {e}")
        return {
            "support_score": 0.0, "contradiction_score": 0.0, "coverage": {},
            "coverage_fraction": 1.0, "n_checked": 0, "n_applicable": 0,
            "verdict": "verified_clean",
            "verdict_reason": "Verification unavailable.", "source_trace": [],
        }
