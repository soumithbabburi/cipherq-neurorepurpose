"""
Global landscape — one consolidated US / EU / India view (data, not link tiles)
═══════════════════════════════════════════════════════════════════════════════
Replaces the row of "go search here" registry buttons with a single fetched,
consolidated table so a global (India-market) client sees approval, generic
availability, IP runway, regulatory pathway and live trial counts per
jurisdiction in one place.

Honesty about sources (this matters in diligence):
  • US     — REAL. FDA Orange Book (patents, exclusivity, ANDA generics) via
             services.acquisition_signal, plus live ClinicalTrials.gov counts.
  • EU     — pathway is factual (Hybrid, Dir 2001/83 Art 10(3)); approval/generic
             status is INFERRED from global approval (EMA has no open API). Each
             row carries an EPAR search link as provenance and an `inferred` flag.
  • India  — pathway is factual (CDSCO new-drug / approved-abroad route); generic
             availability is INFERRED (India is the world's largest generics
             supplier, so a molecule genericised in the US is in practice
             available there). Live India trial count is REAL via ClinicalTrials.gov.

Every cell that is inferred is flagged `inferred: true` so the UI can mark it.
"""
from __future__ import annotations

import logging
import urllib.parse
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)

try:
    from services import http_client
except ImportError:                                   # script context
    import http_client                                # type: ignore

CT_BASE = "https://clinicaltrials.gov/api/v2/studies"

# Essie location expressions for ClinicalTrials.gov country filtering.
_GEO = {
    "US":    "AREA[LocationCountry]United States",
    "EU":    "AREA[LocationCountry](Germany OR France OR Spain OR Italy OR Netherlands "
             "OR Belgium OR Poland OR Sweden OR Austria OR Denmark)",
    "India": "AREA[LocationCountry]India",
}


def _trial_count(drug: str, region: str) -> Optional[int]:
    """Live count of trials for a drug with at least one site in the region."""
    try:
        r = http_client.get(CT_BASE, params={
            "query.term": drug, "filter.advanced": _GEO[region],
            "countTotal": "true", "pageSize": 1, "format": "json"}, timeout=12)
        if r and r.ok:
            return int(r.json().get("totalCount", 0))
    except Exception as e:
        logger.debug("trial count %s/%s failed: %s", drug, region, e)
    return None


def _us_row(drug: str) -> Dict:
    """US: real Orange Book facts via the acquisition-signal engine."""
    row = {"region": "US", "pathway": "505(b)(2) NDA", "authority": "FDA",
           "approved": None, "generics": None, "n_generics": None,
           "ip_until": None, "trials": _trial_count(drug, "US"),
           "inferred": False, "source": "FDA Orange Book + ClinicalTrials.gov",
           "link": f"https://www.accessdata.fda.gov/scripts/cder/daf/index.cfm?event=BasicSearch.process&searchTerm={urllib.parse.quote(drug)}"}
    try:
        from services.acquisition_signal import acquisition_signal
        a = acquisition_signal(drug)
        f = a.get("facts", {}) if a else {}
        if a and a.get("available"):
            row["generics"] = f.get("has_generics")
            row["n_generics"] = f.get("generic_competitors")
            row["ip_until"] = f.get("substance_patent_until") or f.get("latest_any_protection_until")
            # "approved in the US" is implied by having an Orange Book record
            row["approved"] = bool(f.get("originator") or f.get("has_generics"))
            row["note"] = a.get("signal", "")
    except Exception as e:
        logger.debug("US row failed: %s", e)
    return row


def _eu_row(drug: str, approved_mol: bool) -> Dict:
    return {
        "region": "EU", "pathway": "Hybrid (Art 10(3))", "authority": "EMA",
        "approved": (True if approved_mol else None),
        "generics": (True if approved_mol else None),
        "n_generics": None, "ip_until": None,
        "trials": _trial_count(drug, "EU"),
        "inferred": True,
        "note": ("Authorisation/generic status inferred from global approval — verify on the EPAR."
                 if approved_mol else "No global approval on record — status unverified."),
        "source": "Inferred (EMA has no open API) + ClinicalTrials.gov",
        "link": f"https://www.ema.europa.eu/en/medicines?search_api_views_fulltext={urllib.parse.quote(drug)}",
    }


def _india_row(drug: str, us_has_generics: Optional[bool], approved_mol: bool) -> Dict:
    generics = True if (us_has_generics or approved_mol) else None
    return {
        "region": "India", "pathway": "CDSCO New Drug (approved-abroad route)",
        "authority": "CDSCO",
        "approved": (True if approved_mol else None),
        "generics": generics, "n_generics": None, "ip_until": None,
        "trials": _trial_count(drug, "India"),
        "inferred": True,
        "note": ("Generics in practice available — India is the largest global generics supplier; "
                 "a molecule genericised abroad is almost always sourceable locally. CDSCO allows an "
                 "abbreviated route for drugs already approved in well-regulated markets."
                 if generics else "Status unverified — confirm with CDSCO/CTRI."),
        "source": "Inferred (no CDSCO API) + ClinicalTrials.gov (India sites)",
        "link": "https://ctri.nic.in/Clinicaltrials/advsearch.php",
    }


def landscape(drug_name: str, overall_max_phase=0) -> Dict:
    """Consolidated US / EU / India landscape for a molecule."""
    try:
        approved_mol = float(overall_max_phase or 0) >= 4
    except (TypeError, ValueError):
        approved_mol = False

    us = _us_row(drug_name)
    eu = _eu_row(drug_name, approved_mol)
    india = _india_row(drug_name, us.get("generics"), approved_mol)
    return {
        "drug": drug_name,
        "regions": [us, eu, india],
        "note": ("Consolidated from FDA Orange Book (US, authoritative), ClinicalTrials.gov "
                 "(live per-region trial counts), and pathway rules for EU/India. Inferred cells "
                 "are flagged — EU/India lack open regulatory APIs."),
    }


if __name__ == "__main__":
    import json
    for d, ph in [("itraconazole", 4), ("apixaban", 4)]:
        print(f"\n===== {d.upper()} =====")
        L = landscape(d, ph)
        for r in L["regions"]:
            print(f"  {r['region']:<6} pathway={r['pathway']:<34} approved={r['approved']} "
                  f"generics={r['generics']} ip={r['ip_until']} trials={r['trials']} inferred={r['inferred']}")
