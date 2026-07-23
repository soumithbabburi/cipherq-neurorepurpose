"""
DrugPatentWatch connector (config-driven).
═══════════════════════════════════════════════════════════════════════════════
DrugPatentWatch (drugpatentwatch.com) is a PAID, closed enterprise service — it has
no public API and provisions a *custom* REST endpoint per subscribing client. So we
cannot ship live DPW data; this client activates only when credentials are present,
and fails soft to the FDA Orange Book (the authoritative primary source DPW resells)
otherwise.

To activate, set in .env:
    DPW_API_URL   the REST base URL DrugPatentWatch gives your account, e.g.
                  https://<something>.drugpatentwatch.com/api   (must accept a drug
                  name and return JSON patents)
    DPW_API_KEY   your API token
  — OR, if you access DPW through the RapidAPI marketplace listing —
    DPW_RAPIDAPI_KEY    your RapidAPI key
    DPW_RAPIDAPI_HOST   the RapidAPI host (default: drugpatentwatch.p.rapidapi.com)

Response mapping is defensive: DPW's exact field names vary by tier, so we map the
common ones (patent number, expiry, assignee, generic-entry) and pass through the
rest. When you have the account, share one sample response and this maps precisely.
"""
from __future__ import annotations

import logging
import os
from typing import Dict, List

from services import http_client

logger = logging.getLogger(__name__)


def is_configured() -> bool:
    """True only when DPW credentials are present — otherwise the feature is inert."""
    return bool((os.getenv("DPW_API_KEY") and os.getenv("DPW_API_URL"))
                or os.getenv("DPW_RAPIDAPI_KEY"))


def status() -> Dict:
    """Reportable state for the UI / system-status page."""
    if os.getenv("DPW_API_KEY") and os.getenv("DPW_API_URL"):
        return {"configured": True, "mode": "direct"}
    if os.getenv("DPW_RAPIDAPI_KEY"):
        return {"configured": True, "mode": "rapidapi"}
    return {"configured": False, "mode": None,
            "note": "DrugPatentWatch is a paid subscription; set DPW_API_URL + DPW_API_KEY "
                    "(or DPW_RAPIDAPI_KEY) to activate. Falls back to the FDA Orange Book."}


def _normalize(raw: Dict) -> Dict:
    """Map a DPW patent record (field names vary by tier) into our patent dict shape."""
    def pick(*keys):
        for k in keys:
            v = raw.get(k)
            if v not in (None, "", "null"):
                return v
        return ""
    pid = str(pick("patent_number", "patent_no", "patent", "us_patent", "number"))
    expiry = str(pick("expiration_date", "patent_expiry", "expiry", "estimated_expiry",
                      "expiration"))
    return {
        "id": pid,
        "title": str(pick("patent_title", "title", "invention_title")) or (f"Patent {pid}" if pid else "Patent"),
        "assignee": str(pick("applicant", "assignee", "company", "patent_assignee")),
        "expiry": expiry,
        "expiry_iso": expiry[:10] if len(expiry) >= 10 and expiry[4] == "-" else None,
        "generic_entry": str(pick("generic_entry", "earliest_generic", "estimated_generic_entry")),
        "type": str(pick("patent_use_type", "drug_substance_flag", "type")),
        "url": (f"https://patents.google.com/patent/US{pid}/en" if pid[:1].isdigit()
                else (f"https://patents.google.com/patent/{pid}/en" if pid else "")),
        "source": "DrugPatentWatch",
        "authoritative": True,
        "status": str(pick("patent_status", "status")) or "unknown",
    }


def patents_for_drug(drug: str, timeout: float = 12.0) -> Dict:
    """Fetch DrugPatentWatch patents for a drug. {available, configured, patents[], note}.

    Inert (available=False) unless credentials are configured — the caller then uses
    the Orange Book. Never raises: a bad response degrades to available=False.
    """
    drug = (drug or "").strip()
    if not drug:
        return {"available": False, "configured": is_configured(), "patents": []}
    if not is_configured():
        return {"available": False, "configured": False, "patents": [],
                "note": "DrugPatentWatch not configured (paid subscription required)."}

    try:
        if os.getenv("DPW_API_KEY") and os.getenv("DPW_API_URL"):
            base = os.getenv("DPW_API_URL").rstrip("/")
            data = http_client.get_json(
                base, timeout=timeout, default=None,
                params={"drug": drug, "search": drug},
                headers={"Authorization": f"Bearer {os.getenv('DPW_API_KEY')}",
                         "X-Api-Key": os.getenv("DPW_API_KEY")})
        else:
            host = os.getenv("DPW_RAPIDAPI_HOST", "drugpatentwatch.p.rapidapi.com")
            data = http_client.get_json(
                f"https://{host}/", timeout=timeout, default=None,
                params={"drug": drug, "search": drug},
                headers={"X-RapidAPI-Key": os.getenv("DPW_RAPIDAPI_KEY"),
                         "X-RapidAPI-Host": host})
    except Exception as e:
        logger.warning("DrugPatentWatch fetch failed: %s", e)
        return {"available": False, "configured": True, "patents": [],
                "note": "DrugPatentWatch request failed — using Orange Book."}

    if not data:
        return {"available": False, "configured": True, "patents": [],
                "note": "DrugPatentWatch returned no data — using Orange Book."}

    # Accept either a list of patents or an object wrapping them.
    rows = data if isinstance(data, list) else (
        data.get("patents") or data.get("results") or data.get("data") or [])
    patents = [_normalize(r) for r in rows if isinstance(r, dict)]
    patents = [p for p in patents if p["id"]]
    return {"available": bool(patents), "configured": True, "patents": patents,
            "note": "" if patents else "DrugPatentWatch returned an unrecognized shape "
                    "— share a sample response to map the fields."}
