"""
Portfolio Matching Service

Loads a pharma company's generic / ANDA molecule list (data/portfolio.json) and
matches it against repurposing candidates. A candidate that the company already
manufactures is the cleanest 505(b)(2) play — they can cite their own existing
safety / CMC data for the new indication.

Nothing about the portfolio is hardcoded in source: the molecule list lives in
data/portfolio.json (edit it freely). Each molecule name is resolved to a ChEMBL
ID dynamically via the ChEMBL REST API and cached in data/portfolio_resolved.json
so resolution happens once, not on every request.
"""

import json
import logging
import re
import time
from pathlib import Path
from typing import Dict, List, Optional

import requests  # noqa: F401  (kept for type/back-compat)

from services import http_client

logger = logging.getLogger(__name__)

DATA_DIR        = Path(__file__).parent.parent / "data"
PORTFOLIO_FILE  = DATA_DIR / "portfolio.json"
RESOLVED_FILE   = DATA_DIR / "portfolio_resolved.json"
CHEMBL_BASE     = "https://www.ebi.ac.uk/chembl/api/data"

# Boost added to a candidate's composite score when the molecule is already in
# the company portfolio (marketed). Pipeline molecules get a smaller boost.
PORTFOLIO_BOOST_MARKETED = 0.10
PORTFOLIO_BOOST_PIPELINE = 0.05


def _norm(name: str) -> str:
    """Normalise a molecule name for matching: lowercase, strip salts/forms."""
    n = (name or "").lower().strip()
    # Drop common salt / form qualifiers so "Losartan Potassium" == "losartan"
    n = re.sub(
        r"\b(hydrochloride|hcl|sodium|potassium|calcium|maleate|besylate|besilate|"
        r"fumarate|tartrate|sulfate|sulphate|mesylate|succinate|valerate|propionate|"
        r"etabonate|tromethamine|medoxomil|cilexetil|etexilate|acetonide|"
        r"dihydrate|monohydrate|anhydrous)\b",
        "", n,
    )
    n = re.sub(r"[^a-z0-9 ]+", " ", n)
    return re.sub(r"\s+", " ", n).strip()


# ── Load portfolio definition ──────────────────────────────────────────────────

def load_portfolio() -> dict:
    """Read the raw portfolio.json. Returns {} if absent (feature disables cleanly)."""
    try:
        if PORTFOLIO_FILE.exists():
            return json.loads(PORTFOLIO_FILE.read_text(encoding="utf-8"))
    except Exception as e:
        logger.warning(f"Could not read portfolio.json: {e}")
    return {}


def is_available() -> bool:
    return bool(load_portfolio().get("products"))


# ── ChEMBL resolution (dynamic, cached) ─────────────────────────────────────────

def _load_resolved() -> dict:
    try:
        if RESOLVED_FILE.exists():
            return json.loads(RESOLVED_FILE.read_text(encoding="utf-8"))
    except Exception:
        pass
    return {}


def _save_resolved(data: dict):
    try:
        DATA_DIR.mkdir(exist_ok=True)
        RESOLVED_FILE.write_text(json.dumps(data, indent=2), encoding="utf-8")
    except Exception as e:
        logger.debug(f"Could not write portfolio_resolved.json: {e}")


def _chembl_lookup(name: str) -> Optional[Dict]:
    """Resolve a molecule name → {chembl_id, pref_name, max_phase} via ChEMBL."""
    try:
        r = http_client.get(
            f"{CHEMBL_BASE}/molecule/search.json",
            params={"q": name, "limit": 1}, timeout=10,
        )
        if r and r.ok:
            mols = r.json().get("molecules", [])
            if mols:
                m = mols[0]
                return {
                    "chembl_id": m.get("molecule_chembl_id", ""),
                    "pref_name": m.get("pref_name") or name,
                    "max_phase": m.get("max_phase") or 0,
                }
    except Exception as e:
        logger.debug(f"ChEMBL lookup failed for '{name}': {e}")
    return None


def resolve_portfolio(force: bool = False) -> Dict[str, dict]:
    """
    Resolve every portfolio molecule to a ChEMBL ID (cached).
    Returns an index keyed by normalised molecule name:
        { norm_name: {name, chembl_id, pref_name, status, area, form, ...} }
    """
    portfolio = load_portfolio()
    products  = portfolio.get("products", [])
    if not products:
        return {}

    resolved = {} if force else _load_resolved()
    index: Dict[str, dict] = {}
    dirty = False

    for prod in products:
        name = prod.get("name", "")
        key  = _norm(name)
        if not key:
            continue

        cached = resolved.get(key)
        if cached is None:
            looked = _chembl_lookup(name)
            cached = looked or {"chembl_id": "", "pref_name": name, "max_phase": 0}
            resolved[key] = cached
            dirty = True
            time.sleep(0.1)  # be polite to the ChEMBL API

        index[key] = {
            "name":      name,
            "chembl_id": cached.get("chembl_id", ""),
            "pref_name": cached.get("pref_name", name),
            "max_phase": cached.get("max_phase", 0),
            "status":    prod.get("status", "marketed"),
            "area":      prod.get("area", ""),
            "form":      prod.get("form", ""),
            "note":      prod.get("note", ""),
        }

    if dirty:
        _save_resolved(resolved)

    return index


# ── Matching ─────────────────────────────────────────────────────────────────

def _build_lookups(index: Dict[str, dict]):
    by_name   = index  # normalised name → entry
    by_chembl = {e["chembl_id"]: e for e in index.values() if e.get("chembl_id")}
    return by_name, by_chembl


def match(chembl_id: str = "", name: str = "") -> Optional[dict]:
    """Return the portfolio entry for a drug by ChEMBL ID or name, else None."""
    index = resolve_portfolio()
    if not index:
        return None
    by_name, by_chembl = _build_lookups(index)
    if chembl_id and chembl_id in by_chembl:
        return by_chembl[chembl_id]
    if name:
        return by_name.get(_norm(name))
    return None


def annotate_candidates(candidates: List[dict]) -> List[dict]:
    """
    Tag each candidate with portfolio membership and apply a 505(b)(2) boost.
    Each candidate dict may carry 'chembl_id'/'name' and 'composite_score'/'score'.
    Re-sorts in place by the boosted score. Safe no-op if no portfolio is loaded.
    """
    index = resolve_portfolio()
    if not index or not candidates:
        return candidates

    by_name, by_chembl = _build_lookups(index)

    for c in candidates:
        entry = None
        cid = c.get("chembl_id", "")
        if cid and cid in by_chembl:
            entry = by_chembl[cid]
        if entry is None:
            entry = by_name.get(_norm(c.get("name", "")))

        if entry:
            status = entry.get("status", "marketed")
            boost  = PORTFOLIO_BOOST_MARKETED if status == "marketed" else PORTFOLIO_BOOST_PIPELINE
            c["in_portfolio"]    = True
            c["portfolio_status"] = status
            c["portfolio_form"]   = entry.get("form", "")
            c["portfolio_note"]   = entry.get("note", "")
            c["regulatory_ready"] = "505(b)(2)-ready — molecule already in Alembic portfolio"
            base = c.get("composite_score", c.get("score", 0.0)) or 0.0
            c["composite_score"]  = round(min(1.0, base + boost), 4)
            c["score"]            = c["composite_score"]
            c["portfolio_boost"]  = boost
        else:
            c.setdefault("in_portfolio", False)

    candidates.sort(key=lambda x: x.get("composite_score", x.get("score", 0)), reverse=True)
    return candidates
