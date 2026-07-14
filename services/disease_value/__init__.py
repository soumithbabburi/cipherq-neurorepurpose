"""
Disease value — runtime lookup for the precomputed Repurposing Value Score.
═══════════════════════════════════════════════════════════════════════════════
Reads data/disease_reference/disease_value.db (built by build_disease_table.py) and
returns, for a disease (by name / EFO / MONDO id), a 0–1 pharma-repurposing value
score with its pillar breakdown — so the reverse screen can rank/filter by "would a
pharma company actually pursue this", not by raw genetic-association strength.

Fail-soft: no DB, or an unmapped disease → None (caller keeps existing behaviour).
"""
from __future__ import annotations

import logging
import re
import sqlite3
from pathlib import Path
from typing import Dict, Optional

logger = logging.getLogger(__name__)

_DB = Path(__file__).parent.parent.parent / "data" / "disease_reference" / "disease_value.db"

_conn = None
_name_index: Dict[str, str] = {}       # normalized label -> mondo_id
_loaded = False


def _tier(v: float) -> str:
    """How a pharma committee would read the 0–1 value score."""
    if v >= 0.60:
        return "High priority"
    if v >= 0.45:
        return "Attractive"
    if v >= 0.30:
        return "Marginal"
    return "Low value"


def _norm(s: str) -> str:
    return re.sub(r"[^a-z0-9]+", " ", (s or "").lower()).strip()


def _load():
    global _conn, _loaded
    if _loaded:
        return
    _loaded = True
    if not _DB.exists():
        logger.info("disease_value: reference table not built — scores unavailable")
        return
    try:
        _conn = sqlite3.connect(f"file:{_DB}?mode=ro", uri=True, check_same_thread=False)
        cur = _conn.cursor()
        cur.execute("SELECT mondo_id, label FROM diseases "
                    "WHERE is_disease=1 AND value_score IS NOT NULL")
        for mid, label in cur.fetchall():
            key = _norm(label)
            if key:
                _name_index.setdefault(key, mid)
        logger.info("disease_value: loaded %d disease scores", len(_name_index))
    except Exception as e:
        logger.warning("disease_value load failed: %s", e)
        _conn = None


def _row(mondo_id: str) -> Optional[Dict]:
    if not _conn:
        return None
    cur = _conn.cursor()
    cur.execute("SELECT mondo_id,label,category,value_score,burden,unmet,market,"
                "approved_drugs,is_orphan,prevalence_class FROM diseases WHERE mondo_id=?",
                (mondo_id,))
    r = cur.fetchone()
    if not r or r[3] is None:
        return None
    return {
        "mondo_id": r[0], "label": r[1], "category": r[2],
        "value_score": r[3], "tier": _tier(r[3]),
        "pillars": {"burden": r[4], "unmet_need": r[5], "market_fit": r[6]},
        "approved_drugs": r[7], "is_orphan": bool(r[8]) if r[8] is not None else None,
        "prevalence_class": r[9],
    }


def value_for(disease_name: str = "", efo_id: str = "") -> Optional[Dict]:
    """Repurposing Value Score for a disease. Resolves by MONDO/EFO id, then name.
    Returns {value_score, tier, pillars, category, ...} or None (fail-soft)."""
    _load()
    if not _conn:
        return None
    # 1. direct id (Open Targets often returns MONDO_xxxx / MONDO:xxxx)
    eid = (efo_id or "").strip().replace("_", ":")
    if eid.upper().startswith("MONDO:"):
        r = _row(eid.upper())
        if r:
            return r
    # 2. normalized name — exact, then loose containment on the index
    n = _norm(disease_name)
    if not n:
        return None
    if n in _name_index:
        return _row(_name_index[n])
    for k, mid in _name_index.items():
        if (n in k or k in n) and abs(len(n) - len(k)) < 8:
            return _row(mid)
    return None
