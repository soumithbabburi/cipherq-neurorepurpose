"""
MONDO hierarchy — the platform's disease-ontology backbone.
═══════════════════════════════════════════════════════════════════════════════
Reads the MONDO is_a graph (built into disease_value.db by build_mondo_hierarchy.py)
+ the MONDO id ↔ label ↔ MeSH xrefs already in that DB. Gives a complete, orphan-free
disease hierarchy (MONDO links EGPA under vasculitis, which MeSH's orphan row does
not). Used by services/disease_normalize.py as the primary subtype-expansion source,
with MeSH as the fallback when a term isn't in MONDO.

Fail-soft: no DB / not built → available()==False, callers fall back to MeSH.
"""
from __future__ import annotations

import logging
import re
import sqlite3
from functools import lru_cache
from pathlib import Path
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)

_DB = Path(__file__).parent.parent / "data" / "disease_reference" / "disease_value.db"
_conn: Optional[sqlite3.Connection] = None
_MIN_CHILDREN_BROAD = 4


def _c() -> Optional[sqlite3.Connection]:
    global _conn
    if _conn is None:
        try:
            _conn = sqlite3.connect(f"file:{_DB}?mode=ro", uri=True, check_same_thread=False)
        except Exception as e:
            logger.warning("mondo_hierarchy: cannot open %s: %s", _DB, e)
            _conn = None
    return _conn


@lru_cache(maxsize=1)
def available() -> bool:
    c = _c()
    if not c:
        return False
    try:
        c.execute("SELECT 1 FROM mondo_edges LIMIT 1")
        return True
    except Exception:
        return False


def _norm(s: str) -> str:
    return re.sub(r"\s+", " ", (s or "").strip().lower())


@lru_cache(maxsize=4096)
def resolve(name: str, mesh_id: str = "") -> Optional[str]:
    """Resolve a disease name (and/or a MeSH id) to a MONDO id."""
    c = _c()
    if not c:
        return None
    # MeSH xref is the most reliable bridge (we already resolve names → MeSH elsewhere)
    if mesh_id:
        key = mesh_id if mesh_id.startswith("MESH:") else f"MESH:{mesh_id}"
        r = c.execute("SELECT mondo_id FROM diseases WHERE mesh = ? LIMIT 1", (key,)).fetchone()
        if r:
            return r[0]
    n = _norm(name)
    if not n:
        return None
    r = c.execute("SELECT mondo_id FROM diseases WHERE LOWER(label) = ? LIMIT 1", (n,)).fetchone()
    if r:
        return r[0]
    # exact synonym (synonyms stored as " | "-joined text)
    r = c.execute(
        "SELECT mondo_id FROM diseases WHERE LOWER(synonyms) LIKE ? "
        "ORDER BY LENGTH(label) LIMIT 1", (f"%{n}%",)).fetchone()
    return r[0] if r else None


def label(mondo_id: str) -> str:
    c = _c()
    if not c or not mondo_id:
        return ""
    r = c.execute("SELECT COALESCE(d.label, l.label) FROM mondo_labels l "
                  "LEFT JOIN diseases d ON d.mondo_id = l.mondo_id WHERE l.mondo_id = ?",
                  (mondo_id,)).fetchone()
    if r and r[0]:
        return r[0]
    r = c.execute("SELECT label FROM diseases WHERE mondo_id = ?", (mondo_id,)).fetchone()
    return (r[0] if r else "") or ""


def child_count(mondo_id: str) -> int:
    c = _c()
    if not c or not mondo_id:
        return 0
    r = c.execute("SELECT COUNT(*) FROM mondo_edges WHERE parent = ?", (mondo_id,)).fetchone()
    return int(r[0]) if r else 0


def subtypes(mondo_id: str, limit: int = 12, max_depth: int = 10) -> List[Dict]:
    """Clinically-distinct subtypes: disease descendants that are NOT themselves broad
    categories (child count < threshold), ranked by Repurposing Value. Drops the
    intermediate umbrella nodes (e.g. 'systemic vasculitis', 'ANCA-associated vasculitis')
    and keeps the actionable indications (EGPA, GPA, MPA, Takayasu, GCA, Behçet, …)."""
    c = _c()
    if not c or not mondo_id:
        return []
    q = """
    WITH RECURSIVE sub(id, depth) AS (
        SELECT child, 1 FROM mondo_edges WHERE parent = ?
        UNION
        SELECT e.child, s.depth + 1 FROM mondo_edges e JOIN sub s ON e.parent = s.id
        WHERE s.depth < ?
    )
    SELECT sub.id,
           COALESCE(d.label, l.label)                                        AS label,
           d.mesh                                                            AS mesh,
           COALESCE(d.value_score, 0)                                        AS vs,
           (SELECT COUNT(*) FROM mondo_edges ch WHERE ch.parent = sub.id)    AS nkids,
           COALESCE(d.is_disease, 1)                                         AS is_dis
    FROM sub
    LEFT JOIN diseases d      ON d.mondo_id = sub.id
    LEFT JOIN mondo_labels l  ON l.mondo_id = sub.id
    """
    try:
        rows = c.execute(q, (mondo_id, max_depth)).fetchall()
    except Exception as e:
        logger.debug("mondo subtypes query failed: %s", e)
        return []
    seen, out = set(), []
    for mid, lbl, mesh, vs, nkids, is_dis in rows:
        lbl = (lbl or "").strip()
        key = lbl.lower()
        if not lbl or key in seen or lbl.lower().startswith("obsolete"):
            continue
        if not is_dis or nkids >= _MIN_CHILDREN_BROAD:   # skip non-diseases + still-broad nodes
            continue
        seen.add(key)
        out.append({"mondo_id": mid, "label": lbl,
                    "mesh_id": (mesh or "").replace("MESH:", ""),
                    "value_score": round(float(vs or 0.0), 3)})
    out.sort(key=lambda x: (-x["value_score"], x["label"]))
    return out[:limit]


def is_broad(mondo_id: str) -> bool:
    return child_count(mondo_id) >= _MIN_CHILDREN_BROAD


if __name__ == "__main__":
    print("MONDO hierarchy available:", available())
    vid = resolve("vasculitis")
    print("vasculitis →", vid, "| broad:", is_broad(vid))
    for s in subtypes(vid, 12):
        print(f"   {s['value_score']:.2f}  {s['label']}")
