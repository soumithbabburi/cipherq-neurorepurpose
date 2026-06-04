"""
Disease Resolver
Accepts any free-text disease query and resolves it to a canonical MeSH ID
plus all related IDs (parents, children, siblings) using the mesh_diseases table.
Falls back to LIKE-based matching if MeSH lookup fails.
"""

import logging
import os
import re
from functools import lru_cache
from typing import Dict, List, Optional, Tuple

import psycopg2
import psycopg2.extras
import psycopg2.pool

logger = logging.getLogger(__name__)

_DB_PARAMS = dict(
    host=os.getenv("DB_HOST", os.getenv("CHEMBL_DB_HOST", "localhost")),
    port=int(os.getenv("DB_PORT", os.getenv("CHEMBL_DB_PORT", "5433"))),
    user=os.getenv("DB_USER", os.getenv("CHEMBL_DB_USER", "babburisoumith")),
    password=os.getenv("DB_PASSWORD", os.getenv("CHEMBL_DB_PASSWORD", "")),
    dbname=os.getenv("DB_NAME", "neurorepurpose"),
    connect_timeout=3,
)

_pool: Optional[psycopg2.pool.ThreadedConnectionPool] = None
_pool_failed = False


def _get_pool():
    global _pool, _pool_failed
    if _pool is None:
        if _pool_failed:
            return None
        try:
            _pool = psycopg2.pool.ThreadedConnectionPool(1, 5, **_DB_PARAMS)
        except Exception as e:
            _pool_failed = True
            logger.warning(f"DB pool failed: {e}")
    return _pool


def _q(sql: str, params=None) -> List[Dict]:
    pool = _get_pool()
    if pool is None:
        return []
    conn = pool.getconn()
    try:
        with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
            cur.execute(sql, params)
            return [dict(r) for r in cur.fetchall()]
    except Exception as e:
        logger.error(f"Query failed: {e}")
        return []
    finally:
        pool.putconn(conn)


# Common abbreviations and alternate names → canonical form
_ALIASES = {
    "ad": "alzheimer disease",
    "alzheimers": "alzheimer disease",
    "alzheimer's": "alzheimer disease",
    "pd": "parkinson disease",
    "parkinsons": "parkinson disease",
    "parkinson's": "parkinson disease",
    "ms": "multiple sclerosis",
    "als": "amyotrophic lateral sclerosis",
    "lou gehrig": "amyotrophic lateral sclerosis",
    "mnd": "motor neuron disease",
    "hd": "huntington disease",
    "huntingtons": "huntington disease",
    "lbd": "lewy body disease",
    "ftd": "frontotemporal dementia",
    "mci": "mild cognitive impairment",
    "bpd": "bipolar disorder",
    "ocd": "obsessive-compulsive disorder",
    "ptsd": "stress disorders, post-traumatic",
    "adhd": "attention deficit disorder with hyperactivity",
    "asd": "autism spectrum disorder",
    "tbi": "brain injuries, traumatic",
    "cte": "chronic traumatic encephalopathy",
}


def _normalize(query: str) -> str:
    q = query.strip().lower()
    q = re.sub(r"['’]s?\b", "", q)  # remove possessive
    q = re.sub(r"\s+", " ", q).strip()
    return _ALIASES.get(q, q)


def resolve_disease(query: str) -> List[Dict]:
    """
    Resolve a free-text query to one or more MeSH disease records.
    Returns list of dicts with: mesh_id, heading, tree_numbers, entry_terms, parent_ids
    """
    norm = _normalize(query)

    # 1. Exact heading match
    rows = _q(
        "SELECT * FROM mesh_diseases WHERE LOWER(heading) = %s",
        (norm,),
    )
    if rows:
        return rows

    # 2. Exact entry term match
    rows = _q(
        "SELECT * FROM mesh_diseases WHERE %s = ANY(LOWER(entry_terms::text)::text[])",
        (norm,),
    )
    if not rows:
        rows = _q(
            "SELECT * FROM mesh_diseases WHERE EXISTS ("
            "  SELECT 1 FROM unnest(entry_terms) t WHERE LOWER(t) = %s"
            ")",
            (norm,),
        )
    if rows:
        return rows

    # 3. Partial heading match
    rows = _q(
        "SELECT * FROM mesh_diseases WHERE LOWER(heading) LIKE %s ORDER BY heading LIMIT 10",
        (f"%{norm}%",),
    )
    if rows:
        return rows

    # 4. Partial entry term match
    rows = _q(
        "SELECT * FROM mesh_diseases WHERE EXISTS ("
        "  SELECT 1 FROM unnest(entry_terms) t WHERE LOWER(t) LIKE %s"
        ") ORDER BY heading LIMIT 10",
        (f"%{norm}%",),
    )
    if rows:
        return rows

    # 5. Word-by-word fallback — match first significant word
    words = [w for w in norm.split() if len(w) > 3]
    if words:
        pattern = f"%{words[0]}%"
        rows = _q(
            "SELECT * FROM mesh_diseases WHERE LOWER(heading) LIKE %s ORDER BY heading LIMIT 10",
            (pattern,),
        )
        if rows:
            return rows

    return [{
        "mesh_id": norm,
        "heading": query.strip() or norm,
        "tree_numbers": [],
        "entry_terms": [],
        "parent_ids": [],
    }]


def expand_mesh_ids(mesh_ids: List[str], include_children: bool = True) -> List[str]:
    """
    Given a list of MeSH IDs, return all related IDs:
    - The IDs themselves
    - Their parents (broader disease category)
    - Their children (more specific subtypes) — if include_children=True
    """
    if not mesh_ids:
        return []

    if all(str(mid).startswith("LOCAL:") for mid in mesh_ids):
        return mesh_ids

    all_ids = set(mesh_ids)

    # Add parents (one level up)
    for mid in mesh_ids:
        rows = _q("SELECT parent_ids FROM mesh_diseases WHERE mesh_id = %s", (mid,))
        for r in rows:
            for pid in (r.get("parent_ids") or []):
                all_ids.add(pid)

    if include_children:
        # Add children: any disease whose parent_ids contains one of our IDs
        placeholders = ",".join(["%s"] * len(mesh_ids))
        rows = _q(
            f"SELECT mesh_id FROM mesh_diseases "
            f"WHERE parent_ids && ARRAY[{placeholders}]::text[]",
            mesh_ids,
        )
        for r in rows:
            all_ids.add(r["mesh_id"])

    return list(all_ids)


def get_tree_siblings(mesh_id: str) -> List[str]:
    """Return MeSH IDs that share the same immediate parent (siblings in the tree)."""
    rows = _q("SELECT parent_ids FROM mesh_diseases WHERE mesh_id = %s", (mesh_id,))
    if not rows or not rows[0].get("parent_ids"):
        return []
    parent_ids = rows[0]["parent_ids"]
    if not parent_ids:
        return []
    placeholders = ",".join(["%s"] * len(parent_ids))
    sibling_rows = _q(
        f"SELECT mesh_id FROM mesh_diseases "
        f"WHERE parent_ids && ARRAY[{placeholders}]::text[] AND mesh_id != %s",
        parent_ids + [mesh_id],
    )
    return [r["mesh_id"] for r in sibling_rows]


def get_disease_label(mesh_id: str) -> str:
    """Return the canonical heading for a MeSH ID."""
    if str(mesh_id).startswith("LOCAL:"):
        return str(mesh_id).replace("LOCAL:", "").replace("-", " ").title()

    rows = _q("SELECT heading FROM mesh_diseases WHERE mesh_id = %s", (mesh_id,))
    return rows[0]["heading"] if rows else mesh_id


def suggest_diseases(partial: str, limit: int = 10) -> List[str]:
    """Return disease name suggestions for autocomplete."""
    if not partial or len(partial) < 2:
        return []
    norm = partial.strip().lower()
    rows = _q(
        "SELECT DISTINCT heading FROM mesh_diseases "
        "WHERE LOWER(heading) LIKE %s ORDER BY heading LIMIT %s",
        (f"%{norm}%", limit),
    )
    suggestions = [r["heading"] for r in rows]
    if suggestions:
        return suggestions
    fallback = [
        "Alzheimer Disease",
        "Parkinson Disease",
        "Multiple Sclerosis",
        "Epilepsy",
        "Depression",
    ]
    return [s for s in fallback if norm in s.lower()][:limit]


def mesh_available() -> bool:
    """Check whether the mesh_diseases table has data."""
    rows = _q("SELECT COUNT(*) AS n FROM mesh_diseases")
    return bool(rows) and rows[0]["n"] > 0
