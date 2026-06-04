"""
Repurposing Scorer
Scores each compound 0-1 for a given disease using four real signals:
  A. Indication evidence  (ChEMBL drug_indication, matched by MeSH ID or disease name)
  B. Target relevance     (compound targets vs. disease top-target set)
  C. Activity strength    (pChEMBL values; falls back to mechanism count proxy)
  D. Network evidence     (HetioNet CtD / CbG paths when available)
"""

import logging
import os
import re
from functools import lru_cache
from typing import Dict, List, Optional

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

W_INDICATION = 0.40
W_TARGET     = 0.30
W_ACTIVITY   = 0.20
W_NETWORK    = 0.10

_GENERIC_WORDS = {
    'disease', 'disorder', 'syndrome', 'condition', 'illness', 'type', 'with',
    'associated', 'related', 'injury', 'injuries', 'deficit', 'impairment',
    'spectrum', 'traumatic', 'post', 'chronic', 'acute', 'primary',
}


def _get_pool():
    global _pool
    if _pool is None:
        try:
            _pool = psycopg2.pool.ThreadedConnectionPool(1, 5, **_DB_PARAMS)
        except Exception as e:
            logger.warning(f"Scorer DB pool failed: {e}")
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
        logger.debug(f"Scorer query failed: {e}")
        return []
    finally:
        pool.putconn(conn)


def _keyword(text: str) -> str:
    """Extract the most specific word from a disease name for LIKE queries."""
    clean = re.sub(r"['''\-,]", " ", text.lower())
    words = [w for w in clean.split() if len(w) > 3 and w not in _GENERIC_WORDS]
    return words[0] if words else re.sub(r"[^a-z0-9]", "", text.lower())[:12]


def _is_mesh_id(s: str) -> bool:
    """Return True if string looks like a real MeSH descriptor ID."""
    return bool(re.match(r'^D\d{6,}$', s.strip()))


def _disease_like_patterns(mesh_ids: List[str]) -> List[str]:
    """Convert a list of mesh_ids/names to LIKE pattern strings."""
    return [f"%{_keyword(m)}%" for m in mesh_ids if m]


# ─── A: Indication score ──────────────────────────────────────────────────────

def _score_indication(compound_id: int, mesh_ids: List[str]) -> float:
    if not mesh_ids:
        return 0.0

    # Try direct mesh_id match first (works when real MeSH IDs are loaded)
    real_ids = [m for m in mesh_ids if _is_mesh_id(m)]
    if real_ids:
        rows = _q(
            "SELECT MAX(max_phase) AS best "
            "FROM indications WHERE compound_id = %s AND mesh_id = ANY(%s::varchar[])",
            (compound_id, real_ids),
        )
        if rows and rows[0]["best"]:
            return min(1.0, float(rows[0]["best"]) / 4.0)

    # Fallback: disease name LIKE match (handles NULL mesh_id rows)
    best_phase = 0.0
    for pattern in _disease_like_patterns(mesh_ids):
        rows = _q(
            "SELECT MAX(max_phase) AS best "
            "FROM indications WHERE compound_id = %s AND LOWER(disease) LIKE %s",
            (compound_id, pattern),
        )
        if rows and rows[0]["best"]:
            best_phase = max(best_phase, float(rows[0]["best"]))

    return min(1.0, best_phase / 4.0)


# ─── B: Target relevance score ────────────────────────────────────────────────

@lru_cache(maxsize=256)
def _get_disease_top_targets(mesh_ids_key: str) -> List[str]:
    """
    Top internal target IDs for a disease, ranked by how many high-phase
    compounds target them for that indication. Tries MeSH ID then name.
    """
    mesh_ids = [m for m in mesh_ids_key.split("|") if m]
    if not mesh_ids:
        return []

    # Try real MeSH IDs
    real_ids = [m for m in mesh_ids if _is_mesh_id(m)]
    if real_ids:
        rows = _q(
            """
            SELECT m.target_id, COUNT(DISTINCT m.compound_id) AS cnt
            FROM mechanisms m
            JOIN indications i ON i.compound_id = m.compound_id
            WHERE i.mesh_id = ANY(%s::varchar[]) AND i.max_phase >= 2
            GROUP BY m.target_id ORDER BY cnt DESC LIMIT 30
            """,
            (real_ids,),
        )
        if rows:
            return [str(r["target_id"]) for r in rows]

    # Fallback: disease name LIKE
    all_rows: List[Dict] = []
    for pattern in _disease_like_patterns(mesh_ids):
        rows = _q(
            """
            SELECT m.target_id, COUNT(DISTINCT m.compound_id) AS cnt
            FROM mechanisms m
            JOIN indications i ON i.compound_id = m.compound_id
            WHERE LOWER(i.disease) LIKE %s AND i.max_phase >= 2
            GROUP BY m.target_id ORDER BY cnt DESC LIMIT 30
            """,
            (pattern,),
        )
        all_rows.extend(rows)

    # Aggregate across patterns
    tgt_counts: Dict[str, int] = {}
    for r in all_rows:
        tid = str(r["target_id"])
        tgt_counts[tid] = tgt_counts.get(tid, 0) + int(r["cnt"])
    sorted_tgts = sorted(tgt_counts, key=lambda t: tgt_counts[t], reverse=True)
    return sorted_tgts[:30]


def _score_target_relevance(compound_id: int, mesh_ids: List[str]) -> float:
    if not mesh_ids:
        return 0.0
    key = "|".join(sorted(mesh_ids))
    top_targets = _get_disease_top_targets(key)
    if not top_targets:
        return 0.0

    rows = _q("SELECT target_id FROM mechanisms WHERE compound_id = %s", (compound_id,))
    compound_targets = {str(r["target_id"]) for r in rows}
    if not compound_targets:
        return 0.0

    # Fraction of compound's targets that appear in disease top-target set
    overlap = len(compound_targets & set(top_targets))
    # Normalise: matching even 1 of the top-10 disease targets is meaningful
    score = min(1.0, overlap / max(1, min(len(top_targets), 10)))
    return score


# ─── C: Activity strength score ───────────────────────────────────────────────

def _score_activity(compound_id: int, mesh_ids: List[str]) -> float:
    key = "|".join(sorted(mesh_ids)) if mesh_ids else ""
    top_targets = _get_disease_top_targets(key) if key else []

    # Try pChEMBL at disease-relevant targets
    if top_targets:
        rows = _q(
            "SELECT MAX(pchembl_value) AS best "
            "FROM compound_activities "
            "WHERE compound_id = %s AND target_id = ANY(%s::integer[])",
            (compound_id, [int(t) for t in top_targets]),
        )
    else:
        rows = _q(
            "SELECT MAX(pchembl_value) AS best "
            "FROM compound_activities WHERE compound_id = %s",
            (compound_id,),
        )

    if rows and rows[0]["best"] is not None:
        p = float(rows[0]["best"])
        if p >= 8.0:   return 1.0
        if p >= 6.0:   return 0.3 + 0.7 * (p - 6.0) / 2.0
        return max(0.0, 0.1 + 0.2 * (p - 4.0) / 2.0)

    # No pChEMBL data — proxy: number of known mechanisms
    mech_rows = _q(
        "SELECT COUNT(*) AS cnt FROM mechanisms WHERE compound_id = %s",
        (compound_id,),
    )
    cnt = int(mech_rows[0]["cnt"]) if mech_rows else 0
    if cnt == 0:
        return 0.0
    # 1 target → 0.15, 5 targets → 0.55, 10+ targets → 0.75 (cap)
    return min(0.75, 0.10 + cnt * 0.065)


# ─── D: Network evidence score ────────────────────────────────────────────────

def _score_network(compound_name: str, mesh_ids: List[str]) -> float:
    # Find compound node in hetionet by exact then partial name
    rows = _q(
        "SELECT id FROM hetionet_nodes "
        "WHERE kind='Compound' AND LOWER(name) = LOWER(%s) LIMIT 1",
        (compound_name,),
    )
    if not rows:
        rows = _q(
            "SELECT id FROM hetionet_nodes "
            "WHERE kind='Compound' AND LOWER(name) LIKE %s LIMIT 1",
            (f"%{compound_name.lower()[:12]}%",),
        )
    if not rows:
        return 0.0

    c_node_id = rows[0]["id"]

    # Check CtD (Compound treats Disease) edges directly
    ctd_rows = _q(
        "SELECT he.target_id FROM hetionet_edges he "
        "WHERE he.source_id = %s AND he.metaedge = 'CtD'",
        (c_node_id,),
    )
    if ctd_rows:
        # Find disease nodes matching the query
        disease_node_ids = set()
        for mid in mesh_ids:
            kw = f"%{_keyword(mid)}%"
            d_rows = _q(
                "SELECT id FROM hetionet_nodes WHERE kind='Disease' AND LOWER(name) LIKE %s LIMIT 5",
                (kw,),
            )
            disease_node_ids.update(r["id"] for r in d_rows)
        ctd_targets = {r["target_id"] for r in ctd_rows}
        if ctd_targets & disease_node_ids:
            return 1.0  # direct compound-treats-disease edge
        # Has CtD edges but to other diseases — partial credit
        return 0.3

    # Check CbG (Compound binds Gene) — gene count as proxy
    cbg_rows = _q(
        "SELECT COUNT(*) AS cnt FROM hetionet_edges "
        "WHERE source_id = %s AND metaedge = 'CbG'",
        (c_node_id,),
    )
    cbg_cnt = int(cbg_rows[0]["cnt"]) if cbg_rows else 0
    if cbg_cnt > 0:
        return min(0.4, cbg_cnt / 20.0)

    return 0.0


# ─── Main scoring functions ────────────────────────────────────────────────────

def score_compound(
    compound_id: int,
    compound_name: str,
    mesh_ids: List[str],
    compound_phase: float = 0.0,
) -> Dict:
    s_ind = _score_indication(compound_id, mesh_ids)
    s_tgt = _score_target_relevance(compound_id, mesh_ids)
    s_act = _score_activity(compound_id, mesh_ids)
    s_net = _score_network(compound_name, mesh_ids)

    overall = (
        W_INDICATION * s_ind +
        W_TARGET     * s_tgt +
        W_ACTIVITY   * s_act +
        W_NETWORK    * s_net
    )

    phase_bonus = min(0.05, (compound_phase or 0) / 4.0 * 0.05)
    overall = min(1.0, overall + phase_bonus)

    return {
        "overall":           round(overall, 4),
        "indication_score":  round(s_ind, 4),
        "target_score":      round(s_tgt, 4),
        "activity_score":    round(s_act, 4),
        "network_score":     round(s_net, 4),
        "phase_bonus":       round(phase_bonus, 4),
    }


def score_compound_list(compounds: List[Dict], mesh_ids: List[str]) -> List[Dict]:
    for c in compounds:
        breakdown = score_compound(
            compound_id=c["id"],
            compound_name=c.get("name", ""),
            mesh_ids=mesh_ids,
            compound_phase=float(c.get("max_phase") or 0),
        )
        c["score"] = breakdown["overall"]
        c["score_breakdown"] = breakdown

    compounds.sort(key=lambda x: x["score"], reverse=True)
    return compounds
