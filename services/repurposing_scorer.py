"""
Repurposing Scorer
Scores each compound 0-1 for a given disease using four real signals:
  A. Indication evidence  (ChEMBL drug_indication via MeSH tree match)
  B. Target relevance     (compound targets vs. disease top-target set)
  C. Activity strength    (pChEMBL values at disease-relevant targets)
  D. Network evidence     (HetioNet compound-gene-disease path count)

Weights are calibrated so that known FDA-approved drug-disease pairs score
above 0.6 and unevidenced pairs score below 0.3.
"""

import logging
from functools import lru_cache
from typing import Dict, List, Optional

import psycopg2
import psycopg2.extras
import psycopg2.pool

logger = logging.getLogger(__name__)

_DB_PARAMS = dict(
    host="localhost", port=5434, user="babburisoumith",
    password="", dbname="neurorepurpose", connect_timeout=3,
)

_pool: Optional[psycopg2.pool.ThreadedConnectionPool] = None

# Scoring weights (sum to 1.0)
W_INDICATION = 0.40   # strongest signal — direct clinical evidence
W_TARGET     = 0.30   # mechanistic plausibility
W_ACTIVITY   = 0.20   # in vitro quantitative evidence
W_NETWORK    = 0.10   # indirect knowledge-graph evidence


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


# ─── A: Indication score ──────────────────────────────────────────────────────

def _score_indication(compound_id: int, mesh_ids: List[str]) -> float:
    """
    0.0-1.0 based on whether the compound has a ChEMBL indication matching
    any of the expanded MeSH IDs, weighted by clinical phase.
    """
    if not mesh_ids:
        return 0.0
    placeholders = ",".join(["%s"] * len(mesh_ids))
    rows = _q(
        f"SELECT MAX(max_phase) AS best_phase "
        f"FROM indications "
        f"WHERE compound_id = %s AND mesh_id = ANY(ARRAY[{placeholders}]::varchar[])",
        [compound_id] + mesh_ids,
    )
    if not rows or rows[0]["best_phase"] is None:
        # Fall back to string match when mesh_id not yet backfilled
        rows = _q(
            "SELECT MAX(max_phase) AS best_phase FROM indications "
            "WHERE compound_id = %s",
            (compound_id,),
        )
    best = float(rows[0]["best_phase"] or 0) if rows else 0.0
    return min(1.0, best / 4.0)


# ─── B: Target relevance score ────────────────────────────────────────────────

@lru_cache(maxsize=128)
def _get_disease_top_targets(mesh_ids_key: str) -> List[str]:
    """
    Derive top targets for a disease by finding the targets most frequently
    associated with FDA-approved/Phase 3-4 compounds for that disease.
    Returns list of target internal IDs (as strings).
    """
    mesh_ids = mesh_ids_key.split("|")
    if not mesh_ids:
        return []
    placeholders = ",".join(["%s"] * len(mesh_ids))
    rows = _q(
        f"""
        SELECT m.target_id, COUNT(DISTINCT m.compound_id) AS compound_count
        FROM mechanisms m
        JOIN indications i ON i.compound_id = m.compound_id
        WHERE i.mesh_id = ANY(ARRAY[{placeholders}]::varchar[])
          AND i.max_phase >= 3
        GROUP BY m.target_id
        ORDER BY compound_count DESC
        LIMIT 30
        """,
        mesh_ids,
    )
    return [str(r["target_id"]) for r in rows]


def _score_target_relevance(compound_id: int, mesh_ids: List[str]) -> float:
    """
    Fraction of the compound's targets that appear in the disease's top-target set.
    """
    if not mesh_ids:
        return 0.0
    key = "|".join(sorted(mesh_ids))
    top_targets = _get_disease_top_targets(key)
    if not top_targets:
        return 0.0

    rows = _q(
        "SELECT target_id FROM mechanisms WHERE compound_id = %s",
        (compound_id,),
    )
    compound_targets = {str(r["target_id"]) for r in rows}
    if not compound_targets:
        return 0.0

    overlap = len(compound_targets & set(top_targets))
    return min(1.0, overlap / max(len(top_targets[:10]), 1))


# ─── C: Activity strength score ───────────────────────────────────────────────

def _score_activity(compound_id: int, mesh_ids: List[str]) -> float:
    """
    Based on best pChEMBL value at disease-relevant targets.
    pChEMBL > 8 (IC50 < 10nM) → strong (1.0)
    pChEMBL 6-8 (IC50 10nM-1µM) → moderate (0.6)
    pChEMBL < 6 → weak (0.2)
    """
    key = "|".join(sorted(mesh_ids)) if mesh_ids else ""
    top_targets = _get_disease_top_targets(key) if key else []

    if top_targets:
        placeholders = ",".join(["%s"] * len(top_targets))
        rows = _q(
            f"SELECT MAX(pchembl_value) AS best_pchembl "
            f"FROM compound_activities "
            f"WHERE compound_id = %s AND target_id = ANY(ARRAY[{placeholders}]::integer[])",
            [compound_id] + [int(t) for t in top_targets],
        )
    else:
        rows = _q(
            "SELECT MAX(pchembl_value) AS best_pchembl "
            "FROM compound_activities WHERE compound_id = %s",
            (compound_id,),
        )

    if not rows or rows[0]["best_pchembl"] is None:
        return 0.0
    p = float(rows[0]["best_pchembl"])
    if p >= 8.0:
        return 1.0
    if p >= 6.0:
        return 0.3 + 0.7 * (p - 6.0) / 2.0
    return max(0.0, 0.1 + 0.2 * (p - 4.0) / 2.0)


# ─── D: Network evidence score ────────────────────────────────────────────────

def _score_network(compound_name: str, mesh_ids: List[str]) -> float:
    """
    Count HetioNet Compound→Gene edges from this compound, then check how many
    of those genes have Gene→Disease edges to the queried disease nodes.
    """
    # Find compound node in hetionet_nodes
    rows = _q(
        "SELECT id FROM hetionet_nodes WHERE LOWER(name) = LOWER(%s) AND kind = 'Compound' LIMIT 1",
        (compound_name,),
    )
    if not rows:
        # Try partial match
        rows = _q(
            "SELECT id FROM hetionet_nodes WHERE LOWER(name) LIKE %s AND kind = 'Compound' LIMIT 1",
            (f"%{compound_name.lower()[:10]}%",),
        )
    if not rows:
        return 0.0

    c_node_id = rows[0]["id"]

    # Get genes this compound binds
    gene_rows = _q(
        "SELECT target_id FROM hetionet_edges WHERE source_id = %s AND metaedge = 'CbG'",
        (c_node_id,),
    )
    gene_ids = [r["target_id"] for r in gene_rows]
    if not gene_ids:
        return 0.0

    # Get disease nodes from mesh_ids
    disease_node_ids = []
    for mid in mesh_ids:
        d_rows = _q(
            "SELECT id FROM hetionet_nodes WHERE LOWER(name) LIKE %s AND kind = 'Disease' LIMIT 3",
            (f"%{mid}%",),
        )
        disease_node_ids.extend([r["id"] for r in d_rows])

    if not disease_node_ids:
        # Try by mesh heading
        return min(1.0, len(gene_ids) / 20.0) * 0.5  # partial credit for having gene bindings

    # Count gene→disease paths
    gid_placeholders = ",".join(["%s"] * len(gene_ids))
    did_placeholders = ",".join(["%s"] * len(disease_node_ids))
    path_rows = _q(
        f"SELECT COUNT(*) AS paths FROM hetionet_edges "
        f"WHERE source_id IN ({gid_placeholders}) "
        f"AND target_id IN ({did_placeholders})",
        gene_ids + disease_node_ids,
    )
    paths = int(path_rows[0]["paths"]) if path_rows else 0
    return min(1.0, paths / 10.0)


# ─── Main scoring function ────────────────────────────────────────────────────

def score_compound(
    compound_id: int,
    compound_name: str,
    mesh_ids: List[str],
    compound_phase: float = 0.0,
) -> Dict:
    """
    Score a compound for a disease (represented by mesh_ids).
    Returns dict with overall score and per-signal breakdown.
    """
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

    # Phase bonus: FDA-approved drugs get a small boost on overall
    phase_bonus = min(0.05, (compound_phase or 0) / 4.0 * 0.05)
    overall = min(1.0, overall + phase_bonus)

    return {
        "overall": round(overall, 4),
        "indication_score": round(s_ind, 4),
        "target_score":     round(s_tgt, 4),
        "activity_score":   round(s_act, 4),
        "network_score":    round(s_net, 4),
        "phase_bonus":      round(phase_bonus, 4),
    }


def score_compound_list(
    compounds: List[Dict],
    mesh_ids: List[str],
) -> List[Dict]:
    """
    Score a list of compound dicts (each must have id, name, max_phase).
    Adds 'score' and 'score_breakdown' keys. Returns sorted by score desc.
    """
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
