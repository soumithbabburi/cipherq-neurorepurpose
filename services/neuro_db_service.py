"""
NeuroRepurpose Database Service
All queries go to the neurorepurpose PostgreSQL database.
Uses MeSH IDs (not string LIKE) for disease lookups.
Falls back to heading-based LIKE only when mesh_id column is not populated.
"""

import logging
import os
from functools import lru_cache
from typing import Dict, List, Optional

import psycopg2
import psycopg2.extras
import psycopg2.pool

logger = logging.getLogger(__name__)

_DB_PARAMS = dict(
    host=os.getenv("DB_HOST", "localhost"),
    port=int(os.getenv("DB_PORT", "5434")),
    user=os.getenv("DB_USER", "babburisoumith"),
    password=os.getenv("DB_PASSWORD", ""),
    dbname="neurorepurpose",
    connect_timeout=3,
)

_pool: Optional[psycopg2.pool.ThreadedConnectionPool] = None
_db_available: Optional[bool] = None


def _get_pool() -> Optional[psycopg2.pool.ThreadedConnectionPool]:
    global _pool, _db_available
    if _pool is not None:
        return _pool
    if _db_available is False:
        return None
    try:
        _pool = psycopg2.pool.ThreadedConnectionPool(1, 10, **_DB_PARAMS)
        _db_available = True
        return _pool
    except Exception as e:
        _db_available = False
        logger.warning(f"neurorepurpose DB unavailable: {e}")
        return None


def _query(sql: str, params=None) -> List[Dict]:
    pool = _get_pool()
    if pool is None:
        return []
    conn = pool.getconn()
    try:
        with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
            cur.execute(sql, params)
            return [dict(r) for r in cur.fetchall()]
    except Exception as e:
        logger.error(f"DB query error: {e}")
        return []
    finally:
        pool.putconn(conn)


def is_available() -> bool:
    return _get_pool() is not None


# ─── Compound queries ─────────────────────────────────────────────────────────

def get_compounds_for_disease(mesh_ids: List[str], limit: int = 50) -> List[Dict]:
    """
    Return compounds with a ChEMBL indication matching any of the given MeSH IDs.
    Falls back to LIKE on heading when mesh_id column is empty.
    """
    if not mesh_ids:
        return []

    placeholders = ",".join(["%s"] * len(mesh_ids))

    # Primary: use mesh_id column
    rows = _query(
        f"""
        SELECT DISTINCT
            c.id, c.chembl_id, c.name, c.smiles, c.max_phase,
            cp.mw, cp.alogp, cp.psa, cp.hba, cp.hbd, cp.ro5_violations,
            STRING_AGG(DISTINCT i.disease,   '; ') AS indications,
            STRING_AGG(DISTINCT m.mechanism, '; ') AS mechanisms,
            STRING_AGG(DISTINCT t.name,      '; ') AS targets
        FROM compounds c
        JOIN indications i ON i.compound_id = c.id
        LEFT JOIN compound_properties cp ON cp.compound_id = c.id
        LEFT JOIN mechanisms m           ON m.compound_id  = c.id
        LEFT JOIN targets t              ON t.id = m.target_id
        WHERE i.mesh_id = ANY(ARRAY[{placeholders}]::varchar[])
        GROUP BY c.id, c.chembl_id, c.name, c.smiles, c.max_phase,
                 cp.mw, cp.alogp, cp.psa, cp.hba, cp.hbd, cp.ro5_violations
        ORDER BY c.max_phase DESC NULLS LAST
        LIMIT %s
        """,
        mesh_ids + [limit],
    )

    # Fallback: heading LIKE (when mesh backfill not yet run)
    if not rows:
        like_clauses = " OR ".join(
            f"LOWER(i.disease) LIKE %s" for _ in mesh_ids
        )
        patterns = [f"%{mid.lower()}%" for mid in mesh_ids]
        rows = _query(
            f"""
            SELECT DISTINCT
                c.id, c.chembl_id, c.name, c.smiles, c.max_phase,
                cp.mw, cp.alogp, cp.psa, cp.hba, cp.hbd, cp.ro5_violations,
                STRING_AGG(DISTINCT i.disease,   '; ') AS indications,
                STRING_AGG(DISTINCT m.mechanism, '; ') AS mechanisms,
                STRING_AGG(DISTINCT t.name,      '; ') AS targets
            FROM compounds c
            JOIN indications i ON i.compound_id = c.id
            LEFT JOIN compound_properties cp ON cp.compound_id = c.id
            LEFT JOIN mechanisms m           ON m.compound_id  = c.id
            LEFT JOIN targets t              ON t.id = m.target_id
            WHERE {like_clauses}
            GROUP BY c.id, c.chembl_id, c.name, c.smiles, c.max_phase,
                     cp.mw, cp.alogp, cp.psa, cp.hba, cp.hbd, cp.ro5_violations
            ORDER BY c.max_phase DESC NULLS LAST
            LIMIT %s
            """,
            patterns + [limit],
        )

    return rows


def search_compounds(query: str, limit: int = 30) -> List[Dict]:
    """Full-text search across name, mechanism, and target fields."""
    pattern = f"%{query.lower()}%"
    return _query(
        """
        SELECT DISTINCT
            c.id, c.chembl_id, c.name, c.smiles, c.max_phase,
            cp.mw, cp.alogp, cp.psa,
            STRING_AGG(DISTINCT m.mechanism, '; ') AS mechanisms,
            STRING_AGG(DISTINCT t.name,      '; ') AS targets
        FROM compounds c
        LEFT JOIN compound_properties cp ON cp.compound_id = c.id
        LEFT JOIN mechanisms m           ON m.compound_id  = c.id
        LEFT JOIN targets t              ON t.id = m.target_id
        WHERE LOWER(c.name) LIKE %s
           OR LOWER(m.mechanism) LIKE %s
           OR LOWER(t.name) LIKE %s
        GROUP BY c.id, c.chembl_id, c.name, c.smiles, c.max_phase,
                 cp.mw, cp.alogp, cp.psa
        ORDER BY c.max_phase DESC NULLS LAST, c.name
        LIMIT %s
        """,
        (pattern, pattern, pattern, limit),
    )


def get_compound_by_chembl(chembl_id: str) -> Optional[Dict]:
    rows = _query(
        """
        SELECT c.*, cp.mw, cp.alogp, cp.psa, cp.hba, cp.hbd, cp.rtb, cp.ro5_violations
        FROM compounds c
        LEFT JOIN compound_properties cp ON cp.compound_id = c.id
        WHERE c.chembl_id = %s
        """,
        (chembl_id,),
    )
    return rows[0] if rows else None


def get_compound_targets(compound_id: int) -> List[Dict]:
    return _query(
        """
        SELECT t.name, t.target_type, t.gene_symbol, t.organism,
               m.mechanism, m.action_type, m.confidence
        FROM mechanisms m
        JOIN targets t ON t.id = m.target_id
        WHERE m.compound_id = %s
        ORDER BY m.confidence DESC NULLS LAST
        """,
        (compound_id,),
    )


def get_compound_activities(compound_id: int) -> List[Dict]:
    return _query(
        """
        SELECT t.name AS target_name, t.gene_symbol,
               ca.activity_type, ca.pchembl_value, ca.standard_value, ca.standard_units
        FROM compound_activities ca
        JOIN targets t ON t.id = ca.target_id
        WHERE ca.compound_id = %s
        ORDER BY ca.pchembl_value DESC NULLS LAST
        """,
        (compound_id,),
    )


def get_compound_indications(compound_id: int) -> List[Dict]:
    return _query(
        """
        SELECT i.disease, i.mesh_id, i.efo_id, i.max_phase, i.source,
               md.tree_numbers, md.entry_terms
        FROM indications i
        LEFT JOIN mesh_diseases md ON md.mesh_id = i.mesh_id
        WHERE i.compound_id = %s
        ORDER BY i.max_phase DESC NULLS LAST
        """,
        (compound_id,),
    )


def get_hetionet_paths(compound_name: str, mesh_ids: List[str]) -> List[Dict]:
    """Return HetioNet edges linking this compound to disease-relevant nodes."""
    c_rows = _query(
        "SELECT id FROM hetionet_nodes WHERE LOWER(name) = LOWER(%s) AND kind = 'Compound' LIMIT 1",
        (compound_name,),
    )
    if not c_rows:
        c_rows = _query(
            "SELECT id FROM hetionet_nodes WHERE LOWER(name) LIKE %s AND kind = 'Compound' LIMIT 1",
            (f"%{compound_name[:12].lower()}%",),
        )
    if not c_rows:
        return []

    c_id = c_rows[0]["id"]
    return _query(
        """
        SELECT hn_s.name AS source_name, hn_s.kind AS source_kind,
               he.metaedge,
               hn_t.name AS target_name, hn_t.kind AS target_kind
        FROM hetionet_edges he
        JOIN hetionet_nodes hn_s ON hn_s.id = he.source_id
        JOIN hetionet_nodes hn_t ON hn_t.id = he.target_id
        WHERE he.source_id = %s
        ORDER BY he.metaedge
        LIMIT 50
        """,
        (c_id,),
    )


def get_disease_top_compounds(mesh_ids: List[str], limit: int = 3) -> List[Dict]:
    """
    Return the top compounds for a disease — Phase 4 first, then by mechanism count.
    Used for the top-N recommendation feature.
    """
    if not mesh_ids:
        return []
    placeholders = ",".join(["%s"] * len(mesh_ids))
    return _query(
        f"""
        SELECT DISTINCT
            c.id, c.chembl_id, c.name, c.smiles, c.max_phase,
            cp.mw, cp.alogp, cp.psa,
            COUNT(DISTINCT m.id) AS mechanism_count,
            STRING_AGG(DISTINCT i.disease,   '; ') AS indications,
            STRING_AGG(DISTINCT m.mechanism, '; ') AS mechanisms,
            STRING_AGG(DISTINCT t.name,      '; ') AS targets
        FROM compounds c
        JOIN indications i ON i.compound_id = c.id
        LEFT JOIN compound_properties cp ON cp.compound_id = c.id
        LEFT JOIN mechanisms m           ON m.compound_id  = c.id
        LEFT JOIN targets t              ON t.id = m.target_id
        WHERE i.mesh_id = ANY(ARRAY[{placeholders}]::varchar[])
          AND c.max_phase >= 3
        GROUP BY c.id, c.chembl_id, c.name, c.smiles, c.max_phase,
                 cp.mw, cp.alogp, cp.psa
        ORDER BY c.max_phase DESC NULLS LAST, mechanism_count DESC
        LIMIT %s
        """,
        mesh_ids + [limit],
    )


def get_stats() -> Dict:
    results = {}
    for table in ["compounds", "targets", "mechanisms", "indications", "hetionet_nodes", "mesh_diseases"]:
        rows = _query(f"SELECT COUNT(*) AS n FROM {table}")
        results[table] = rows[0]["n"] if rows else 0
    return results


def get_available_diseases(limit: int = 200) -> List[str]:
    """Return disease headings that have at least one compound in the DB."""
    rows = _query(
        """
        SELECT DISTINCT i.disease
        FROM indications i
        JOIN compounds c ON c.id = i.compound_id
        WHERE i.disease IS NOT NULL
        ORDER BY i.disease
        LIMIT %s
        """,
        (limit,),
    )
    return [r["disease"] for r in rows]
