"""
NeuroRepurpose Database Service
Queries the neurorepurpose PostgreSQL database (backed by ChEMBL 33 + HetioNet).
Falls back gracefully if the database isn't available.
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
        _pool = psycopg2.pool.ThreadedConnectionPool(1, 5, **_DB_PARAMS)
        _db_available = True
        logger.info("neurorepurpose DB pool ready")
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
        logger.error(f"DB query failed: {e}")
        return []
    finally:
        pool.putconn(conn)


def is_available() -> bool:
    return _get_pool() is not None


# ─── Public API ───────────────────────────────────────────────────────────────

def get_neuro_compounds(disease: str, limit: int = 20) -> List[Dict]:
    """Return compounds indicated for a neurological disease."""
    pattern = f"%{disease.lower()}%"
    return _query(
        """
        SELECT DISTINCT c.id, c.chembl_id, c.name, c.smiles, c.max_phase,
               cp.mw, cp.alogp, cp.psa, cp.hba, cp.hbd,
               STRING_AGG(DISTINCT i.disease, '; ') AS indications,
               STRING_AGG(DISTINCT m.mechanism, '; ') AS mechanisms,
               STRING_AGG(DISTINCT t.name, '; ') AS targets
        FROM compounds c
        JOIN indications i ON i.compound_id = c.id
        LEFT JOIN compound_properties cp ON cp.compound_id = c.id
        LEFT JOIN mechanisms m ON m.compound_id = c.id
        LEFT JOIN targets t ON t.id = m.target_id
        WHERE LOWER(i.disease) LIKE %s
        GROUP BY c.id, c.chembl_id, c.name, c.smiles, c.max_phase,
                 cp.mw, cp.alogp, cp.psa, cp.hba, cp.hbd
        ORDER BY c.max_phase DESC NULLS LAST, c.name
        LIMIT %s
        """,
        (pattern, limit),
    )


def search_compounds(query: str, limit: int = 20) -> List[Dict]:
    """Full-text search across compound names, mechanisms, and targets."""
    pattern = f"%{query.lower()}%"
    return _query(
        """
        SELECT DISTINCT c.id, c.chembl_id, c.name, c.smiles, c.max_phase,
               cp.mw, cp.alogp, cp.psa,
               STRING_AGG(DISTINCT m.mechanism, '; ') AS mechanisms,
               STRING_AGG(DISTINCT t.name, '; ') AS targets
        FROM compounds c
        LEFT JOIN compound_properties cp ON cp.compound_id = c.id
        LEFT JOIN mechanisms m ON m.compound_id = c.id
        LEFT JOIN targets t ON t.id = m.target_id
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


def get_compound_targets(chembl_id: str) -> List[Dict]:
    """Return protein targets for a compound."""
    return _query(
        """
        SELECT t.name, t.target_type, t.gene_symbol, t.organism,
               m.mechanism, m.action_type, m.confidence
        FROM compounds c
        JOIN mechanisms m ON m.compound_id = c.id
        JOIN targets t ON t.id = m.target_id
        WHERE c.chembl_id = %s
        ORDER BY m.confidence DESC NULLS LAST
        """,
        (chembl_id,),
    )


def get_compound_mechanisms(chembl_id: str) -> List[Dict]:
    """Return mechanisms of action for a compound."""
    return _query(
        """
        SELECT m.mechanism, m.action_type, m.confidence,
               t.name AS target_name, t.gene_symbol
        FROM compounds c
        JOIN mechanisms m ON m.compound_id = c.id
        LEFT JOIN targets t ON t.id = m.target_id
        WHERE c.chembl_id = %s
        ORDER BY m.confidence DESC NULLS LAST
        """,
        (chembl_id,),
    )


def get_compound_indications(chembl_id: str) -> List[Dict]:
    """Return disease indications for a compound."""
    return _query(
        """
        SELECT i.disease, i.efo_id, i.max_phase, i.source
        FROM compounds c
        JOIN indications i ON i.compound_id = c.id
        WHERE c.chembl_id = %s
        ORDER BY i.max_phase DESC NULLS LAST
        """,
        (chembl_id,),
    )


def get_hetionet_paths(compound_name: str, disease: str) -> List[Dict]:
    """Return HetioNet edges linking a compound to a disease."""
    return _query(
        """
        SELECT hn_s.name AS source_name, hn_s.kind AS source_kind,
               he.metaedge,
               hn_t.name AS target_name, hn_t.kind AS target_kind
        FROM hetionet_edges he
        JOIN hetionet_nodes hn_s ON hn_s.id = he.source_id
        JOIN hetionet_nodes hn_t ON hn_t.id = he.target_id
        WHERE (LOWER(hn_s.name) LIKE %s AND LOWER(hn_t.name) LIKE %s)
           OR (LOWER(hn_s.name) LIKE %s AND LOWER(hn_t.name) LIKE %s)
        LIMIT 50
        """,
        (
            f"%{compound_name.lower()}%", f"%{disease.lower()}%",
            f"%{disease.lower()}%", f"%{compound_name.lower()}%",
        ),
    )


def get_stats() -> Dict:
    """Return row counts for the main tables."""
    results = {}
    for table in ["compounds", "targets", "mechanisms", "indications", "hetionet_nodes"]:
        rows = _query(f"SELECT COUNT(*) AS n FROM {table}")
        results[table] = rows[0]["n"] if rows else 0
    return results
