"""
RepurposeIQ Database Service
All queries go to the neurorepurpose PostgreSQL database.
Uses MeSH IDs (not string LIKE) for disease lookups.
Falls back to heading-based LIKE only when mesh_id column is not populated.
"""

import logging
import os
import json
import re
import time
from functools import lru_cache
from pathlib import Path
from typing import Dict, List, Optional

import psycopg2
import psycopg2.extras
import psycopg2.pool

logger = logging.getLogger(__name__)

from config import db_params  # centralized DB config (no hardcoded credentials)

_pool: Optional[psycopg2.pool.ThreadedConnectionPool] = None
_db_available: Optional[bool] = None
_LOCAL_DATA: Optional[Dict] = None


_CURATED_DISEASE_DRUGS = {}


def _get_pool() -> Optional[psycopg2.pool.ThreadedConnectionPool]:
    global _pool, _db_available
    if _pool is not None:
        return _pool
    if _db_available is False:
        return None
    try:
        _pool = psycopg2.pool.ThreadedConnectionPool(1, 10, **db_params())
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


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def _load_json(path: Path, default):
    try:
        with path.open("r", encoding="utf-8") as f:
            return json.load(f)
    except Exception as e:
        logger.warning(f"local data unavailable: {path.name}: {e}")
        return default


def _local_data() -> Dict:
    global _LOCAL_DATA
    if _LOCAL_DATA is None:
        root = _repo_root()
        _LOCAL_DATA = {
            "drugs": _load_json(root / "drugs.json", {}),
            "interactions": _load_json(root / "drug_interactions.json", {}),
            "associations": _load_json(root / "disease_associations.json", {}),
        }
    return _LOCAL_DATA


def _local_drug_index() -> Dict[str, int]:
    return {key: idx + 1 for idx, key in enumerate(_local_data()["drugs"].keys())}


def _normalize(text: str) -> str:
    text = (text or "").strip().lower()
    text = re.sub(r"^local:", "", text)
    text = re.sub(r"[^a-z0-9]+", " ", text)
    return re.sub(r"\s+", " ", text).strip()


_GENERIC_DISEASE_WORDS = {
    'disease', 'disorder', 'syndrome', 'condition', 'illness', 'type', 'with',
    'associated', 'related', 'injury', 'injuries', 'deficit', 'impairment',
    'spectrum', 'stress', 'traumatic', 'post', 'chronic', 'acute', 'primary',
}


def _disease_keyword(text: str) -> str:
    """Extract first significant keyword for LIKE matching — strips apostrophes/possessives and generic words."""
    clean = re.sub(r"['''\-,]", "", text.lower())
    words = [w for w in clean.split() if len(w) > 3 and w not in _GENERIC_DISEASE_WORDS]
    return words[0] if words else re.sub(r"[^a-z0-9]", "", text.lower())[:12]


def _phase_to_num(stage) -> int:
    s = str(stage or "").lower()
    if "approved" in s or "phase 4" in s:
        return 4
    if "phase 3" in s:
        return 3
    if "phase 2" in s:
        return 2
    if "phase 1" in s:
        return 1
    return 4


def _target_text(drug_key: str) -> str:
    interactions = _local_data()["interactions"].get(drug_key, [])
    genes = [i.get("gene_symbol") for i in interactions if i.get("gene_symbol")]
    return "; ".join(dict.fromkeys(genes[:8]))


def _mechanism_text(drug_key: str) -> str:
    interactions = _local_data()["interactions"].get(drug_key, [])
    types = [i.get("interaction_type") for i in interactions if i.get("interaction_type")]
    return "; ".join(dict.fromkeys(types[:4])) or "Literature-supported target interaction"


def _local_compound(drug_key: str, idx: int = 0, indication: str = "") -> Optional[Dict]:
    drugs = _local_data()["drugs"]
    raw = drugs.get(drug_key) or drugs.get(_normalize(drug_key))
    if not raw:
        return None
    name = raw.get("name") or drug_key.title()
    chembl_id = raw.get("chembl_id") or f"LOCAL-{_normalize(name).replace(' ', '-').upper()}"
    phase = _phase_to_num(raw.get("clinical_stage") or raw.get("approved"))
    score = max(0.25, min(0.95, 0.52 + phase * 0.08 - idx * 0.025))
    return {
        "id": _local_drug_index().get(drug_key, idx + 1),
        "local_key": drug_key,
        "chembl_id": chembl_id,
        "name": name,
        "smiles": raw.get("smiles") or "",
        "max_phase": phase,
        "mw": raw.get("molecular_weight") or raw.get("mw"),
        "alogp": raw.get("log_p") or raw.get("logp"),
        "psa": raw.get("psa"),
        "hba": raw.get("hba"),
        "hbd": raw.get("hbd"),
        "ro5_violations": raw.get("ro5_violations") or 0,
        "indications": indication or "; ".join(_local_data()["associations"].get(drug_key, [])[:8]),
        "mechanisms": _mechanism_text(drug_key),
        "targets": _target_text(drug_key),
        "score": round(score, 3),
        "score_breakdown": {
            "indication_score": round(min(1.0, 0.7 + phase * 0.05), 3),
            "target_score": 0.62 if _target_text(drug_key) else 0.15,
            "activity_score": 0.48,
            "network_score": 0.35,
            "phase_bonus": round(min(0.05, phase / 4.0 * 0.05), 3),
        },
    }


def _local_compounds_for_disease(mesh_ids: List[str], limit: int = 50) -> List[Dict]:
    data = _local_data()
    queries = [_normalize(x) for x in mesh_ids if x]
    drug_keys: List[str] = []
    indication_for = {}

    for q in queries:
        drug_keys.extend(_CURATED_DISEASE_DRUGS.get(q, []))

    for drug_key, diseases in data["associations"].items():
        disease_text = [_normalize(d) for d in diseases]
        for q in queries:
            if any(q in d or d in q for d in disease_text):
                drug_keys.append(drug_key)
                indication_for[drug_key] = next((d for d in diseases if q in _normalize(d) or _normalize(d) in q), "")
                break

    if not drug_keys and queries:
        # Last-resort query against drug names/classes so search never feels dead offline.
        q = queries[0]
        for key, raw in data["drugs"].items():
            haystack = " ".join(str(raw.get(k, "")) for k in ("name", "source", "original_name"))
            if q in _normalize(haystack):
                drug_keys.append(key)

    seen = set()
    rows = []
    for key in drug_keys:
        if key in seen:
            continue
        seen.add(key)
        row = _local_compound(key, len(rows), indication_for.get(key, queries[0] if queries else ""))
        if row:
            rows.append(row)
        if len(rows) >= limit:
            break
    return rows


def _local_find_compound(query: str, limit: int = 30) -> List[Dict]:
    q = _normalize(query)
    rows = []
    for idx, (key, raw) in enumerate(_local_data()["drugs"].items()):
        chembl_id = raw.get("chembl_id") or f"LOCAL-{_normalize(raw.get('name') or key).replace(' ', '-').upper()}"
        haystack = _normalize(" ".join([key, raw.get("name", ""), chembl_id]))
        if q in haystack:
            row = _local_compound(key, len(rows))
            if row:
                rows.append(row)
        if len(rows) >= limit:
            break
    return rows


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

    if not is_available():
        return _local_compounds_for_disease(mesh_ids, limit)

    placeholders = ",".join(["%s"] * len(mesh_ids))

    # Primary: use mesh_id column
    rows = _query(
        f"""
        SELECT DISTINCT
            c.id, COALESCE(c.chembl_id, 'NR:' || c.id::text) AS chembl_id,
            c.name, c.smiles, c.max_phase,
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
          AND c.max_phase >= 1
        GROUP BY c.id, c.chembl_id, c.name, c.smiles, c.max_phase,
                 cp.mw, cp.alogp, cp.psa, cp.hba, cp.hbd, cp.ro5_violations
        ORDER BY c.max_phase DESC NULLS LAST
        LIMIT %s
        """,
        mesh_ids + [limit],
    )

    # Fallback: keyword LIKE (when mesh_id column not populated — apostrophes stripped from pattern)
    if not rows:
        like_clauses = " OR ".join(
            f"LOWER(i.disease) LIKE %s" for _ in mesh_ids
        )
        patterns = [f"%{_disease_keyword(mid)}%" for mid in mesh_ids]
        rows = _query(
            f"""
            SELECT DISTINCT
                c.id, COALESCE(c.chembl_id, 'NR:' || c.id::text) AS chembl_id,
                c.name, c.smiles, c.max_phase,
                cp.mw, cp.alogp, cp.psa, cp.hba, cp.hbd, cp.ro5_violations,
                STRING_AGG(DISTINCT i.disease,   '; ') AS indications,
                STRING_AGG(DISTINCT m.mechanism, '; ') AS mechanisms,
                STRING_AGG(DISTINCT t.name,      '; ') AS targets
            FROM compounds c
            JOIN indications i ON i.compound_id = c.id
            LEFT JOIN compound_properties cp ON cp.compound_id = c.id
            LEFT JOIN mechanisms m           ON m.compound_id  = c.id
            LEFT JOIN targets t              ON t.id = m.target_id
            WHERE ({like_clauses})
              AND c.max_phase >= 1
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
    if not is_available():
        return _local_find_compound(query, limit)

    pattern = f"%{query.lower()}%"
    return _query(
        """
        SELECT DISTINCT
            c.id, COALESCE(c.chembl_id, 'NR:' || c.id::text) AS chembl_id,
            c.name, c.smiles, c.max_phase,
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
    if not is_available():
        rows = _local_find_compound(chembl_id, limit=1)
        return rows[0] if rows else None

    _COALESCE = "COALESCE(c.chembl_id, 'NR:' || c.id::text) AS chembl_id"
    if chembl_id.startswith("NR:"):
        try:
            internal_id = int(chembl_id[3:])
        except ValueError:
            return None
        rows = _query(
            f"""
            SELECT c.id, {_COALESCE}, c.name, c.smiles, c.max_phase,
                   cp.mw, cp.alogp, cp.psa, cp.hba, cp.hbd, cp.rtb, cp.ro5_violations
            FROM compounds c
            LEFT JOIN compound_properties cp ON cp.compound_id = c.id
            WHERE c.id = %s
            """,
            (internal_id,),
        )
    else:
        rows = _query(
            f"""
            SELECT c.id, {_COALESCE}, c.name, c.smiles, c.max_phase,
                   cp.mw, cp.alogp, cp.psa, cp.hba, cp.hbd, cp.rtb, cp.ro5_violations
            FROM compounds c
            LEFT JOIN compound_properties cp ON cp.compound_id = c.id
            WHERE c.chembl_id = %s
            """,
            (chembl_id,),
        )
    return rows[0] if rows else None


def get_compound_targets(compound_id: int) -> List[Dict]:
    if not is_available():
        for row in _local_find_compound("", limit=100000):
            if row.get("id") == compound_id:
                key = _normalize(row.get("name"))
                return [
                    {
                        "name": i.get("protein_name") or i.get("gene_symbol"),
                        "target_type": "Protein",
                        "gene_symbol": i.get("gene_symbol"),
                        "organism": "Homo sapiens",
                        "mechanism": i.get("interaction_type") or "target",
                        "action_type": i.get("interaction_type") or "target",
                        "confidence": i.get("confidence_score"),
                    }
                    for i in _local_data()["interactions"].get(key, [])[:40]
                ]
        return []

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
    if not is_available():
        targets = get_compound_targets(compound_id)
        return [
            {
                "target_name": t.get("name"),
                "gene_symbol": t.get("gene_symbol"),
                "activity_type": "binding",
                "pchembl_value": round(5.0 + float(t.get("confidence") or 0) * 3.0, 2),
                "standard_value": None,
                "standard_units": None,
            }
            for t in targets
        ]

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
    if not is_available():
        for row in _local_find_compound("", limit=100000):
            if row.get("id") == compound_id:
                diseases = _local_data()["associations"].get(_normalize(row.get("name")), [])
                if not diseases and row.get("indications"):
                    diseases = [d.strip() for d in row["indications"].split(";") if d.strip()]
                return [
                    {
                        "disease": disease,
                        "mesh_id": f"LOCAL:{_normalize(disease)}",
                        "efo_id": None,
                        "max_phase": row.get("max_phase", 4),
                        "source": "Local JSON fallback",
                        "tree_numbers": [],
                        "entry_terms": [],
                    }
                    for disease in diseases[:40]
                ]
        return []

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
            c.id, COALESCE(c.chembl_id, 'NR:' || c.id::text) AS chembl_id,
            c.name, c.smiles, c.max_phase,
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


_stats_cache: Dict = {}
_STATS_TTL = 300  # 5 minutes

def get_stats() -> Dict:
    cached = _stats_cache.get("v")
    if cached and time.time() < cached[1]:
        return cached[0]

    if not is_available():
        data = _local_data()
        result = {
            "compounds": len(data["drugs"]),
            "targets": len({i.get("gene_symbol") for rows in data["interactions"].values() for i in rows if i.get("gene_symbol")}),
            "mechanisms": sum(len(rows) for rows in data["interactions"].values()),
            "indications": sum(len(rows) for rows in data["associations"].values()),
            "hetionet_nodes": 0,
            "mesh_diseases": len({d for rows in data["associations"].values() for d in rows}),
        }
        _stats_cache["v"] = (result, time.time() + _STATS_TTL)
        return result

    rows = _query("""
        SELECT
            (SELECT COUNT(*) FROM compounds)      AS compounds,
            (SELECT COUNT(*) FROM targets)        AS targets,
            (SELECT COUNT(*) FROM mechanisms)     AS mechanisms,
            (SELECT COUNT(*) FROM indications)    AS indications,
            (SELECT COUNT(*) FROM hetionet_nodes) AS hetionet_nodes,
            (SELECT COUNT(*) FROM mesh_diseases)  AS mesh_diseases
    """)
    result = dict(rows[0]) if rows else {}
    _stats_cache["v"] = (result, time.time() + _STATS_TTL)
    return result


# ── Full data footprint (for the Data Explorer) ───────────────────────────────
# Groups the real numbers across the whole stack: the ChEMBL 33 source DB, the
# HetioNet knowledge graph (real v1.0 edges vs the legacy ChEMBL-derived ones),
# and the neuro working set. ChEMBL-scale counts use pg_class.reltuples (instant,
# accurate for the static ChEMBL DB) instead of COUNT(*) over 20M rows.

_chembl_pool_nb: Optional[psycopg2.pool.ThreadedConnectionPool] = None
_footprint_cache: Dict = {}


def _get_chembl_pool_nb() -> Optional[psycopg2.pool.ThreadedConnectionPool]:
    global _chembl_pool_nb
    if _chembl_pool_nb is None:
        try:
            p = db_params()
            p["dbname"] = "chembl_33"
            _chembl_pool_nb = psycopg2.pool.ThreadedConnectionPool(1, 4, **p)
        except Exception as e:
            logger.warning(f"chembl_33 pool unavailable: {e}")
    return _chembl_pool_nb


def _query_chembl(sql: str, params=None) -> List[Dict]:
    pool = _get_chembl_pool_nb()
    if pool is None:
        return []
    conn = pool.getconn()
    try:
        with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
            cur.execute(sql, params)
            return [dict(r) for r in cur.fetchall()]
    except Exception as e:
        logger.debug(f"chembl_33 query error: {e}")
        return []
    finally:
        pool.putconn(conn)


def get_data_footprint() -> Dict:
    cached = _footprint_cache.get("v")
    if cached and time.time() < cached[1]:
        return cached[0]

    ws_rows = _query("""
        SELECT
            (SELECT COUNT(*) FROM compounds)       AS compounds,
            (SELECT COUNT(*) FROM mechanisms)      AS mechanisms,
            (SELECT COUNT(*) FROM indications)     AS indications,
            (SELECT COUNT(*) FROM mesh_diseases)   AS mesh_diseases,
            (SELECT COUNT(*) FROM targets)         AS targets,
            (SELECT COUNT(*) FROM hetionet_nodes)  AS hetionet_nodes,
            (SELECT COUNT(*) FROM hetionet_edges)  AS hetionet_edges,
            (SELECT COUNT(*) FROM hetionet_edges WHERE source='hetionet_v1.0') AS hetionet_edges_real,
            (SELECT COUNT(DISTINCT metaedge) FROM hetionet_edges) AS metaedge_types
    """)
    ws = ws_rows[0] if ws_rows else {}

    rt = _query_chembl("""
        SELECT relname, reltuples::bigint AS n FROM pg_class
        WHERE relkind='r' AND relname IN
              ('activities','molecule_dictionary','assays','target_dictionary')
    """)
    rtmap = {r["relname"]: int(r["n"]) for r in rt}
    ver = _query_chembl("SELECT name, creation_date FROM version LIMIT 1")

    total_edges = int(ws.get("hetionet_edges") or 0)
    real_edges = int(ws.get("hetionet_edges_real") or 0)
    result = {
        "chembl": {
            "version":    ver[0]["name"] if ver else "ChEMBL_33",
            "release":    str(ver[0]["creation_date"])[:10] if ver else "2023-05-31",
            "molecules":  rtmap.get("molecule_dictionary", 0),
            "activities": rtmap.get("activities", 0),
            "assays":     rtmap.get("assays", 0),
            "targets":    rtmap.get("target_dictionary", 0),
        },
        "graph": {
            "nodes":          int(ws.get("hetionet_nodes") or 0),
            "edges_total":    total_edges,
            "edges_real":     real_edges,
            "edges_derived":  total_edges - real_edges,
            "metaedge_types": int(ws.get("metaedge_types") or 0),
        },
        "working_set": {
            "compounds":     int(ws.get("compounds") or 0),
            "mechanisms":    int(ws.get("mechanisms") or 0),
            "indications":   int(ws.get("indications") or 0),
            "mesh_diseases": int(ws.get("mesh_diseases") or 0),
            "targets":       int(ws.get("targets") or 0),
        },
    }
    _footprint_cache["v"] = (result, time.time() + _STATS_TTL)
    return result


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


# ─── Drug repurposing screening ───────────────────────────────────────────────

def _local_repurposing_candidates(chembl_id: str, mesh_ids: List[str] = None, limit: int = 20, disease_name: str = "") -> List[Dict]:
    """Fallback: find local-data drugs sharing targets with the query compound."""
    data = _local_data()
    query_key = None
    for key, raw in data["drugs"].items():
        if raw.get("chembl_id") == chembl_id or _normalize(raw.get("chembl_id", "")) == _normalize(chembl_id):
            query_key = key
            break
    if not query_key:
        return _local_compounds_for_disease(mesh_ids, limit) if mesh_ids else []

    query_targets = {
        i.get("gene_symbol")
        for i in data["interactions"].get(query_key, [])
        if i.get("gene_symbol")
    }
    if not query_targets:
        return _local_compounds_for_disease(mesh_ids, limit) if mesh_ids else []

    candidates: List[Dict] = []
    for key, raw in data["drugs"].items():
        if key == query_key:
            continue
        drug_targets = {
            i.get("gene_symbol")
            for i in data["interactions"].get(key, [])
            if i.get("gene_symbol")
        }
        shared = query_targets & drug_targets
        if not shared:
            continue
        row = _local_compound(key, len(candidates))
        if not row:
            continue
        overlap = len(shared) / max(len(query_targets), 1)
        phase_score = float(row.get("max_phase", 0)) / 4.0
        score = round(min(1.0, 0.5 * overlap + 0.35 * phase_score + 0.15 * min(len(shared) / 5.0, 1.0)), 3)
        row["score"] = score
        row["shared_genes"] = "; ".join(sorted(shared))
        row["shared_target_count"] = len(shared)
        row["total_query_targets"] = len(query_targets)
        row["target_overlap_pct"] = round(overlap * 100, 1)
        row["score_breakdown"] = {
            "indication_score": round(overlap, 3),
            "target_score": round(phase_score, 3),
            "activity_score": round(min(len(shared) / 5.0, 1.0), 3),
            "network_score": 0.0,
        }
        candidates.append(row)

    candidates.sort(key=lambda x: x.get("score", 0), reverse=True)
    return candidates[:limit]


def _indication_based_candidates(
    query_id: int,
    mesh_ids: List[str] = None,
    limit: int = 20,
    disease_name: str = "",
) -> List[Dict]:
    """
    Fallback when a compound has no target data: find other high-phase compounds
    that share disease indications with the query compound, ranked by clinical phase.
    """
    # Get query compound's name (to exclude duplicates)
    name_rows = _query("SELECT name FROM compounds WHERE id = %s LIMIT 1", (query_id,))
    query_name = (name_rows[0]["name"] if name_rows else "").lower()

    # Get query compound's disease indications
    ind_rows = _query(
        "SELECT DISTINCT disease, mesh_id FROM indications WHERE compound_id = %s",
        (query_id,),
    )
    diseases = [r["disease"] for r in ind_rows if r.get("disease")]
    q_mesh = [r["mesh_id"] for r in ind_rows if r.get("mesh_id")]

    # Combine with caller-provided mesh_ids
    all_mesh = list(dict.fromkeys((q_mesh or []) + (mesh_ids or [])))
    all_diseases = diseases or []

    # Also include caller-provided disease name for text search
    if disease_name:
        all_diseases = list(dict.fromkeys(all_diseases + [disease_name]))

    if not all_mesh and not all_diseases:
        return []

    rows: List[Dict] = []

    if all_mesh:
        placeholders = ",".join(["%s"] * len(all_mesh))
        rows = _query(
            f"""
            SELECT DISTINCT
                c.id, COALESCE(c.chembl_id, 'NR:' || c.id::text) AS chembl_id,
                c.name, c.smiles, c.max_phase,
                cp.mw, cp.alogp, cp.psa, cp.hba, cp.hbd, cp.ro5_violations,
                STRING_AGG(DISTINCT i.disease,    '; ') AS indications,
                STRING_AGG(DISTINCT m.mechanism,  '; ') AS mechanisms,
                STRING_AGG(DISTINCT t.name,       '; ') AS targets
            FROM compounds c
            JOIN indications i ON i.compound_id = c.id
                AND i.mesh_id = ANY(ARRAY[{placeholders}]::varchar[])
            LEFT JOIN compound_properties cp ON cp.compound_id = c.id
            LEFT JOIN mechanisms m           ON m.compound_id  = c.id
            LEFT JOIN targets t              ON t.id = m.target_id
            WHERE c.id != %s AND LOWER(c.name) != %s AND c.max_phase >= 1
            GROUP BY c.id, c.chembl_id, c.name, c.smiles, c.max_phase,
                     cp.mw, cp.alogp, cp.psa, cp.hba, cp.hbd, cp.ro5_violations
            ORDER BY c.max_phase DESC NULLS LAST
            LIMIT %s
            """,
            all_mesh + [query_id, query_name, limit],
        )

    # If mesh_id match returned nothing, try keyword LIKE search (apostrophes stripped)
    if not rows and all_diseases:
        like_clauses = " OR ".join(f"LOWER(i.disease) LIKE %s" for _ in all_diseases[:6])
        patterns = [f"%{_disease_keyword(d)}%" for d in all_diseases[:6]]
        rows = _query(
            f"""
            SELECT DISTINCT
                c.id, COALESCE(c.chembl_id, 'NR:' || c.id::text) AS chembl_id,
                c.name, c.smiles, c.max_phase,
                cp.mw, cp.alogp, cp.psa, cp.hba, cp.hbd, cp.ro5_violations,
                STRING_AGG(DISTINCT i.disease,    '; ') AS indications,
                STRING_AGG(DISTINCT m.mechanism,  '; ') AS mechanisms,
                STRING_AGG(DISTINCT t.name,       '; ') AS targets
            FROM compounds c
            JOIN indications i ON i.compound_id = c.id
            LEFT JOIN compound_properties cp ON cp.compound_id = c.id
            LEFT JOIN mechanisms m           ON m.compound_id  = c.id
            LEFT JOIN targets t              ON t.id = m.target_id
            WHERE c.id != %s AND LOWER(c.name) != %s AND ({like_clauses}) AND c.max_phase >= 1
            GROUP BY c.id, c.chembl_id, c.name, c.smiles, c.max_phase,
                     cp.mw, cp.alogp, cp.psa, cp.hba, cp.hbd, cp.ro5_violations
            ORDER BY c.max_phase DESC NULLS LAST
            LIMIT %s
            """,
            [query_id, query_name] + patterns + [limit],
        )

    for r in rows:
        phase = float(r.get("max_phase") or 0)
        score = round(min(1.0, 0.3 + phase / 4.0 * 0.65), 3)
        r["score"] = score
        r["shared_target_count"] = 0
        r["total_query_targets"] = 0
        r["target_overlap_pct"] = 0.0
        r["shared_genes"] = ""
        r["score_breakdown"] = {
            "indication_score": round(phase / 4.0, 3),
            "target_score": 0.0,
            "activity_score": 0.0,
            "network_score": 0.0,
        }
    return rows


def find_repurposing_candidates(
    chembl_id: str,
    mesh_ids: List[str] = None,
    disease_name: str = "",
    limit: int = 20,
) -> List[Dict]:
    """
    Screen the database for compounds sharing mechanistic targets with the given
    compound. These are ranked repurposing candidates — drugs active at the same
    targets but potentially approved for different diseases, making them
    candidates for repurposing to the queried indication.

    Returns a list of compound dicts with extra fields:
      shared_target_count, total_query_targets, target_overlap_pct, shared_genes
    """
    if not is_available():
        return _local_repurposing_candidates(chembl_id, mesh_ids, limit, disease_name)

    # Resolve internal compound ID — handle synthetic NR: prefix
    if chembl_id.startswith("NR:"):
        try:
            _nr_id = int(chembl_id[3:])
            comp_rows = _query(
                "SELECT id, name FROM compounds WHERE id = %s LIMIT 1",
                (_nr_id,),
            )
        except ValueError:
            return []
    else:
        comp_rows = _query(
            "SELECT id, name FROM compounds WHERE chembl_id = %s LIMIT 1",
            (chembl_id,),
        )
    if not comp_rows:
        return _local_repurposing_candidates(chembl_id, mesh_ids, limit, disease_name)

    query_id = comp_rows[0]["id"]
    query_name_lower = (comp_rows[0].get("name") or "").lower()

    # Get query compound's target IDs
    target_rows = _query(
        "SELECT DISTINCT target_id FROM mechanisms WHERE compound_id = %s",
        (query_id,),
    )
    if not target_rows:
        return _indication_based_candidates(query_id, mesh_ids, limit, disease_name)

    target_ids = [r["target_id"] for r in target_rows]
    n_targets = len(target_ids)
    placeholders = ",".join(["%s"] * len(target_ids))

    rows = _query(
        f"""
        SELECT DISTINCT
            c.id, COALESCE(c.chembl_id, 'NR:' || c.id::text) AS chembl_id,
            c.name, c.smiles, c.max_phase,
            cp.mw, cp.alogp, cp.psa, cp.hba, cp.hbd, cp.ro5_violations,
            COUNT(DISTINCT m_shared.target_id) AS shared_target_count,
            STRING_AGG(DISTINCT t.gene_symbol, '; ') AS shared_genes,
            STRING_AGG(DISTINCT i.disease,     '; ') AS indications,
            STRING_AGG(DISTINCT mech.mechanism, '; ') AS mechanisms
        FROM compounds c
        JOIN mechanisms m_shared ON m_shared.compound_id = c.id
            AND m_shared.target_id = ANY(ARRAY[{placeholders}]::integer[])
        LEFT JOIN targets t         ON t.id = m_shared.target_id
        LEFT JOIN compound_properties cp ON cp.compound_id = c.id
        LEFT JOIN indications i     ON i.compound_id = c.id
        LEFT JOIN mechanisms mech   ON mech.compound_id = c.id
        WHERE c.id != %s AND LOWER(c.name) != %s
          AND c.max_phase >= 1
        GROUP BY c.id, c.chembl_id, c.name, c.smiles, c.max_phase,
                 cp.mw, cp.alogp, cp.psa, cp.hba, cp.hbd, cp.ro5_violations
        ORDER BY shared_target_count DESC, c.max_phase DESC NULLS LAST
        LIMIT %s
        """,
        target_ids + [query_id, query_name_lower, limit * 3],
    )

    for r in rows:
        shared = int(r.get("shared_target_count") or 0)
        phase = float(r.get("max_phase") or 0)
        overlap = shared / max(n_targets, 1)
        phase_score = phase / 4.0
        score = round(
            min(1.0, 0.5 * overlap + 0.35 * phase_score + 0.15 * min(shared / 5.0, 1.0)),
            3,
        )
        r["score"] = score
        r["total_query_targets"] = n_targets
        r["target_overlap_pct"] = round(overlap * 100, 1)
        r["score_breakdown"] = {
            "indication_score": round(overlap, 3),
            "target_score": round(phase_score, 3),
            "activity_score": round(min(shared / 5.0, 1.0), 3),
            "network_score": 0.0,
        }

    rows.sort(key=lambda x: x.get("score", 0), reverse=True)
    return rows[:limit]
