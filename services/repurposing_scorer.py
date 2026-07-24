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

from config import db_params  # centralized DB config (no hardcoded credentials)

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
            _pool = psycopg2.pool.ThreadedConnectionPool(1, 8, **db_params())   # headroom for the parallel reverse-screen scoring loop
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


# ─── D: Network evidence score (HetioNet) ─────────────────────────────────────
# Uses the REAL HetioNet v1.0 graph (source='hetionet_v1.0'), never the legacy
# ChEMBL-derived fallback edges, so this signal stays ORTHOGONAL to A/B/C.
#
# Evidence tiers, strongest first:
#   1.0        direct    Compound -treats-> Disease (CtD to THIS disease)
#   0.47–0.85  indirect  Compound -> Gene <- Disease: the drug perturbs genes
#              that the disease implicates (Compound binds/up/down-regulates a
#              Gene that the Disease associates/up/down-regulates). Scaled by the
#              number of shared genes. This is the multi-hop biology that direct
#              ChEMBL evidence cannot capture.
#   0.3        generic   Compound has CtD edges (treats some other disease)
#   ≤0.4       weak      Compound -binds-> Gene count proxy
HETIONET_SRC = "hetionet_v1.0"


def _resolve_hetionet_compound(compound_name: str) -> Optional[str]:
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
    return rows[0]["id"] if rows else None


@lru_cache(maxsize=256)
def _resolve_hetionet_diseases_cached(mesh_ids_key: str) -> tuple:
    mesh_ids = [m for m in mesh_ids_key.split("|") if m]
    ids: set = set()
    for pattern in _disease_like_patterns(mesh_ids):
        d_rows = _q(
            "SELECT id FROM hetionet_nodes "
            "WHERE kind='Disease' AND LOWER(name) LIKE %s LIMIT 5",
            (pattern,),
        )
        ids.update(r["id"] for r in d_rows)
    return tuple(ids)


def _resolve_hetionet_diseases(mesh_ids: List[str]) -> List[str]:
    # Cached on the disease key — a screen scores many compounds against the
    # SAME disease, so resolve its HetioNet nodes once.
    return list(_resolve_hetionet_diseases_cached("|".join(sorted(mesh_ids))))


def _network_evidence(compound_name: str, mesh_ids: List[str]) -> Dict:
    """Return {score, basis, genes} for the HetioNet network signal.

    Follows the indirect Compound→Gene←Disease path, not just direct CtD edges.
    """
    out: Dict = {"score": 0.0, "basis": "no compound node", "genes": []}

    c_node_id = _resolve_hetionet_compound(compound_name)
    if not c_node_id:
        return out

    disease_ids = _resolve_hetionet_diseases(mesh_ids) if mesh_ids else []

    if disease_ids:
        # 1. Direct: Compound -treats-> THIS disease (strongest)
        ctd = _q(
            "SELECT 1 FROM hetionet_edges "
            "WHERE source = %s AND metaedge = 'CtD' "
            "  AND source_id = %s AND target_id = ANY(%s) LIMIT 1",
            (HETIONET_SRC, c_node_id, disease_ids),
        )
        if ctd:
            out.update(score=1.0, basis="direct: treats this disease")
            return out

        # 2. Indirect: Compound → Gene ← Disease (shared disease genes)
        shared = _q(
            "SELECT DISTINCT hn.name AS gene "
            "FROM hetionet_edges cg "
            "JOIN hetionet_edges dg ON dg.target_id = cg.target_id "
            "JOIN hetionet_nodes hn ON hn.id = cg.target_id "
            "WHERE cg.source = %s AND cg.source_id = %s "
            "  AND cg.metaedge IN ('CbG','CuG','CdG') "
            "  AND dg.source = %s AND dg.source_id = ANY(%s) "
            "  AND dg.metaedge IN ('DaG','DuG','DdG') "
            "LIMIT 25",
            (HETIONET_SRC, c_node_id, HETIONET_SRC, disease_ids),
        )
        genes = [r["gene"] for r in shared if r.get("gene")]
        if genes:
            n = len(genes)
            out.update(
                score=min(0.85, 0.35 + 0.12 * n),
                basis=f"indirect: shares {n} disease gene(s) via the knowledge graph",
                genes=genes[:8],
            )
            return out

    # 3. Generic fallback: treats SOME disease
    ctd_any = _q(
        "SELECT 1 FROM hetionet_edges "
        "WHERE source = %s AND metaedge = 'CtD' AND source_id = %s LIMIT 1",
        (HETIONET_SRC, c_node_id),
    )
    if ctd_any:
        out.update(score=0.3, basis="generic: treats other disease(s)")
        return out

    # 4. Weak: binds genes (no disease link resolved)
    cbg = _q(
        "SELECT COUNT(*) AS cnt FROM hetionet_edges "
        "WHERE source = %s AND metaedge = 'CbG' AND source_id = %s",
        (HETIONET_SRC, c_node_id),
    )
    cbg_cnt = int(cbg[0]["cnt"]) if cbg else 0
    if cbg_cnt > 0:
        out.update(score=min(0.4, cbg_cnt / 20.0), basis="weak: binds genes")
    return out


def _score_network(compound_name: str, mesh_ids: List[str]) -> float:
    return _network_evidence(compound_name, mesh_ids)["score"]


def hetionet_novel_compounds(disease_terms: List[str], limit: int = 20) -> List[Dict]:
    """Compounds connected to the disease through shared genes in the REAL HetioNet
    graph (Compound→Gene←Disease), ranked by the number of shared disease genes.

    These are candidates for NOVEL repurposing — drugs the biology links to the
    disease even when they were never indicated for it. Returns [{name, shared_genes}].
    """
    disease_ids = _resolve_hetionet_diseases(disease_terms) if disease_terms else []
    if not disease_ids:
        return []
    rows = _q(
        """
        WITH dgenes AS (
            SELECT DISTINCT he.target_id AS gene
            FROM hetionet_edges he
            WHERE he.source = %s AND he.metaedge IN ('DaG','DuG','DdG')
              AND he.source_id = ANY(%s)
        ), cmpd AS (
            SELECT he.source_id AS cnode, COUNT(DISTINCT he.target_id) AS shared
            FROM hetionet_edges he JOIN dgenes ON he.target_id = dgenes.gene
            WHERE he.source = %s AND he.metaedge IN ('CbG','CuG','CdG')
            GROUP BY he.source_id
        )
        SELECT hn.name AS name, c.shared AS shared
        FROM cmpd c JOIN hetionet_nodes hn ON hn.id = c.cnode
        WHERE hn.kind = 'Compound' AND c.shared >= 1
        ORDER BY c.shared DESC, hn.name LIMIT %s
        """,
        (HETIONET_SRC, disease_ids, HETIONET_SRC, limit),
    )
    return [{"name": r["name"], "shared_genes": int(r["shared"])} for r in rows]


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
    net   = _network_evidence(compound_name, mesh_ids)
    s_net = net["score"]

    overall = (
        W_INDICATION * s_ind +
        W_TARGET     * s_tgt +
        W_ACTIVITY   * s_act +
        W_NETWORK    * s_net
    )

    phase_bonus = min(0.05, (compound_phase or 0) / 4.0 * 0.05)
    overall = min(1.0, overall + phase_bonus)

    breakdown = {
        "overall":           round(overall, 4),
        "indication_score":  round(s_ind, 4),
        "target_score":      round(s_tgt, 4),
        "activity_score":    round(s_act, 4),
        "network_score":     round(s_net, 4),
        # How the network signal was earned + the genes linking drug→disease
        # (empty unless an indirect Compound→Gene←Disease path was found).
        "network_basis":     net["basis"],
        "network_genes":     net["genes"],
        "phase_bonus":       round(phase_bonus, 4),
        # True only when this molecule is APPROVED (phase 4) for THIS disease —
        # the one genuine novelty disqualifier for a 505(b)(2) repurposing claim.
        # Prior development at a lower phase is NOT flagged (it is de-risking, not
        # a blocker). s_ind reaches 1.0 only when a phase-4 indication matched.
        "approved_here":     bool(s_ind >= 0.999),
    }

    # Asset-Scoring step: predicted likelihood of progressing through each
    # remaining development phase. The overall relevance score is the efficacy
    # evidence for the new indication; an already-approved molecule (phase 4) is
    # treated as a repurposing candidate (Phase-1 safety de-risked). The molecule
    # feature profile feeds the fitted calibrator when it is present. Fail-soft.
    try:
        from services.pos_model import predict_progression
        breakdown["pos"] = predict_progression(
            current_phase=compound_phase,
            evidence_score=overall,
            is_repurposing=(compound_phase or 0) >= 4,
            features=_molecule_features(compound_id),
        )
    except Exception as e:
        logger.debug(f"PoS prediction skipped: {e}")

    return breakdown


@lru_cache(maxsize=4096)
def _molecule_features(compound_id: int) -> Dict:
    """Disease-independent biological-profile features for the PoS calibrator —
    the same columns services/train_pos_model.py fits on. One indexed lookup,
    cached. Returns {} if the DB is unavailable so PoS falls back to analytic."""
    rows = _q(
        """
        SELECT COUNT(DISTINCT m.target_id)        AS n_targets,
               COUNT(DISTINCT m.id)               AS n_mechanisms,
               COUNT(DISTINCT i.id)               AS n_indications,
               COALESCE(MAX(a.pchembl_value), 0)  AS max_pchembl
        FROM compounds c
        LEFT JOIN mechanisms  m         ON m.compound_id = c.id
        LEFT JOIN indications i         ON i.compound_id = c.id
        LEFT JOIN compound_activities a ON a.compound_id = c.id
        WHERE c.id = %s
        """,
        (compound_id,),
    )
    if not rows:
        return {}
    r = rows[0]
    return {
        "n_targets":     float(r.get("n_targets") or 0),
        "n_mechanisms":  float(r.get("n_mechanisms") or 0),
        "n_indications": float(r.get("n_indications") or 0),
        "max_pchembl":   float(r.get("max_pchembl") or 0),
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
        c["approved_here"] = breakdown.get("approved_here", False)

    compounds.sort(key=lambda x: x["score"], reverse=True)
    return compounds


def approved_chembls_for_disease(chembl_ids: List[str], mesh_ids: List[str],
                                 disease_name: str = "", min_phase: int = 4) -> set:
    """Of the given molecules, which have a development record at >= `min_phase` for
    THIS disease? One batched query keyed by chembl_id, so it works regardless of which
    scoring path produced the candidates (DB scorer or the repurposing engine).

    min_phase=4 (default) → APPROVED, the single genuine novelty disqualifier for a
    505(b)(2) repurposing claim. min_phase=1 → any development footprint, used as the
    confounding-by-indication guard (a drug DEVELOPED for a disease will have that
    disease in its FAERS, so its adverse-event overlap must not be read as toxicity)."""
    chembl_ids = [c for c in (chembl_ids or []) if c]
    if not chembl_ids:
        return set()
    real_ids = [m for m in (mesh_ids or []) if _is_mesh_id(m)]
    patterns = set()
    kw = _keyword(disease_name or "")
    if kw:
        patterns.add(f"%{kw}%")
    for m in (mesh_ids or []):
        if not _is_mesh_id(m):
            k = _keyword(m)
            if k:
                patterns.add(f"%{k}%")

    where, params = [], [chembl_ids]
    if real_ids:
        where.append("i.mesh_id = ANY(%s::varchar[])"); params.append(real_ids)
    if patterns:
        where.append("LOWER(i.disease) LIKE ANY(%s::text[])"); params.append(list(patterns))
    if not where:
        return set()
    sql = ("SELECT DISTINCT c.chembl_id FROM indications i "
           "JOIN compounds c ON c.id = i.compound_id "
           "WHERE c.chembl_id = ANY(%s::varchar[]) AND i.max_phase >= %s "
           "AND (" + " OR ".join(where) + ")")
    params.insert(1, int(min_phase))   # placeholder order: chembl_ids, min_phase, where…
    rows = _q(sql, tuple(params))
    return {r["chembl_id"] for r in rows if r.get("chembl_id")}
