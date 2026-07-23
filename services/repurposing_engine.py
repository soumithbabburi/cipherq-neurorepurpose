"""
5-Screen Drug Repurposing Engine — Solix EAI Framework

Scoring dimensions (weights):
  1. Target-Based     25%  Jaccard(drug targets, disease gene set)
  2. Pathway-Based    20%  Shared Reactome pathways (via Open Targets)
  3. PPI Network      20%  STRING network proximity
  4. Clinical Signal  15%  Phase / approval / indication evidence
  5. Indication Sem.  10%  Existing indication ∩ disease keywords
  6. 505(b)(2) Reg.   10%  Regulatory feasibility (orphan + safety data)

Results cached in data/repurposing_cache.json (TTL 6h).
"""

import json
import logging
import re
import time
from pathlib import Path
from typing import Dict, List, Optional

import requests  # noqa: F401

import psycopg2
import psycopg2.pool

from services import http_client
from config import db_params

logger = logging.getLogger(__name__)
CACHE_FILE = Path(__file__).parent.parent / "data" / "repurposing_cache.json"

WEIGHTS = {
    "target":     0.25,
    "pathway":    0.20,
    "ppi":        0.20,
    "clinical":   0.15,
    "indication": 0.10,
    "regulatory": 0.10,
}

CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"
CACHE_TTL   = 21600  # 6 hours


# ── Cache helpers ─────────────────────────────────────────────────────────────

def _load_cache() -> dict:
    try:
        if CACHE_FILE.exists():
            data = json.loads(CACHE_FILE.read_text(encoding="utf-8"))
            now = time.time()
            return {k: v for k, v in data.items() if now - v.get("_ts", 0) < CACHE_TTL}
    except Exception:
        pass
    return {}


def _save_cache(cache: dict):
    try:
        CACHE_FILE.parent.mkdir(exist_ok=True)
        CACHE_FILE.write_text(json.dumps(cache), encoding="utf-8")
    except Exception:
        pass


# ── ChEMBL REST API helpers ───────────────────────────────────────────────────

_chembl_pool: Optional[psycopg2.pool.ThreadedConnectionPool] = None


def _get_chembl_pool() -> Optional[psycopg2.pool.ThreadedConnectionPool]:
    global _chembl_pool
    if _chembl_pool is None:
        try:
            params = db_params()
            params["dbname"] = "chembl_33"   # indications live in chembl_33, not neurorepurpose
            _chembl_pool = psycopg2.pool.ThreadedConnectionPool(1, 4, **params)
        except Exception as e:
            logger.warning(f"chembl_33 pool unavailable (will fall back to REST): {e}")
    return _chembl_pool


def _local_drug_indications(disease_name: str, limit: int) -> Optional[List[Dict]]:
    """Fast local indication lookup from the validated chembl_33 DB (~0.2s vs the
    REST API's ~80s). Returns None (not []) when the local DB is unavailable, so
    the caller falls back to the REST API."""
    pool = _get_chembl_pool()
    if pool is None:
        return None
    conn = None
    try:
        conn = pool.getconn()
        with conn.cursor() as cur:
            cur.execute(
                """
                SELECT md.chembl_id,
                       COALESCE(MAX(di.mesh_heading), MAX(di.efo_term)) AS indication,
                       MAX(di.max_phase_for_ind) AS phase
                FROM drug_indication di
                JOIN molecule_dictionary md ON md.molregno = di.molregno
                WHERE di.mesh_heading ILIKE %s OR di.efo_term ILIKE %s
                GROUP BY md.chembl_id
                ORDER BY phase DESC NULLS LAST
                LIMIT %s
                """,
                (f"%{disease_name}%", f"%{disease_name}%", limit),
            )
            rows = cur.fetchall()
        return [
            {"chembl_id": r[0],
             "indication": r[1] or disease_name,
             "max_phase_for_ind": float(r[2]) if r[2] is not None else 0}
            for r in rows
        ]
    except Exception as e:
        logger.debug(f"local drug_indication query failed: {e}")
        return None
    finally:
        if conn is not None:
            pool.putconn(conn)


def _chembl_drug_indications(disease_name: str, limit: int = 50) -> List[Dict]:
    """Compounds indicated for a disease.

    Fast path: the local chembl_33 `drug_indication` table (~0.2s). Falls back to
    the ChEMBL REST API only when the local DB is unreachable. The two REST field
    queries (MeSH heading + EFO term) run in PARALLEL with a tight timeout and a
    single retry; returns [] if ChEMBL is unreachable so the caller degrades to a
    clean "no candidates" result rather than hanging.
    """
    # Fast path — local validated DB (avoids the slow ChEMBL REST __icontains call).
    local = _local_drug_indications(disease_name, limit)
    if local is not None:
        return local

    from concurrent.futures import ThreadPoolExecutor

    def _fetch(param: str) -> List[Dict]:
        try:
            r = http_client.get(
                f"{CHEMBL_BASE}/drug_indication.json",
                params={param: disease_name, "limit": limit, "format": "json"},
                timeout=8, retries=1,        # fail fast — don't stack 12s x 3 retries
            )
            if r and r.ok:
                return r.json().get("drug_indications", [])
        except Exception as e:
            logger.debug(f"ChEMBL {param} error: {e}")
        return []

    params = ("mesh_heading__icontains", "efo_term__icontains")
    seen: set = set()
    results: List[Dict] = []
    with ThreadPoolExecutor(max_workers=len(params)) as pool:
        for inds in pool.map(_fetch, params):
            for ind in inds:
                mol_id = ind.get("molecule_chembl_id", "")
                if mol_id and mol_id not in seen:
                    seen.add(mol_id)
                    results.append({
                        "chembl_id":         mol_id,
                        "indication":        ind.get("mesh_heading") or ind.get("efo_term") or disease_name,
                        "max_phase_for_ind": ind.get("max_phase_for_ind") or 0,
                    })
            if len(results) >= limit:
                break

    return results


def _local_molecule_details(chembl_ids: List[str]) -> Optional[Dict[str, Dict]]:
    """Molecule details from local chembl_33. None if DB unavailable."""
    pool = _get_chembl_pool()
    if pool is None:
        return None
    ids = [c for c in chembl_ids if c]
    if not ids:
        return {}
    conn = None
    try:
        conn = pool.getconn()
        with conn.cursor() as cur:
            cur.execute(
                """
                SELECT md.chembl_id, md.pref_name, md.max_phase,
                       cs.canonical_smiles, cp.mw_freebase, cp.alogp, cp.psa,
                       cp.hba, cp.hbd, cp.num_ro5_violations
                FROM molecule_dictionary md
                LEFT JOIN compound_structures cs ON cs.molregno = md.molregno
                LEFT JOIN compound_properties cp ON cp.molregno = md.molregno
                WHERE md.chembl_id = ANY(%s)
                """,
                (ids,),
            )
            rows = cur.fetchall()
        out: Dict[str, Dict] = {}
        for r in rows:
            out[r[0]] = {
                "chembl_id": r[0],
                "name": r[1] or r[0],
                "smiles": r[3] or "",
                "max_phase": float(r[2]) if r[2] is not None else 0,
                "mw":    float(r[4]) if r[4] is not None else None,
                "alogp": float(r[5]) if r[5] is not None else None,
                "psa":   float(r[6]) if r[6] is not None else None,
                "hba":   int(r[7]) if r[7] is not None else None,
                "hbd":   int(r[8]) if r[8] is not None else None,
                "ro5_violations": int(r[9]) if r[9] is not None else None,
                "indications": "", "mechanisms": "", "targets": "",
                "source": "chembl_33 (local)",
            }
        return out
    except Exception as e:
        logger.debug(f"local molecule details failed: {e}")
        return None
    finally:
        if conn is not None:
            pool.putconn(conn)


def _local_targets_for_molecules(chembl_ids: List[str]) -> Optional[Dict[str, List[str]]]:
    """Target gene symbols per molecule from local chembl_33. None if DB unavailable.

    Coverage strategy (validated on repoDB — fallback mode lifts drug coverage
    46%→68% and AUROC 0.71→0.73 vs mechanism-only; see validate_predictions.py):
      1. Curated mechanism targets (drug_mechanism), folding salt→parent so a salt
         inherits the parent's targets.
      2. ONLY for molecules still with no target, fall back to high-confidence
         single-protein bioactivity (pchembl≥6, confidence≥8). Gap-fill only —
         well-annotated drugs keep their clean curated targets (broadening every
         drug with activity targets diluted precision in testing)."""
    pool = _get_chembl_pool()
    if pool is None:
        return None
    ids = [c for c in chembl_ids if c]
    if not ids:
        return {}
    drug_genes: Dict[str, set] = {cid: set() for cid in ids}
    conn = None
    try:
        conn = pool.getconn()
        with conn.cursor() as cur:
            # chembl_id → own + parent molregno (salt→parent)
            cur.execute(
                """
                SELECT md.chembl_id, md.molregno, COALESCE(mh.parent_molregno, md.molregno)
                FROM molecule_dictionary md
                LEFT JOIN molecule_hierarchy mh ON mh.molregno = md.molregno
                WHERE md.chembl_id = ANY(%s)
                """,
                (ids,),
            )
            mol_to_cids: Dict[int, set] = {}
            all_mols: set = set()
            for chembl, mol, parent in cur.fetchall():
                for m in (mol, parent):
                    if m is not None:
                        mol_to_cids.setdefault(m, set()).add(chembl)
                        all_mols.add(m)
            if not all_mols:
                return {cid: [] for cid in ids}
            mol_list = list(all_mols)

            # 1. curated mechanism targets
            cur.execute(
                """
                SELECT dm.molregno, csyn.component_synonym
                FROM drug_mechanism dm
                JOIN target_components tc ON tc.tid = dm.tid
                JOIN component_synonyms csyn ON csyn.component_id = tc.component_id
                                            AND csyn.syn_type = 'GENE_SYMBOL'
                WHERE dm.molregno = ANY(%s)
                """,
                (mol_list,),
            )
            for mol, gene in cur.fetchall():
                if gene:
                    for cid in mol_to_cids.get(mol, ()):
                        drug_genes[cid].add(gene)

            # 2. gap-fill from high-confidence bioactivity — but SELECTIVITY-FILTERED.
            #    A promiscuous kinase inhibitor binds dozens of kinases at ≤1µM
            #    (pchembl≥6) that it never engages therapeutically; crediting those
            #    off-targets lets any of them force-match a disease (Erlotinib's
            #    AURKB/JAK3/PDGFRB → phantom leads). So we take the per-gene MAX
            #    potency and keep only targets within ~10× (1.0 log) of the drug's
            #    MOST-potent target — its genuinely selective, therapeutic set.
            gap_cids = {cid for cid in ids if not drug_genes[cid]}
            if gap_cids:
                gap_mols = [m for m, cids in mol_to_cids.items() if cids & gap_cids]
                if gap_mols:
                    cur.execute(
                        """
                        SELECT a.molregno, csyn.component_synonym, MAX(a.pchembl_value)
                        FROM activities a
                        JOIN assays ass            ON ass.assay_id = a.assay_id
                        JOIN target_dictionary td  ON td.tid = ass.tid
                                                  AND td.target_type = 'SINGLE PROTEIN'
                        JOIN target_components tc  ON tc.tid = ass.tid
                        JOIN component_synonyms csyn ON csyn.component_id = tc.component_id
                                                    AND csyn.syn_type = 'GENE_SYMBOL'
                        WHERE a.molregno = ANY(%s)
                          AND a.pchembl_value >= 6
                          AND a.standard_relation = '='
                          AND ass.confidence_score >= 8
                        GROUP BY a.molregno, csyn.component_synonym
                        """,
                        (gap_mols,),
                    )
                    # collect per-molecule (gene, max_pchembl), then keep the
                    # selective window relative to each molecule's best target.
                    by_mol: Dict[int, list] = {}
                    for mol, gene, pot in cur.fetchall():
                        if gene and pot is not None:
                            by_mol.setdefault(mol, []).append((gene, float(pot)))
                    for mol, rows in by_mol.items():
                        best = max(p for _, p in rows)
                        floor = max(6.5, best - 1.0)          # within 10× of the most potent
                        for gene, pot in rows:
                            if pot >= floor:
                                for cid in (mol_to_cids.get(mol, set()) & gap_cids):
                                    drug_genes[cid].add(gene)
        return {cid: sorted(v) for cid, v in drug_genes.items()}
    except Exception as e:
        logger.debug(f"local targets failed: {e}")
        return None
    finally:
        if conn is not None:
            pool.putconn(conn)


def _actions_for_molecules(chembl_ids: List[str]) -> Dict[str, Dict[str, str]]:
    """Map each molecule → {GENE_SYMBOL: action_type} from local chembl_33
    drug_mechanism. Powers the direction-aware pathway score. Returns {} when the
    local DB is unavailable, so scoring falls back to direction-blind coverage."""
    pool = _get_chembl_pool()
    ids = [c for c in chembl_ids if c]
    if pool is None or not ids:
        return {}
    out: Dict[str, Dict[str, str]] = {cid: {} for cid in ids}
    conn = None
    try:
        conn = pool.getconn()
        with conn.cursor() as cur:
            cur.execute(
                """
                SELECT md.chembl_id, UPPER(csyn.component_synonym), dm.action_type
                FROM molecule_dictionary md
                JOIN drug_mechanism dm ON dm.molregno = md.molregno
                JOIN target_components tc ON tc.tid = dm.tid
                JOIN component_synonyms csyn ON csyn.component_id = tc.component_id
                                            AND csyn.syn_type = 'GENE_SYMBOL'
                WHERE md.chembl_id = ANY(%s) AND dm.action_type IS NOT NULL
                """,
                (ids,),
            )
            for cid, gene, action in cur.fetchall():
                if cid and gene and action:
                    out.setdefault(cid, {}).setdefault(gene, action)
        return out
    except Exception as e:
        logger.debug(f"local actions failed: {e}")
        return {}
    finally:
        if conn is not None:
            pool.putconn(conn)


def _local_molecules_by_name(names: List[str]) -> Dict[str, Dict]:
    """Map lowercased compound names -> ChEMBL molecule dicts (local chembl_33).
    Used to turn HetioNet novelty compound names into scorable candidates."""
    pool = _get_chembl_pool()
    names = [n for n in names if n]
    if pool is None or not names:
        return {}
    conn = None
    try:
        conn = pool.getconn()
        with conn.cursor() as cur:
            cur.execute(
                """
                SELECT LOWER(md.pref_name), md.chembl_id, md.pref_name, md.max_phase,
                       cs.canonical_smiles, cp.mw_freebase, cp.alogp, cp.psa
                FROM molecule_dictionary md
                LEFT JOIN compound_structures cs ON cs.molregno = md.molregno
                LEFT JOIN compound_properties cp ON cp.molregno = md.molregno
                WHERE LOWER(md.pref_name) = ANY(%s)
                """,
                ([n.lower() for n in names],),
            )
            out: Dict[str, Dict] = {}
            for r in cur.fetchall():
                out[r[0]] = {
                    "chembl_id": r[1], "name": r[2] or r[1],
                    "max_phase": float(r[3]) if r[3] is not None else 0,
                    "smiles": r[4] or "",
                    "mw": float(r[5]) if r[5] is not None else None,
                    "alogp": float(r[6]) if r[6] is not None else None,
                    "psa": float(r[7]) if r[7] is not None else None,
                    "indications": "", "mechanisms": "", "targets": "",
                }
            return out
    except Exception as e:
        logger.debug(f"local molecules-by-name failed: {e}")
        return {}
    finally:
        if conn is not None:
            pool.putconn(conn)


def _chembl_molecule_details(chembl_ids: List[str]) -> Dict[str, Dict]:
    """Molecule details — local chembl_33 first, ChEMBL REST API as fallback."""
    local = _local_molecule_details(chembl_ids)
    if local is not None:
        return local
    results: Dict[str, Dict] = {}
    for i in range(0, min(len(chembl_ids), 60), 20):
        chunk = chembl_ids[i : i + 20]
        try:
            r = http_client.get(
                f"{CHEMBL_BASE}/molecule.json",
                params={"molecule_chembl_id__in": ",".join(chunk), "limit": 20, "format": "json"},
                timeout=12,
            )
            if r and r.ok:
                for mol in r.json().get("molecules", []):
                    mid = mol.get("molecule_chembl_id", "")
                    props = mol.get("molecule_properties") or {}
                    struct = mol.get("molecule_structures") or {}
                    results[mid] = {
                        "chembl_id": mid,
                        "name":  mol.get("pref_name") or mid,
                        "smiles": struct.get("canonical_smiles", ""),
                        "max_phase": mol.get("max_phase") or 0,
                        "mw":    props.get("mw_freebase"),
                        "alogp": props.get("alogp"),
                        "psa":   props.get("psa"),
                        "hba":   props.get("hba"),
                        "hbd":   props.get("hbd"),
                        "ro5_violations": props.get("num_ro5_violations"),
                        "indications": "",
                        "mechanisms": "",
                        "targets": "",
                        "source": "ChEMBL API",
                    }
        except Exception as e:
            logger.debug(f"ChEMBL molecule detail error: {e}")
    return results


def _chembl_targets_for_molecules(chembl_ids: List[str]) -> Dict[str, List[str]]:
    """Map each ChEMBL molecule ID to its target gene symbols.

    Batched in two passes so the whole candidate set is covered with a handful of
    requests, instead of the old 1 + (mechanisms) sequential per-compound calls
    that timed out and left most candidates with no genes (→ target/pathway = 0):
      1. mechanism.json?molecule_chembl_id__in=…  → {molecule: {target_chembl_id}}
      2. target.json?target_chembl_id__in=…       → {target_chembl_id: [gene_symbol]}
    Then join the two maps. Uses __in batching (the same filter the molecule
    fetch already relies on), so it stays robust for 40-50 candidates.
    """
    # Fast path — local validated DB.
    local = _local_targets_for_molecules(chembl_ids)
    if local is not None:
        return local

    ids = [c for c in chembl_ids if c]
    drug_genes: Dict[str, List[str]] = {cid: [] for cid in ids}
    if not ids:
        return drug_genes

    # ── Pass 1: molecules → target ChEMBL ids ─────────────────────────────────
    mol_targets: Dict[str, set] = {}
    all_target_ids: set = set()
    for i in range(0, len(ids), 20):
        chunk = ids[i : i + 20]
        try:
            r = http_client.get(
                f"{CHEMBL_BASE}/mechanism.json",
                params={"molecule_chembl_id__in": ",".join(chunk),
                        "limit": 1000, "format": "json"},
                timeout=15,
            )
            if not r or not r.ok:
                continue
            for mech in r.json().get("mechanisms", []):
                mid = mech.get("molecule_chembl_id", "")
                tid = mech.get("target_chembl_id", "")
                if mid and tid:
                    mol_targets.setdefault(mid, set()).add(tid)
                    all_target_ids.add(tid)
        except Exception as e:
            logger.debug(f"mechanism batch error: {e}")

    if not all_target_ids:
        return drug_genes

    # ── Pass 2: target ChEMBL ids → gene symbols ──────────────────────────────
    target_genes: Dict[str, List[str]] = {}
    tids = list(all_target_ids)
    for i in range(0, len(tids), 20):
        chunk = tids[i : i + 20]
        try:
            r = http_client.get(
                f"{CHEMBL_BASE}/target.json",
                params={"target_chembl_id__in": ",".join(chunk),
                        "limit": 1000, "format": "json"},
                timeout=15,
            )
            if not r or not r.ok:
                continue
            for t in r.json().get("targets", []):
                tid = t.get("target_chembl_id", "")
                genes: List[str] = []
                for comp in (t.get("target_components") or []):
                    for syn in (comp.get("target_component_synonyms") or []):
                        if syn.get("syn_type") == "GENE_SYMBOL" and syn.get("component_synonym"):
                            genes.append(syn["component_synonym"])
                if tid:
                    target_genes[tid] = genes
        except Exception as e:
            logger.debug(f"target batch error: {e}")

    # ── Join ──────────────────────────────────────────────────────────────────
    for cid in ids:
        joined: set = set()
        for tid in mol_targets.get(cid, ()):
            joined.update(target_genes.get(tid, []))
        drug_genes[cid] = sorted(joined)

    # ── Fallback for molecules ChEMBL has no curated mechanism row for ────────
    # Drugs like ropinirole / bromocriptine have NO `mechanism` record under their
    # ChEMBL id, yet well-known targets (DRD2/3/4). resolve_drug unions ChEMBL
    # molecule forms (salts/parents) + activity to recover the gene symbols the
    # mechanism endpoint misses, so they no longer score target = 0. Only the gaps
    # are resolved, in parallel, so the cost stays bounded.
    missing = [cid for cid in ids if not drug_genes.get(cid)]
    if missing:
        try:
            from concurrent.futures import ThreadPoolExecutor
            from services.reverse_repurposing import resolve_drug

            def _fallback(cid: str):
                try:
                    info = resolve_drug(cid) or {}
                    return cid, [g for g in (info.get("targets") or []) if g]
                except Exception:
                    return cid, []

            with ThreadPoolExecutor(max_workers=6) as pool:
                for cid, genes in pool.map(_fallback, missing):
                    if genes:
                        drug_genes[cid] = sorted({g.upper() for g in genes})
        except Exception as e:
            logger.debug(f"resolve_drug target fallback failed: {e}")

    return drug_genes


# ── Scoring functions ─────────────────────────────────────────────────────────

def _score_target_overlap(drug_genes: List[str], disease_genes: List[str]) -> float:
    """Drug-centric target overlap: what fraction of the DRUG's targets are disease
    genes (recall), with a bonus when those targets rank among the disease's top genes.
    This replaces a Jaccard index, which structurally collapsed toward zero because a
    drug's handful of targets was divided by the disease's large gene set."""
    if not drug_genes or not disease_genes:
        return 0.0
    a = {g.upper() for g in drug_genes}
    ranked = [g.upper() for g in disease_genes]
    bset = set(ranked)
    overlap = a & bset
    if not overlap:
        return 0.0
    recall = len(overlap) / len(a)                       # fraction of drug targets that hit
    top = set(ranked[:10])
    in_top = bool(a & top)
    rank_bonus = 0.3 * (len(a & top) / len(a))            # weight high-ranked disease genes
    hit_floor = 0.25 if in_top else 0.0                  # a top-gene hit is meaningful on its own
    return min(1.0, 0.6 * recall + rank_bonus + hit_floor)


def _score_target_overlap_weighted(drug_genes: List[str],
                                   disease_weights: Dict[str, float]) -> float:
    """Genetics-weighted target overlap. disease_weights: {GENE_UPPER:
    genetic_association_score 0..1}. Rewards a drug for hitting GENETICALLY
    supported disease genes (strength), for hitting disease genes at all (recall),
    and for hitting several (breadth). Reduces to a sensible overlap score; used
    when genetics weights are available, else the engine falls back to the flat
    _score_target_overlap."""
    if not drug_genes or not disease_weights:
        return 0.0
    a = {g.upper() for g in drug_genes}
    hits = [disease_weights[g] for g in a if g in disease_weights]
    if not hits:
        return 0.0
    # STRENGTH-dominated. The old form was recall-dominated, so a SELECTIVE drug
    # hitting ONE weak-association gene got recall=1.0 → ~0.55 regardless of how
    # weak the link was (Erlotinib→EGFR@0.35 for Cowden read like a real hit). Now
    # the association STRENGTH gates the score: a weak best-hit caps it low, while
    # hitting several genetically-supported drivers (weighted breadth) lifts it.
    strength = max(hits)                                  # strongest disease-gene association hit
    weighted_breadth = 1.0 - pow(2.71828, -sum(hits) / 1.0)   # saturating Σ of hit strengths
    recall = len(hits) / len(a)                           # selectivity toward the disease (minor)
    # strength sets the ceiling; breadth fills it in; recall a small bonus.
    return min(1.0, round(strength * (0.65 + 0.35 * weighted_breadth) + 0.10 * recall, 4))


def _score_pathway_overlap(drug_genes: List[str], disease_pathways: List[Dict],
                           drug_signature: Optional[Dict[str, float]] = None,
                           gene_damp: Optional[Dict[str, float]] = None) -> float:
    """Fraction of the disease's pathways the drug touches.

    When `drug_signature` (gene→signed action, from signature_engine.drug_signature)
    is supplied, the raw coverage is modulated by DIRECTION: disease pathway genes
    are treated as elevated in disease (the same documented proxy the signature
    engine uses), so a drug that INHIBITS them reverses the disease (favourable)
    while one that ACTIVATES them mimics it (unfavourable). Without a signature the
    behaviour is identical to the original direction-blind coverage score."""
    if not drug_genes or not disease_pathways:
        return 0.0
    dg = {g.upper() for g in drug_genes}
    # Degree-damped hit count: a pathway touched only through a HUB gene (EGFR…)
    # counts for less than one touched through a specific gene (DWPC, Himmelstein).
    hits = 0.0
    for pw in disease_pathways[:20]:
        touching = dg & {g.upper() for g in pw.get("genes", [])}
        if touching:
            hits += (max((gene_damp.get(g, 1.0) for g in touching), default=1.0)
                     if gene_damp else 1.0)
    base = min(1.0, hits / max(1, min(10, len(disease_pathways))))
    if not base or not drug_signature:
        return base
    try:
        from services.signature_engine import disease_signature, connectivity
        pw_genes: set = set()
        for pw in disease_pathways[:20]:
            pw_genes.update(g.upper() for g in pw.get("genes", []))
        conn = connectivity(disease_signature(list(pw_genes)), drug_signature)
        if conn["n_shared"] == 0:
            return base
        c = conn["connectivity"]
        if c > 0.1:                       # mimics disease direction → unfavourable
            return round(base * max(0.2, 1.0 - c), 4)
        # reversing / neutral → favourable; small boost for strong reversal
        return round(min(1.0, base * (1.0 + 0.3 * conn["reversal_score"])), 4)
    except Exception as e:
        logger.debug(f"pathway direction modifier skipped: {e}")
        return base


def _score_ppi_network(drug_genes: List[str], ppi_adj: Dict[str, List[str]],
                       gene_damp: Optional[Dict[str, float]] = None) -> float:
    if not drug_genes or not ppi_adj:
        return 0.0
    dg = {g.upper() for g in drug_genes}
    ppi_nodes = {g.upper() for g in ppi_adj}
    ppi_neighbors: set = set()
    for neighbors in ppi_adj.values():
        ppi_neighbors.update(g.upper() for g in neighbors)
    # Each drug gene's contribution is degree-damped: a HUB target that "interacts"
    # with a disease gene is far less specific than a low-degree target doing so.
    total = 0.0
    for g in dg:
        d = gene_damp.get(g, 1.0) if gene_damp else 1.0
        if g in ppi_nodes:
            total += d * 1.0
        elif g in ppi_neighbors:
            total += d * 0.5
    return min(1.0, total / len(dg))


def _to_phase(max_phase) -> int:
    try:
        return int(float(max_phase or 0))
    except (TypeError, ValueError):
        return 0


def _indication_matches_disease(existing_ind: str, disease_name: str) -> bool:
    """True when the drug's existing indication text specifically covers THIS
    disease. Requires ALL of the disease's SPECIFIC tokens (generic words like
    'syndrome', 'disease', 'carcinoma' stripped) to appear — so 'Cowden syndrome'
    is NOT matched to a drug approved for 'Myelodysplastic syndromes' merely
    because both contain 'syndrome' (which was inflating the clinical dimension)."""
    ind = (existing_ind or "").lower()
    if not ind:
        return False
    try:
        from services.disease_id import canonical_key
        toks = canonical_key(disease_name)
    except Exception:
        toks = frozenset(w for w in disease_name.lower().split() if len(w) > 4)
    # Drop oncology/generic descriptors that are not disease-specific.
    _generic = {"carcinoma", "neoplasm", "neoplasms", "cancer", "tumor", "tumour",
                "deficiency", "familial", "congenital", "hereditary", "type"}
    specific = [t for t in toks if t not in _generic and len(t) > 3]
    if not specific:
        specific = list(toks)
    if not specific:
        return False
    return all(t in ind for t in specific)


def _score_clinical(max_phase, existing_ind: str, disease_name: str,
                    trial_count: int = 0, max_trial_phase: int = 0,
                    trial_outcome: float = 0.0, indication_phase=None) -> float:
    """Clinical evidence SPECIFIC to this drug-disease pair.

    A drug's global development phase is clinical evidence only for the indications
    it was developed FOR (an approved oncology drug has NO clinical evidence for
    diabetic nephropathy — crediting its global Phase 4 there is reverse-causation).
    So global-phase credit is gated on the disease matching an existing indication.

    CONFOUNDING GUARD (audited 2026-07): when the disease matches an existing
    indication, credit the PER-INDICATION development phase (`indication_phase`), NOT
    the drug's global max phase. A globally-approved drug that only reached Phase 2 in
    THIS disease (imatinib→pulmonary hypertension, mepolizumab→eosinophilic
    esophagitis — studied, not approved) must NOT be credited a Phase-4 "established
    therapy" clinical score. When per-indication phase is unknown we cap conservatively
    at Phase 2 rather than assuming global approval applies here.

    Real trials of THIS drug in THIS indication ARE disease-specific evidence — a human
    decided it was worth testing — so they earn direct (phase-scaled) credit."""
    if _indication_matches_disease(existing_ind, disease_name):
        eff_phase = (int(indication_phase) if indication_phase is not None
                     else min(2, _to_phase(max_phase)))   # unknown → conservative Phase-2 cap
        phase_map = {4: 1.0, 3: 0.75, 2: 0.50, 1: 0.25}
        return phase_map.get(eff_phase, 0.05)
    if trial_count > 0 or max_trial_phase > 0:
        base = min(0.7, 0.25 + 0.12 * max_trial_phase + 0.05 * min(trial_count, 5))
        # Outcome direction (P3): a trial that FAILED for efficacy/safety here is
        # negative evidence — the drug was tested and did not work — so the credit
        # shrinks toward 0 (outcome ∈ [-1,0]); a completed / with-results trial
        # earns a small boost (outcome > 0).
        if trial_outcome < 0:
            return round(max(0.0, base * (1.0 + trial_outcome)), 4)
        return round(min(0.7, base * (1.0 + 0.3 * trial_outcome)), 4)
    return 0.05                            # untested repurposing hypothesis


def _score_indication(existing_ind: str, disease_name: str, indication_phase=None) -> float:
    """Credit for the disease already being an indication of the drug. Phase-scaled
    (audited 2026-07): full 'this IS an established use' credit only for an APPROVED
    indication (Phase 4). A merely-STUDIED indication is prior clinical art, not an
    established use, so its credit is scaled down by development phase — otherwise a
    Phase-2 studied pair reads as an established therapy (the confounding that ranked
    tried/failed indications at the top)."""
    stop = {"disease", "disorder", "syndrome", "with", "associated"}
    dw = {w.lower() for w in re.split(r'\W+', disease_name) if len(w) > 3 and w.lower() not in stop}
    iw = {w.lower() for w in re.split(r'\W+', existing_ind or "") if len(w) > 3 and w.lower() not in stop}
    if not dw:
        return 0.0
    frac = min(1.0, len(dw & iw) / len(dw))
    if frac <= 0:
        return 0.0
    ph = None if indication_phase is None else int(indication_phase)
    # 4=approved→full; 3/2/1=studied→partial; unknown→0.5 (can't confirm approval here)
    phase_factor = {4: 1.0, 3: 0.6, 2: 0.4, 1: 0.25}.get(ph, 0.5)
    return round(frac * phase_factor, 4)


def _score_regulatory(max_phase, disease_name: str, existing_ind: str) -> float:
    score = 0.0
    phase = _to_phase(max_phase)
    if phase >= 4:   score += 0.50
    elif phase == 3: score += 0.35
    elif phase == 2: score += 0.20
    elif phase == 1: score += 0.10
    try:
        from services.disease_ontology import is_orphan_candidate
        if is_orphan_candidate(disease_name):
            score += 0.30
    except Exception:
        pass
    disease_kws = [w for w in disease_name.lower().split() if len(w) > 4]
    if any(kw in (existing_ind or "").lower() for kw in disease_kws):
        score += 0.20
    return min(1.0, score)


# ── Compound scorer ───────────────────────────────────────────────────────────

def score_compound_for_disease(
    compound: Dict,
    disease_name: str,
    disease_genes: List[str],
    disease_pathways: List[Dict],
    ppi_adj: Dict[str, List[str]],
    drug_genes: Optional[List[str]] = None,
    drug_actions: Optional[Dict[str, str]] = None,
    disease_gene_weights: Optional[Dict[str, float]] = None,
    trial_count: int = 0,
    max_trial_phase: int = 0,
    trial_outcome: float = 0.0,
    mechanistic_prior: Optional[float] = None,
    studied_for_disease: bool = False,
    indication_phase: Optional[int] = None,
) -> Dict:
    existing_ind = compound.get("indications", "") or ""
    raw_phase    = compound.get("max_phase") or compound.get("max_phase_for_ind") or 0
    try:
        max_phase = int(float(raw_phase))
    except (TypeError, ValueError):
        max_phase = 0
    if drug_genes is None:
        drug_genes = []

    # Build a signed drug signature for the direction-aware pathway score when we
    # have the drug's mechanism action types (gene→INHIBITOR/AGONIST/…). Optional:
    # absent → pathway score stays direction-blind, exactly as before.
    drug_sig = None
    if drug_actions:
        try:
            from services.signature_engine import drug_signature
            drug_sig = drug_signature(
                [{"gene": g, "action": a} for g, a in drug_actions.items()])
        except Exception as e:
            logger.debug(f"drug signature build skipped: {e}")

    # Genetics-weighted target overlap when weights are available (validated to
    # improve recovery + reduce leakage); flat overlap otherwise.
    target_score = (_score_target_overlap_weighted(drug_genes, disease_gene_weights)
                    if disease_gene_weights
                    else _score_target_overlap(drug_genes, disease_genes))
    # Hub-degree damping (DWPC): weight each drug target's pathway/PPI contribution
    # by its global connectivity, so a network hub (EGFR, TP53…) cannot trivially
    # max out cohesion for every disease it touches. Fail-soft → no damping.
    gene_damp = {}
    try:
        from services.hub_degree import damping_map
        gene_damp = damping_map(drug_genes)
    except Exception as e:
        logger.debug(f"hub-degree damping skipped: {e}")

    scores = {
        "target":     target_score,
        "pathway":    _score_pathway_overlap(drug_genes, disease_pathways,
                                             drug_signature=drug_sig, gene_damp=gene_damp),
        "ppi":        _score_ppi_network(drug_genes, ppi_adj, gene_damp=gene_damp),
        "clinical":   _score_clinical(max_phase, existing_ind, disease_name,
                                       trial_count, max_trial_phase, trial_outcome,
                                       indication_phase=indication_phase),
        "indication": _score_indication(existing_ind, disease_name,
                                        indication_phase=indication_phase),
        "regulatory": _score_regulatory(max_phase, disease_name, existing_ind),
    }

    # ── Mechanism direction (fix the direction-blindness bug) ────────────────────
    # Target overlap alone can't tell an INHIBITOR of a disease-driving gene (treats)
    # from an ACTIVATOR of it (exacerbates). Fold the drug's action × the disease's
    # up/down-regulated genes (HetioNet) into the target dimension, and raise a
    # contraindication flag when the net effect would push the disease the wrong way.
    direction = {"factor": 1.0, "net": "neutral", "covered": False, "flag": ""}
    try:
        from services.mechanism_direction import mechanism_direction
        direction = mechanism_direction(drug_genes, drug_actions, disease_name)
        scores["target"] = round(scores["target"] * direction["factor"], 4)
    except Exception as e:
        logger.debug(f"mechanism direction skipped: {e}")

    _overlap_n = len({g.upper() for g in drug_genes} & {g.upper() for g in disease_genes})

    # ── Target-driven pathogenic proliferation ───────────────────────────────────
    # A mechanism the genetic-overlap sub-scores are STRUCTURALLY blind to: the drug
    # inhibits a positive proliferation driver (CDK4/6, cyclins…) and the disease is
    # built on pathogenic hyperproliferation (RA synovial hyperplasia, fibrosis,
    # restenosis…). Data-driven + direction-aware (services/proliferation.py). When it
    # fires it becomes the mechanistic prior (un-phantoms the pair + credits it) and
    # adds a bounded bonus — so palbociclib→RA is no longer buried at ~0.01.
    prolif = {"match": False, "score": 0.0}
    try:
        from services.proliferation import proliferation_match
        _pname = compound.get("name") or compound.get("pref_name") or ""
        prolif = proliferation_match(drug_genes, drug_actions, disease_name)
        if prolif.get("match"):
            mechanistic_prior = max(float(mechanistic_prior or 0.0), float(prolif["score"]))
    except Exception as e:
        logger.debug(f"proliferation match skipped: {e}")

    # ── Signature reversal (CMap/LINCS, gap #3) ──────────────────────────────────
    # Real transcriptional connectivity: does the drug REVERSE the disease's expression
    # signature? Orthogonal to target/KG/genetics — catches mechanisms none of them see.
    # A reversing signal feeds the mechanistic prior (un-phantoms) + a bonus; a MIMICKING
    # signal (drug amplifies the disease signature) is a contraindication-like penalty.
    reversal = {"covered": False, "score": 0.0, "direction": "none"}
    try:
        from services.signature_reversal import reversal_score
        _cname = compound.get("name") or compound.get("pref_name") or ""
        if _cname:
            reversal = reversal_score(_cname, disease_name)
            if reversal.get("direction") == "reversing":
                mechanistic_prior = max(float(mechanistic_prior or 0.0),
                                        0.5 * float(reversal.get("score", 0.0)))
    except Exception as e:
        logger.debug(f"signature reversal skipped: {e}")

    # ── Externally-established mechanistic linkage (surfacing evidence) ───────────
    # When a candidate was surfaced by an engine that has ALREADY established a
    # mechanistic drug→disease link — the pathway screen (drug direction-matches a
    # disease-driving pathway) or novel-target discovery (drug hits an inferred
    # disease target) — that linkage is real evidence, yet it is invisible to the
    # OT-gene-overlap sub-scores (the surfacing genes/pathways are by construction
    # NOT in the disease's top-40 OT set). Discarding it is why strong pathway
    # modulators read as "~5". We fold the prior into the PATHWAY dimension (bounded,
    # so it credits a real-but-unproven hypothesis without flattening every
    # direction-matched drug to the top). This also lets the pair clear the phantom-
    # cohesion gate. Disease-specific evidence (real target overlap, trials, an
    # existing indication) still differentiates candidates above this floor.
    if mechanistic_prior is not None:
        try:
            scores["pathway"] = max(scores["pathway"],
                                    round(0.5 * max(0.0, min(1.0, float(mechanistic_prior))), 4))
        except (TypeError, ValueError):
            pass

    # ── Applicable-weight renormalization ────────────────────────────────────────
    # For a GENUINE repurposing pair (drug not yet used for this disease), two of the
    # six dimensions are structurally near-zero by construction: `indication` is 0
    # (if it matched, it wouldn't be repurposing) and `clinical` sits at its 0.05
    # "untested-hypothesis" floor unless a real trial exists. Summing them at 0 into a
    # weighted mean silently caps the achievable score at ~0.3 no matter how strong
    # the mechanism is — the exact reason good leads read as "~10". So we renormalize:
    # the disease-specific mechanistic mass (target/pathway/ppi/clinical/indication,
    # 0.90 of the weight) is spread over only the dimensions that actually carry
    # disease-specific signal, with a denominator floor so a single weak signal is not
    # over-amplified. `regulatory` is a GLOBAL de-risking prior (drug's overall phase),
    # not disease-specific evidence, so it is kept as a small bounded term and is NOT
    # eligible for amplification — otherwise an approved drug with zero mechanism for
    # the disease (reverse causation) would inflate. A pair with NO disease-specific
    # evidence therefore still scores low: this separates leads, it does not inflate.
    # Audited 2026-07: the previous version renormalized clinical+indication (prior
    # clinical art) INTO the "mechanistic mass", so a ZERO-mechanism pair that had merely
    # been trialed/approved here renormalized up to "Strong" (imatinib→pulmonary
    # hypertension 0.66 with target=pathway=ppi=0). Prior clinical art is real evidence
    # but it is NOT mechanism. So we (a) renormalize ONLY the mechanism dims
    # (target/pathway/ppi) — a strong mechanism across few active dims still reads strong —
    # and (b) add clinical+indication as a bounded ADDITIVE term at nominal weight, never
    # amplified. _score_clinical / _score_indication already phase-scale that art, so an
    # APPROVED own-indication keeps full credit while a merely-studied pair gets only
    # partial. Net: mechanism separates leads; prior art can lift a score but can never
    # manufacture a "Strong" out of no mechanism.
    _MECH_DIMS = ("target", "pathway", "ppi")
    _MECH_MASS = sum(WEIGHTS[k] for k in _MECH_DIMS)              # 0.65
    _DENOM_FLOOR = 0.45          # never divide the mechanistic mass by less than this
    _active = [k for k in _MECH_DIMS if scores[k] > 1e-9]
    if _active:
        _w_active = sum(WEIGHTS[k] for k in _active)
        _mech_raw = sum(scores[k] * WEIGHTS[k] for k in _active)
        mech_component = (_mech_raw / max(_DENOM_FLOOR, _w_active)) * _MECH_MASS
    else:
        mech_component = 0.0
    prior_art = (scores["clinical"] * WEIGHTS["clinical"]
                 + scores["indication"] * WEIGHTS["indication"])
    composite = min(1.0, mech_component + prior_art
                    + scores["regulatory"] * WEIGHTS["regulatory"])

    # Real HetioNet network evidence (Compound→Gene←Disease path). Independent of
    # ChEMBL, so it can legitimately surface a candidate as a lead even when its
    # targets don't overlap the Open Targets disease-gene set. Folded in as a
    # bounded bonus and exposed in score_breakdown for the lead-mechanism filter
    # and the dossier. Fail-soft — never block scoring if the graph is unavailable.
    net = {"score": 0.0, "genes": [], "basis": ""}
    try:
        from services.repurposing_scorer import _network_evidence
        cname = compound.get("name") or compound.get("pref_name") or ""
        if cname:
            net = _network_evidence(cname, [disease_name])
    except Exception as e:
        logger.debug(f"network evidence skipped: {e}")
    net_score = float(net.get("score") or 0.0)
    # Specificity damping (audited 2026-07): the HetioNet Compound→Gene←Disease signal
    # is promiscuous — a single shared gene, or a "generic/weak" basis, lights up at
    # 0.3–0.85 for pairs as unrelated as loratadine→Parkinson, so it cannot separate a
    # plausible-wrong pair from a random one. A DIRECT "treats this disease" edge is
    # kept at full credit; an indirect path is only credited in full when it rests on
    # ≥2 shared disease genes, otherwise it is damped to a hint. Keeps the signal
    # (it never breaks true-positive ranking) without letting one promiscuous gene inflate.
    _net_basis = net.get("basis", "") or ""
    _net_genes = net.get("genes", []) or []
    if _net_basis.startswith("direct"):
        net_eff = net_score
    elif len(_net_genes) >= 2:
        net_eff = net_score
    else:
        net_eff = net_score * 0.35          # single-gene / generic / weak: promiscuous → hint only
    composite = min(1.0, composite + 0.10 * net_eff)

    # Target-driven proliferation bonus (bounded) — the mechanism is credited on top of
    # the renormalised pathway prior it already set, so a genuine proliferation-arrest
    # lead (palbociclib→RA) reaches a real Moderate score instead of ~0.01.
    if prolif.get("match"):
        composite = min(1.0, composite + 0.18 * float(prolif["score"]))

    # Signature-reversal bonus / mimic penalty (bounded, orthogonal).
    if reversal.get("direction") == "reversing":
        composite = min(1.0, composite + 0.15 * float(reversal.get("score", 0.0)))
    elif reversal.get("direction") == "mimicking":
        composite *= 0.85          # the drug amplifies the disease signature (wrong way)

    # Directional literature evidence (P3): a typed drug→gene→disease path or a
    # direct drug-treats-disease edge from DRKG (GNBR/DGIDB/DrugBank triples) —
    # relational support, not co-occurrence. Bounded bonus; fail-soft (0) off-coverage.
    directional = {"covered": False, "signal": 0.0, "direct_treat": False,
                   "n_paths": 0, "path_genes": [], "note": ""}
    try:
        from services.directional_evidence import directional_evidence
        dname = compound.get("name") or compound.get("pref_name") or ""
        if dname and disease_name:
            directional = directional_evidence(dname, disease_name)
            composite = min(1.0, composite + 0.12 * float(directional.get("signal", 0.0)))
    except Exception as e:
        logger.debug(f"directional evidence skipped: {e}")

    # ── Learned composite (gap #4) ───────────────────────────────────────────────
    # If a supervised model has been fitted AND validated to BEAT the hand weights
    # (services/train_composite.py only saves it in that case), use its evidence score
    # in place of the hand-weighted sum + bonuses. Fail-soft: no model → keep the hand
    # composite, which is itself validated (AUC 0.97 CtD-vs-random). The clinical
    # PENALTIES below still apply on top either way. DWPC plausibility stays its own
    # axis (excluded from the features to avoid circularity).
    try:
        from services.composite_model import learned_composite
        _lc = learned_composite([
            scores["target"], scores["pathway"], scores["ppi"], scores["clinical"],
            scores["indication"], scores["regulatory"], net_score,
            float(directional.get("signal", 0.0) or 0.0),
            float(prolif.get("score", 0.0) or 0.0),
            float(reversal.get("connectivity", 0.0) or 0.0),
            float(direction.get("factor", 1.0) or 1.0),
        ])
        if _lc is not None:
            composite = float(_lc)
    except Exception as e:
        logger.debug(f"learned composite skipped: {e}")

    # ── Track 1 snapshot: MECHANISTIC PLAUSIBILITY ────────────────────────────────
    # The composite at THIS point is pure "does the biology connect?" — mechanism +
    # KG/network + directional + learned, before ANY clinical/feasibility penalty.
    # Captured separately so a strong hypothesis is never dressed up as a ready lead
    # (principle: separate mechanistic plausibility from clinically actionable repurposing).
    mech_plausibility = round(min(1.0, composite), 4)

    # Negative-safety cross-filter (generalized): down-rank a match when the drug
    # carries a SERIOUS toxicity signal for the organ system the disease lives in
    # (e.g. a pulmonary-embolism-risk drug for a pulmonary indication). Applied as
    # a bounded multiplier on the composite. Fail-soft — multiplier 1.0 when no
    # adverse-event data is available, so we never invent a penalty.
    # ── Confounding-by-indication guard (shared by the penalty blocks) ────────────
    # A drug that is an ESTABLISHED/DEVELOPED therapy for THIS disease (approved, in
    # trials, or with a matching indication record) has a real-world benefit/risk that
    # is already settled by its development — so the class-heuristic penalties below do
    # not apply to it. Worse, its adverse-event / organ-toxicity signals are CONFOUNDED
    # by the indication itself (imatinib's haematologic AEs arise BECAUSE it treats a
    # blood cancer; sildenafil is not "high-tox for a low-severity indication" — it is
    # THE therapy). Audited 2026-07: without this guard the penalty stack halved the
    # scores of drugs' own approved indications (imatinib→CML 0.85→0.42, sildenafil→ED
    # 0.64→0.31). For a genuine NOVEL pair scores["indication"]≈0, so the guard is off
    # and the constraints still apply — this fixes false penalties, it does not inflate.
    # Confounding-by-indication guard — suppress the safety/clinical/coverage penalties
    # ONLY for a genuinely APPROVED own-indication (Phase 4 in THIS disease). Audited
    # 2026-07: the old trigger fired on ANY trialed pair (trial_count>0 / phase>=2 /
    # name-match / indication>=0.20), so a merely-STUDIED or FAILED pair (imatinib→
    # pulmonary hypertension, mepolizumab→eosinophilic esophagitis) had its guardrails
    # switched off as if it were established therapy — inflating it. Now suppression
    # requires real approval here; studied and novel pairs keep every gate.
    _own_therapy = bool(
        (indication_phase is not None and int(indication_phase) >= 4)
        or studied_for_disease                          # forward screen asserts drug is used for disease
        or (indication_phase is None
            and _indication_matches_disease(existing_ind, disease_name)
            and _to_phase(max_phase) >= 4)              # legacy fallback only when per-indication phase unknown
    )

    safety = {"multiplier": 1.0, "penalized": False, "flags": []}
    try:
        from services.safety_filter import assess as _assess_safety, therapeutic_organs as _ther_organs
        cname = compound.get("name") or compound.get("pref_name") or ""
        if cname and disease_name and not _own_therapy:
            # Layer 3 — pass the organs the drug is DEVELOPED to act on, so a toxicity
            # signal in one of them (a lung drug's respiratory FAERS) is treated as
            # therapeutic-organ overlap, not a false contraindication (mepolizumab→ABPA).
            _tho = None
            try:
                _tho = _ther_organs(existing_ind)
            except Exception:
                _tho = None
            safety = _assess_safety(cname, disease_name, therapeutic_organs=_tho)
            # Recorded, not applied inline — soft penalties are consolidated below so
            # two moderate ones don't compound multiplicatively into annihilation.
        elif _own_therapy:
            safety = {"multiplier": 1.0, "penalized": False, "flags": [],
                      "suppressed": "established therapy for this disease — organ-toxicity "
                                    "overlap is confounded by the indication, not a contraindication"}
    except Exception as e:
        logger.debug(f"safety filter skipped: {e}")

    # Clinical Constraint Harmonizer (CCH): multiply by clinical-reality fit —
    # severity/therapeutic-window and chronicity/tolerability match. Distinct from
    # the safety filter above (organ-toxicity↔disease-organ); this is a
    # severity/duration MATCH check. Fail-soft (1.0 without data).
    clinical = {"multiplier": 1.0, "penalized": False, "flags": [], "factors": {}}
    try:
        from services.clinical_constraints import harmonize as _harmonize
        cname = compound.get("name") or compound.get("pref_name") or ""
        if cname and disease_name and not _own_therapy:
            clinical = _harmonize(cname, disease_name, smiles=compound.get("smiles", ""),
                                  has_trials=(trial_count > 0 or max_trial_phase > 0),
                                  indications=existing_ind)
            # Recorded, not applied inline (consolidated with the other soft penalties).
        elif _own_therapy:
            clinical = {"multiplier": 1.0, "penalized": False, "flags": [], "factors": {},
                        "suppressed": "established therapy for this disease — the clinical "
                                      "therapeutic window is set by its approval/development"}
    except Exception as e:
        logger.debug(f"clinical constraints skipped: {e}")

    # Demonstrated-failure penalty (P3 trial outcomes): if this drug was actually
    # TRIED and FAILED for efficacy/safety in this indication (trial_outcome < 0),
    # down-rank it — a real clinical failure is stronger evidence than any
    # inferred association. Bounded; positive/absent outcomes do nothing here.
    trial_failure = {"multiplier": 1.0, "penalized": False, "outcome": round(trial_outcome, 3)}
    if trial_outcome < -0.3:
        tf_mult = max(0.35, 1.0 + 0.55 * trial_outcome)   # outcome -1 → 0.45, -0.6 → 0.67
        # Hard, evidence-based penalty (a trial actually failed here) — applied below.
        trial_failure = {"multiplier": round(tf_mult, 3), "penalized": True,
                         "outcome": round(trial_outcome, 3),
                         "flag": "Down-ranked: a trial of this drug in this indication "
                                 "was terminated for efficacy/safety."}

    # Target-coverage / mandatory-intersection gate (generalized): for a polygenic
    # disease, penalize a drug that covers only part of the driver-target set. Uses
    # the genetic-association weights; fail-soft (multiplier 1.0) without them.
    coverage = {"multiplier": 1.0, "penalized": False, "coverage": None, "flags": []}
    try:
        if disease_gene_weights and not _own_therapy:
            from services.target_coverage import assess_coverage
            coverage = assess_coverage(drug_genes, disease_gene_weights)
            # Recorded, not applied inline (consolidated with the other soft penalties).
        elif _own_therapy and disease_gene_weights:
            # An established therapy for this disease has DEMONSTRATED sufficiency —
            # "incomplete polygenic-driver coverage" is moot when the drug already works
            # (sildenafil→ED via PDE5 alone, imatinib→CML via BCR-ABL alone). Guard off
            # for genuine novel pairs, where partial coverage is a real concern.
            coverage = {"multiplier": 1.0, "penalized": False, "coverage": None, "flags": [],
                        "suppressed": "established therapy — demonstrated sufficiency overrides "
                                      "the polygenic-coverage heuristic"}
    except Exception as e:
        logger.debug(f"target-coverage gate skipped: {e}")

    # ── Appropriateness gate (adverse-event / phenotype-causation) — PLATFORM-WIDE ──
    # A disease that is itself a SERIOUS FAERS reaction of the drug is a toxicity the
    # drug CAUSES, not a target it can treat — adverse-event overlap is NEGATIVE
    # evidence, never target engagement. Also catches inhibitor↔loss-of-function
    # mismatch. Bounded factor (0.3–1.0), demote-not-drop. Applied in the reverse
    # screen already; wired here so forward discovery inherits the same gate.
    # CONFOUNDING-BY-INDICATION GUARD: a drug DEVELOPED for this disease (approved,
    # in trials, or with a matching indication record) will list the disease among its
    # FAERS reports — that is the indication co-reported, NOT a toxicity the drug causes.
    # Feeding FAERS to the AE check for such a drug would flag its OWN therapy as harmful
    # (e.g. methotrexate → RA). So we withhold FAERS from the gate when the drug is an
    # established/developed therapy for the disease; the LoF/direction checks still run.
    appropriate = {"factor": 1.0, "appropriate": True, "flags": [], "reasons": []}
    try:
        from services.disease_appropriateness import appropriateness as _appr, infer_drug_action
        from services.safety_filter import _faers_serious_reactions
        _cn = compound.get("name") or compound.get("pref_name") or ""
        if _cn and disease_name:
            # _own_therapy computed once above (confounding-by-indication guard)
            _rx = [] if _own_therapy else (_faers_serious_reactions(_cn) or [])
            _rx_total = sum(c for _, c in _rx)
            _act = infer_drug_action(_cn, drug_genes)
            appropriate = _appr(_cn, disease_name, existing_ind or [], _act,
                                faers_reactions=_rx, faers_total=_rx_total)
    except Exception as e:
        logger.debug(f"appropriateness gate skipped: {e}")

    # ── Delivery feasibility for LOCALIZED indications — PLATFORM-WIDE ─────────────
    # For a compartmentalized disease (CNS / eye / skin / joint / local-airway) a
    # mechanistically plausible pair is only ACTIONABLE if the drug can physically
    # reach the tissue. Bounded multiplier, demote-not-drop; fail-soft (1.0) when the
    # indication is systemic or delivery is unassessable.
    delivery = {"localized": False, "compartment": "systemic", "deliverable": True,
                "multiplier": 1.0, "flag": "", "assessed": False}
    try:
        from services.delivery_feasibility import assess_delivery
        _cn = compound.get("name") or compound.get("pref_name") or ""
        _ind_blob = (" ".join(existing_ind or [])).lower()
        _cns_ind = any(k in _ind_blob for k in
                       ("cns", "alzheimer", "parkinson", "epilep", "depress", "schizophren",
                        "migraine", "neuropath", "dementia", "seizure", "psychos"))
        delivery = assess_delivery(disease_name, smiles=compound.get("smiles", ""),
                                   drug_name=_cn, cns_indicated=_cns_ind)
    except Exception as e:
        logger.debug(f"delivery feasibility skipped: {e}")

    # Mechanism scope (generalized): the composite is target/genetics-centric, so it
    # is only meaningful for target-mediated diseases. A drug with no resolvable
    # HUMAN protein targets — an anti-infective (bacterial/viral targets), a
    # symptomatic/nutritional agent, or simply an unmapped molecule — will score
    # near-zero on the mechanistic dimensions for reasons of SCOPE, not weak
    # evidence. Flag it so the UI can say so instead of showing a silent ~0.
    mechanism_applicable = bool(drug_genes)
    mechanism_scope = {
        "target_mediated": mechanism_applicable,
        "note": ("" if mechanism_applicable else
                 "No human protein target resolved — the mechanistic score does "
                 "not apply (e.g. an anti-infective or symptomatic drug). Ranking "
                 "reflects clinical / literature evidence only."),
    }

    # ── CTPA Rule 1 — harmonic affinity / cohesion gate (PLATFORM-WIDE) ──────────
    # In the canonical scorer so EVERY surface (forward discovery, pathway, novel
    # targets, reverse) inherits it. A pair resting on a single target/text match
    # with NO functional cohesion (≤1 shared gene AND no pathway/PPI/directional/
    # trial support) is a phantom — a low-affinity off-target string match. Crush it.
    #
    # A single driver hit is NOT a phantom when that driver is a STRONG genetic
    # association: the weighted target score already gates on association strength, so
    # target ≥ 0.40 means at least one high-confidence disease-gene hit — real (if
    # narrow) cohesion, not a coincidental string match. Without this, a selective
    # drug hitting one strongly-associated driver (e.g. CDK4 for Alzheimer, assoc
    # 0.54 → target 0.52) got both a high target sub-score AND the phantom penalty.
    _cname = compound.get("name") or compound.get("pref_name") or ""
    _overlap_n = len({g.upper() for g in drug_genes} & {g.upper() for g in disease_genes})
    _cohesion = (_overlap_n >= 2 or scores["pathway"] > 0.05 or scores["ppi"] > 0.05
                 or scores["target"] >= 0.40 or float(net_score) >= 0.35
                 or float(directional.get("signal", 0.0)) > 0 or trial_count > 0)
    ctpa = {"cohesion": bool(_cohesion), "phantom": False, "multiplier": 1.0}
    if drug_genes and not _cohesion:
        # Phantom gate: a pair resting on a single target/text match with no pathway/
        # PPI/directional/trial cohesion is likely a low-affinity off-target string
        # match. Still a strong down-rank (0.35), but no longer an annihilation to
        # ~0.15 that a genuinely selective single-target lead couldn't recover from.
        # Applied as a hard gate in the consolidation block below.
        ctpa = {"cohesion": False, "phantom": True, "multiplier": 0.35,
                "flag": ("No functional cohesion beyond a single target/text "
                         "association - likely a low-affinity off-target match.")}

    # ── CTPA Rule 2 — registry ghost audit (PLATFORM-WIDE) ───────────────────────
    # Zero a drug whose late-phase program in this indication failed / was abandoned
    # (recommending it means repeating a dead trial). Only queried when there is a
    # late-phase clinical footprint (bounded, cached); fail-soft.
    registry = {"ghost": False, "multiplier": 1.0}
    if _cname and disease_name and (trial_count > 0 or max_trial_phase >= 2):
        try:
            from services.registry_audit import audit as _reg_audit
            registry = _reg_audit(_cname, disease_name)
            # Hard gate (a late-phase program here failed/was abandoned) — applied below.
        except Exception as e:
            logger.debug(f"registry audit skipped: {e}")

    # ── Penalty consolidation ────────────────────────────────────────────────────
    # SOFT penalties (organ-toxicity fit, clinical-reality fit, incomplete polygenic
    # coverage) are environmental-fit signals of similar kind. Compounding them
    # multiplicatively double-counts the same "imperfect fit" and was a major driver
    # of single-digit scores, so we take the SINGLE worst soft penalty instead of the
    # product. HARD, evidence-based gates (phantom off-target, a trial that actually
    # failed here, a dead registry program) remain multiplicative — each is a distinct
    # real-world kill signal that should stack.
    # Delivery feasibility joins the SOFT group (environmental-fit, single worst wins) —
    # it should not compound multiplicatively with organ-toxicity/coverage into annihilation.
    _soft_mult = min(float(safety.get("multiplier", 1.0)),
                     float(clinical.get("multiplier", 1.0)),
                     float(coverage.get("multiplier", 1.0)),
                     float(delivery.get("multiplier", 1.0)))
    composite *= _soft_mult
    # A net-HARMFUL mechanism direction is a contraindication signal — down-rank the
    # whole pair (the drug would push the disease the wrong way), not just its target dim.
    if direction.get("net") == "harmful":
        composite *= float(direction.get("factor", 1.0))
    composite *= float(ctpa.get("multiplier", 1.0))
    composite *= float(trial_failure.get("multiplier", 1.0))
    # Appropriateness is a HARD, evidence-based gate (the disease is a toxicity the drug
    # CAUSES, or an inhibitor↔loss-of-function mismatch) — a distinct real-world kill
    # signal, so it stacks with the other hard gates rather than joining the soft min.
    composite *= float(appropriate.get("factor", 1.0))
    if registry.get("ghost"):
        composite *= float(registry.get("multiplier", 1.0))

    final_score = round(min(1.0, composite), 4)

    # Calibration (P1): make the raw score interpretable by ranking it against a
    # null distribution of random drug-disease pairs → percentile + tier, so a
    # strong lead reads as "top 3%" rather than a weak-looking absolute number.
    calibration = {"percentile": None, "enrichment": None, "tier": None, "basis": "none"}
    try:
        from services.score_calibration import calibrate as _calibrate
        calibration = _calibrate(final_score)
    except Exception as e:
        logger.debug(f"score calibration skipped: {e}")

    # ── Evidence-quality gate (anti-overselling + hub de-biasing) ────────────────
    # Calibration ranks a score against RANDOM pairs, so any mechanistically-linked
    # candidate looks top-percentile — that is not the same as a strong LEAD. A
    # single-gene genetic association with no corroboration cannot be "Strong",
    # however high its percentile. CORROBORATION deliberately EXCLUDES pathway/PPI,
    # because a network HUB (EGFR, TP53…) trivially maxes those for every disease it
    # touches — the very inflation that made a weak Cowden/Erlotinib hit read #1.
    # Real corroboration = breadth of target overlap, real clinical/trial evidence,
    # directional literature, or genuine driver coverage.
    # Audited 2026-07: corroboration must be MECHANISTIC. The old set counted
    # `clinical > 0.10` and `trial_count > 0`, so PRIOR CLINICAL ART alone (a drug merely
    # trialed here) satisfied "corroboration" and could keep a Strong tier — the same
    # confounding fixed in the composite. A "Strong" claim now requires real mechanism:
    # breadth of target overlap, directional KG evidence, genuine driver coverage, a real
    # network path, or pathway mechanism. Trial existence is prior art, not mechanism, so
    # it no longer buys a tier.
    corroboration = sum([
        _overlap_n >= 2,                                    # breadth of target overlap
        scores["target"] >= 0.40,                          # a genuine on-target hit (single strong target counts)
        float(directional.get("signal", 0.0)) > 0,         # directional KG evidence
        float(coverage.get("coverage") or 0.0) >= 0.5,     # genuine driver coverage
        net_eff >= 0.30,                                    # real (non-promiscuous) network path
        scores["pathway"] >= 0.30,                          # pathway mechanism
    ])
    if corroboration < 1 and calibration.get("tier") in ("Strong", "Promising"):
        calibration = dict(calibration, tier="Moderate", evidence_capped=True,
                           evidence_note=("No mechanistic corroboration (no target breadth, "
                                          "directional, coverage, network-path or pathway "
                                          "support) - capped below Strong. Prior clinical "
                                          "art and hub-target PPI signals are excluded as "
                                          "non-discriminating."))

    # ── Two-track separation: mechanistic PLAUSIBILITY vs clinical ACTIONABILITY ──
    # Rank a repurposing pair as ACTIONABLE only when mechanistic, clinical, and
    # feasibility evidence ALIGN. Otherwise it stays a labelled hypothesis — demoted,
    # never dropped — so a strong mechanism is never mistaken for a ready-to-trial lead.
    mech_present = mech_plausibility >= 0.20 or corroboration >= 1
    clinical_present = (scores["clinical"] > 0.10 or trial_count > 0
                        or max_trial_phase >= 2 or bool(existing_ind))
    feasible = bool(delivery.get("deliverable", True))
    ae_caused = "reported adverse event" in appropriate.get("flags", [])
    harmful = direction.get("net") == "harmful"

    reasons: List[str] = []
    if not mech_present:     reasons.append("no mechanistic linkage to the disease")
    if not clinical_present: reasons.append("no clinical or trial support")
    if not feasible:         reasons.append(delivery.get("flag") or "delivery-limited to the target tissue")
    if ae_caused:            reasons.append("disease is a serious adverse event the drug causes")
    if harmful:              reasons.append("mechanism direction would worsen the disease")

    aligned = mech_present and clinical_present and feasible and not ae_caused and not harmful
    if aligned:
        act_tier = "Actionable"
    elif ae_caused or harmful:
        act_tier = "Contraindicated by evidence"
    elif not feasible:
        act_tier = "Delivery-limited"
    elif mech_present and not clinical_present:
        act_tier = "Mechanistic hypothesis"
    else:
        act_tier = "Insufficiently supported"

    actionability = {
        "aligned": bool(aligned), "tier": act_tier, "reasons": reasons,
        "mechanistic_plausibility": mech_plausibility,
        "clinical_present": bool(clinical_present), "feasible": feasible,
        # ranking key — aligned leads sort above hypotheses of equal raw score.
        "rank_score": round(final_score * (1.0 if aligned else 0.6), 4),
    }

    return {
        "composite_score": final_score,
        "mechanistic_plausibility": mech_plausibility,
        "actionability":   actionability,
        "appropriateness": appropriate,
        "delivery":        delivery,
        "scores":          {k: round(v, 4) for k, v in scores.items()},
        "score_breakdown": {
            "network_score": round(net_score, 4),
            "network_genes": net.get("genes", []),
            "network_basis": net.get("basis", ""),
        },
        "proliferation":   prolif,     # target-driven proliferation match (RA-type mechanism)
        "signature_reversal": reversal, # CMap/LINCS transcriptional connectivity (gap #3)
        "direction":       direction,  # mechanism direction (corrective/harmful) + contraindication flag
        "safety":          safety,
        "coverage":        coverage,
        "clinical_constraints": clinical,
        "trial_failure":   trial_failure,
        "directional_evidence": directional,
        "ctpa":            ctpa,
        "registry":        registry,
        "mechanism_scope": mechanism_scope,
        "calibration":     calibration,
        "drug_genes":      drug_genes[:20],
    }


# ── Main engine ───────────────────────────────────────────────────────────────

def run_repurposing_screen(
    disease_name: str,
    max_candidates: int = 50,
    db_compounds: Optional[List[Dict]] = None,
) -> Dict:
    """
    5-screen repurposing analysis for a disease.
    db_compounds: pre-fetched compounds from local DB (may be empty).
    Returns dict with ranked candidates + disease context.
    """
    cache_key = f"screen:{disease_name.lower().strip()}"
    cache = _load_cache()
    if cache_key in cache:
        return cache[cache_key]

    # ── Disease context from Open Targets ─────────────────────────────────────
    disease_info: Dict = {}
    disease_genes: List[str] = []
    disease_pathways: List[Dict] = []
    ppi_adj: Dict = {}

    try:
        from services.disease_ontology import resolve_disease as ot_resolve, get_ppi_network
        disease_info    = ot_resolve(disease_name)
        disease_genes   = [t["gene_symbol"] for t in disease_info.get("targets", [])[:40]]
        disease_pathways = disease_info.get("pathways", [])
        if disease_genes:
            ppi_adj = get_ppi_network(disease_genes[:15])
    except Exception as e:
        logger.warning(f"Open Targets lookup failed for '{disease_name}': {e}")

    # ── Build candidate pool ──────────────────────────────────────────────────
    seen_ids: set = set()
    candidates: List[Dict] = []

    if db_compounds:
        for c in db_compounds:
            cid = c.get("chembl_id") or str(c.get("id", ""))
            if cid and cid not in seen_ids:
                seen_ids.add(cid)
                candidates.append(c)

    # ChEMBL REST API fallback when DB has < 10 results
    if len(candidates) < 10:
        try:
            api_inds = _chembl_drug_indications(disease_name, limit=40)
            new_ids  = [x["chembl_id"] for x in api_inds if x["chembl_id"] not in seen_ids]
            mol_map  = _chembl_molecule_details(new_ids[:30])
            for ind in api_inds:
                mid = ind["chembl_id"]
                if mid in seen_ids or mid not in mol_map:
                    continue
                seen_ids.add(mid)
                mol = {**mol_map[mid]}
                mol["max_phase_for_ind"] = ind.get("max_phase_for_ind", 0)
                mol["indications"] = ind.get("indication", disease_name)
                candidates.append(mol)
        except Exception as e:
            logger.debug(f"ChEMBL API pool-build error: {e}")

    # ── Augment with NOVEL HetioNet-connected candidates ──────────────────────
    # Drugs the biology links to the disease (Compound→Gene←Disease) that are NOT
    # already indicated for it — genuine repurposing leads, not confirmations.
    # This is what makes the engine discover rather than merely confirm.
    try:
        from services.repurposing_scorer import hetionet_novel_compounds
        novel = hetionet_novel_compounds([disease_name], limit=25)
        if novel:
            by_name = _local_molecules_by_name([n["name"] for n in novel])
            added = 0
            for n in novel:
                mol = by_name.get(n["name"].lower())
                if not mol or not mol.get("smiles"):
                    continue
                cid = mol["chembl_id"]
                if cid in seen_ids:
                    continue
                seen_ids.add(cid)
                candidates.append({**mol, "source": "Knowledge-graph novelty",
                                   "novelty": True, "shared_disease_genes": n["shared_genes"]})
                added += 1
                if added >= 15:
                    break
    except Exception as e:
        logger.debug(f"novelty augmentation skipped: {e}")

    if not candidates:
        return {
            "disease": disease_name, "disease_info": disease_info,
            "candidates": [], "screens": WEIGHTS,
            "error": "No candidates found", "_ts": time.time(),
        }

    # ── Fetch drug genes for EVERY candidate that lacks targets ───────────────
    # (DB-sourced rows often have an empty `targets` field too — previously only
    # ChEMBL-API compounds, capped at 20, were resolved, which is why most
    # candidates scored target/pathway = 0). Batched, so covering the full set is cheap.
    need_genes = [
        c.get("chembl_id", "") for c in candidates[:max_candidates]
        if (c.get("chembl_id", "") or "").startswith("CHEMBL") and not (c.get("targets") or "").strip()
    ]
    api_gene_map = _chembl_targets_for_molecules(need_genes) if need_genes else {}

    # Mechanism action types (gene→INHIBITOR/AGONIST/…) for the candidate set —
    # feeds the direction-aware pathway score. Local chembl_33; {} if unavailable.
    all_cids = [c.get("chembl_id", "") for c in candidates[:max_candidates]
                if (c.get("chembl_id", "") or "").startswith("CHEMBL")]
    actions_map = _actions_for_molecules(all_cids) if all_cids else {}

    # Confounding-by-indication set: molecules with ANY development record (phase >= 1)
    # for THIS disease. A drug developed for the disease has it in its FAERS, so its
    # adverse-event overlap must NOT be read as drug-caused toxicity. One batched query.
    studied_set: set = set()
    try:
        from services.repurposing_scorer import approved_chembls_for_disease
        if all_cids:
            studied_set = approved_chembls_for_disease(all_cids, [], disease_name, min_phase=1)
    except Exception as e:
        logger.debug(f"studied-for-disease set skipped: {e}")

    # ── Score each candidate ──────────────────────────────────────────────────
    scored: List[Dict] = []
    for comp in candidates[:max_candidates]:
        cid = comp.get("chembl_id", "")
        targets_raw = comp.get("targets", "") or ""
        drug_genes_list = [g.strip() for g in re.split(r"[;,]", targets_raw) if g.strip()]
        if not drug_genes_list and cid in api_gene_map:
            drug_genes_list = api_gene_map[cid]

        # KG-novelty candidates were surfaced by a Compound→Gene←Disease path (shared
        # disease genes in HetioNet) whose genes are, by construction, outside the
        # disease's OT top-40 — so the OT-overlap sub-scores miss them. Pass that shared
        # -gene evidence as the mechanistic prior (same saturating scale the network
        # signal uses) so a genuine KG lead isn't scored as zero-mechanism / phantom.
        _prior = None
        if comp.get("novelty") and comp.get("shared_disease_genes"):
            _prior = min(0.85, 0.35 + 0.12 * int(comp["shared_disease_genes"]))
        sr = score_compound_for_disease(
            comp, disease_name, disease_genes, disease_pathways, ppi_adj, drug_genes_list,
            drug_actions=actions_map.get(cid), mechanistic_prior=_prior,
            studied_for_disease=(cid in studied_set),
        )
        sc = {**comp, **sr, "score": sr["composite_score"]}
        # Same quality filters as every other surface (plausibility + lead-viability;
        # disease is the fixed input here so its value is attached once on the result).
        try:
            from services.quality_overlay import overlay
            sc.update(overlay(comp.get("name", ""), cid, disease_name, drug_genes_list,
                              smiles=comp.get("smiles", ""), with_disease_value=False))
        except Exception as e:
            logger.debug(f"discover overlay skipped: {e}")
        scored.append(sc)

    # Rank by ACTIONABILITY: aligned (mechanistic + clinical + feasibility) leads first,
    # then by the actionability-weighted score — so a labelled hypothesis or a delivery-
    # limited pair sinks below a genuinely actionable one of equal raw mechanism. Nothing
    # is dropped; the demoted candidates stay visible with their reason.
    scored.sort(key=lambda x: (
        1 if x.get("actionability", {}).get("aligned") else 0,
        x.get("actionability", {}).get("rank_score", x.get("composite_score", 0)),
    ), reverse=True)

    _dv = None
    try:
        from services.disease_value import value_for
        _dv = value_for(disease_name)
    except Exception:
        pass
    result = {
        "disease":               disease_name,
        "disease_value":         _dv,      # is THIS disease worth a repurposing program?
        "disease_info":          disease_info,
        "disease_gene_count":    len(disease_genes),
        "disease_pathway_count": len(disease_pathways),
        "top_disease_genes":     disease_genes[:15],
        "candidates":            scored[:max_candidates],
        "screens":               WEIGHTS,
        "_ts":                   time.time(),
    }
    cache[cache_key] = result
    _save_cache(cache)
    return result
