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

import requests

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

def _chembl_drug_indications(disease_name: str, limit: int = 50) -> List[Dict]:
    """Query ChEMBL for compounds with matching disease indication."""
    seen = set()
    results = []

    for param in ("mesh_heading__icontains", "efo_term__icontains"):
        if len(results) >= limit:
            break
        try:
            r = requests.get(
                f"{CHEMBL_BASE}/drug_indication.json",
                params={param: disease_name, "limit": limit, "format": "json"},
                timeout=12,
            )
            if r.ok:
                for ind in r.json().get("drug_indications", []):
                    mol_id = ind.get("molecule_chembl_id", "")
                    if mol_id and mol_id not in seen:
                        seen.add(mol_id)
                        results.append({
                            "chembl_id":         mol_id,
                            "indication":        ind.get("mesh_heading") or ind.get("efo_term") or disease_name,
                            "max_phase_for_ind": ind.get("max_phase_for_ind") or 0,
                        })
        except Exception as e:
            logger.debug(f"ChEMBL {param} error: {e}")

    return results


def _chembl_molecule_details(chembl_ids: List[str]) -> Dict[str, Dict]:
    """Batch-fetch molecule details from ChEMBL REST API."""
    results: Dict[str, Dict] = {}
    for i in range(0, min(len(chembl_ids), 60), 20):
        chunk = chembl_ids[i : i + 20]
        try:
            r = requests.get(
                f"{CHEMBL_BASE}/molecule.json",
                params={"molecule_chembl_id__in": ",".join(chunk), "limit": 20, "format": "json"},
                timeout=12,
            )
            if r.ok:
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
    """Batch-fetch gene targets for a list of ChEMBL molecule IDs."""
    drug_genes: Dict[str, List[str]] = {cid: [] for cid in chembl_ids}
    for cid in chembl_ids[:20]:  # limit API calls
        try:
            r = requests.get(
                f"{CHEMBL_BASE}/mechanism.json",
                params={"molecule_chembl_id": cid, "limit": 20, "format": "json"},
                timeout=8,
            )
            if not r.ok:
                continue
            for mech in r.json().get("mechanisms", []):
                target_id = mech.get("target_chembl_id", "")
                if not target_id:
                    continue
                tr = requests.get(f"{CHEMBL_BASE}/target/{target_id}.json", timeout=6)
                if not tr.ok:
                    continue
                for comp in (tr.json().get("target_components") or []):
                    for syn in (comp.get("target_component_synonyms") or []):
                        if syn.get("syn_type") == "GENE_SYMBOL":
                            g = syn.get("component_synonym", "")
                            if g:
                                drug_genes[cid].append(g)
        except Exception:
            continue
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


def _score_pathway_overlap(drug_genes: List[str], disease_pathways: List[Dict]) -> float:
    if not drug_genes or not disease_pathways:
        return 0.0
    dg = {g.upper() for g in drug_genes}
    hits = sum(1 for pw in disease_pathways[:20] if dg & {g.upper() for g in pw.get("genes", [])})
    return min(1.0, hits / max(1, min(10, len(disease_pathways))))


def _score_ppi_network(drug_genes: List[str], ppi_adj: Dict[str, List[str]]) -> float:
    if not drug_genes or not ppi_adj:
        return 0.0
    dg = {g.upper() for g in drug_genes}
    ppi_nodes = {g.upper() for g in ppi_adj}
    ppi_neighbors: set = set()
    for neighbors in ppi_adj.values():
        ppi_neighbors.update(g.upper() for g in neighbors)
    direct   = len(dg & ppi_nodes)
    neighbor = len(dg & ppi_neighbors)
    return min(1.0, (direct * 1.0 + neighbor * 0.5) / len(dg))


def _to_phase(max_phase) -> int:
    try:
        return int(float(max_phase or 0))
    except (TypeError, ValueError):
        return 0


def _score_clinical(max_phase, existing_ind: str, disease_name: str) -> float:
    phase_map = {4: 1.0, 3: 0.75, 2: 0.50, 1: 0.25}
    score = phase_map.get(_to_phase(max_phase), 0.05)
    disease_kws = [w for w in disease_name.lower().split() if len(w) > 4]
    if any(kw in (existing_ind or "").lower() for kw in disease_kws):
        score = min(1.0, score + 0.15)
    return score


def _score_indication(existing_ind: str, disease_name: str) -> float:
    stop = {"disease", "disorder", "syndrome", "with", "associated"}
    dw = {w.lower() for w in re.split(r'\W+', disease_name) if len(w) > 3 and w.lower() not in stop}
    iw = {w.lower() for w in re.split(r'\W+', existing_ind or "") if len(w) > 3 and w.lower() not in stop}
    if not dw:
        return 0.0
    return min(1.0, len(dw & iw) / len(dw))


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
) -> Dict:
    existing_ind = compound.get("indications", "") or ""
    raw_phase    = compound.get("max_phase") or compound.get("max_phase_for_ind") or 0
    try:
        max_phase = int(float(raw_phase))
    except (TypeError, ValueError):
        max_phase = 0
    if drug_genes is None:
        drug_genes = []

    scores = {
        "target":     _score_target_overlap(drug_genes, disease_genes),
        "pathway":    _score_pathway_overlap(drug_genes, disease_pathways),
        "ppi":        _score_ppi_network(drug_genes, ppi_adj),
        "clinical":   _score_clinical(max_phase, existing_ind, disease_name),
        "indication": _score_indication(existing_ind, disease_name),
        "regulatory": _score_regulatory(max_phase, disease_name, existing_ind),
    }
    composite = sum(scores[k] * WEIGHTS[k] for k in WEIGHTS)
    return {
        "composite_score": round(composite, 4),
        "scores":          {k: round(v, 4) for k, v in scores.items()},
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

    if not candidates:
        return {
            "disease": disease_name, "disease_info": disease_info,
            "candidates": [], "screens": WEIGHTS,
            "error": "No candidates found", "_ts": time.time(),
        }

    # ── Fetch drug genes for ChEMBL API compounds ─────────────────────────────
    api_cids = [
        c.get("chembl_id", "") for c in candidates
        if c.get("source") == "ChEMBL API" and not c.get("targets")
    ]
    api_gene_map = _chembl_targets_for_molecules(api_cids[:20]) if api_cids else {}

    # ── Score each candidate ──────────────────────────────────────────────────
    scored: List[Dict] = []
    for comp in candidates[:max_candidates]:
        cid = comp.get("chembl_id", "")
        targets_raw = comp.get("targets", "") or ""
        drug_genes_list = [g.strip() for g in re.split(r"[;,]", targets_raw) if g.strip()]
        if not drug_genes_list and cid in api_gene_map:
            drug_genes_list = api_gene_map[cid]

        sr = score_compound_for_disease(
            comp, disease_name, disease_genes, disease_pathways, ppi_adj, drug_genes_list,
        )
        sc = {**comp, **sr, "score": sr["composite_score"]}
        scored.append(sc)

    scored.sort(key=lambda x: x.get("composite_score", 0), reverse=True)

    result = {
        "disease":               disease_name,
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
