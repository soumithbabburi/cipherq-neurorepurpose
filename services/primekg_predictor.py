"""
PrimeKG predictor — broad-coverage link prediction + contraindication lookup.
════════════════════════════════════════════════════════════════════════════════
Serves queries against the DistMult model trained on PrimeKG (17,080 diseases, the graph
TxGNN is built on) and its labeled drug-disease edges. Two things our older KGs could not do:
  • broad `treat` plausibility + top-k over 17k diseases (vs Hetionet 137 / DRKG ~5k)
  • a first-class CONTRAINDICATION output, from PrimeKG's 61k real contraindication edges
    (labeled) augmented by the model where a pair is unlabeled (predicted).

Query resolution uses PrimeKG's own node space: drug by DrugBank id / name, disease by name
(exact, then normalized, then token-subset), gene by symbol. Fully fail-soft: if the model
has not been trained yet (no pkg_embeddings.npz) every call returns None/empty.
"""
from __future__ import annotations

import json
import logging
import re
from functools import lru_cache
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np

logger = logging.getLogger(__name__)

_D = Path(__file__).parent.parent / "data" / "primekg"
_state: Optional[Dict] = None
_loaded = False


def _norm(s: str) -> str:
    return re.sub(r"[^a-z0-9 ]", " ", (s or "").lower()).strip()


def _load() -> Optional[Dict]:
    global _state, _loaded
    if _loaded:
        return _state
    _loaded = True
    try:
        emb = np.load(_D / "pkg_embeddings.npz")
        rels = json.load(open(_D / "pkg_rels.json"))
        idx = json.load(open(_D / "pkg_index.json"))
        nodes = json.load(open(_D / "pkg_nodes.json"))
        labels = json.load(open(_D / "pkg_labels.json")) if (_D / "pkg_labels.json").exists() else {}
        disease_ids = np.array([int(i) for i, v in nodes.items() if v[0] == "disease"], dtype=np.int64)
        drug_ids = np.array([int(i) for i, v in nodes.items() if v[0] == "drug"], dtype=np.int64)
        # Direction-aware TREATS classifier (the precise generator). The raw DistMult is kept
        # only for the embedding features it feeds this model; it is NOT used for ranking.
        treats = None
        try:
            import joblib
            if (_D / "pkg_treats.pkl").exists():
                treats = joblib.load(_D / "pkg_treats.pkl")["model"]
        except Exception as e:
            logger.debug("treats classifier unavailable: %s", e)
        # R-GCN refined embeddings (the message-passing universe ranker). Fused with the
        # treats classifier in top_diseases_for_drug; validated to beat the cascade on the
        # held-out ranking harness. Fail-soft: needs torch + pkg_gnn.pt, else we fall back
        # to the classifier-only cascade. Kept as numpy for the same vectorised path.
        gnn_H = gnn_r = None
        try:
            if (_D / "pkg_gnn.pt").exists():
                import torch
                _blob = torch.load(_D / "pkg_gnn.pt", weights_only=False, map_location="cpu")
                gnn_H = _blob["H"].numpy(); gnn_r = _blob["r_ind"].numpy()
        except Exception as e:
            logger.debug("GNN ranker unavailable: %s", e)
        _state = {
            "E": emb["E"], "R": emb["R"], "rels": rels, "treats": treats,
            "gnn_H": gnn_H, "gnn_r": gnn_r,
            "drug_by_name": idx.get("drug_by_name", {}), "drug_by_dbid": idx.get("drug_by_dbid", {}),
            "disease_by_name": idx.get("disease_by_name", {}),
            "nodes": nodes, "labels": labels,
            "disease_ids": disease_ids, "drug_ids": drug_ids,
            # normalized disease-name index for fuzzy resolution
            "disease_norm": {_norm(k): v for k, v in idx.get("disease_by_name", {}).items()},
        }
        logger.info("PrimeKG predictor loaded (%s diseases, %s drugs)",
                    len(disease_ids), len(drug_ids))
    except Exception as e:
        logger.info("PrimeKG predictor unavailable (train it first): %s", e)
        _state = None
    return _state


def available() -> bool:
    return _load() is not None


def _name(idx: int) -> str:
    s = _load()
    return s["nodes"].get(str(idx), ["", "", "?"])[2] if s else "?"


def resolve_drug(drug: str, chembl_id: str = "") -> Optional[int]:
    s = _load()
    if not s or not drug:
        return None
    d = drug.strip()
    # DrugBank id direct
    if d.upper() in s["drug_by_dbid"]:
        return s["drug_by_dbid"][d.upper()]
    low = d.lower()
    if low in s["drug_by_name"]:
        return s["drug_by_name"][low]
    # DrugBank id via UniChem (reuse the DRKG resolver's mapping if present)
    return None


@lru_cache(maxsize=4096)
def resolve_disease(disease: str) -> Optional[int]:
    s = _load()
    if not s or not disease:
        return None
    low = disease.strip().lower()
    if low in s["disease_by_name"]:
        return s["disease_by_name"][low]
    n = _norm(disease)
    if n in s["disease_norm"]:
        return s["disease_norm"][n]
    # synonym bridge via our disease resolver (e.g. "chronic myeloid leukemia" ->
    # "chronic myelogenous leukemia")
    try:
        from services.disease_ontology import resolve_disease as _rd
        info = _rd(disease) or {}
        for cand in [info.get("name", "")] + (info.get("synonyms", []) or []):
            cl = (cand or "").strip().lower()
            if cl and cl in s["disease_by_name"]:
                return s["disease_by_name"][cl]
            cn = _norm(cand)
            if cn and cn in s["disease_norm"]:
                return s["disease_norm"][cn]
    except Exception:
        pass
    # last resort: token-subset containment against normalized names
    toks = set(n.split())
    if len(toks) >= 2:
        for kn, vi in s["disease_norm"].items():
            kt = set(kn.split())
            if toks <= kt or kt <= toks:
                return vi
    return None


def _rel_vec(name: str):
    s = _load()
    rid = s["rels"].get(name) if s else None
    return s["R"][rid] if (s and rid is not None) else None


def _treats_scores(d_idx: int, target_ids: np.ndarray):
    """Direction-aware P(treats) for (d_idx, indication, each target) via the supervised
    classifier over node-embedding pair features [Ed, Ez, Ed*Ez, |Ed-Ez|]. Vectorised."""
    s = _load()
    clf = s.get("treats") if s else None
    if clf is None:
        return None
    E = s["E"]; ed = E[d_idx]; Ez = E[target_ids]
    n = len(target_ids)
    feats = np.concatenate(
        [np.broadcast_to(ed, (n, ed.shape[0])), Ez, ed * Ez, np.abs(ed - Ez)], axis=1)
    return clf.predict_proba(feats)[:, 1]


def _treats_scores_drugs(z_idx: int, drug_ids: np.ndarray):
    """P(treats) for (each drug, indication, z_idx)."""
    s = _load()
    clf = s.get("treats") if s else None
    if clf is None:
        return None
    E = s["E"]; ez = E[z_idx]; Ed = E[drug_ids]
    n = len(drug_ids)
    feats = np.concatenate(
        [Ed, np.broadcast_to(ez, (n, ez.shape[0])), Ed * ez, np.abs(Ed - ez)], axis=1)
    return clf.predict_proba(feats)[:, 1]


def _score_vs_all(node_idx: int, rel: str, target_ids: np.ndarray):
    """DistMult scores of (node_idx, rel, target) for every target id (vectorised)."""
    s = _load()
    rv = _rel_vec(rel)
    if s is None or rv is None:
        return None
    h = s["E"][node_idx] * rv                       # elementwise (dim,)
    return h @ s["E"][target_ids].T                 # (n_targets,)


def treat_plausibility(drug: str, disease: str, chembl_id: str = "") -> Optional[float]:
    """Direction-aware P(treats), 0..1, from the supervised treats classifier (compound-
    disjoint AUC 0.98 vs real contraindications). Unlike the raw DistMult it distinguishes
    treat from contraindicate. Still a plausibility / candidate-generation signal that the
    composite + guardrails rank and gate, NOT a probability of clinical success."""
    s = _load()
    di = resolve_drug(drug, chembl_id)
    zi = resolve_disease(disease)
    if s is None or di is None or zi is None or s.get("treats") is None:
        return None
    sc = _treats_scores(di, np.array([zi], dtype=np.int64))
    return round(float(sc[0]), 4) if sc is not None else None


_POOL = 400   # recall pool size before direction-aware re-ranking


def _rrf_rerank(scores_list, weights, k_rrf: int = 60):
    """Reciprocal Rank Fusion of several score vectors over the SAME candidate pool.
    Returns a fused score per candidate (higher = better)."""
    fused = np.zeros(len(scores_list[0]))
    for sc, w in zip(scores_list, weights):
        order = np.argsort(-sc)
        rank = np.empty(len(sc)); rank[order] = np.arange(1, len(sc) + 1)
        fused += w * (1.0 / (k_rrf + rank))
    return fused


def top_diseases_for_drug(drug: str, k: int = 20, chembl_id: str = "",
                          min_score: float = 0.0) -> List[Dict]:
    """Universe generator: rank candidate indications for a drug over PrimeKG's 17k diseases.

    Method (validated on the held-out ranking harness, validation/validate_primekg_*):
    a DistMult relatedness pool (recall) is re-ranked by fusing the direction-aware treats
    classifier (pair features) with the R-GCN (graph structure) via Reciprocal Rank Fusion,
    weighted 3:1. This FUSION beats the classifier-only cascade on the compound-disjoint
    held-out set (R@20 0.392 vs 0.378, R@100 0.763 vs 0.754, median rank 34 vs 36) while
    tying it at R@10. Fail-soft: with no GNN artifact it degrades to the classifier-only
    cascade (the prior behaviour).

    CAVEAT (honest): even the fusion recovers only ~28% of a held-out drug's true indications
    in its top 10 — a strong CANDIDATE-GENERATION / triage signal, not a calibrated
    probability of success. The trustworthy use for a decision remains treat_plausibility()
    on a SPECIFIC pair the biology engine already surfaced."""
    s = _load()
    di = resolve_drug(drug, chembl_id)
    if s is None or di is None or s.get("treats") is None:
        return []
    rel = _score_vs_all(di, "indication", s["disease_ids"])   # DistMult relatedness (recall)
    if rel is None:
        return []
    pool = np.argsort(-rel)[:_POOL]
    pool_ids = s["disease_ids"][pool]
    prob = _treats_scores(di, pool_ids)                        # classifier (precision/direction)
    # Fuse with the R-GCN structural ranker when available (validated to beat the cascade).
    gnn_H, gnn_r = s.get("gnn_H"), s.get("gnn_r")
    if gnn_H is not None and gnn_r is not None:
        gnn = (gnn_H[di] * gnn_r * gnn_H[pool_ids]).sum(-1)   # R-GCN score over the pool
        fused = _rrf_rerank([prob, gnn], [3.0, 1.0])
        order = [j for j in np.argsort(-fused) if prob[j] >= min_score][:k]
        return [{"disease": _name(int(pool_ids[j])),
                 "node": int(pool_ids[j]),
                 "score": round(float(prob[j]), 4),   # displayed direction-aware P(treats)
                 "fused_rank": int(r + 1)} for r, j in enumerate(order)]
    order = [j for j in np.argsort(-prob) if prob[j] >= min_score][:k]  # fail-soft: classifier-only
    return [{"disease": _name(int(pool_ids[j])),
             "node": int(pool_ids[j]),
             "score": round(float(prob[j]), 4)} for j in order]


def top_drugs_for_disease(disease: str, k: int = 20) -> List[Dict]:
    """Drugs for a disease, same recall-then-precision cascade."""
    s = _load()
    zi = resolve_disease(disease)
    if s is None or zi is None or s.get("treats") is None:
        return []
    rv = _rel_vec("indication")
    rel = s["E"][s["drug_ids"]] @ (s["E"][zi] * rv)
    pool = np.argsort(-rel)[:_POOL]
    pool_ids = s["drug_ids"][pool]
    prob = _treats_scores_drugs(zi, pool_ids)
    order = np.argsort(-prob)[:k]
    return [{"drug": _name(int(pool_ids[j])),
             "node": int(pool_ids[j]),
             "score": round(float(prob[j]), 4)} for j in order]


def contraindications_for_drug(drug: str, k: int = 25, chembl_id: str = "") -> Dict:
    """Diseases the drug should NOT be repurposed into, from PrimeKG's REAL labeled
    `contraindication` edges (ground truth). Off-label edges are returned separately (they
    feed the commercial-friction layer).

    NOTE on method: validation showed the DistMult model recovers indications vs random
    diseases at AUC 0.99 but only 0.57 vs real contraindications — i.e. it captures
    drug-disease RELATEDNESS, not DIRECTION. So this view does NOT surface model-predicted
    contraindications (they would be near-random). Direction-aware contraindication reasoning
    for UNLABELED pairs is the job of services.mechanism_direction (drug action x disease gene
    direction), not this KG embedding."""
    s = _load()
    di = resolve_drug(drug, chembl_id)
    if s is None or di is None:
        return {"labeled": [], "off_label": [], "resolved": False}
    lab = s["labels"].get(str(di), {})
    labeled = [{"disease": _name(z), "node": z} for z in lab.get("contraindication", [])]
    off = [{"disease": _name(z), "node": z} for z in lab.get("off-label use", [])]
    return {"labeled": labeled[:k], "off_label": off[:k], "resolved": True, "drug": _name(di),
            "note": ("PrimeKG labeled contraindications (ground truth). For unlabeled pairs, "
                     "direction-aware contraindication risk comes from the mechanism-direction "
                     "engine, not this relatedness model.")}


def labeled_relation(drug: str, disease: str, chembl_id: str = "") -> Optional[str]:
    """'indication' | 'contraindication' | 'off-label use' | None for a (drug, disease) pair
    from PrimeKG's real edges — a ground-truth cross-check for any candidate."""
    s = _load()
    di = resolve_drug(drug, chembl_id); zi = resolve_disease(disease)
    if s is None or di is None or zi is None:
        return None
    lab = s["labels"].get(str(di), {})
    for r in ("indication", "contraindication", "off-label use"):
        if zi in lab.get(r, []):
            return r
    return None
