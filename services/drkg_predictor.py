"""
DRKG link-prediction predictor  —  modern-KG plausibility + candidate generation.
═══════════════════════════════════════════════════════════════════════════════════
Runtime over the DistMult embeddings trained on DRKG (services/drkg_train.py). Scores a
(Compound, treats, Disease) link the way the Rephetio/TxGNN lineage does, but on DRKG's
~5.1k diseases / ~23k compounds instead of Hetionet's 213 / 8894 — so newer drugs and
far more diseases get a knowledge-graph signal.

Two uses, mirroring services/repurposing_predictor so wiring is a drop-in:
  • treat_probability(drug, disease)         -> plausibility axis for a specific pair
  • top_diseases_for_drug / top_drugs_for_disease -> KG-as-candidate-generator

ID resolution (our platform IDs -> DRKG namespaces):
  drug   : chembl_id -(UniChem)-> DrugBank -> Compound::DB####   (name fallback)
  disease: mesh_id -> Disease::MESH:D####  |  doid -> Disease::DOID:######  (name fallback)

Fail-soft everywhere: if the embeddings aren't built or an ID can't be mapped, callers
get None / [] and fall back to the existing Hetionet path.
"""
from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Optional

import numpy as np

logger = logging.getLogger(__name__)

_DRKG = Path(__file__).parent.parent / "data" / "drkg"
_EMB = _DRKG / "drkg_embeddings.npz"
_NAMES = _DRKG / "drkg_names.json"
_UNICHEM = _DRKG / "unichem_chembl_drugbank.txt"

# "Treatment"-type Compound->Disease relations across DRKG's sources. We score against the
# mean of these so the signal is "does the graph place this pair like a known treatment".
_TREAT_RELS = (
    "DRUGBANK::treats::Compound:Disease",
    "Hetionet::CtD::Compound:Disease",
    "GNBR::T::Compound:Disease",
)

_MIN_DIS_DEG = 20   # min subgraph degree to rank a disease (drops noisy rare-concept nodes)
_MIN_CMP_DEG = 10   # min subgraph degree to rank a compound

_S: dict = {}          # lazy singleton state
_LOADED = False


def _load() -> bool:
    global _LOADED
    if _LOADED:
        return bool(_S)
    _LOADED = True
    if not _EMB.exists():
        logger.info("drkg predictor: embeddings not built — using Hetionet path")
        return False
    try:
        d = np.load(_EMB, allow_pickle=True)
        E, R = d["E"].astype(np.float32), d["R"].astype(np.float32)
        nodes = [str(x) for x in d["nodes"]]
        rels = [str(x) for x in d["rels"]]
        node2id = {n: i for i, n in enumerate(nodes)}
        rel2id = {r: i for i, r in enumerate(rels)}

        treat_ids = [rel2id[r] for r in _TREAT_RELS if r in rel2id]
        if not treat_ids:
            logger.warning("drkg predictor: no treat relation in embeddings")
            return False
        treat_vec = R[treat_ids].mean(axis=0)

        dis_idx = np.array([i for i, n in enumerate(nodes) if n.startswith("Disease::")], dtype=np.int64)
        cmp_idx = np.array([i for i, n in enumerate(nodes) if n.startswith("Compound::")], dtype=np.int64)

        # Restrict the RANKING universe to sufficiently-connected nodes: low-degree DRKG
        # nodes (rare MeSH C-record supplementary concepts, one-edge compounds) have
        # under-trained embeddings that otherwise dominate the top of every list with
        # noise. Well-connected diseases keep the real signal (RA, the proliferative
        # cancers) while dropping the junk. Degrees come from the built subgraph.
        deg = None
        if "deg" in d.files:
            deg = d["deg"]
        else:  # older embeddings without baked degrees — fall back to the build graph
            try:
                g = np.load(_DRKG / "drkg_graph.npz", allow_pickle=True)["triples"]
                deg = np.bincount(np.concatenate([g[:, 0], g[:, 2]]), minlength=len(nodes))
            except Exception as e:
                logger.debug("drkg predictor: degree load failed (%s); no degree filter", e)
        if deg is not None:
            dis_scored = dis_idx[deg[dis_idx] >= _MIN_DIS_DEG]
            cmp_scored = cmp_idx[deg[cmp_idx] >= _MIN_CMP_DEG]
        else:
            dis_scored, cmp_scored = dis_idx, cmp_idx

        names = {}
        if _NAMES.exists():
            names = json.loads(_NAMES.read_text(encoding="utf-8"))
        # reverse name -> node (lowercased) for the name fallback
        name2node = {}
        for n in nodes:
            nm = names.get(n)
            if nm:
                name2node.setdefault(nm.strip().lower(), n)
        # Enrich disease name resolution with MeSH synonyms (entry_terms) so a query like
        # "rheumatoid arthritis" resolves to Disease::MESH:D001172 ("Arthritis, Rheumatoid").
        try:
            import psycopg2
            from config import db_params
            conn = psycopg2.connect(**db_params())
            cur = conn.cursor()
            cur.execute("SELECT mesh_id, heading, entry_terms FROM mesh_diseases")
            for mid, heading, terms in cur.fetchall():
                node = f"Disease::MESH:{mid}"
                if node not in node2id:
                    continue
                # OVERWRITE (not setdefault): MeSH disease nodes are far better connected
                # than DRKG's duplicate DOID nodes (177k vs 1.7k edges), so a name should
                # resolve to the MeSH node — otherwise 'rheumatoid arthritis' hits the
                # sparse DOID node and scores near-zero instead of its true rank.
                for t in ([heading] + list(terms or [])):
                    if t:
                        name2node[t.strip().lower()] = node
            conn.close()
        except Exception as e:  # pragma: no cover — synonym enrichment is best-effort
            logger.debug("drkg predictor: mesh synonym enrichment skipped: %s", e)

        chembl2db = {}
        if _UNICHEM.exists():
            for ln in _UNICHEM.read_text(encoding="utf-8").splitlines()[1:]:
                p = ln.split("\t")
                if len(p) == 2:
                    chembl2db[p[0].strip().upper()] = p[1].strip()

        _S.update(E=E, R=R, treat_vec=treat_vec, nodes=nodes, node2id=node2id,
                  dis_idx=dis_idx, cmp_idx=cmp_idx, dis_scored=dis_scored,
                  cmp_scored=cmp_scored, names=names, name2node=name2node, chembl2db=chembl2db)
        logger.info("drkg predictor: loaded (%d/%d scored diseases, %d/%d scored compounds)",
                    len(dis_scored), len(dis_idx), len(cmp_scored), len(cmp_idx))
        return True
    except Exception as e:  # pragma: no cover
        logger.warning("drkg predictor load failed: %s", e)
        _S.clear()
        return False


def available() -> bool:
    return _load()


# ── ID resolution ────────────────────────────────────────────────────────────────

def _resolve_drug(name: Optional[str], chembl_id: Optional[str]) -> Optional[int]:
    n2i = _S["node2id"]
    if chembl_id:
        db = _S["chembl2db"].get(str(chembl_id).strip().upper())
        if db and f"Compound::{db}" in n2i:
            return n2i[f"Compound::{db}"]
    if name:
        node = _S["name2node"].get(name.strip().lower())
        if node and _S["nodes"][n2i[node]].startswith("Compound::"):
            return n2i[node]
    return None


def _resolve_disease(name: Optional[str], mesh_id: Optional[str],
                     doid: Optional[str]) -> Optional[int]:
    n2i = _S["node2id"]
    if mesh_id:
        m = str(mesh_id).strip()
        m = m.split(":")[-1]  # 'MESH:D000544' or 'D000544' -> 'D000544'
        if f"Disease::MESH:{m}" in n2i:
            return n2i[f"Disease::MESH:{m}"]
    if doid:
        d = str(doid).strip().replace("DOID:", "")
        if f"Disease::DOID:{d}" in n2i:
            return n2i[f"Disease::DOID:{d}"]
    if name:
        node = _S["name2node"].get(name.strip().lower())
        if node and _S["nodes"][n2i[node]].startswith("Disease::"):
            return n2i[node]
    return None


def _label(node: str) -> str:
    nm = _S["names"].get(node) or node.split("::", 1)[-1]
    # MeSH headings are comma-inverted ("Arthritis, Rheumatoid"). De-invert disease names
    # to the natural form so the downstream composite resolver (which keys on common names
    # for target/pathway/mechanism signals) recognises them — otherwise a real candidate
    # like RA collapses to a near-zero score purely on a name-format mismatch.
    if node.startswith("Disease::") and ", " in nm and "(" not in nm:
        parts = [p.strip() for p in nm.split(",") if p.strip()]
        nm = " ".join(reversed(parts))   # "Leukemia, Myeloid, Acute" -> "Acute Myeloid Leukemia"
    return nm


# ── scoring ──────────────────────────────────────────────────────────────────────

def _treat_scores(query_idx: int, tgt_idx: np.ndarray) -> np.ndarray:
    """DistMult treat-link scores for one node against a set of targets: (E[q]⊙treat)·E[t]."""
    return _S["E"][tgt_idx] @ (_S["E"][query_idx] * _S["treat_vec"])


def treat_probability(drug: Optional[str] = None, disease: Optional[str] = None,
                      chembl_id: Optional[str] = None, mesh_id: Optional[str] = None,
                      doid: Optional[str] = None) -> Optional[float]:
    """KG plausibility that `drug` treats `disease`: the percentile of the pair's treat
    score within the drug's score distribution over the connected disease universe.
    None if either side can't be mapped to DRKG."""
    if not _load():
        return None
    di = _resolve_drug(drug, chembl_id)
    si = _resolve_disease(disease, mesh_id, doid)
    if di is None or si is None:
        return None
    scores = _treat_scores(di, _S["dis_scored"])
    pair = float(_treat_scores(di, np.array([si]))[0])
    return round(float((scores < pair).mean()), 4)


def top_diseases_for_drug(drug: Optional[str] = None, chembl_id: Optional[str] = None,
                          k: int = 30, min_p: float = 0.90) -> list[dict]:
    """Rank the connected DRKG diseases for a drug (the Rephetio/TxGNN inversion).
    Returns [{disease, node, identifier, probability}], probability = within-drug percentile."""
    if not _load():
        return []
    di = _resolve_drug(drug, chembl_id)
    if di is None:
        return []
    dis_scored = _S["dis_scored"]
    scores = _treat_scores(di, dis_scored)
    order = np.argsort(scores)[::-1]
    n = len(scores)
    out = []
    for rank, j in enumerate(order[: k * 3]):
        p = 1.0 - rank / n
        if p < min_p and len(out) >= k:
            break
        node = _S["nodes"][dis_scored[j]]
        out.append({"disease": _label(node), "node": node,
                    "identifier": node.split("::", 1)[-1], "probability": round(float(p), 4)})
        if len(out) >= k:
            break
    return out


def top_drugs_for_disease(disease: Optional[str] = None, mesh_id: Optional[str] = None,
                          doid: Optional[str] = None, k: int = 30,
                          min_p: float = 0.90) -> list[dict]:
    """Rank the connected DRKG compounds for a disease. [{drug, node, identifier, probability}]."""
    if not _load():
        return []
    si = _resolve_disease(disease, mesh_id, doid)
    if si is None:
        return []
    cmp_idx = _S["cmp_scored"]
    scores = _treat_scores(si, cmp_idx)
    order = np.argsort(scores)[::-1]
    n = len(scores)
    out = []
    for rank, j in enumerate(order[:k]):
        node = _S["nodes"][cmp_idx[j]]
        out.append({"drug": _label(node), "node": node,
                    "identifier": node.split("::", 1)[-1],
                    "probability": round(1.0 - rank / n, 4)})
    return out
