"""
KG-embedding scoring + Hetionet crosswalk
═══════════════════════════════════════════════════════════════════════════
Turns the trained DistMult embedding (data/kg_embeddings.npz) into a usable
(drug, disease) → "treats" score for the live engine, by crosswalking the
engine's identifiers to Hetionet nodes:

  • Compound : by name (8.9k Hetionet compound names) or DrugBank id.
  • Disease  : by normalised name against Hetionet's 137 disease names
               (de-pluralising token match + a small alias map for the platform's
               headline diseases that name differently, e.g. glioblastoma →
               "malignant glioma").

Coverage is intentionally bounded — Hetionet has only 137 diseases — so the score
is returned as None when either side does not map. The engine then ensembles the
KG score where available and falls back to the mechanistic score otherwise
(never worse). Validated as a *complementary* signal: gene-overlap + KG ensemble
beats either alone by +0.042 AUROC on held-out treat edges (KGE-03).
"""
from __future__ import annotations

import re
import logging
from pathlib import Path
from typing import Dict, Optional

import numpy as np

logger = logging.getLogger(__name__)

_EMB = Path(__file__).parent.parent / "data" / "kg_embeddings.npz"
TREATS = "CtD"

# Platform diseases whose Hetionet name differs enough to miss fuzzy matching.
_DISEASE_ALIASES = {
    "glioblastoma": "malignant glioma",
    "glioblastoma multiforme": "malignant glioma",
    "glioma": "malignant glioma",
    "als": "amyotrophic lateral sclerosis",
    "lou gehrig disease": "amyotrophic lateral sclerosis",
    "epilepsy": "epilepsy syndrome",
}

_model = None
_node2id: Dict[str, int] = {}
_rel2id: Dict[str, int] = {}
_cmp_by_name: Dict[str, str] = {}     # normalised name -> full node id
_cmp_by_db: Dict[str, str] = {}       # DrugBank id -> full node id
_dis_by_norm: Dict[str, str] = {}     # normalised disease name -> full node id
_loaded = False
_available = False


_COMPOUND_ALIASES = {
    "aspirin": "acetylsalicylic acid",
    "acetaminophen": "paracetamol",
    "vitamin c": "ascorbic acid",
}


def _norm(s: str) -> str:
    s = re.sub(r"[^a-z0-9 ]", " ", (s or "").lower())
    # drop single-char tokens (e.g. the stray "s" from "Parkinson's") and de-pluralise
    toks = [t[:-1] if (len(t) > 3 and t.endswith("s")) else t
            for t in s.split() if len(t) > 1]
    return " ".join(toks).strip()


def _ensure() -> bool:
    """Lazy-load embeddings + build the crosswalk indices. Fail-soft."""
    global _model, _node2id, _rel2id, _loaded, _available
    if _loaded:
        return _available
    _loaded = True
    try:
        if not _EMB.exists():
            logger.info("kg_score: no embeddings file; KG signal disabled.")
            return False
        from services.kg_embedding import DistMultKGE
        _model, _node2id, _rel2id = DistMultKGE.load(_EMB)
        if TREATS not in _rel2id:
            return False
        import psycopg2
        from config import db_params
        conn = psycopg2.connect(**db_params())
        cur = conn.cursor()
        cur.execute("SELECT id, name, kind FROM hetionet_nodes WHERE kind IN ('Compound','Disease')")
        for nid, name, kind in cur.fetchall():
            if nid not in _node2id:          # only embedded nodes are scorable
                continue
            if kind == "Compound":
                if name:
                    _cmp_by_name.setdefault(_norm(name), nid)
                if "::" in nid:
                    _cmp_by_db[nid.split("::", 1)[1]] = nid
            elif kind == "Disease" and name:
                _dis_by_norm.setdefault(_norm(name), nid)
        conn.close()
        _available = bool(_cmp_by_name and _dis_by_norm)
        logger.info(f"kg_score: loaded {len(_cmp_by_name)} compounds, "
                    f"{len(_dis_by_norm)} diseases.")
        return _available
    except Exception as e:
        logger.warning(f"kg_score load failed: {e}")
        return False


def is_available() -> bool:
    return _ensure()


def resolve_compound(name: str = "", drugbank_id: str = "") -> Optional[str]:
    if not _ensure():
        return None
    if drugbank_id and drugbank_id in _cmp_by_db:
        return _cmp_by_db[drugbank_id]
    n = _norm(name)
    if n in _cmp_by_name:
        return _cmp_by_name[n]
    alias = _COMPOUND_ALIASES.get(n) or _COMPOUND_ALIASES.get((name or "").lower().strip())
    if alias:
        return _cmp_by_name.get(_norm(alias))
    return None


def resolve_disease(name: str) -> Optional[str]:
    if not _ensure():
        return None
    n = _norm(name)
    if n in _dis_by_norm:
        return _dis_by_norm[n]
    # alias map (normalise the alias target too)
    if name and name.lower().strip() in _DISEASE_ALIASES:
        return _dis_by_norm.get(_norm(_DISEASE_ALIASES[name.lower().strip()]))
    a = _DISEASE_ALIASES.get(n)
    if a:
        hit = _dis_by_norm.get(_norm(a))
        if hit:
            return hit
    # de-pluralised token-set match (handles "parkinson disease" vs "parkinsons disease")
    qtok = set(n.split())
    if qtok:
        best, best_j = None, 0.0
        for k, v in _dis_by_norm.items():
            kt = set(k.split())
            inter = len(qtok & kt)
            if not inter:
                continue
            j = inter / len(qtok | kt)
            if j > best_j:
                best, best_j = v, j
        if best_j >= 0.67:
            return best
    # substring fallback
    for k, v in _dis_by_norm.items():
        if n and (n in k or k in n):
            return v
    return None


def treats_score(drug_name: str = "", disease: str = "",
                 drugbank_id: str = "") -> Optional[float]:
    """Raw DistMult treats-score for (drug, disease), or None if either side
    does not map to an embedded Hetionet node. Higher = more likely a treatment."""
    if not _ensure():
        return None
    c = resolve_compound(drug_name, drugbank_id)
    d = resolve_disease(disease)
    if not c or not d:
        return None
    try:
        ci, di, ti = _node2id[c], _node2id[d], _rel2id[TREATS]
        return float(_model.score(np.array([ci]), np.array([ti]), np.array([di]))[0])
    except Exception:
        return None


if __name__ == "__main__":
    print("available:", is_available())
    for drug, dis in [("Ropinirole", "Parkinson disease"),
                      ("Memantine", "Alzheimer disease"),
                      ("Imatinib", "glioblastoma"),
                      ("Aspirin", "asthma")]:
        print(f"  {drug:12} / {dis:22} -> {treats_score(drug, dis)}")
