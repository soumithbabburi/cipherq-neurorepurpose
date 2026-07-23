"""
Directional literature evidence (P3)  —  typed triples, not co-occurrence.
═══════════════════════════════════════════════════════════════════════════════
Plain literature co-mention cannot tell whether a drug TREATS or CAUSES a disease,
nor trace a mechanism. This module uses DRKG's directional, literature-derived
typed triples (GNBR gene-disease relations from PubMed, DGIDB drug-target actions,
DrugBank drug-treats-disease) to answer a relational question:

    is there a DIRECTIONAL path  Compound —targets→ Gene ←associated— Disease ,
    or a direct  Compound —treats→ Disease  edge?

That is real mechanistic support ("drug modulates gene X, gene X is implicated in
the disease"), distinct from two entities appearing in the same abstract.

COVERAGE: the bundled DRKG is the neurological subset (46 diseases / 400
compounds), so this fires for neuro drugs/diseases and is FAIL-SOFT (covered=False,
signal 0) elsewhere — never a fabricated signal.
"""
from __future__ import annotations

import logging
from typing import Dict, Optional

logger = logging.getLogger(__name__)

_idx: Optional[Dict] = None


def _index() -> Dict:
    global _idx
    if _idx is not None:
        return _idx
    _idx = {}
    try:
        from services.biocypher_graph import _get_drkg
        drkg = _get_drkg()
        if not drkg or "edges" not in drkg:
            return _idx
        name2comp: Dict[str, str] = {}
        for cid, meta in (drkg.get("compounds") or {}).items():
            nm = (meta.get("name") or "").upper()
            if nm:
                name2comp[nm] = cid
            ch = (meta.get("chembl_id") or "").upper()
            if ch:
                name2comp[ch] = cid
        edges = drkg["edges"]
        comp_genes: Dict[str, Dict[str, str]] = {}
        for e in edges.get("target", []):
            comp_genes.setdefault(e["h"], {})[e["t"]] = e["r"]
        dis_genes: Dict[str, Dict[str, str]] = {}
        for e in edges.get("dis_gene", []):
            dis_genes.setdefault(e["h"], {})[e["t"]] = e["r"]
        treat = {(e["h"], e["t"]) for e in edges.get("treat", [])}
        _idx = {"name2comp": name2comp, "dni": drkg.get("disease_name_index", {}),
                "comp_genes": comp_genes, "dis_genes": dis_genes, "treat": treat,
                "genes": drkg.get("genes", {})}
    except Exception as e:
        logger.debug(f"DRKG index build failed: {e}")
        _idx = {}
    return _idx


def _resolve_disease(idx: Dict, disease_name: str) -> Optional[str]:
    dni = idx.get("dni", {})
    d = (disease_name or "").lower().strip()
    if d in dni:
        return dni[d]
    try:
        from services.disease_id import same_disease
        for name, did in dni.items():
            if same_disease(disease_name, name):
                return did
    except Exception:
        pass
    return None


def _gene_symbol(idx: Dict, gene_id: str) -> str:
    meta = (idx.get("genes") or {}).get(gene_id) or {}
    return meta.get("symbol") or meta.get("name") or gene_id.replace("Gene::", "gene ")


def directional_evidence(drug_name: str, disease_name: str) -> Dict:
    """Directional DRKG evidence for a drug-disease pair. Returns
    {covered, direct_treat, n_paths, path_genes, signal 0..1, note}. Fail-soft."""
    out = {"covered": False, "direct_treat": False, "n_paths": 0,
           "path_genes": [], "signal": 0.0, "note": ""}
    idx = _index()
    if not idx:
        return out
    cid = idx["name2comp"].get((drug_name or "").upper())
    did = _resolve_disease(idx, disease_name)
    if not cid or not did:
        # HONESTY (audited 2026-07): the loaded DRKG is the NEUROLOGY subset, so most
        # non-neuro drugs/diseases (mepolizumab, imatinib…) fall outside it. Say so
        # explicitly — n_paths=0 here means "not evaluated / out of coverage", NOT
        # "no mechanism exists". The broader KG-plausibility model (DRKG DistMult,
        # ~5k diseases) is the coverage-wide KG signal for these pairs.
        out["out_of_coverage"] = True
        out["note"] = ("Not evaluated — this pair is outside the loaded directional-KG "
                       "subset (neurology). Absence of a path here is not evidence of no "
                       "mechanism; see the KG-plausibility signal for wider coverage.")
        return out
    out["covered"] = True
    out["direct_treat"] = (cid, did) in idx["treat"]
    cg = idx["comp_genes"].get(cid, {})
    dg = idx["dis_genes"].get(did, {})
    shared = list(set(cg) & set(dg))
    out["n_paths"] = len(shared)
    out["path_genes"] = [_gene_symbol(idx, g) for g in shared[:8]]
    # Signal: a curated drug→treats→disease edge is the strongest; otherwise a
    # mechanistic drug→gene→disease path, saturating with the number of paths.
    if out["direct_treat"]:
        out["signal"] = 1.0
        out["note"] = "Directional DRKG edge: drug treats disease (DrugBank)."
    elif shared:
        out["signal"] = round(min(1.0, 0.35 + 0.18 * len(shared)), 3)
        out["note"] = (f"Directional path via {', '.join(out['path_genes'][:3])} "
                       f"(drug targets the gene; gene linked to disease in literature).")
    else:
        out["note"] = "In directional-KG coverage, but no drug→gene→disease path found here."
    return out
