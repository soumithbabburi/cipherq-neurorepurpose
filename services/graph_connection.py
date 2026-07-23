"""
Drug ⇄ Indication connection (Knowledge Graph rationale).
════════════════════════════════════════════════════════════════════════════════
Answers the question the graph should answer but did not: given a molecule and a
candidate indication, HOW are they connected? We compute the real mechanistic bridge
and return three things:

  1. SHARED TARGETS — the intersection of the drug's protein targets and the disease's
     evidence-associated genes. This is the mechanistic path: drug → (inhibits/binds)
     → shared gene → (associated with) → disease.
  2. PATHWAYS — disease pathways that the shared targets participate in, so the
     connection is placed in biological context.
  3. A grounded NARRATIVE — a plain-language rationale and description built ONLY from
     the facts above (shared genes, pathways, association scores). The local LLM
     (llama3.1:8b) phrases it when available; a deterministic template is used
     otherwise. The model is never given latitude to introduce a mechanism we did not
     compute.

Plus a small cytoscape-style subgraph (drug, shared-target genes, disease) so the UI
can draw the exact path instead of the drug's whole neighbourhood.

Grounded, fail-soft: an unreachable source degrades the rationale, never fabricates it.
"""
from __future__ import annotations

import logging
from typing import Dict, List

logger = logging.getLogger(__name__)


def _norm(g: str) -> str:
    return (g or "").strip().upper()


def _drug_targets(drug: str) -> Dict:
    """Resolve a drug (name or ChEMBL id) to {name, chembl_id, genes[]}."""
    try:
        from services.reverse_repurposing import resolve_drug
        info = resolve_drug(drug) or {}
        genes = []
        for t in info.get("targets", []) or []:
            g = t.get("gene_symbol") if isinstance(t, dict) else t
            if g:
                genes.append(str(g))
        nm = info.get("name") or drug
        # DB names come back upper-cased; title-case for prose but leave any short
        # all-caps token (an acronym) alone.
        nm = " ".join(w if (w.isupper() and len(w) <= 5) else w.capitalize()
                      for w in str(nm).split())
        return {"name": nm,
                "chembl_id": info.get("chembl_id") or "",
                "genes": genes}
    except Exception as e:
        logger.debug("drug target resolve failed: %s", e)
        return {"name": drug, "chembl_id": "", "genes": []}


def connect(drug: str, disease: str, *, use_llm: bool = True) -> Dict:
    """Compute the mechanistic connection between a drug and a candidate indication."""
    d = _drug_targets(drug)
    drug_name = d["name"]
    drug_genes = d["genes"]
    drug_up = {_norm(g) for g in drug_genes}

    disease_genes: List[str] = []
    weights: Dict[str, float] = {}
    pathways: List = []
    try:
        from services.disease_ontology import (
            get_disease_gene_set, get_disease_gene_weights, get_disease_pathways,
        )
        disease_genes = get_disease_gene_set(disease, top_n=40) or []
        weights = {k: v for k, v in (get_disease_gene_weights(disease) or {}).items()}
        pathways = get_disease_pathways(disease) or []
    except Exception as e:
        logger.debug("disease evidence lookup failed: %s", e)

    dis_up = {_norm(g) for g in disease_genes}

    # The mechanistic bridge: targets the drug hits that the disease's genetics implicate.
    shared = []
    for g in drug_genes:
        gu = _norm(g)
        if gu in dis_up:
            shared.append({"gene": gu, "weight": round(float(weights.get(gu, 0.0)), 3)})
    shared.sort(key=lambda s: s["weight"], reverse=True)

    # Pathways the shared targets sit in (place the connection in biological context).
    shared_up = {s["gene"] for s in shared}
    ctx_pathways = []
    for p in pathways:
        pname = p.get("name") if isinstance(p, dict) else str(p)
        pgenes = {_norm(x) for x in (p.get("genes", []) if isinstance(p, dict) else [])}
        hit = sorted(shared_up & pgenes) if pgenes else []
        if pname and (hit or not ctx_pathways):
            ctx_pathways.append({"name": pname, "genes": hit})
    ctx_pathways = [p for p in ctx_pathways if p["genes"]][:6] or ctx_pathways[:4]

    connected = bool(shared)
    facts = {
        "drug": drug_name,
        "disease": disease,
        "shared_targets": shared,
        "n_drug_targets": len(drug_genes),
        "n_disease_genes": len(disease_genes),
        "pathways": ctx_pathways,
        "connected": connected,
    }
    facts["rationale"] = _deterministic_rationale(facts)
    facts["description"] = _describe(facts, use_llm=use_llm)
    facts["elements"] = _subgraph(facts, d["chembl_id"])
    facts["drug_targets"] = drug_genes[:20]
    facts["disease_genes"] = disease_genes[:20]
    return facts


def _deterministic_rationale(f: Dict) -> str:
    """A grounded one-line rationale built purely from the computed facts."""
    drug, disease = f["drug"], f["disease"]
    shared = f["shared_targets"]
    if not shared:
        if f["n_drug_targets"] and f["n_disease_genes"]:
            return (f"{drug} and {disease} share no direct protein target in the current "
                    f"evidence. The drug acts on {f['n_drug_targets']} target(s); none is "
                    f"among the {f['n_disease_genes']} genes this indication implicates, so any "
                    f"link would be indirect (shared pathway or downstream signalling).")
        return (f"Insufficient target evidence to connect {drug} to {disease} from the "
                f"current sources.")
    genes = ", ".join(s["gene"] for s in shared[:5])
    lead = shared[0]["gene"]
    pw = f["pathways"][0]["name"] if f["pathways"] else ""
    base = (f"{drug} acts on {genes}, which the evidence associates with {disease}. "
            f"{lead} is the strongest shared target")
    base += f", and it participates in {pw}." if pw else "."
    return base


def _describe(f: Dict, *, use_llm: bool) -> str:
    """A short plain-language description of the mechanistic connection. Phrased by the
    local LLM from the grounded facts when available; deterministic template otherwise."""
    det = _deterministic_description(f)
    if not use_llm or not f["connected"]:
        return det
    try:
        from services import llm
        if not llm.available():
            return det
        shared = ", ".join(s["gene"] for s in f["shared_targets"][:6]) or "none"
        pws = "; ".join(p["name"] for p in f["pathways"][:4]) or "none identified"
        system = (
            "You are a pharmacology analyst. Write two or three plain sentences for a drug "
            "repurposing dossier, explaining ONLY the mechanistic connection given. Do not "
            "introduce any target, pathway, or claim that is not in the facts. Do not claim "
            "clinical efficacy. Do not use emojis, symbols, markdown, hyphens, or dashes. "
            "Write in plain professional prose."
        )
        prompt = (
            f"Drug: {f['drug']}\n"
            f"Candidate indication: {f['disease']}\n"
            f"Shared protein targets (drug targets that are also disease associated genes): {shared}\n"
            f"Relevant pathways containing those shared targets: {pws}\n\n"
            f"Explain how {f['drug']} could be mechanistically relevant to {f['disease']} "
            f"through these shared targets and pathways. Keep it to two or three sentences."
        )
        out = llm.generate(prompt, system=system, max_tokens=220, temperature=0.2, timeout=75)
        return out or det
    except Exception as e:
        logger.debug("LLM description failed: %s", e)
        return det


def _deterministic_description(f: Dict) -> str:
    drug, disease = f["drug"], f["disease"]
    shared = f["shared_targets"]
    if not shared:
        return _deterministic_rationale(f)
    genes = [s["gene"] for s in shared[:6]]
    glist = ", ".join(genes[:-1]) + (" and " + genes[-1] if len(genes) > 1 else genes[0])
    txt = (f"{drug} engages {glist}, "
           f"{'a target' if len(genes) == 1 else 'targets'} that current evidence links to "
           f"{disease}. Because the same "
           f"{'protein sits' if len(genes) == 1 else 'proteins sit'} on both the drug's "
           f"mechanism and the disease's biology, modulating "
           f"{'it' if len(genes) == 1 else 'them'} is a plausible route to affect the disease.")
    if f["pathways"]:
        pw = f["pathways"][0]
        txt += (f" These targets act within {pw['name']}, "
                f"placing the connection in a defined biological pathway.")
    return txt


def _subgraph(f: Dict, chembl_id: str) -> List[Dict]:
    """Cytoscape-style elements for the exact drug → shared target → disease path."""
    els: List[Dict] = []
    drug_id = "drug::" + (chembl_id or f["drug"])
    dis_id = "disease::" + f["disease"]
    els.append({"data": {"id": drug_id, "label": f["drug"], "full_label": f["drug"],
                         "kind": "Compound", "size": 42, "anchor": True}})
    els.append({"data": {"id": dis_id, "label": f["disease"], "full_label": f["disease"],
                         "kind": "Disease", "size": 40, "anchor": True}})
    for s in f["shared_targets"][:12]:
        gid = "gene::" + s["gene"]
        w = 0.5 + 0.5 * min(1.0, s["weight"])
        els.append({"data": {"id": gid, "label": s["gene"], "full_label": s["gene"],
                             "kind": "Gene", "size": 22}})
        els.append({"data": {"source": drug_id, "target": gid, "label": "targets", "weight": 0.8}})
        els.append({"data": {"source": gid, "target": dis_id, "label": "associated with",
                             "weight": round(w, 3)}})
    return els
