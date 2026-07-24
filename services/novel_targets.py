"""
Novel-target discovery (P2)  —  from validation to discovery.
═══════════════════════════════════════════════════════════════════════════════
The platform otherwise only scores drugs against KNOWN Open Targets disease
associations — it CONFIRMS biology, it does not DISCOVER it. This module infers
targets NOT yet annotated to a disease, by network guilt-by-association:

    disease → known OT seed genes (weighted by genetic association)
            → STRING PPI partners of those seeds (incl. genes outside the seed set)
            → a non-seed gene wired to many / high-weight seeds is a candidate
              NOVEL target (a hidden node one hop off the known disease module)

Each inferred target carries its evidence (which seeds it connects to, at what
confidence) and a normalized confidence, so it is clearly labelled as INFERRED
(network-based), not an experimentally annotated association.

`drugs_via_novel_targets` then closes the loop: drugs that hit those inferred
targets become repurposing hypotheses reachable ONLY through the novel target —
the "give me a novel target AND results" ask.

FAIL-SOFT: no seeds / no PPI → empty result (never fabricates a target).
"""
from __future__ import annotations

import logging
from typing import Dict, List

import requests

from services import http_client

logger = logging.getLogger(__name__)

STRING_PARTNERS = "https://string-db.org/api/json/interaction_partners"
_STRING_MIN = 700           # STRING confidence floor (0-1000); 700 = high confidence
_SEED_CAP = 15              # seeds queried against STRING
_PARTNERS_PER_SEED = 15


def _string_partners(genes: List[str]) -> List[tuple]:
    """(seed, partner, score 0..1) — top PPI partners of each seed from STRING
    interaction_partners (partners MAY be outside the seed set = discovery)."""
    if not genes:
        return []
    try:
        r = requests.get(STRING_PARTNERS, params={
            "identifiers": "%0d".join(genes[:_SEED_CAP]),
            "species": 9606, "required_score": _STRING_MIN,
            "limit": _PARTNERS_PER_SEED,
            "caller_identity": "neurorepurpose_platform"}, timeout=15)
        if not r.ok:
            return []
        out = []
        for link in r.json():
            seed = link.get("preferredName_A", "")
            partner = link.get("preferredName_B", "")
            score = float(link.get("score", 0.0))       # already 0..1
            if seed and partner:
                out.append((seed, partner, score))
        return out
    except Exception as e:
        logger.debug(f"STRING interaction_partners failed: {e}")
        return []


def infer_novel_targets(disease_name: str, max_targets: int = 15) -> Dict:
    """Infer NOVEL (network-based) candidate targets for a disease via PPI
    guilt-by-association off the known OT seed genes. Returns ranked inferred
    targets with evidence + confidence. Fail-soft → empty."""
    result = {"disease": disease_name, "seed_count": 0, "novel_targets": []}
    try:
        from services.disease_ontology import get_disease_gene_weights
        seeds = get_disease_gene_weights(disease_name, top_n=25) or {}
    except Exception as e:
        logger.debug(f"seed fetch failed: {e}")
        seeds = {}
    seeds = {g.upper(): float(w) for g, w in seeds.items() if w and w > 0}
    if not seeds:
        return result
    result["seed_count"] = len(seeds)
    seed_set = set(seeds)

    partners = _string_partners(list(seeds))
    if not partners:
        return result

    # Aggregate: a candidate's support = Σ (seed genetic weight × PPI confidence)
    # over the seeds it connects to. Non-seed genes only (seeds are already known).
    cand: Dict[str, Dict] = {}
    for seed, partner, sc in partners:
        p = partner.upper()
        if p in seed_set:
            continue
        w = seeds.get(seed.upper(), 0.0)
        e = cand.setdefault(p, {"support": 0.0, "via": [], "best_ppi": 0.0})
        e["support"] += w * sc
        e["via"].append(seed.upper())
        e["best_ppi"] = max(e["best_ppi"], sc)

    if not cand:
        return result
    top = sorted(cand.items(), key=lambda kv: -kv[1]["support"])[:max_targets]
    max_support = top[0][1]["support"] or 1.0

    result["novel_targets"] = [{
        "gene":        g,
        "confidence":  round(min(1.0, d["support"] / max_support), 3),
        "seed_links":  len(set(d["via"])),
        "via_seeds":   sorted(set(d["via"]))[:6],
        "best_ppi":    round(d["best_ppi"], 3),
        "basis":       "inferred (PPI guilt-by-association)",
    } for g, d in top]
    return result


def drugs_via_novel_targets(disease_name: str, max_targets: int = 8,
                            max_drugs_per_target: int = 5) -> Dict:
    """Close the loop: drugs that hit the inferred novel targets → repurposing
    hypotheses reachable ONLY through a novel target. Each is clearly flagged as
    resting on an INFERRED (network-based) target, not a known association."""
    inf = infer_novel_targets(disease_name, max_targets=max_targets)
    _dv = None
    try:
        from services.disease_value import value_for
        _dv = value_for(disease_name)
    except Exception:
        pass
    out = {"disease": disease_name, "disease_value": _dv,
           "novel_targets": inf["novel_targets"], "drugs": []}
    seen: set = set()
    for nt in inf["novel_targets"][:max_targets]:
        gene = nt["gene"]
        for d in _drugs_for_gene(gene, limit=max_drugs_per_target):
            key = d["chembl_id"] or d["name"].lower()
            if key in seen:
                continue
            seen.add(key)
            drug = {
                **d,
                "via_novel_target": gene,
                "target_confidence": nt["confidence"],
                "basis": "inferred novel target (network-based, not a known association)",
            }
            # CCH — a toxic drug reached via a novel target is flagged/penalised here
            # too, not just on the repurpose page.
            try:
                from services.clinical_constraints import harmonize
                cc = harmonize(d.get("name", ""), disease_name)
                drug["clinical_constraints"] = cc
                drug["clinical_multiplier"] = cc.get("multiplier", 1.0)
            except Exception:
                drug["clinical_multiplier"] = 1.0
            # CTPA Rule 2 — registry ghost audit: drop drugs whose program in this
            # indication definitively failed / was abandoned (e.g. Tarenflurbil for
            # Alzheimer's) so we never surface a "repeat the failed trial" lead.
            try:
                from services.registry_audit import audit
                ga = audit(d.get("name", ""), disease_name)
                drug["registry"] = ga
                if ga.get("ghost"):
                    drug["clinical_multiplier"] = min(drug["clinical_multiplier"], ga["multiplier"])
            except Exception:
                pass
            # Same platform quality filters as the other surfaces (plausibility +
            # lead-viability); the novel target gene is the mechanistic hypothesis.
            try:
                from services.quality_overlay import overlay
                drug.update(overlay(d.get("name", ""), d.get("chembl_id", ""),
                                    disease_name, [gene], with_disease_value=False))
            except Exception:
                pass
            out["drugs"].append(drug)
    # Development-reality guardrails (same as the /discover path): a drug already APPROVED
    # for this disease is not a "novel-target discovery" even if it happens to hit an
    # inferred target, and a drug withdrawn for safety is not a viable lead. Exclude those
    # (surfaced separately) and flag me-too mechanisms; keep the inferred-target framing.
    try:
        from services import forward_guardrails as _fg
        part = _fg.apply(out["drugs"], disease_name)
        out["drugs"] = part["leads"]
        out["excluded"] = [{"name": c.get("name"), "chembl_id": c.get("chembl_id"),
                            "reason": c.get("removed_reason", ""),
                            "market_status": c.get("market_status")}
                           for c in part["removed"]]
    except Exception:
        out.setdefault("excluded", [])
    # Rank so clinically-viable drugs surface above crushed/ghost ones.
    out["drugs"].sort(key=lambda x: x.get("clinical_multiplier", 1.0), reverse=True)
    return out


def _drugs_for_gene(gene: str, limit: int = 5) -> List[Dict]:
    """Approved/clinical drugs whose mechanism targets a gene, via Open Targets."""
    q = """
    query($sym: String!) {
      search(queryString: $sym, entityNames: ["target"], page:{index:0,size:1}) {
        hits { id }
      }
    }"""
    try:
        from services.disease_ontology import _ot_graphql
    except Exception:
        return []
    try:
        hits = (_ot_graphql(q, {"sym": gene}).get("search", {}) or {}).get("hits", [])
        if not hits:
            return []
        ensembl = hits[0]["id"]
        # OT v4: drugs/clinical candidates for a target live under
        # drugAndClinicalCandidates → rows{ maxClinicalStage drug{id name} }.
        q2 = """
        query($id: String!) {
          target(ensemblId: $id) {
            drugAndClinicalCandidates {
              rows { maxClinicalStage drug { id name } }
            }
          }
        }"""
        rows = (((_ot_graphql(q2, {"id": ensembl}).get("target") or {})
                 .get("drugAndClinicalCandidates") or {}).get("rows") or [])
        # Highest clinical stage first — the most advanced, actionable candidates.
        rows.sort(key=lambda r: (r.get("maxClinicalStage") or 0), reverse=True)
        out, seen = [], set()
        for r in rows:
            d = r.get("drug") or {}
            cid = d.get("id", "")
            if not cid or cid in seen:
                continue
            seen.add(cid)
            out.append({"chembl_id": cid, "name": d.get("name", ""),
                        "max_phase": r.get("maxClinicalStage")})
            if len(out) >= limit:
                break
        return out
    except Exception as e:
        logger.debug(f"drugs_for_gene({gene}) failed: {e}")
        return []
