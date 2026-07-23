"""
Positive Pivot engine  —  turn a clinical MISMATCH into a discovery lead.
═══════════════════════════════════════════════════════════════════════════════
The Clinical Constraint Harmonizer (CCH) crushes a biologically strong hit when it
is clinically unviable (a cytotoxic drug for a benign chronic disease, a
non-penetrant drug for a CNS indication…). A pure filter stops there. This engine
instead uses the SAME CCH verdict to PIVOT the hypothesis toward a form that could
work — the difference between "what won't work" and "here is the viable angle":

  • Severity / chronicity crush → PRECISION pivot: a high-severity / orphan variant
    of the disease where the toxicity is justified (e.g. fulminant vs general MS).
  • CNS-barrier crush          → ANALOGUE pivot: the limiting physicochemical
    property to modify (scaffold-hopping direction), no fabricated molecule.
  • Chronic-tolerability crush  → COMBINATION pivot: a low-toxicity co-drug on the
    same target for a dose-sparing combination.

Every pivot is explicitly a HYPOTHESIS TO VALIDATE, not a de-risked asset — the
platform generates leads, it does not certify them. FAIL-SOFT throughout.
"""
from __future__ import annotations

import logging
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)

_SEVERE_MARKERS = ("acute", "malignant", "fulminant", "aggressive", "severe",
                   "rapidly progressive", "blast", "metasta", "advanced",
                   "refractory", "crisis", "fatal")


def _severe_variants(disease_name: str, max_n: int = 4) -> List[Dict]:
    """High-severity / orphan variants of a disease from Open Targets search, where
    a toxic drug's risk could be clinically justified. Fail-soft → []."""
    try:
        from services.disease_ontology import _ot_graphql
    except Exception:
        return []
    q = """
    query($q: String!) {
      search(queryString: $q, entityNames: ["disease"], page: {index: 0, size: 12}) {
        hits { id name object { __typename ... on Disease { therapeuticAreas { name } } } }
      }
    }"""
    seen, out = set(), []
    for term in (f"acute {disease_name}", f"malignant {disease_name}", disease_name):
        try:
            hits = _ot_graphql(q, {"q": term}).get("search", {}).get("hits", [])
        except Exception:
            hits = []
        for h in hits:
            obj = h.get("object") or {}
            if obj.get("__typename") != "Disease":
                continue
            nm = h.get("name", "")
            low = nm.lower()
            if h.get("id") in seen or low == (disease_name or "").lower():
                continue
            # a genuine severe variant: name carries a severity marker AND shares a
            # word with the base disease (so "acute leukemia" for "leukemia", not noise)
            base_tokens = {w for w in (disease_name or "").lower().split() if len(w) > 3}
            if any(m in low for m in _SEVERE_MARKERS) and (base_tokens & set(low.split())):
                seen.add(h.get("id"))
                out.append({"name": nm, "efo_id": h.get("id", "")})
            if len(out) >= max_n:
                return out
    return out


def _limiting_cns_property(smiles: str) -> Optional[Dict]:
    """Which physicochemical property most limits CNS penetration, and the fix
    direction — the scaffold-modification target for an analogue. None if fine."""
    if not smiles:
        return None
    try:
        from services.cns_mpo import cns_mpo
        r = cns_mpo(smiles=smiles) or {}
    except Exception:
        return None
    comps = r.get("components") or r.get("desirability") or {}
    # Prefer explicit weak components if the module exposes them; else derive from props.
    props = r.get("props") or r.get("properties") or {}
    tpsa = props.get("tpsa") or props.get("psa")
    hbd = props.get("hbd") or props.get("h_donors")
    mw = props.get("mw") or props.get("mol_weight")
    limiters = []
    if tpsa and tpsa > 90:
        limiters.append(("polar surface area", f"TPSA {round(tpsa)} > 90 Å²",
                         "reduce polar surface area / mask H-bond donors"))
    if hbd and hbd > 3:
        limiters.append(("H-bond donors", f"{hbd} donors > 3",
                         "cap or bioisosterically replace H-bond donors"))
    if mw and mw > 450:
        limiters.append(("molecular weight", f"MW {round(mw)} > 450",
                         "trim the scaffold to lower molecular weight"))
    if not limiters:
        return {"summary": "CNS-MPO borderline; consider reducing efflux (P-gp/BCRP) liability",
                "targets": []}
    return {"summary": "; ".join(x[1] for x in limiters),
            "targets": [{"property": p, "issue": i, "fix": f} for p, i, f in limiters]}


def _safe_codrug(disease_genes: List[str], exclude_drug: str = "") -> Optional[Dict]:
    """A LOW-toxicity drug hitting a disease driver target — a dose-sparing
    combination partner. Fail-soft → None."""
    try:
        from services.novel_targets import _drugs_for_gene
        from services.clinical_constraints import drug_risk_profile
    except Exception:
        return None
    for gene in (disease_genes or [])[:5]:
        for d in _drugs_for_gene(gene, limit=6):
            nm = d.get("name", "")
            if not nm or nm.lower() == (exclude_drug or "").lower():
                continue
            rp = drug_risk_profile(nm)
            if not rp["cytotoxic"] and rp["risk"] < 0.3:      # genuinely tolerable
                return {"name": nm, "chembl_id": d.get("chembl_id", ""),
                        "on_target": gene, "risk": rp["risk"]}
    return None


def generate_pivots(drug_name: str, disease_name: str, *, cch: Optional[Dict] = None,
                    smiles: str = "", disease_genes: Optional[List[str]] = None) -> List[Dict]:
    """From a CCH verdict, generate positive pivot strategies (hypotheses). []
    when the candidate is not clinically penalised (nothing to pivot)."""
    if not cch or not cch.get("penalized"):
        return []
    factors = cch.get("factors", {})
    pivots: List[Dict] = []

    # PRECISION — severity / chronicity mismatch → severe/orphan variant.
    if factors.get("severity", 1.0) < 0.8 or factors.get("duration", 1.0) < 0.8:
        variants = _severe_variants(disease_name)
        pivots.append({
            "type": "precision",
            "title": "Pivot to a high-severity / orphan variant",
            "rationale": (f"{drug_name} is too toxic for general {disease_name}, but its "
                          "risk is justified where the disease is life-threatening."),
            "suggestion": (f"Re-target the aggressive / rare variant of {disease_name} "
                           "where the therapeutic window opens."),
            "candidates": variants,
            "hypothesis": True,
        })
        # COMBINATION — chronic tolerability → low-dose rescue with a safe co-drug.
        if factors.get("duration", 1.0) < 0.8:
            co = _safe_codrug(disease_genes or [], exclude_drug=drug_name)
            pivots.append({
                "type": "combination",
                "title": "Pivot to a dose-sparing combination",
                "rationale": (f"{drug_name}'s mechanism is on-target but its full dose is "
                              "intolerable for chronic use."),
                "suggestion": (f"Combine an ultra-low dose with {co['name']} (low-toxicity, "
                               f"acts on {co['on_target']}) for a synergistic, dose-sparing "
                               "blockade." if co else
                               "Screen the pathway for a low-toxicity partner enabling a "
                               "dose-sparing combination."),
                "candidates": [co] if co else [],
                "hypothesis": True,
            })

    # ANALOGUE — CNS-barrier mismatch → scaffold modification direction.
    if factors.get("tissue", 1.0) < 0.85:
        lim = _limiting_cns_property(smiles)
        pivots.append({
            "type": "analogue",
            "title": "Pivot to a brain-penetrant analogue",
            "rationale": (f"{drug_name} engages the target but does not reach the CNS at "
                          "therapeutic exposure."),
            "suggestion": ("Scaffold-hop toward an analogue with the same pharmacophore but "
                           + (lim["summary"] if lim else "improved CNS exposure") + "."),
            "modify": (lim.get("targets") if lim else []),
            "hypothesis": True,
        })
    return pivots
