"""
Trial participant-selection assistant  —  turn a repurposing candidate into a
population-selection recommendation.
═══════════════════════════════════════════════════════════════════════════════════
For a (drug, disease) candidate, propose WHO to enrol so a trial has the best chance
of showing a signal — the piece KG-plausibility models never give you:

  • enrichment biomarker(s): the drug's OWN targets, prioritised by whether they also
    DRIVE this disease → the subgroup most likely to respond. This is exactly why an
    UNSELECTED trial can fail while a marker-defined subgroup works (nintedanib failed
    in non-selected soft-tissue sarcoma).
  • safety exclusions: from the drug's serious FAERS reactions → who to keep out.
  • dose / exposure note: from the lead-viability potency × window read.
  • enrolment feasibility: from the disease Value Score / burden.
  • a one-line suggested design.

Fail-soft everywhere: any missing input drops its section, never errors. This is a
decision aid, not a protocol — every field is sourced.
"""
from __future__ import annotations

import logging
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)

# Serious-AE keyword → exclusion criterion (curated, standard oncology/PK exclusions).
_EXCLUSION_MAP = [
    (("hypertension", "blood pressure"),                         "uncontrolled hypertension"),
    (("haemorrhage", "hemorrhage", "bleeding", "epistaxis", "haematemesis"), "recent major bleeding / therapeutic anticoagulation"),
    (("thromboemb", "embolism", "thrombosis", "infarction", "stroke", "ischaemic"), "recent arterial/venous thromboembolic event"),
    (("perforation", "fistula"),                                  "GI perforation / fistula history"),
    (("hepat", "liver", "transaminase", "hepatotox", "bilirubin", "hepatic failure"), "hepatic impairment (Child-Pugh B/C)"),
    (("proteinuria", "nephro", "renal failure", "creatinine"),    "significant renal impairment / proteinuria"),
    (("neutropenia", "thrombocytopenia", "anaemia", "anemia", "cytopenia", "febrile"), "baseline cytopenia below protocol threshold"),
    (("cardiac", "qt", "ejection fraction", "cardiomyopathy", "cardiac failure"), "significant cardiac disease / QTc prolongation"),
    (("interstitial", "pneumonitis", "ild", "respiratory failure"), "pre-existing interstitial lung disease"),
    (("wound", "healing"),                                        "major surgery within 4 weeks / impaired wound healing"),
    (("hypothyroid", "thyroid"),                                  "uncontrolled thyroid dysfunction (monitor)"),
]


def _drug_targets(drug_name: str, chembl_id: str, drug_genes: Optional[List[str]]) -> List[str]:
    if drug_genes:
        return [g.upper() for g in drug_genes]
    try:
        from services.reverse_repurposing import resolve_drug
        info = resolve_drug(chembl_id or drug_name) or {}
        return [g.upper() for g in (info.get("targets") or [])]
    except Exception:
        return []


def _disease_drivers(disease: str) -> Dict[str, float]:
    out: Dict[str, float] = {}
    try:
        from services.disease_ontology import resolve_disease
        d = resolve_disease(disease) or {}
        for t in d.get("targets", []):
            g = (t.get("gene_symbol") or "").upper()
            if g:
                out[g] = max(out.get(g, 0.0),
                             float(t.get("score") or t.get("quality_score") or 0.0))
    except Exception as e:
        logger.debug(f"disease drivers failed: {e}")
    return out


def participant_selection(drug_name: str, disease: str, chembl_id: str = "",
                          drug_genes: Optional[List[str]] = None, efo_id: str = "",
                          evidence_tier: Optional[Dict] = None) -> Dict:
    """Population-selection recommendation for a (drug, disease) candidate. Sourced,
    fail-soft. Returns {enrichment, exclusions, dose_note, feasibility, design, rationale}."""
    out: Dict = {"drug": drug_name, "disease": disease, "enrichment": [],
                 "exclusions": [], "dose_note": "", "feasibility": {}, "design": "",
                 "rationale": ""}

    genes = _drug_targets(drug_name, chembl_id, drug_genes)
    drivers = _disease_drivers(disease)
    try:
        from services.disease_appropriateness import infer_drug_action
        _action = infer_drug_action(drug_name, genes)
    except Exception:
        _action = ""
    inhibitor = _action == "inhibitor"

    # ── enrichment: score each of the drug's targets as a candidate biomarker on a
    #    3-tier support model — an association is a HYPOTHESIS, not a validated marker.
    #    support = biology (target↔disease assoc) + directional fit; the top tier is gated
    #    on real clinical-response evidence, which we do not have a source for, so it never
    #    fires (deliberately — we never claim a marker is clinically validated on biology alone).
    enr = []
    for g in genes:
        assoc = drivers.get(g, 0.0)
        biology = min(1.0, assoc)
        if assoc >= 0.30 and inhibitor:
            direction = 1.0          # inhibiting a disease driver is directionally corrective
        elif assoc >= 0.30:
            direction = 0.6          # relevant target, action direction unclear
        else:
            direction = 0.4
        evidence = 0.0               # no marker-level clinical-response source available
        support = round(0.55 * biology + 0.30 * direction + 0.15 * evidence, 2)
        if support >= 0.70 and evidence >= 0.5:
            tier, label = "clinically-supported", "Clinically validated subgroup"
        elif support >= 0.52 and assoc >= 0.30:
            tier, label = "enrichment-candidate", "Plausible enrichment · mechanistic"
        else:
            tier, label = "hypothesis", "Hypothesis-generating"
        enr.append({"gene": g, "assoc": round(assoc, 2), "support": support,
                    "tier": tier, "label": label,
                    "criterion": f"{g}-altered / high-{g} (amplification, overexpression, or activating alteration)"})
    enr.sort(key=lambda x: -x["support"])
    out["enrichment"] = enr[:5]

    # ── safety exclusions from the drug's serious FAERS reactions
    try:
        from services.safety_filter import _faers_serious_reactions
        rx = _faers_serious_reactions(drug_name) or []
        total = sum(c for _, c in rx) or 1
        seen = set()
        for term, cnt in rx:
            tl = term.lower()
            for keys, excl in _EXCLUSION_MAP:
                if excl in seen:
                    continue
                if any(k in tl for k in keys) and cnt / total >= 0.01:
                    out["exclusions"].append({"criterion": excl, "faers_signal": term,
                                              "share": round(cnt / total, 3)})
                    seen.add(excl)
    except Exception as e:
        logger.debug(f"exclusions failed: {e}")

    # ── dose / exposure read from the lead-viability funnel
    try:
        from services.lead_viability import assess as _lv
        v = _lv(drug_name, chembl_id, disease, list(drivers.keys())[:10] or genes)
        out["dose_note"] = v.get("verdict") or v.get("summary") or ""
    except Exception as e:
        logger.debug(f"dose note failed: {e}")

    # ── enrolment feasibility from the disease Value Score / burden
    try:
        from services.disease_value import value_for
        dv = value_for(disease, efo_id) or {}
        if dv:
            burden = dv.get("burden", dv.get("burden_score"))
            market = (dv.get("market") or dv.get("market_fit") or "")
            rare = ("orphan" in str(market).lower()) or (isinstance(burden, (int, float)) and burden < 0.30)
            out["feasibility"] = {"value_score": dv.get("value_score"), "burden": burden,
                                  "note": "orphan / low-prevalence — expect slow accrual, plan multi-site"
                                          if rare else "prevalent — single-region accrual feasible"}
    except Exception as e:
        logger.debug(f"feasibility failed: {e}")

    # ── safety gate: a biomarker can be biologically valid but the drug still a poor fit
    try:
        from services.safety_filter import assess as _safety
        sg = _safety(drug_name, disease) or {}
        if sg.get("penalized"):
            out["safety_gate"] = {"compatible": False, "flags": (sg.get("flags") or [])[:2],
                                  "note": "drug toxicity may conflict with this indication/population — verify before enrichment"}
        else:
            out["safety_gate"] = {"compatible": True, "note": "no blocking toxicity–indication conflict detected"}
    except Exception:
        out["safety_gate"] = {}

    # ── a HYPOTHESIS does not trigger enrichment; only an enrichment-candidate+ does.
    strong = [e for e in out["enrichment"] if e["tier"] in ("enrichment-candidate", "clinically-supported")]
    marker = strong[0] if strong else None
    tier = (evidence_tier or {}).get("tier", "")
    phase = "Phase 2 (proof-of-concept)" if tier in ("mechanistic", "preclinical", "literature", "promising") else "Phase 2/3"
    parts = [phase]
    parts.append(f"enrich for {marker['gene']}+ (mechanistic hypothesis)" if marker
                 else "unselected — no marker rises above hypothesis-generating")
    if out["exclusions"]:
        parts.append(f"exclude {out['exclusions'][0]['criterion']}")
    out["design"] = " · ".join(parts)

    # ── plain-language rationale — honest, no subgroup-response claim
    r = []
    if marker:
        r.append(f"{marker['gene']} is biologically relevant to {disease} (association {marker['assoc']}) and the drug's "
                 f"action is directionally appropriate, so a {marker['gene']}-enriched arm is a defensible mechanistic "
                 f"hypothesis — NOT a validated biomarker. There is no subgroup-response data here, so treat it as "
                 f"hypothesis-generating, not a guarantee of benefit.")
    elif out["enrichment"]:
        r.append(f"The drug's targets appear in {disease} biology but none rises above hypothesis-generating "
                 f"(association only, no directional or clinical support). An unselected design or translational "
                 f"biomarker work should precede any enrichment claim.")
    if out.get("safety_gate", {}).get("compatible") is False:
        r.append("Safety gate: " + out["safety_gate"]["note"] + ".")
    if out["exclusions"]:
        r.append("Exclude patients with " + ", ".join(e["criterion"] for e in out["exclusions"][:2]) + ".")
    out["rationale"] = " ".join(r)
    return out
