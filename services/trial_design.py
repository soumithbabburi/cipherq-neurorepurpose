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

    # ── enrichment: the drug's targets as candidate biomarkers, ranked by whether they
    #    also drive THIS disease (a driver = strong enrichment; else a subset marker).
    enr = []
    for g in genes:
        assoc = drivers.get(g, 0.0)
        enr.append({"gene": g, "assoc": round(assoc, 2),
                    "basis": "disease driver" if assoc >= 0.30 else "target-expressing subset",
                    "criterion": f"{g}-altered / high-{g} disease (amplification, overexpression, or activating alteration)"})
    enr.sort(key=lambda x: -x["assoc"])
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

    # ── one-line design sketch
    marker = out["enrichment"][0]["gene"] if out["enrichment"] else None
    tier = (evidence_tier or {}).get("tier", "")
    phase = "Phase 2 (proof-of-concept)" if tier in ("mechanistic", "preclinical", "literature", "promising") else "Phase 2/3"
    parts = [phase]
    if marker:
        parts.append(f"biomarker-enriched ({marker}+)")
    else:
        parts.append("unselected — no enrichment marker found (higher risk)")
    if out["exclusions"]:
        parts.append(f"excluding {out['exclusions'][0]['criterion']}")
    out["design"] = " · ".join(parts)

    # ── plain-language rationale
    r = []
    if marker and out["enrichment"][0]["assoc"] >= 0.30:
        r.append(f"The drug's target {marker} is a driver of {disease} (association {out['enrichment'][0]['assoc']}), "
                 f"so enriching for {marker}-altered patients raises the chance of a signal — this is exactly why an "
                 f"unselected trial can fail while a marker-defined subgroup responds.")
    elif marker:
        r.append(f"No target is a top disease driver, but testing in the {marker}-expressing subset (the drug's "
                 f"primary target) is a more informative first trial than an unselected population.")
    if out["exclusions"]:
        r.append("Its FAERS profile flags " + ", ".join(e["criterion"] for e in out["exclusions"][:2])
                 + " — exclude those patients up front.")
    out["rationale"] = " ".join(r)
    return out
