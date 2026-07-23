"""
Delivery-feasibility gate for LOCALIZED indications.
═══════════════════════════════════════════════════════════════════════════════════
A mechanistically-plausible pair is only clinically actionable if the drug can PHYSICALLY
REACH the tissue where the disease lives. For a compartmentalized indication — brain (BBB),
eye, skin, joint, local airway — a systemically-dosed molecule that engages the right
target still does nothing if it never arrives at an effective concentration.

This runs a physical-access check BEFORE a pair is allowed to rank as actionable:

  • CNS / brain / CNS-tumor  → the BBB is the barrier. Score BBB penetrance from CNS-MPO
    (physchem). No easy local route exists (you cannot topically dose the brain), so a
    non-penetrant molecule is genuinely delivery-limited.
  • eye / skin / joint / local-airway → a LOCAL route EXISTS (drops, topical, intra-
    articular injection, inhaler). The molecule is deliverable, but usually only after
    REFORMULATION — a real development step, so a mild flag, not a physical block.
  • solid tumor (non-CNS) / systemic disease → reached by the circulation; no barrier.

Bounded, fail-soft, DEMOTE-not-drop: returns a multiplier in [0.6, 1.0] plus a
human-readable flag. Missing SMILES → cannot assess → neutral (1.0) with an honest note,
never an invented penalty.
"""
from __future__ import annotations

import logging
import re
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)

# compartment → (name-match keywords). Order matters: CNS-tumor must beat generic "tumor".
_CNS = ("brain", "cerebr", "cerebell", "cortical", "glioma", "glioblast", "cns ",
        "central nervous", "alzheimer", "parkinson", "huntington", "dementia",
        "epилeps", "epilep", "multiple sclerosis", "amyotrophic", "als ", "meningi",
        "neurodegener", "encephal", "intracranial", "medulloblast", "astrocytoma")
_OCULAR = ("ocular", "retina", "retinal", "macular", "glaucoma", "uveitis", "corneal",
           "conjunctiv", "eye ", "ophthalm", "diabetic macular")
_DERMAL = ("skin", "dermal", "dermat", "psoriasis", "atopic", "eczema", "acne",
           "cutaneous", "epidermal", "vitiligo", "rosacea")
_JOINT = ("intra-articular", "synovial", "joint ", "osteoarthritis", "cartilage")
_AIRWAY = ("asthma", "copd", "chronic obstructive", "bronch", "airway", "cystic fibrosis",
           "idiopathic pulmonary fibrosis", "pulmonary fibrosis")

# a local route exists for these — deliverable after reformulation, not a physical wall.
_LOCAL_ROUTE = {
    "ocular":       "topical / intravitreal (eye drops or injection)",
    "dermal":       "topical (cream / gel)",
    "intra-articular": "intra-articular injection",
    "local-airway": "inhaled (nebuliser / DPI)",
}


def _compartment(disease: str) -> str:
    d = " " + (disease or "").lower() + " "
    if any(k in d for k in _CNS):
        return "cns"
    if any(k in d for k in _OCULAR):
        return "ocular"
    if any(k in d for k in _DERMAL):
        return "dermal"
    if any(k in d for k in _JOINT):
        return "intra-articular"
    if any(k in d for k in _AIRWAY):
        return "local-airway"
    return "systemic"


def assess_delivery(disease_name: str, smiles: str = "", drug_name: str = "",
                    props: Optional[Dict] = None,
                    cns_indicated: bool = False) -> Dict:
    """Delivery feasibility for a (drug, disease) pair. Returns
    {localized, compartment, deliverable, route, multiplier, flag, assessed}."""
    out = {"localized": False, "compartment": "systemic", "deliverable": True,
           "route": "systemic", "multiplier": 1.0, "flag": "", "assessed": False}

    comp = _compartment(disease_name)
    out["compartment"] = comp
    if comp == "systemic":
        return out                      # reached by the circulation — no barrier to check
    out["localized"] = True

    # ── compartments with a real LOCAL delivery route: deliverable, mild dev-step flag ──
    if comp in _LOCAL_ROUTE:
        out["assessed"] = True
        out["route"] = _LOCAL_ROUTE[comp]
        out["deliverable"] = True
        out["multiplier"] = 0.9
        out["flag"] = (f"Localized ({comp}) indication — reachable via {_LOCAL_ROUTE[comp]}; "
                       f"systemic dosing may under-expose the site, so a local formulation "
                       f"is likely required.")
        return out

    # ── CNS: the BBB is the barrier and there is no easy local route ──────────────────
    if comp == "cns":
        out["route"] = "systemic → blood-brain barrier"
        if cns_indicated:               # already an approved CNS drug — it crosses.
            out["assessed"] = True
            out["deliverable"] = True
            out["flag"] = "CNS indication — drug already has CNS exposure precedent."
            return out
        verdict, mpo = _bbb(smiles, props)
        if verdict is None:
            out["flag"] = ("CNS indication — BBB penetrance UNASSESSED (no structure); "
                           "confirm brain exposure before treating as actionable.")
            return out                  # fail-soft: no invented penalty
        out["assessed"] = True
        if verdict == "low":
            out["deliverable"] = False
            out["multiplier"] = 0.6
            out["flag"] = (f"CNS indication but the molecule is unlikely to cross the BBB "
                           f"(CNS-MPO {mpo:.1f}/6) — delivery-limited regardless of target fit.")
        elif verdict == "borderline":
            out["deliverable"] = True
            out["multiplier"] = 0.8
            out["flag"] = (f"CNS indication with BORDERLINE BBB penetrance (CNS-MPO {mpo:.1f}/6) "
                           f"— brain exposure uncertain; verify before ranking as actionable.")
        else:                           # good / high
            out["deliverable"] = True
            out["flag"] = f"CNS indication — favourable BBB penetrance (CNS-MPO {mpo:.1f}/6)."
        return out

    return out


def _bbb(smiles: str, props: Optional[Dict]):
    """→ (verdict in {'low','borderline','good'}|None, cns_mpo score). None = unassessable."""
    try:
        from services.cns_mpo import cns_mpo
        r = cns_mpo(props=props, smiles=smiles) or {}
        score = r.get("score")
        if score is None:
            return None, 0.0
        if score < 3.0:
            return "low", float(score)
        if score < 4.0:
            return "borderline", float(score)
        return "good", float(score)
    except Exception as e:
        logger.debug(f"BBB assess failed: {e}")
        return None, 0.0
