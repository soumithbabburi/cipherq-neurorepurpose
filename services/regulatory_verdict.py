"""
Regulatory & novelty verdict — single source of truth
═══════════════════════════════════════════════════════════════════════════════
Reconciles two subsystems that used to contradict each other in the dossier:
  • evidence_dossier._regulatory_assessment  → "505(b)(2) STRONGLY ELIGIBLE"
  • class_comparison._novelty                → "Novelty failure — already indicated"

The bug those produced: a drug could be recommended as a 505(b)(2) candidate AND
declared non-novel in the same dossier, because novelty was computed from ChEMBL
`known_indications`, which lists every indication studied at ANY phase (not just
APPROVED ones). A low-phase / drug-interaction co-mention (e.g. itraconazole
appearing in BACE-inhibitor PK studies that are tagged "Alzheimer Disease") was
read as "already approved for Alzheimer's." It is not.

This module computes novelty PHASE-AWARE and emits one coherent verdict:

  approved_here     (best phase ≥ 4)  → genuinely not novel; same-use claim blocked
  investigated_here (best phase 1-3)  → NEW for approval; 505(b)(2) viable; prior
                                        clinical art → secure a method-of-use claim
  prior_signal      (phase 0 / weak)  → early/PK/preclinical co-mention only; novel
  novel             (no match)        → clean repurposing hypothesis

Two facts are kept strictly separate (their conflation was the original sin):
  1. Is the MOLECULE approved (for any indication)?  → enables citing its safety
     package via 505(b)(2)/hybrid pathways. Good thing.
  2. Is it approved FOR THIS indication?             → novelty for the new claim.
     Bad thing only if (1)-style "approved" is mistaken for (2).

Jurisdiction-aware: maps the US 505(b)(2) route to its EU and India equivalents.
"""
from __future__ import annotations

import re
from typing import Dict, List, Optional

# Generic disease words that must NOT, on their own, count as an indication match
# (so "Cancer" doesn't match "Lung Cancer", but "Alzheimer" matches "Alzheimer Disease").
_GENERIC = {
    "disease", "disorder", "syndrome", "cancer", "tumor", "tumour", "carcinoma",
    "chronic", "acute", "primary", "secondary", "type", "stage", "condition",
    "neoplasm", "infection", "deficiency", "failure", "injury", "of", "the", "and",
}

# US route → its established analogues elsewhere.
_PATHWAYS = [
    ("US",    "505(b)(2) NDA",                                "FDA"),
    ("EU",    "Hybrid application (Dir. 2001/83, Art 10(3))", "EMA"),
    ("India", "New Drug application (CDSCO, approved-abroad route)", "CDSCO"),
]


def _tokens(name: str) -> set:
    return {w for w in re.split(r"[^a-z0-9]+", (name or "").lower())
            if len(w) > 3 and w not in _GENERIC}


def _matches(disease: str, indication: str) -> bool:
    """True if an indication name refers to the same disease — by specific-token
    overlap, falling back to exact string equality for token-less names."""
    a, b = _tokens(disease), _tokens(indication)
    if not a or not b:
        return (disease or "").lower().strip() == (indication or "").lower().strip()
    return bool(a & b)


def novelty(known_indications: List, disease_name: str) -> Dict:
    """Phase-aware novelty of a (molecule, indication) pair."""
    best = -1.0
    matched: Optional[str] = None
    for k in known_indications or []:
        name = (k.get("name") if isinstance(k, dict) else str(k)) or ""
        try:
            ph = float(k.get("max_phase") or 0) if isinstance(k, dict) else 0.0
        except (TypeError, ValueError):
            ph = 0.0
        if name and _matches(disease_name, name) and ph > best:
            best, matched = ph, name

    dz = disease_name or "this indication"
    if matched is None:
        return {
            "state": "novel", "best_phase": None, "matched": None,
            "is_known_indication": False, "blocks_b2": False,
            "headline": f"Novel indication — no prior development of this molecule in {dz}.",
            "note": (f"No recorded development of this molecule for \"{dz}\" — a genuine "
                     "repurposing hypothesis (novelty preserved)."),
            "patentability": ("Composition-of-matter is typically expired for a repurposing "
                              "candidate, but a new method-of-use claim on this indication is open."),
        }
    if best >= 4:
        return {
            "state": "approved_here", "best_phase": best, "matched": matched,
            "is_known_indication": True, "blocks_b2": True,
            "headline": "Already approved for this indication — not a novel repurposing.",
            "note": (f"\"{matched}\" is an APPROVED indication for this molecule (Phase 4). A "
                     "same-use 505(b)(2) or method-of-use patent FAILS the novelty requirement; "
                     "differentiate via a new formulation/route/combination, or choose another indication."),
            "patentability": "Same-use claim blocked; only a new formulation/route/combination may be patentable.",
        }
    if best >= 1:
        return {
            "state": "investigated_here", "best_phase": best, "matched": matched,
            "is_known_indication": False, "blocks_b2": False,
            "headline": (f"New indication with prior Phase-{int(best)} investigation — novel for "
                         "approval, prior art to address."),
            "note": (f"This molecule reached Phase {int(best)} for \"{dz}\" but is NOT approved for it. "
                     "First-approval via 505(b)(2)/hybrid remains viable; secure a method-of-use or "
                     "formulation claim and address the prior clinical art in the patent strategy."),
            "patentability": ("Method-of-use claim viable but must distinguish prior disclosures; a new "
                              "formulation/dose/route strengthens the position."),
        }
    return {
        "state": "prior_signal", "best_phase": best, "matched": matched,
        "is_known_indication": False, "blocks_b2": False,
        "headline": "New indication with an early/exploratory prior signal — novelty preserved.",
        "note": (f"An early record links this molecule to \"{dz}\" (preclinical, or a "
                 "pharmacokinetic/drug-interaction study rather than a therapeutic trial). There is no "
                 "clinical development for the indication, so 505(b)(2)/hybrid and a method-of-use claim "
                 "remain open — confirm the signal is a genuine therapeutic investigation, not a "
                 "drug-interaction co-mention."),
        "patentability": "Method-of-use claim open; verify the prior-art record is not an enabling disclosure.",
    }


def assess(drug_name: str, overall_max_phase, known_indications: List,
           disease_name: str) -> Dict:
    """One reconciled verdict: molecule-approval status (enables data citation) +
    indication novelty (claimability) + the regulatory route per jurisdiction."""
    nv = novelty(known_indications, disease_name)
    try:
        mol_phase = int(float(overall_max_phase or 0))
    except (TypeError, ValueError):
        mol_phase = 0
    approved_mol = mol_phase >= 4

    if approved_mol:
        mol_line = (f"{drug_name} is an approved molecule, so its established safety/PK package "
                    "can be cited rather than regenerated")
    else:
        mol_line = (f"{drug_name} is investigational (Phase {mol_phase}), so the reference safety "
                    "package is not yet citable")

    pathways: List[Dict] = []
    for region, route, ref in _PATHWAYS:
        if nv["blocks_b2"]:
            eligible, note = False, ("Same-use claim blocked — differentiate the product or pick "
                                     "another indication.")
        elif approved_mol:
            eligible, note = True, "Cite the approved molecule's data; run a bridging study for the new indication."
        else:
            eligible, note = False, "Opens once the reference molecule is approved; full data package required for now."
        pathways.append({"region": region, "pathway": route, "reference": ref,
                         "eligible": eligible, "note": note})

    return {
        "headline": f"{mol_line}. {nv['headline']}",
        "molecule_approved": approved_mol,
        "molecule_phase": mol_phase,
        "novelty": nv,
        "pathways": pathways,
        "patentability": nv["patentability"],
    }


if __name__ == "__main__":
    import json
    # itraconazole-style: approved molecule, AD only as a low-phase/PK co-mention
    ki = [{"name": "Tinea Pedis", "max_phase": 4}, {"name": "Aspergillosis", "max_phase": 4},
          {"name": "Alzheimer Disease", "max_phase": 1}]
    print("ITRACONAZOLE → Alzheimer (prior Phase-1 record):")
    print(json.dumps(assess("Itraconazole", 4, ki, "Alzheimer Disease"), indent=2))
    print("\nSame molecule → a truly novel indication:")
    print(json.dumps(novelty(ki, "Diabetic Nephropathy"), indent=2))
    print("\nApproved-here case (blocks):")
    print(json.dumps(novelty(ki, "Aspergillosis"), indent=2)["state"] if False else
          novelty(ki, "Aspergillosis")["state"])
