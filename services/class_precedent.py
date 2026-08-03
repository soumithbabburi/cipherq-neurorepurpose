"""
Class-level mechanism -> APPROVED-PRECEDENT engine (the positive mirror of class_safety).
═══════════════════════════════════════════════════════════════════════════════════════
When a congener that shares the drug's primary target CLASS is ALREADY FDA-approved for the
SAME indication, that is the strongest possible validation of a mechanistic repurposing
hypothesis — far stronger than any single Open Targets association score. A pure target-
overlap ranker will present such a pair as merely "novel", hiding the fact that the class is
already clinically proven for that disease. This module surfaces that signal.

Canonical examples (the reason this module exists):
    PDE3 inhibition x intermittent claudication  -> cilostazol, FDA-approved 1999
    PDE3 inhibition x essential thrombocythemia   -> anagrelide, FDA-approved 1997
    PDE4 inhibition x atopic dermatitis/eczema     -> crisaborole, FDA-approved 2016 (topical)

Curated, evidence-linked, DETERMINISTIC (not LLM output). Matches on TARGET-FAMILY
membership + inhibitor direction + disease, exactly like class_safety, so ANY PDE3/PDE4
inhibitor triggers — not a hardcoded drug name. This is display-only: it NEVER changes the
composite score or rank; it attaches a "Class-validated" badge the card can render.
Designed for extension: add an entry to CLASS_PRECEDENTS to cover a new class precedent.
"""
from __future__ import annotations

import re
from typing import Dict, List, Optional

_INHIBITORY = ("inhib", "antagon", "block", "negative", "suppress")
_ACTIVATING = ("agonist", "activator", "positive modulator", "positive allosteric")


def _norm(s: str) -> str:
    return re.sub(r"[^a-z0-9 ]", " ", (s or "").lower()).strip()


def _norm_id(s: str) -> str:
    return (s or "").strip().upper().replace(":", "_")


# ── The curated class-precedent table ────────────────────────────────────────
# Each mechanism entry lists the target CLASS and, per indication, the already-approved
# congener that anchors the precedent.
CLASS_PRECEDENTS: List[Dict] = [
    {
        "mechanism": "PDE3 inhibition",
        "target_family": {"PDE3A", "PDE3B"},
        "action": "inhibitor",
        "indications": [
            {
                "name_substrings": ("intermittent claudication",
                                    "intermittent vascular claudication", "claudication"),
                "mondo_ids": {"MONDO_0005295"},
                "congener": "cilostazol",
                "congener_class": "PDE3 inhibitor",
                "year": 1999,
            },
            {
                "name_substrings": ("essential thrombocythemia", "essential thrombocythaemia"),
                "mondo_ids": {"MONDO_0005029"},
                "congener": "anagrelide",
                "congener_class": "PDE3A inhibitor",
                "year": 1997,
            },
        ],
    },
    {
        "mechanism": "PDE4 inhibition",
        "target_family": {"PDE4A", "PDE4B", "PDE4C", "PDE4D"},
        "action": "inhibitor",
        "indications": [
            {
                "name_substrings": ("atopic dermatitis", "atopic eczema", "eczema"),
                "mondo_ids": {"MONDO_0004980"},
                "congener": "crisaborole",
                "congener_class": "PDE4 inhibitor",
                "year": 2016,
                "route_note": "topical",
            },
        ],
    },
]


def _is_inhibitor(drug_action: str, drug_action_map: Optional[Dict],
                  matched_targets: List[str]) -> bool:
    """The precedent applies to INHIBITORS of the class (the approved congeners are all
    inhibitors). Fire on positive inhibition evidence, or on an unknown action (a
    phosphodiesterase 'activator' is not a real drug class); never on a known activator."""
    blobs: List[str] = []
    if drug_action:
        blobs.append(drug_action.lower())
    if drug_action_map:
        for g in matched_targets:
            v = drug_action_map.get(g) or drug_action_map.get(g.upper())
            if v:
                blobs.append(str(v).lower())
    blob = " ".join(blobs)
    if any(w in blob for w in _INHIBITORY):
        return True
    if any(w in blob for w in _ACTIVATING):
        return False
    return not blob   # unknown action -> allow (these are inhibitor classes in practice)


def _indication_matches(ind: Dict, disease_name: str, disease_efo: str) -> bool:
    if disease_efo and _norm_id(disease_efo) in ind.get("mondo_ids", set()):
        return True
    dn = _norm(disease_name)
    if not dn:
        return False
    return any(sub in dn for sub in ind.get("name_substrings", ()))


def class_precedent(drug_genes: Optional[List[str]], disease_name: str,
                    disease_efo: str = "", drug_action: str = "",
                    drug_action_map: Optional[Dict] = None) -> Optional[Dict]:
    """Return a class-validated-precedent verdict if a congener sharing the drug's target
    CLASS is already FDA-approved for this indication, else None.

    Verdict fields: mechanism, matched_targets, congener, congener_class, year, indication,
    label (the badge text), note (the tooltip)."""
    genes = {g.upper() for g in (drug_genes or []) if g}
    if not genes or not disease_name:
        return None
    for rule in CLASS_PRECEDENTS:
        matched = sorted(genes & rule["target_family"])
        if not matched:
            continue
        for ind in rule["indications"]:
            if not _indication_matches(ind, disease_name, disease_efo):
                continue
            if not _is_inhibitor(drug_action, drug_action_map, matched):
                continue
            route = f" ({ind['route_note']})" if ind.get("route_note") else ""
            label = (f"Class-validated: {ind['congener']} ({ind['congener_class']}) is "
                     f"FDA-approved{route} for this indication since {ind['year']}")
            note = (label + ". A drug sharing this target class is already clinically proven "
                    "for this disease — the strongest mechanistic validation. Note a 505(b)(2) "
                    "would face this approved competitor. Source: DrugBank / FDA label.")
            return {
                "mechanism": rule["mechanism"],
                "matched_targets": matched,
                "congener": ind["congener"],
                "congener_class": ind["congener_class"],
                "year": ind["year"],
                "indication": disease_name,
                "label": label,
                "note": note,
            }
    return None
