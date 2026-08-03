"""
Class-level mechanism -> contraindication engine.
═══════════════════════════════════════════════════════════════════════════════════════
Some drug MECHANISMS carry a disease-specific NEGATIVE safety signal that is a property of
the target CLASS, not of any single molecule — so a target-overlap / relatedness ranker,
which sees only "this gene is associated with this disease", will happily surface a
CONTRAINDICATED pair as a lead. The canonical example this module was built for:

    PDE3 inhibition x heart failure
    ----------------------------------
    Inhibiting PDE3 raises cardiac cAMP -> positive inotropy/chronotropy, which acutely
    looks helpful in heart failure but chronically INCREASES MORTALITY:
      - Milrinone (PDE3 inhibitor), PROMISE trial, NEJM 1991: +28% all-cause and +34%
        cardiovascular mortality vs placebo in severe chronic heart failure.
      - Cilostazol (PDE3 inhibitor) carries an FDA BOXED WARNING: contraindicated in heart
        failure of ANY severity, because "several drugs with this pharmacologic effect have
        caused decreased survival compared to placebo in patients with class III-IV heart
        failure."
      - Anagrelide (PDE3 inhibitor) label warns "other drugs that inhibit PDE3 have caused
        decreased survival" in heart failure.

The rule is deliberately DISEASE-SPECIFIC: cilostazol is FDA-approved for intermittent
CLAUDICATION (a PDE3 inhibitor indication), so the contraindication must be "heart failure",
NEVER "any cardiovascular disease". Over-firing on claudication / peripheral vascular
disease would be a false positive.

This is a CURATED, evidence-linked, deterministic table (NOT LLM output). It matches on
TARGET-FAMILY membership + action DIRECTION, so ANY drug that inhibits a PDE3 isoform
(PDE3A or PDE3B) triggers — not a hardcoded drug name. Designed for extension: add another
entry to CLASS_RISKS to cover a new class risk.
"""
from __future__ import annotations

import re
from typing import Dict, List, Optional

# Words in a curated/inferred mechanism action that mark INHIBITION of the target.
_INHIBITORY = ("inhib", "antagon", "block", "negative", "suppress")
# Words that mark the OPPOSITE (activation) — a positive signal that the class risk,
# which is driven by INHIBITION, does not apply. Used to avoid over-firing.
_ACTIVATING = ("agonist", "activator", "positive modulator", "positive allosteric")


def _norm(s: str) -> str:
    return re.sub(r"[^a-z0-9 ]", " ", (s or "").lower()).strip()


def _norm_id(s: str) -> str:
    """Normalize an ontology id to the underscore form used by Open Targets candidates
    (MONDO:0005252 / MONDO_0005252 -> MONDO_0005252)."""
    return (s or "").strip().upper().replace(":", "_")


# ── The curated class-risk table ─────────────────────────────────────────────
# Each entry:
#   mechanism        human-readable mechanism string (e.g. "PDE3 inhibition")
#   target_family    the target-CLASS gene symbols (uppercase) any of which, when the
#                    drug's ACTION opposes them per `action`, means the drug is a member
#   action           the drug action that triggers the risk ("inhibitor" | "agonist")
#   diseases         list of {name substrings, mondo ids} that ARE contraindicated
#   severity         worst documented consequence ("mortality" here)
#   citations        the evidence anchoring the rule
#   flag             a one-line directional flag for the "why not" (Wrong direction: ...)
#   tier_label       the Contradicted badge label
CLASS_RISKS: List[Dict] = [
    {
        "mechanism": "PDE3 inhibition",
        "target_family": {"PDE3A", "PDE3B"},
        "action": "inhibitor",
        # Matched disease-specifically: name-substring OR exact MONDO id. Substrings are
        # specific enough ("heart failure" / "cardiac failure") that claudication, peripheral
        # vascular disease and pulmonary hypertension do NOT match.
        "disease_name_substrings": (
            "heart failure", "cardiac failure", "ventricular failure",
        ),
        "disease_mondo_ids": {
            "MONDO_0005252",   # heart failure
            "MONDO_0005009",   # congestive heart failure
        },
        "severity": "mortality",
        "citations": [
            "Milrinone, PROMISE trial (Packer et al., NEJM 1991): +28% all-cause / "
            "+34% cardiovascular mortality in severe chronic heart failure.",
            "Cilostazol FDA label BOXED WARNING: contraindicated in heart failure of any "
            "severity (PDE3 inhibitors have caused decreased survival in class III-IV HF).",
            "Anagrelide FDA label: 'other drugs that inhibit PDE3 have caused decreased "
            "survival' in heart failure.",
        ],
        "flag": ("PDE3 inhibition increases mortality in heart failure "
                 "(PROMISE, NEJM 1991; cilostazol FDA boxed warning; anagrelide label)"),
        "tier_label": "Contradicted: class-level mortality signal",
    },
]


def _action_is(action_word: str, drug_action: str, drug_action_map: Optional[Dict],
               matched_targets: List[str]) -> bool:
    """Does the drug's action match the risk's triggering action for the matched targets?

    For an inhibition-driven risk we require positive evidence of inhibition (curated or
    inferred). Because these risks are safety-critical and a phosphodiesterase 'activator'
    is not a therapeutic drug class, an UNKNOWN action defaults to firing UNLESS the drug is
    positively known to ACTIVATE the target (which would make the risk inapplicable)."""
    want_inhibitor = action_word == "inhibitor"
    # Gather any action words we know for the matched targets + the drug-level action.
    blobs: List[str] = []
    if drug_action:
        blobs.append(drug_action.lower())
    if drug_action_map:
        for g in matched_targets:
            v = drug_action_map.get(g) or drug_action_map.get(g.upper())
            if v:
                blobs.append(str(v).lower())
    blob = " ".join(blobs)
    inhibitory = any(w in blob for w in _INHIBITORY)
    activating = any(w in blob for w in _ACTIVATING)
    if want_inhibitor:
        if inhibitory:
            return True
        if activating:
            return False            # positively an activator -> risk does not apply
        return not blob             # unknown action -> conservative fire (safety gate)
    # agonist-driven risk (symmetric; none defined yet)
    if activating:
        return True
    if inhibitory:
        return False
    return not blob


def _disease_matches(rule: Dict, disease_name: str, disease_efo: str) -> bool:
    if disease_efo and _norm_id(disease_efo) in rule.get("disease_mondo_ids", set()):
        return True
    dn = _norm(disease_name)
    if not dn:
        return False
    return any(sub in dn for sub in rule.get("disease_name_substrings", ()))


def class_contraindication(drug_genes: Optional[List[str]], disease_name: str,
                           disease_efo: str = "", drug_action: str = "",
                           drug_action_map: Optional[Dict] = None) -> Optional[Dict]:
    """Return a contraindication verdict if the (drug mechanism, disease) pair hits a curated
    class risk, else None. Matches on target-family membership + action direction, so ANY
    inhibitor of the target class triggers — not a per-drug name.

    Verdict fields: mechanism, matched_targets, disease, mondo, severity, citations, flag,
    tier_label, note (flag + citations, for the Contradicted evidence-tier note)."""
    genes = {g.upper() for g in (drug_genes or []) if g}
    if not genes or not disease_name:
        return None
    for rule in CLASS_RISKS:
        matched = sorted(genes & rule["target_family"])
        if not matched:
            continue
        if not _disease_matches(rule, disease_name, disease_efo):
            continue
        if not _action_is(rule["action"], drug_action, drug_action_map, matched):
            continue
        note = rule["flag"] + ". " + " ".join(rule["citations"])
        return {
            "mechanism": rule["mechanism"],
            "matched_targets": matched,
            "disease": disease_name,
            "mondo": _norm_id(disease_efo),
            "severity": rule["severity"],
            "citations": list(rule["citations"]),
            "flag": rule["flag"],
            "tier_label": rule["tier_label"],
            "note": note,
        }
    return None
