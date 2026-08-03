"""
Disease-appropriateness  —  could the drug TREAT this indication, or would it worsen /
does it already CAUSE it?
═══════════════════════════════════════════════════════════════════════════════════════
An Open Targets gene↔disease association means the gene is INVOLVED, not that inhibiting
it TREATS the disease — often the reverse. Three complementary checks catch the backwards
candidates that a direction-blind target-overlap score lets through:

 A. DIRECTION (handled by services.mechanism_direction, which we now feed an INFERRED
    action when ChEMBL has no curated mechanism — see infer_drug_action). An inhibitor of
    a gene a disease has turned DOWN is harmful.
 B. LOSS-OF-FUNCTION / DEVELOPMENTAL: an INHIBITOR proposed for a congenital/developmental
    genetic disease is usually backwards — those are driven by loss of the gene and need
    restoration, not further inhibition (RET → Hirschsprung: RET loss causes it, so an RET
    inhibitor would worsen it, even though the RET↔Hirschsprung association is very strong).
 C. ADVERSE EVENT: the disease is a known SERIOUS adverse reaction of the drug (FAERS) — a
    toxicity the drug CAUSES, not a target (angiokinase inhibitors → hypothyroidism).

All checks are fail-soft and return a bounded multiplier + a human-readable flag; they
GATE the ranking, they do not silently rewrite the mechanistic score.
"""
from __future__ import annotations

import logging
import re
from typing import Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

# Suffixes / words that reliably mark a small molecule or antibody as an INHIBITOR of its
# target (kinase -ib, protease/enzyme -stat, blocking mAb -mab). Deliberately high-precision
# so we only infer an action when confident; anything else stays unknown (neutral).
_INHIBITOR_WORDS = ("inhibitor", "antagonist", "blocker", "negative")
_AGONIST_WORDS = ("agonist", "activator", "positive modulator")


def infer_drug_action(name: str, drug_genes: Optional[List[str]] = None) -> str:
    """Best-guess pharmacological action ('inhibitor' | 'agonist' | '') from the drug name,
    for when ChEMBL drug_mechanism is empty. Used only to power the direction / appropriateness
    checks — never to assert mechanism. High-precision: unknown → '' (stays neutral)."""
    n = (name or "").strip().lower()
    if not n:
        return ""
    if any(w in n for w in _INHIBITOR_WORDS):
        return "inhibitor"
    if any(w in n for w in _AGONIST_WORDS):
        return "agonist"
    # INN stems: -ib (kinase/enzyme inhibitors: -nib/-tinib/-ciclib/-parib/-lisib/-zomib…),
    # -stat (enzyme inhibitors), -mab (therapeutic antibodies are overwhelmingly blocking).
    if n.endswith("ib") or n.endswith("stat") or n.endswith("mab"):
        return "inhibitor"
    return ""


# Open Targets therapeutic-area label for congenital/developmental genetic disorders —
# the ones most likely to be loss-of-function (an inhibitor is backwards for them).
_DEVELOPMENTAL_AREA = "genetic, familial or congenital"


# FAERS terms that are OUTCOMES or the drug's own INDICATION, not a drug-caused disease —
# they must never trigger the adverse-event flag (confounding-by-indication guard).
_NON_AE_TERMS = frozenset({
    "death", "disease progression", "malignant neoplasm progression", "drug ineffective",
    "off label use", "condition aggravated", "general physical health deterioration",
    "drug interaction", "product use issue", "therapeutic response decreased",
    "no adverse event", "drug intolerance", "treatment failure",
})


def _norm(s: str) -> str:
    return re.sub(r"[^a-z0-9 ]", " ", (s or "").lower()).strip()


def _matches_adverse_event(disease: str, reactions: List[Tuple[str, int]],
                           total: int) -> Optional[Tuple[str, int, float]]:
    """Is `disease` essentially one of the drug's serious FAERS reactions, with enough
    reporting mass to be a real signal? Returns (term, count, share) or None. Match is
    whole-phrase containment either way (so 'hypothyroidism' ↔ 'Hypothyroidism'), guarded
    by a count/share floor so an incidental low-count term never penalises a candidate."""
    dn = _norm(disease)
    if len(dn) < 4:
        return None
    for term, cnt in reactions:
        tn = _norm(term)
        if not tn or tn in _NON_AE_TERMS:
            continue
        share = cnt / max(1, total)
        if (dn == tn or f" {tn} " in f" {dn} " or f" {dn} " in f" {tn} "):
            if cnt >= 15 or share >= 0.005:
                return (term, cnt, round(share, 4))
    return None


def appropriateness(drug_name: str, disease_name: str, therapeutic_areas: List[str],
                    action: str, faers_reactions: Optional[List[Tuple[str, int]]] = None,
                    faers_total: int = 0, drug_genes: Optional[List[str]] = None,
                    drug_action_map: Optional[Dict] = None, disease_efo: str = "") -> Dict:
    """Bounded appropriateness verdict for a (drug, candidate-disease) pair.
    Returns {factor 0.3–1.0, appropriate, flags, reasons}. factor GATES ranking /
    actionability; it does not overwrite the mechanistic composite."""
    out = {"factor": 1.0, "appropriate": True, "flags": [], "reasons": []}
    areas_blob = " ".join(therapeutic_areas or []).lower()

    # D — CLASS-LEVEL mechanism contraindication (curated, evidence-linked). A drug whose
    # target CLASS + action direction carries a disease-SPECIFIC mortality/serious signal
    # (e.g. PDE3 inhibition x heart failure — PROMISE / cilostazol boxed warning) is
    # CONTRAINDICATED, not merely weak. This is the strongest negative: it drives the
    # candidate to the "Contradicted" tier and out of the actionable bucket (see the
    # reverse screen wiring). Disease-specific by design, so cilostazol's approved
    # claudication use is NOT flagged. Fail-soft.
    try:
        from services.class_safety import class_contraindication
        cs = class_contraindication(drug_genes, disease_name, disease_efo=disease_efo,
                                    drug_action=action, drug_action_map=drug_action_map)
        if cs:
            out["factor"] = min(out["factor"], 0.2)
            out["flags"].append(cs["flag"])
            out["reasons"].append(cs["note"])
            out["class_safety"] = cs
    except Exception as e:                                   # pragma: no cover
        logger.debug("class-safety check skipped: %s", e)

    # B — inhibitor proposed for a congenital/developmental (loss-of-function) disorder.
    if action == "inhibitor" and _DEVELOPMENTAL_AREA in areas_blob:
        out["factor"] = min(out["factor"], 0.4)
        out["flags"].append("loss-of-function mismatch")
        out["reasons"].append(
            f"{disease_name} is a congenital/developmental (typically loss-of-function) "
            f"disorder; an inhibitor would suppress the target further. Appropriate only "
            f"if the disease is gain-of-function.")

    # C — the candidate disease is a serious adverse event the drug is reported to CAUSE.
    if faers_reactions:
        ae = _matches_adverse_event(disease_name, faers_reactions, faers_total)
        if ae:
            term, cnt, share = ae
            out["factor"] = min(out["factor"], 0.3)
            out["flags"].append("reported adverse event")
            out["reasons"].append(
                f"{disease_name} is a serious adverse reaction of {drug_name} in FAERS "
                f"({cnt:,} reports, {share*100:.1f}% of its serious mass) — a toxicity the "
                f"drug causes, not a target it could treat.")

    out["appropriate"] = out["factor"] >= 0.75
    return out
