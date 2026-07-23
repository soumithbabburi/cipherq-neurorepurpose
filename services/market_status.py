"""
Drug market / regulatory status (forward-path development-reality filter).
════════════════════════════════════════════════════════════════════════════════
A computational target match does not know that a drug was pulled from the market for
safety, discontinued when its propellant was banned, or rendered obsolete by a safer
in-class successor. Surfacing such a drug as a repurposing "lead" is misleading.

This module answers "is this molecule still a viable development asset?" for the
FORWARD discovery path (disease → drugs). The local compound records carry no
availability/withdrawal field (their ids are internal, not real ChEMBL ids, so the
live ChEMBL withdrawn_flag cannot be keyed reliably), so we maintain a curated,
CITED registry of the well-documented cases. Each entry is a regulatory fact with a
reason and a source; the registry is intentionally conservative (only clear cases)
and easy to extend.

Status meanings:
  withdrawn    — removed from market for safety/efficacy (a hard disqualifier).
  discontinued — no longer marketed for commercial/technical reasons, no safety flag
                 (e.g. a propellant ban) — a soft disqualifier: a reformulation could
                 in principle revive it, but it is not an off-the-shelf asset.
  obsolete     — superseded by a safer/more selective in-class drug; still available
                 in places but not a rational development target (soft disqualifier).
"""
from __future__ import annotations

from typing import Dict, Optional

# Curated regulatory-status registry. Keys are normalised (lower-case) drug names and
# common synonyms. Every entry cites the reason so the flag is auditable, never a bare
# assertion. Extend as new cases are verified.
_REGISTRY: Dict[str, Dict] = {
    # β2-agonists — the asthma false-positive class that motivated this filter.
    "orciprenaline": {"status": "obsolete", "reason": (
        "Non-selective beta agonist largely phased out in favour of selective beta 2 "
        "agonists (e.g. albuterol) because its lack of receptor selectivity carries a "
        "higher rate of cardiac effects (tachycardia, arrhythmia)."),
        "source": "Clinical pharmacology consensus"},
    "metaproterenol": {"status": "obsolete", "reason": (
        "Metaproterenol is orciprenaline; a non-selective beta agonist superseded by "
        "selective beta 2 agonists on cardiac-safety grounds."),
        "source": "Clinical pharmacology consensus"},
    "pirbuterol": {"status": "discontinued", "reason": (
        "Maxair (pirbuterol) was discontinued after the FDA phase-out of "
        "chlorofluorocarbon inhaler propellants; it was not reformulated as an HFA "
        "product because selective beta 2 agonists already dominated the market."),
        "source": "FDA CFC inhaler phase-out"},
    "isoprenaline": {"status": "obsolete", "reason": (
        "Non-selective beta agonist superseded by selective agents for airway disease."),
        "source": "Clinical pharmacology consensus"},
    "isoproterenol": {"status": "obsolete", "reason": (
        "Isoproterenol is isoprenaline; non-selective beta agonist, superseded for "
        "airway indications."),
        "source": "Clinical pharmacology consensus"},

    # Classic safety withdrawals — kept as hard disqualifiers wherever they surface.
    "cerivastatin": {"status": "withdrawn", "reason": (
        "Withdrawn worldwide in 2001 for fatal rhabdomyolysis risk."),
        "source": "Safety withdrawal, 2001"},
    "rofecoxib": {"status": "withdrawn", "reason": (
        "Withdrawn in 2004 for increased cardiovascular thrombotic risk."),
        "source": "Safety withdrawal, 2004"},
    "troglitazone": {"status": "withdrawn", "reason": (
        "Withdrawn for hepatotoxicity."),
        "source": "Safety withdrawal, 2000"},
    "sibutramine": {"status": "withdrawn", "reason": (
        "Withdrawn for increased cardiovascular event risk."),
        "source": "Safety withdrawal, 2010"},
    "cisapride": {"status": "withdrawn", "reason": (
        "Withdrawn/restricted for QT prolongation and fatal arrhythmia."),
        "source": "Safety withdrawal, 2000"},
    "terfenadine": {"status": "withdrawn", "reason": (
        "Withdrawn for QT prolongation; superseded by fexofenadine."),
        "source": "Safety withdrawal, 1998"},
    "astemizole": {"status": "withdrawn", "reason": (
        "Withdrawn for cardiac arrhythmia (QT prolongation)."),
        "source": "Safety withdrawal, 1999"},
    "phenylpropanolamine": {"status": "withdrawn", "reason": (
        "Withdrawn for haemorrhagic stroke risk."),
        "source": "Safety withdrawal, 2000"},
    "thioridazine": {"status": "obsolete", "reason": (
        "Restricted to refractory use for QT prolongation / arrhythmia risk; "
        "superseded by safer antipsychotics."),
        "source": "Safety restriction"},
}

# Statuses that should REMOVE a candidate from the lead list entirely vs merely demote.
_HARD = {"withdrawn"}
_SOFT = {"discontinued", "obsolete"}

# How hard to demote a soft-status candidate's rank (multiplicative).
_SOFT_MULT = 0.35


def _norm(name: str) -> str:
    return (name or "").strip().lower()


def status(drug_name: str) -> Optional[Dict]:
    """Return {status, reason, source, disqualifier, multiplier} for a drug, or None if
    it is not in the registry (treated as an active, viable asset)."""
    rec = _REGISTRY.get(_norm(drug_name))
    if not rec:
        return None
    st = rec["status"]
    return {
        "status": st,
        "reason": rec["reason"],
        "source": rec.get("source", ""),
        "disqualifier": "hard" if st in _HARD else "soft",
        "multiplier": 0.0 if st in _HARD else _SOFT_MULT,
        "label": {"withdrawn": "Withdrawn for safety",
                  "discontinued": "Discontinued",
                  "obsolete": "Obsolete / superseded"}.get(st, st.title()),
    }
