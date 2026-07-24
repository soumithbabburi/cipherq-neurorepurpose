"""
Mechanism direction  —  does the drug push the disease the RIGHT way?
═══════════════════════════════════════════════════════════════════════════════
The genetic-overlap target score is DIRECTION-BLIND: a drug that ACTIVATES a gene
over-expressed in a disease scores identically to one that INHIBITS it — even though
the first would EXACERBATE the disease and the second would treat it. That is a
reasoning bug, not a missing feature.

This module makes target scoring direction-aware, and doubles as a lightweight
CONTRAINDICATION signal. For each gene shared by the drug's targets and the disease:
  • disease direction  d(g) = +1 if UP in disease (HetioNet DuG), −1 if DOWN (DdG)
  • drug effect        e(g) = −1 if the drug inhibits g, +1 if it activates g (ChEMBL)
  • CORRECTIVE when e·d < 0 (drug pushes the gene back toward normal — inhibit an
    up-gene, or activate a down-gene) → keep/boost credit.
  • HARMFUL when e·d > 0 (drug amplifies the disease direction) → penalise + flag
    "may exacerbate".
Fail-soft: no direction / no action data → neutral factor 1.0 (never invent a penalty).
"""
from __future__ import annotations

import logging
from functools import lru_cache
from typing import Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

_INHIBIT = ("inhibitor", "antagonist", "blocker", "negative modulator",
            "inverse agonist", "degrader", "disruptor", "suppressor")
_ACTIVATE = ("agonist", "activator", "positive modulator", "opener",
             "positive allosteric", "partial agonist", "inducer")


def _drug_effect(action: str) -> int:
    a = (action or "").lower()
    if any(t in a for t in _INHIBIT):
        return -1
    if any(t in a for t in _ACTIVATE):
        return +1
    return 0


@lru_cache(maxsize=512)
def _disease_direction(disease_key: str) -> Tuple[frozenset, frozenset]:
    """(up_genes, down_genes) for a disease from HetioNet DuG/DdG. Cached per disease."""
    from services.repurposing_scorer import _resolve_hetionet_diseases, _q
    dids = _resolve_hetionet_diseases([disease_key])
    if not dids:
        return frozenset(), frozenset()
    up = _q("SELECT DISTINCT hn.name AS g FROM hetionet_edges e "
            "JOIN hetionet_nodes hn ON hn.id=e.target_id "
            "WHERE e.metaedge='DuG' AND e.source='hetionet_v1.0' AND e.source_id = ANY(%s)",
            (dids,))
    dn = _q("SELECT DISTINCT hn.name AS g FROM hetionet_edges e "
            "JOIN hetionet_nodes hn ON hn.id=e.target_id "
            "WHERE e.metaedge='DdG' AND e.source='hetionet_v1.0' AND e.source_id = ANY(%s)",
            (dids,))
    return (frozenset(r["g"].upper() for r in up if r.get("g")),
            frozenset(r["g"].upper() for r in dn if r.get("g")))


def mechanism_direction(drug_genes: List[str], drug_actions: Optional[Dict[str, str]],
                        disease: str) -> Dict:
    """{factor, net, corrective, harmful, flag} — a bounded multiplier on the target
    score and a contraindication flag. factor ∈ [0.35, 1.15]; 1.0 when undetermined."""
    out = {"factor": 1.0, "net": "neutral", "corrective": [], "harmful": [],
           "covered": False, "assessed": False, "flag": ""}
    actions = {k.upper(): _drug_effect(v) for k, v in (drug_actions or {}).items()}
    if not actions or not drug_genes:
        return out
    up, down = _disease_direction(disease)
    if not up and not down:
        # We know HOW the drug acts (agonist/antagonist) but have no directional
        # (up/down-regulation) data for this disease, so we cannot confirm the drug pushes
        # it the right way. Say so honestly rather than defaulting to a silent neutral pass,
        # which would let an activator of a disease-driving target read as a clean lead.
        out["flag"] = ("Direction not assessed: no directional expression data for this "
                       "disease, so the platform cannot verify the drug corrects rather than "
                       "worsens it. Treat the mechanistic direction as unconfirmed.")
        return out
    out["covered"] = True
    out["assessed"] = True
    corrective, harmful = [], []
    for g in {x.upper() for x in drug_genes}:
        e = actions.get(g, 0)
        if e == 0:
            continue
        d = +1 if g in up else -1 if g in down else 0
        if d == 0:
            continue
        (harmful if e * d > 0 else corrective).append(g)
    out["corrective"], out["harmful"] = sorted(corrective), sorted(harmful)
    nc, nh = len(corrective), len(harmful)
    if nc == 0 and nh == 0:
        return out
    # net-corrective boosts the mechanistic score modestly; net-harmful crushes it and
    # raises a contraindication flag (the drug would push the disease the wrong way).
    if nh > nc:
        out["net"] = "harmful"
        out["factor"] = round(max(0.35, 1.0 - 0.35 * (nh - nc)), 3)
        out["flag"] = (f"May EXACERBATE: the drug amplifies disease-direction "
                       f"gene(s) {', '.join(out['harmful'][:4])} (activates an up-regulated "
                       "or inhibits a down-regulated disease gene).")
    elif nc > nh:
        out["net"] = "corrective"
        out["factor"] = round(min(1.15, 1.0 + 0.08 * (nc - nh)), 3)
    return out
