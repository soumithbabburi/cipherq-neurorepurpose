"""
Quality overlay — one place that attaches the platform's validated quality signals
to a (drug, disease) candidate, so EVERY repurposing surface (repurpose-disease /
discover, repurpose-molecule / reverse, pathway screen, novel targets) applies the
same filters consistently:

  • plausibility   — validated DWPC P(treats), where the pair maps into Hetionet
  • lead_viability — potency (IC50) gate; PBPK window deferred to the dossier
  • disease_value  — Repurposing Value Score for the indication (burden × unmet-need
                     × market-fit) — "would a pharma company pursue this disease?"

All three are fail-soft: a missing signal is simply absent, never fabricated.
"""
from __future__ import annotations

import logging
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)


def overlay(drug_name: str, chembl_id: str, disease: str,
            genes: Optional[List[str]] = None, *, smiles: str = "",
            efo_id: str = "", with_disease_value: bool = True) -> Dict:
    """Return {plausibility, lead_viability, disease_value} for a candidate pair.
    Cheap + cached under the hood; the slow PBPK window is deferred to the dossier."""
    out: Dict = {}
    try:
        from services.repurposing_predictor import plausibility
        out["plausibility"] = plausibility(drug_name, disease)
    except Exception as e:
        logger.debug(f"overlay plausibility skipped: {e}")
    try:
        from services.lead_viability import assess
        pl = (out.get("plausibility") or {}).get("probability")
        out["lead_viability"] = assess(drug_name, chembl_id, disease, list(genes or []),
                                       smiles=smiles, plausibility_p=pl, run_window=False)
    except Exception as e:
        logger.debug(f"overlay lead_viability skipped: {e}")
    if with_disease_value:
        try:
            from services.disease_value import value_for
            out["disease_value"] = value_for(disease, efo_id)
        except Exception as e:
            logger.debug(f"overlay disease_value skipped: {e}")
    return out
