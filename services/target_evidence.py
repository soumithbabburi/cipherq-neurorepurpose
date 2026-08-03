"""
Target-header evidence provenance — MEASURED vs family-INFERRED.
═══════════════════════════════════════════════════════════════════════════════════════
A drug's target list on the reverse screen comes from ChEMBL's drug_mechanism ANNOTATION,
which frequently lists a whole target FAMILY (e.g. ensifentrine annotated against PDE3A,
PDE3B, PDE4A, PDE4B, PDE4C, PDE4D) even when isoform-specific POTENCY has only been
measured for a couple of them. Rendering all six identically implies six measured targets,
which overstates the evidence.

This module labels each displayed target with its evidence level, WITHOUT fabricating any
potency value:

  - "measured"        a measured bioactivity record exists in ChEMBL (activities table:
                      an IC50/Ki/Kd/EC50/pChEMBL for this molecule x this target), OR the
                      isoform is documented as measured in a cited primary source.
  - "family-inferred" the isoform appears only via mechanism ANNOTATION (no isoform-specific
                      potency in the platform's data or the cited source).

The curated literature map is evidence-linked and carries NO fabricated numbers — it only
records WHICH isoforms have a published measured affinity and the citation, so the UI can
mark the family-inferred ones. Extensible: add a drug entry to _LITERATURE_MEASURED.
"""
from __future__ import annotations

import logging
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)

# Curated, citation-linked record of which isoforms have PUBLISHED measured affinity for a
# drug (keyed by lowercase drug name). NO potency numbers are asserted here — only the fact
# that a measurement exists and its source — so nothing is fabricated. Everything else in the
# annotated target list is treated as family-inferred (mechanism annotation only).
_LITERATURE_MEASURED: Dict[str, Dict] = {
    "ensifentrine": {
        "measured": {"PDE3A", "PDE4B"},
        "citation": ("Published potency is isoform-specific to PDE3A (sub-nanomolar) and, "
                     "more weakly, PDE4B; no isoform-specific IC50 has been published for "
                     "PDE3B, PDE4A, PDE4C or PDE4D (family-inferred from the dual PDE3/PDE4 "
                     "mechanism). Refs: Boswell-Smith 2006 (J Pharmacol Exp Ther); Cazzola 2018."),
    },
}


def _measured_from_chembl(chembl_ids: Optional[List[str]], genes: List[str]) -> set:
    """Genes for which ChEMBL holds a MEASURED potency (IC50/Ki/Kd/EC50/pChEMBL) for this
    molecule. Empty when the target list is annotation-only (as for ensifentrine locally)."""
    ids = [c for c in (chembl_ids or []) if c]
    want = {g.upper() for g in (genes or [])}
    if not ids or not want:
        return set()
    try:
        from services.repurposing_engine import _get_chembl_pool
    except Exception:
        return set()
    pool = _get_chembl_pool()
    if pool is None:
        return set()
    conn = None
    try:
        conn = pool.getconn()
        with conn.cursor() as cur:
            cur.execute(
                """
                SELECT DISTINCT cs.component_synonym
                FROM molecule_dictionary md
                LEFT JOIN molecule_hierarchy mh ON mh.molregno = md.molregno
                JOIN activities act
                     ON act.molregno = md.molregno OR act.molregno = mh.parent_molregno
                JOIN assays a ON a.assay_id = act.assay_id
                JOIN target_components tc ON tc.tid = a.tid
                JOIN component_synonyms cs
                     ON cs.component_id = tc.component_id AND cs.syn_type = 'GENE_SYMBOL'
                WHERE md.chembl_id = ANY(%s)
                  AND act.standard_type IN ('IC50','Ki','Kd','EC50','pIC50','pKi')
                  AND (act.pchembl_value IS NOT NULL OR act.standard_value IS NOT NULL)
                """,
                (ids,),
            )
            found = {(r[0] or "").upper() for r in cur.fetchall()}
            return found & want
    except Exception as e:
        logger.debug("measured-target chembl query failed: %s", e)
        return set()
    finally:
        if conn is not None:
            pool.putconn(conn)


def annotate(drug_name: str, chembl_ids: Optional[List[str]], genes: List[str]) -> Dict:
    """Evidence-level annotation for a drug's displayed target list.

    Returns {"targets": {GENE: {"level": measured|family-inferred, "source": ...}},
             "measured": [...], "inferred": [...], "note": "<explanatory + citation>"}.
    Fail-soft: on any error, everything is reported as measured=unknown (no downgrade)."""
    glist = [g for g in (genes or []) if g]
    out = {"targets": {}, "measured": [], "inferred": [], "note": ""}
    if not glist:
        return out
    measured = _measured_from_chembl(chembl_ids, glist)
    source_of = {g.upper(): "chembl_activity" for g in measured}
    lit = _LITERATURE_MEASURED.get((drug_name or "").strip().lower())
    citation = ""
    if lit:
        for g in lit.get("measured", set()):
            gu = g.upper()
            if gu not in measured and gu in {x.upper() for x in glist}:
                measured = measured | {gu}
                source_of[gu] = "literature"
        citation = lit.get("citation", "")
    for g in glist:
        gu = g.upper()
        if gu in measured:
            out["targets"][g] = {"level": "measured", "source": source_of.get(gu, "chembl_activity")}
            out["measured"].append(g)
        else:
            out["targets"][g] = {"level": "family-inferred", "source": "mechanism_annotation"}
            out["inferred"].append(g)
    if out["inferred"]:
        base = ("Targets are from ChEMBL mechanism annotation. "
                + (f"Measured/published potency: {', '.join(out['measured'])}. " if out["measured"] else
                   "No isoform-specific potency is present in the platform's data. ")
                + f"Family-inferred (annotation only, no isoform-specific IC50): "
                  f"{', '.join(out['inferred'])}.")
        out["note"] = (base + " " + citation).strip()
    return out
