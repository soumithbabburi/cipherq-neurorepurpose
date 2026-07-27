"""
Clinical-evidence miner — structured, provenance-tagged, non-keyword.
═══════════════════════════════════════════════════════════════════════════════
REPORTS trial-registry and literature FACTS with provenance and caveats. It never
asserts an efficacy/safety determination of its own (non-device / CDS lane, see
validation/REGULATORY_POSITIONING.md). Extraction is DETERMINISTIC from official
STRUCTURED fields, not free-text keyword matching:

  * Trial outcomes  -> delegated to endpoint_parser.aggregate (already typed by
    ClinicalTrials.gov v2 unitOfMeasure / paramType, not title keywords).
  * Adverse events  -> resultsSection.adverseEventsModule organ-system tables WITH
    numAtRisk denominators (so rates are real, unlike FAERS, which has none).
  * Literature      -> PubMed publication types + MeSH indexing (study-design
    weighting), not naive [tiab] co-mention.

The pure parsers (parse_adverse_events, literature_tier) are unit-tested with
fixtures; mine_clinical_evidence adds the fetch layer and assembles the record.
Fail-soft: any missing structured data yields 'not reported', never a fabricated
value.
"""
from __future__ import annotations

import logging
from typing import Dict, List, Optional
from xml.etree import ElementTree as ET

logger = logging.getLogger(__name__)

CT_STUDIES = "https://clinicaltrials.gov/api/v2/studies"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# MeSH PublicationType descriptor UIs (also matched by name, case-insensitive)
_PT = {"meta": ("D017418", "meta-analysis"),
       "sysrev": ("D000078182", "systematic review"),
       "rct": ("D016449", "randomized controlled trial"),
       "case_report": ("D002363", "case reports")}


# ── Pure parsers (deterministic, unit-tested) ────────────────────────────────

def parse_adverse_events(study: Dict) -> Optional[Dict]:
    """Structured AE summary from resultsSection.adverseEventsModule, WITH the
    numAtRisk denominators that make trial rates real. None if not posted."""
    aem = ((study.get("resultsSection") or {}).get("adverseEventsModule")) or {}
    groups = aem.get("eventGroups") or []
    if not aem or not groups:
        return None
    def _isum(key):
        return sum(int(g.get(key) or 0) for g in groups)
    at_risk = _isum("seriousNumAtRisk")
    affected = _isum("seriousNumAffected")
    organs: Dict[str, int] = {}
    for ev in (aem.get("seriousEvents") or []):
        soc = ev.get("organSystem") or "Unspecified"
        n = sum(int(s.get("numAffected") or 0) for s in (ev.get("stats") or []))
        organs[soc] = organs.get(soc, 0) + n
    top = sorted(organs.items(), key=lambda x: -x[1])[:5]
    return {
        "serious_num_affected": affected,
        "serious_num_at_risk": at_risk,
        "serious_rate": (round(affected / at_risk, 4) if at_risk else None),
        "deaths": _isum("deathsNumAffected"),
        "top_serious_organ_systems": [{"organ_system": k, "events": v} for k, v in top],
        "no_denominator": at_risk == 0,
        "caveat": ("Trial-reported serious-AE counts; denominator = participants at risk in THIS "
                   "trial, so rates are trial-specific, not population incidence. Do NOT compare or "
                   "sum with FAERS (which has no denominator)."),
    }


def literature_tier(records: List[Dict]) -> Dict:
    """Evidence tier from PubMed STUDY DESIGN + MeSH, not co-mention counts.
    records: [{pub_types:[...], mesh:[{descriptor, qualifiers:[...]}], mesh_confirmed:bool}]."""
    counts = {"rct": 0, "meta": 0, "sysrev": 0, "case_report": 0, "other": 0}
    direction = {"treats": 0, "causes": 0}
    mesh_confirmed = False
    for r in records or []:
        pts = {str(p).lower() for p in (r.get("pub_types") or [])}
        if any(k in pts for k in _PT["meta"]):
            counts["meta"] += 1
        elif any(k in pts for k in _PT["sysrev"]):
            counts["sysrev"] += 1
        elif any(k in pts for k in _PT["rct"]):
            counts["rct"] += 1
        elif any(k in pts for k in _PT["case_report"]):
            counts["case_report"] += 1
        else:
            counts["other"] += 1
        for m in (r.get("mesh") or []):
            quals = " ".join(m.get("qualifiers", [])).lower()
            if "drug therapy" in quals or "therapeutic use" in quals:
                direction["treats"] += 1
            if "chemically induced" in quals or "adverse effects" in quals:
                direction["causes"] += 1
        if r.get("mesh_confirmed"):
            mesh_confirmed = True
    total = sum(counts.values())
    if counts["meta"] >= 1 or counts["sysrev"] >= 1 or counts["rct"] >= 2:
        tier = "high"
    elif counts["rct"] >= 1 or (mesh_confirmed and counts["other"] >= 3):
        tier = "moderate"
    elif mesh_confirmed or total > 0:
        tier = "low"
    else:
        tier = "none"
    d = ("treats" if direction["treats"] > direction["causes"]
         else "causes" if direction["causes"] > direction["treats"] else "unknown")
    return {"tier": tier, "counts": counts, "mesh_confirmed": mesh_confirmed,
            "directional": d, "n": total}


def _parse_efetch_xml(xml_text: str) -> List[Dict]:
    """PubMed efetch XML -> [{pmid, pub_types, mesh, mesh_confirmed}]. Deterministic."""
    out: List[Dict] = []
    try:
        root = ET.fromstring(xml_text)
    except Exception:
        return out
    for art in root.findall(".//PubmedArticle"):
        pts = [pt.text for pt in art.findall(".//PublicationTypeList/PublicationType") if pt.text]
        mesh = []
        for mh in art.findall(".//MeshHeadingList/MeshHeading"):
            dn = mh.find("DescriptorName")
            quals = [q.text for q in mh.findall("QualifierName") if q.text]
            if dn is not None and dn.text:
                mesh.append({"descriptor": dn.text, "qualifiers": quals})
        pmid_el = art.find(".//PMID")
        out.append({"pmid": pmid_el.text if pmid_el is not None else None,
                    "pub_types": pts, "mesh": mesh, "mesh_confirmed": bool(mesh)})
    return out


# ── Fetch + assemble (network; fail-soft) ────────────────────────────────────

def mine_clinical_evidence(drug_name: str, disease_name: str,
                           max_trials: int = 100, include_literature: bool = True) -> Dict:
    """Provenance-tagged clinical-evidence record for a (drug, disease) pair.
    Reports facts with caveats; never asserts efficacy/safety. Fail-soft."""
    from services import http_client
    try:
        from services import provenance as prov
    except Exception:
        prov = None

    rec: Dict = {"drug": drug_name, "disease": disease_name,
                 "trials": None, "adverse_events": None, "literature": None, "caveats": []}

    # 1) Trials + outcomes (delegate typing to endpoint_parser)
    try:
        data = http_client.get_json(
            CT_STUDIES, default={},
            params={"query.intr": drug_name, "query.cond": disease_name,
                    "pageSize": max_trials, "format": "json"}) or {}
        studies = data.get("studies", []) or []
        from services.endpoint_parser import aggregate
        agg = aggregate(studies) if studies else {"verdict": "not_reported", "outcome_signal": 0.0}
        # Report the CLASS DISTRIBUTION across matched trials, not just the single most-
        # informative verdict — a broad drug-condition query returns many trials, and one
        # terminated arm should not read as "the drug failed".
        rec["trials"] = {"n_matched": len(studies), "verdict": agg.get("verdict"),
                         "class_counts": agg.get("counts", {}),
                         "outcome_signal": agg.get("outcome_signal"), "note": agg.get("note", ""),
                         "caveat": ("verdict = the single most-informative trial; class_counts shows "
                                    "the full mix across matched trials.")}
        if prov:
            rec["trials"]["provenance"] = prov.stamp("clinicaltrials")
        # AE: the study with the largest serious denominator
        aes = [ae for ae in (parse_adverse_events(s) for s in studies) if ae]
        if aes:
            rec["adverse_events"] = max(aes, key=lambda a: a.get("serious_num_at_risk") or 0)
            if prov:
                rec["adverse_events"]["provenance"] = prov.stamp("clinicaltrials")
    except Exception as e:
        logger.debug("clinical_evidence trials fetch failed: %s", e)

    # 2) Literature tier (MeSH + publication type)
    if include_literature:
        try:
            def _esearch(term):
                return (http_client.get_json(
                    f"{EUTILS}/esearch.fcgi", default={},
                    params={"db": "pubmed", "term": term, "retmode": "json", "retmax": 50})
                    or {}).get("esearchresult", {}).get("idlist", [])
            # Prefer MeSH co-indexing (directional, study-design aware); fall back to
            # title/abstract co-mention ONLY if MeSH matches nothing, flagged honestly.
            mesh_used = True
            ids = _esearch(f'"{drug_name}"[MeSH Terms] AND "{disease_name}"[MeSH Terms]')
            if not ids:
                mesh_used = False
                ids = _esearch(f'"{drug_name}"[tiab] AND "{disease_name}"[tiab]')
            records = []
            if ids:
                xml = http_client.get_text(
                    f"{EUTILS}/efetch.fcgi",
                    params={"db": "pubmed", "id": ",".join(ids), "retmode": "xml"}) or ""
                records = _parse_efetch_xml(xml)
            lt = literature_tier(records)
            lt["query"] = "mesh" if mesh_used else "tiab_fallback"
            if prov:
                lt["provenance"] = prov.stamp("pubmed")
            rec["literature"] = lt
        except Exception as e:
            logger.debug("clinical_evidence literature fetch failed: %s", e)

    if prov:
        try:
            rec["lineage"] = prov.lineage(["clinicaltrials", "openfda_faers", "pubmed"])
        except Exception:
            pass
    rec["caveats"] = [
        "Registry completeness: only a subset of trials post structured results; absence of a "
        "result is 'not reported', not 'no effect'.",
        "This record reports what sponsors/indexers posted; it is not RepurposeIQ's own clinical "
        "determination.",
    ]
    return rec
