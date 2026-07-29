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


# ── Card display columns (Safety / Efficacy) ─────────────────────────────────
# These SURFACE facts already computed by parse_adverse_events (serious-AE rate with
# a real denominator) and endpoint_parser.aggregate (primary-endpoint verdict) as two
# human-readable card columns. They add NOTHING to any composite/rank — display only,
# so the clinical signal the engine already consumes internally is not double counted.
#
# HARD HONESTY RULE: a pair with NO trials reads "No trial data" for BOTH columns
# (muted, never green, never blank). Absence of a trial is not evidence of safety or
# efficacy. A pair WITH trials but no posted AE table / no parseable primary reads a
# distinct, still-muted "not reported", never a default that implies safe/effective.

_EFFICACY_LABELS = {
    "met_primary": "Met primary endpoint",
    "failed_primary": "Failed primary",
    "biomarker_only": "Biomarker only",
    "terminated_efficacy": "Terminated (efficacy)",
    "terminated_safety": "Terminated (safety)",
}


def _phase_phrase(phase) -> str:
    ph = int(phase or 0)
    return f", Phase {ph}" if ph else ""


def efficacy_column(verdict, *, n_trials: int = 0, phase: int = 0, note: str = "") -> Dict:
    """Plain-language Efficacy column from an endpoint_parser verdict.
    No trials -> 'No trial data' (muted). Trials present but no parseable primary
    endpoint -> 'Outcome not reported' (muted). A met primary is the only positive."""
    n = int(n_trials or 0)
    ph = int(phase or 0)
    if n <= 0:
        return {"status": "no_data", "verdict": None, "label": "No trial data", "muted": True,
                "tooltip": "No registered trials for this drug and disease pair."}
    label = _EFFICACY_LABELS.get(verdict)
    if not label:
        return {"status": "not_reported", "verdict": verdict or "unknown",
                "label": "Outcome not reported", "muted": True,
                "tooltip": (f"{n} trial(s)" + (f" up to Phase {ph}" if ph else "")
                            + "; no primary endpoint result posted or parseable in platform.")}
    tip = (note + " " if note else "") + f"Based on {n} trial(s)" + (f" up to Phase {ph}" if ph else "") + "."
    return {"status": "computed", "verdict": verdict, "label": label,
            "muted": (verdict != "met_primary"), "phase": ph, "n_trials": n,
            "tooltip": tip.strip()}


def safety_column(ae, *, n_trials: int = 0, phase: int = 0) -> Dict:
    """Serious-AE Safety column from parse_adverse_events output.
    No trials -> 'No trial data'. Trials but no AE table -> 'Serious AE not reported'.
    AE table without an at-risk denominator -> 'Rate not computable'. Denominator present
    -> 'Serious AE X% (n=N, Phase P)'. All non-computed states are muted, never green."""
    n = int(n_trials or 0)
    ph = int(phase or (ae or {}).get("phase") or 0)
    if n <= 0:
        return {"status": "no_data", "label": "No trial data", "muted": True,
                "tooltip": "No registered trials for this drug and disease pair."}
    if not ae:
        return {"status": "not_reported", "label": "Serious AE not reported", "muted": True,
                "tooltip": (f"{n} trial(s)" + (f" up to Phase {ph}" if ph else "")
                            + "; no structured adverse event table posted.")}
    at_risk = int(ae.get("serious_num_at_risk") or 0)
    if ae.get("no_denominator") or not at_risk:
        return {"status": "no_denominator", "label": "Rate not computable", "muted": True,
                "tooltip": ("Serious adverse events were posted without a denominator for "
                            "participants at risk, so a rate cannot be computed.")}
    rate = ae.get("serious_rate") or 0.0
    affected = int(ae.get("serious_num_affected") or 0)
    pct = round(rate * 100)
    label = f"Serious AE {pct}% (n={at_risk}" + _phase_phrase(ph) + ")"
    tip = (f"{affected} of {at_risk} participants had a serious adverse event in the largest "
           f"trial that posted results" + (f" (Phase {ph})" if ph else "")
           + ". This is a trial specific rate, not population incidence; do not compare with FAERS.")
    return {"status": "computed", "label": label, "muted": False, "serious_rate": rate,
            "serious_num_affected": affected, "serious_num_at_risk": at_risk,
            "phase": ph, "n_trials": n, "tooltip": tip}


def safety_approval_line(max_phase, approved_here: bool = False,
                         market_status: str = "") -> Dict:
    """DISPLAY-ONLY established-safety header for the Safety column.

    An already-approved drug has cleared human safety trials — a DERISKING POSITIVE
    for repurposing, not a risk. So the Safety column LEADS with what human safety
    exposure the molecule already has, from data already on the candidate, and the
    trial serious-AE rate is demoted to secondary context (see safety_column).

    Reads ONLY: global max_phase, approved_here, market_status. Changes no score.
    Returns {label, muted}. Never green (established safety is stated, not scored)."""
    try:
        mp = int(float(max_phase or 0))
    except (TypeError, ValueError):
        mp = 0
    ms = (market_status or "").strip().lower()
    if mp >= 4 or approved_here:
        label = "Approved drug, established human safety profile"
        muted = False
    elif mp >= 1:
        label = f"Phase {mp}, human safety data available"
        muted = False
    else:
        label = "No approved human safety data"
        muted = True
    # A withdrawn/obsolete market status is a caution even for an approved drug —
    # append it honestly rather than implying an unqualified safe profile.
    if ms and ("withdrawn" in ms or "obsolete" in ms or "discontinued" in ms):
        label += " (withdrawn from market)"
        muted = True
    return {"label": label, "muted": muted}


def severity_qualifier(clinical_constraints: Optional[Dict]) -> str:
    """DISPLAY-ONLY qualifier that reads the ALREADY-computed clinical_constraints
    severity of the TARGET indication and frames the acceptable safety bar for it.
    Returns "" when the field is absent (never fabricated)."""
    cc = clinical_constraints or {}
    sev = cc.get("disease_severity")
    if sev == "high":
        return "acceptable bar for a life threatening indication"
    if sev == "low":
        return "high safety bar for a benign indication"
    return ""


def _study_max_phase(study: Dict) -> int:
    """Highest numeric trial phase from a CT.gov v2 study (0 if none / N/A)."""
    _MAP = {"PHASE4": 4, "PHASE3": 3, "PHASE2": 2, "PHASE1": 1, "EARLY_PHASE1": 1}
    phases = ((study.get("protocolSection") or {}).get("designModule") or {}).get("phases") or []
    return max((_MAP.get((p or "").upper(), 0) for p in phases), default=0)


# Per-process memo so repeated forward-screen renders of the same (drug, disease) pair
# do not re-hit ClinicalTrials.gov. Bounded; display-only, so staleness is harmless.
_PAIR_COL_CACHE: Dict[str, Dict] = {}


def pair_clinical_columns(drug_name: str, disease_name: str, max_trials: int = 50) -> Dict:
    """ONE ClinicalTrials.gov call for a (drug, disease) pair -> {safety_ct, efficacy_ct}
    display columns. Intended for the FORWARD screen, where the engine does not otherwise
    fetch per-candidate trials — so the CALLER must cap this to a few top leads to avoid a
    live-call storm. Fail-soft: any error or no trials yields honest 'No trial data' for both."""
    key = f"{(drug_name or '').strip().lower()}|{(disease_name or '').strip().lower()}"
    if key in _PAIR_COL_CACHE:
        return _PAIR_COL_CACHE[key]
    result = {"safety_ct": safety_column(None, n_trials=0),
              "efficacy_ct": efficacy_column(None, n_trials=0)}
    if not (drug_name or "").strip() or not (disease_name or "").strip():
        return result
    try:
        from services import http_client
        data = http_client.get_json(
            CT_STUDIES, default={},
            params={"query.intr": drug_name, "query.cond": disease_name,
                    "pageSize": max_trials, "format": "json"}) or {}
        studies = data.get("studies", []) or []
        n = len(studies)
        if studies:
            from services.endpoint_parser import aggregate
            agg = aggregate(studies)
            max_ph = max((_study_max_phase(s) for s in studies), default=0)
            result["efficacy_ct"] = efficacy_column(
                agg.get("verdict"), n_trials=n, phase=max_ph, note=agg.get("note", ""))
            best_ae, best_n = None, -1
            for s in studies:
                ae = parse_adverse_events(s)
                if ae is None:
                    continue
                m = ae.get("serious_num_at_risk") or 0
                if m > best_n:
                    best_ae, best_n = {**ae, "phase": _study_max_phase(s)}, m
            result["safety_ct"] = safety_column(best_ae, n_trials=n, phase=max_ph)
    except Exception as e:
        logger.debug("pair_clinical_columns failed for %s / %s: %s", drug_name, disease_name, e)
    _PAIR_COL_CACHE[key] = result
    return result
