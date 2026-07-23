"""
Evidence Dossier Service
Generates structured evidence reports for drug-disease pairs.
Supports 505(b)(2) regulatory documentation and clinical rationale.
"""

import logging
import re
from typing import Dict, List, Optional

import requests  # noqa: F401

from services import http_client

logger = logging.getLogger(__name__)

CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"
NCBI_BASE   = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
CT_BASE     = "https://clinicaltrials.gov/api/v2/studies"

# 505(b)(2) pathway descriptions
_505B2_LEVELS = {
    "strong":   "505(b)(2) STRONGLY ELIGIBLE — Approved safety/efficacy data fully established for citation.",
    "moderate": "505(b)(2) ELIGIBLE — Phase 3 data available; can cite existing clinical evidence.",
    "possible": "505(b)(2) POSSIBLE — Phase 2 data available; partial safety data established.",
    "weak":     "505(b)(2) LIMITED — Early-stage; limited existing data for citation.",
}


def _fetch_pubmed(drug: str, disease: str, n: int = 6) -> List[Dict]:
    papers = []
    try:
        sr = http_client.get(f"{NCBI_BASE}/esearch.fcgi",
            params={"db": "pubmed", "term": f"{drug}[tiab] AND {disease}[tiab]",
                    "retmax": n, "retmode": "json"}, timeout=8)
        ids = (sr.json() if sr and sr.ok else {}).get("esearchresult", {}).get("idlist", [])
        if ids:
            smr = http_client.get(f"{NCBI_BASE}/esummary.fcgi",
                params={"db": "pubmed", "id": ",".join(ids), "retmode": "json"}, timeout=8)
            res = (smr.json() if smr and smr.ok else {}).get("result", {})
            for pid in ids:
                if pid in res:
                    papers.append({
                        "pmid":    pid,
                        "title":   res[pid].get("title", "")[:100],
                        "journal": res[pid].get("source", ""),
                        "year":    res[pid].get("pubdate", "")[:4],
                        "url":     f"https://pubmed.ncbi.nlm.nih.gov/{pid}/",
                    })
    except Exception as e:
        logger.debug(f"PubMed fetch failed for {drug}/{disease}: {e}")
    return papers


def _fetch_trials(drug: str, disease: str, n: int = 5) -> List[Dict]:
    trials = []
    try:
        r = http_client.get(CT_BASE,
            params={"query.cond": disease, "query.term": drug,
                    "pageSize": n, "format": "json"}, timeout=10)
        if r and r.ok:
            for s in r.json().get("studies", []):
                pm  = s.get("protocolSection", {})
                nct = pm.get("identificationModule", {}).get("nctId", "")
                trials.append({
                    "nct":    nct,
                    "title":  pm.get("identificationModule", {}).get("briefTitle", "")[:80],
                    "status": pm.get("statusModule", {}).get("overallStatus", ""),
                    "phase":  ", ".join(pm.get("designModule", {}).get("phases", [])) or "N/A",
                    "url":    f"https://clinicaltrials.gov/study/{nct}",
                })
    except Exception as e:
        logger.debug(f"ClinicalTrials fetch failed for {drug}/{disease}: {e}")
    return trials


def _regulatory_assessment(max_phase, disease_name: str, existing_ind: str) -> Dict:
    try:
        phase = int(float(max_phase or 0))
    except (TypeError, ValueError):
        phase = 0
    is_approved = phase >= 4

    orphan = False
    try:
        from services.disease_ontology import is_orphan_candidate
        orphan = is_orphan_candidate(disease_name)
    except Exception:
        pass

    disease_kws = [w for w in disease_name.lower().split() if len(w) > 4]
    direct_match = any(kw in (existing_ind or "").lower() for kw in disease_kws)

    if phase >= 4:
        level = "strong"
        pathway = "NDA via 505(b)(2)"
        timeline = "~2–4 years (cites approved NDA/BLA data)"
    elif phase == 3:
        level = "moderate"
        pathway = "NDA via 505(b)(2) or 505(b)(1) depending on data package"
        timeline = "~3–5 years"
    elif phase == 2:
        level = "possible"
        pathway = "IND → Phase 3 → NDA"
        timeline = "~5–8 years"
    else:
        level = "weak"
        pathway = "IND → Full development"
        timeline = "~8–12 years"

    exclusivity = []
    if orphan:
        exclusivity.append("7-year Orphan Drug Exclusivity (< 200,000 US patients)")
    if is_approved:
        exclusivity.append("3-year NCE exclusivity if new indication + clinical data")
    if not exclusivity:
        exclusivity.append("Standard 5-year data exclusivity")

    return {
        "level":          level,
        "label":          _505B2_LEVELS[level],
        "pathway":        pathway,
        "timeline":       timeline,
        "orphan_drug":    orphan,
        "direct_match":   direct_match,
        "is_approved":    is_approved,
        "exclusivity":    exclusivity,
        "phase":          phase,
    }


def _mechanistic_rationale(
    drug_genes: List[str],
    disease_genes: List[str],
    disease_pathways: List[Dict],
    drug_name: str,
    disease_name: str,
) -> Dict:
    overlap = list({g.upper() for g in drug_genes} & {g.upper() for g in disease_genes})
    shared_pathways = [
        pw for pw in disease_pathways[:10]
        if {g.upper() for g in drug_genes} & {g.upper() for g in pw.get("genes", [])}
    ]

    if overlap:
        rationale = (
            f"{drug_name} directly modulates {', '.join(overlap[:5])}, "
            f"{'which are' if len(overlap) > 1 else 'which is'} among the top-ranked "
            f"target genes for {disease_name} by Open Targets evidence score."
        )
    elif shared_pathways:
        pw_name = shared_pathways[0].get("pathway_name", "a key disease pathway")
        rationale = (
            f"{drug_name} engages {', '.join(drug_genes[:3]) or 'molecular targets'} "
            f"that participate in {pw_name}, a pathway dysregulated in {disease_name}."
        )
    else:
        rationale = (
            f"{drug_name} has established clinical activity relevant to {disease_name} "
            f"based on indication-level evidence and phase data, though direct target "
            f"overlap data is limited."
        )

    return {
        "text":            rationale,
        "overlapping_targets": overlap[:10],
        "shared_pathways": [
            {"name": pw.get("pathway_name", ""), "id": pw.get("pathway_id", "")}
            for pw in shared_pathways[:5]
        ],
    }


OPENFDA_EVENT = "https://api.fda.gov/drug/event.json"


def _fetch_trials_molecule(drug: str, max_pages: int = 3, cap: int = 25) -> Dict:
    """ALL ClinicalTrials.gov trials for a drug across EVERY condition (not just the
    selected indication) — the molecule-wide clinical picture."""
    from collections import Counter
    out = {"total": 0, "by_phase": {}, "by_status": {}, "top_conditions": [],
           "top_sponsors": [], "trials": [], "registry": "ClinicalTrials.gov"}
    cond_c, sp_c, ph_c, st_c = Counter(), Counter(), Counter(), Counter()
    fetched, token = [], None
    try:
        for _ in range(max_pages):
            params = {"query.term": drug, "pageSize": 100, "format": "json", "countTotal": "true"}
            if token:
                params["pageToken"] = token
            r = http_client.get(CT_BASE, params=params, timeout=15)
            if not (r and r.ok):
                break
            j = r.json()
            out["total"] = j.get("totalCount", out["total"])
            for s in j.get("studies", []):
                pm = s.get("protocolSection", {})
                idm = pm.get("identificationModule", {})
                status = pm.get("statusModule", {}).get("overallStatus", "")
                phases = pm.get("designModule", {}).get("phases", []) or ["N/A"]
                conds = pm.get("conditionsModule", {}).get("conditions", []) or []
                sponsor = (pm.get("sponsorCollaboratorsModule", {}).get("leadSponsor", {}) or {}).get("name", "")
                st_c[status] += 1
                for p in phases:
                    ph_c[p] += 1
                for cc in conds:
                    cond_c[cc] += 1
                if sponsor:
                    sp_c[sponsor] += 1
                fetched.append({"nct": idm.get("nctId", ""), "title": idm.get("briefTitle", "")[:90],
                                "status": status, "phase": ", ".join(phases),
                                "conditions": conds[:3], "sponsor": sponsor,
                                "url": f"https://clinicaltrials.gov/study/{idm.get('nctId','')}"})
            token = j.get("nextPageToken")
            if not token:
                break
        out["by_phase"] = dict(ph_c.most_common())
        out["by_status"] = dict(st_c.most_common())
        out["top_conditions"] = [{"name": k, "count": v} for k, v in cond_c.most_common(10)]
        out["top_sponsors"] = [{"name": k, "count": v} for k, v in sp_c.most_common(6)]
        order = {"RECRUITING": 0, "ACTIVE_NOT_RECRUITING": 1, "ENROLLING_BY_INVITATION": 2,
                 "NOT_YET_RECRUITING": 3, "COMPLETED": 4}
        fetched.sort(key=lambda t: order.get(t["status"], 9))
        out["trials"] = fetched[:cap]
    except Exception as e:
        logger.debug(f"molecule trials fetch failed for {drug}: {e}")
    return out


def _fetch_pubmed_molecule(drug: str, n: int = 12) -> Dict:
    """All PubMed literature for a drug across indications (molecule-wide)."""
    from collections import Counter
    out = {"total": 0, "recent": [], "top_journals": []}
    try:
        sr = http_client.get(f"{NCBI_BASE}/esearch.fcgi",
            params={"db": "pubmed", "term": f"{drug}[tiab]", "retmax": n,
                    "retmode": "json", "sort": "relevance"}, timeout=10)
        j = (sr.json() if sr and sr.ok else {}).get("esearchresult", {})
        out["total"] = int(j.get("count", 0) or 0)
        ids = j.get("idlist", [])
        if ids:
            smr = http_client.get(f"{NCBI_BASE}/esummary.fcgi",
                params={"db": "pubmed", "id": ",".join(ids), "retmode": "json"}, timeout=10)
            res = (smr.json() if smr and smr.ok else {}).get("result", {})
            jc = Counter()
            for pid in ids:
                if pid in res:
                    src = res[pid].get("source", "")
                    jc[src] += 1
                    out["recent"].append({"pmid": pid, "title": res[pid].get("title", "")[:110],
                                           "journal": src, "year": res[pid].get("pubdate", "")[:4],
                                           "url": f"https://pubmed.ncbi.nlm.nih.gov/{pid}/"})
            out["top_journals"] = [{"name": k, "count": v} for k, v in jc.most_common(5)]
    except Exception as e:
        logger.debug(f"molecule pubmed fetch failed for {drug}: {e}")
    return out


def _fetch_adverse_events(drug: str, n: int = 12) -> Dict:
    """Top adverse reactions for a drug from FDA FAERS (openFDA). NOTE: US
    spontaneous reports — reporting bias, counts are not incidence."""
    out = {"total": 0, "serious": 0, "top_reactions": [],
           "source": "FDA FAERS via openFDA", "caveat": "Spontaneous US reports — reporting bias; counts ≠ incidence."}
    d = (drug or "").strip().lower()
    if not d:
        return out
    search = (f'(patient.drug.openfda.generic_name:"{d}"'
              f'+patient.drug.openfda.brand_name:"{d}"'
              f'+patient.drug.medicinalproduct:"{d}")')
    try:
        r = http_client.get(OPENFDA_EVENT, params={
            "search": search, "count": "patient.reaction.reactionmeddrapt.exact"}, timeout=12)
        if r and r.ok:
            out["top_reactions"] = [{"term": x["term"].title(), "count": x["count"]}
                                    for x in r.json().get("results", [])[:n]]
        rt = http_client.get(OPENFDA_EVENT, params={"search": search, "limit": 1}, timeout=10)
        if rt and rt.ok:
            out["total"] = rt.json().get("meta", {}).get("results", {}).get("total", 0)
        rs = http_client.get(OPENFDA_EVENT, params={"search": search, "count": "serious"}, timeout=10)
        if rs and rs.ok:
            for x in rs.json().get("results", []):
                if str(x.get("term")) == "1":   # 1 = serious, 2 = non-serious
                    out["serious"] = x.get("count", 0)
    except Exception as e:
        logger.debug(f"FAERS fetch failed for {drug}: {e}")
    return out


def _registry_links(drug: str) -> List[Dict]:
    """Multi-jurisdiction trial registries + regulators (deep-link searches — these
    sources have no clean public API)."""
    import urllib.parse as _up
    q = _up.quote((drug or "").strip())
    return [
        {"name": "ClinicalTrials.gov", "region": "US / global", "kind": "Trials",
         "url": f"https://clinicaltrials.gov/search?term={q}"},
        {"name": "WHO ICTRP", "region": "Global", "kind": "Trials",
         "desc": "Aggregates CTRI, EU-CTR, ISRCTN & more",
         "url": f"https://trialsearch.who.int/Default.aspx?SearchTermStat={q}"},
        {"name": "CTRI (India)", "region": "India", "kind": "Trials",
         "desc": "Clinical Trials Registry – India",
         "url": "https://ctri.nic.in/Clinicaltrials/advsearch.php"},
        {"name": "Drugs@FDA", "region": "US", "kind": "Regulatory",
         "url": f"https://www.accessdata.fda.gov/scripts/cder/daf/index.cfm?event=BasicSearch.process&searchTerm={q}"},
        {"name": "EMA", "region": "EU", "kind": "Regulatory",
         "url": f"https://www.ema.europa.eu/en/search?search_api_fulltext={q}"},
        {"name": "MHRA (UK)", "region": "UK", "kind": "Regulatory",
         "desc": "UK products & Public Assessment Reports",
         "url": f"https://products.mhra.gov.uk/search/?search={q}&page=1"},
        {"name": "EudraGMDP", "region": "EU", "kind": "Manufacturing / GMP",
         "desc": "GMP/GDP authorisations & inspections — CMC / supply-chain due-diligence",
         "url": "https://eudragmdp.ema.europa.eu/inspections/displayWelcome.do"},
    ]


def generate_dossier(
    chembl_id: str,
    disease_name: str,
    compound: Optional[Dict] = None,
    screen_result: Optional[Dict] = None,
) -> Dict:
    """
    Generate evidence dossier for a drug–disease pair.
    compound: pre-loaded compound dict (from DB / ChEMBL API)
    screen_result: output from repurposing_engine.run_repurposing_screen (optional)
    """
    drug_name = (compound or {}).get("name", chembl_id)
    max_phase = (compound or {}).get("max_phase", 0) or 0
    existing_ind = (compound or {}).get("indications", "") or ""
    smiles = (compound or {}).get("smiles", "")

    # ── Disease context ───────────────────────────────────────────────────────
    disease_info: Dict = {}
    disease_genes: List[str] = []
    disease_pathways: List[Dict] = []

    try:
        from services.disease_ontology import resolve_disease as ot_resolve
        disease_info    = ot_resolve(disease_name)
        disease_genes   = [t["gene_symbol"] for t in disease_info.get("targets", [])[:30]]
        disease_pathways = disease_info.get("pathways", [])
    except Exception as e:
        logger.debug(f"Disease ontology error: {e}")

    # ── Drug targets ──────────────────────────────────────────────────────────
    targets_raw = (compound or {}).get("targets", "") or ""
    drug_genes = [g.strip() for g in re.split(r"[;,]", targets_raw) if g.strip()]

    # ── Score (reuse from screen if available, else compute) ──────────────────
    scores: Dict = {}
    composite_score = 0.0
    if screen_result:
        for cand in screen_result.get("candidates", []):
            cid = cand.get("chembl_id", "")
            if cid == chembl_id or cand.get("name", "").lower() == drug_name.lower():
                scores = cand.get("scores", {})
                composite_score = cand.get("composite_score", 0.0)
                if not drug_genes and cand.get("drug_genes"):
                    drug_genes = cand["drug_genes"]
                break

    if not scores:
        # Canonical score — the SAME function the Repurpose card and analysis header use,
        # so the dossier headline matches the rest of the platform for this pair.
        try:
            from services.reverse_repurposing import canonical_pair_score
            ps = canonical_pair_score(
                chembl_id=chembl_id, disease=disease_name,
                drug_genes=(drug_genes or None),
                max_phase=int(float(max_phase or 0)),
                indications=existing_ind, drug_name=drug_name,
            )
            scores = ps.get("scores", {})
            composite_score = ps.get("score", 0.0)
            if not drug_genes and ps.get("drug_genes"):
                drug_genes = ps["drug_genes"]
        except Exception as e:
            logger.debug(f"Dossier canonical scoring error: {e}")

    # ── Sub-components ────────────────────────────────────────────────────────
    regulatory    = _regulatory_assessment(max_phase, disease_name, existing_ind)
    mechanistic   = _mechanistic_rationale(drug_genes, disease_genes, disease_pathways, drug_name, disease_name)
    literature    = _fetch_pubmed(drug_name, disease_name)
    trials        = _fetch_trials(drug_name, disease_name)
    # ── Molecule-wide evidence (all indications) + safety + registries ─────────
    molecule_trials     = _fetch_trials_molecule(drug_name)
    molecule_literature = _fetch_pubmed_molecule(drug_name)
    adverse_events      = _fetch_adverse_events(drug_name)
    registry_links      = _registry_links(drug_name)

    # ── Physicochemical profile ───────────────────────────────────────────────
    physchem: Dict = {}
    if smiles:
        try:
            from rdkit import Chem
            from rdkit.Chem import Descriptors, rdMolDescriptors
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                mw   = round(Descriptors.MolWt(mol), 1)
                logp = round(Descriptors.MolLogP(mol), 2)
                psa  = round(Descriptors.TPSA(mol), 1)
                hbd  = rdMolDescriptors.CalcNumHBD(mol)
                # BBB score
                bbb = 0
                if psa < 60:  bbb += 30
                elif psa < 90: bbb += 15
                if mw < 400:  bbb += 25
                elif mw < 500: bbb += 10
                if 1 <= logp <= 3:     bbb += 25
                elif 0 <= logp <= 4.5: bbb += 10
                if hbd <= 2:  bbb += 20
                elif hbd <= 3: bbb += 10
                physchem = {
                    "mw": mw, "logp": logp, "psa": psa, "hbd": hbd,
                    "bbb_score": bbb,
                    "bbb_level": "High" if bbb >= 75 else "Moderate" if bbb >= 45 else "Low",
                    "qed": round(Descriptors.qed(mol), 3),
                }
        except Exception:
            pass

    # ── Area-aware developability (route-appropriate, not CNS-only) ────────────
    developability: Dict = {}
    if smiles:
        try:
            from services import developability as _dev
            therapeutic_areas: List[str] = []
            efo = disease_info.get("disease_id", "")
            if efo:
                ta_q = """query($id:String!){ disease(efoId:$id){ therapeuticAreas { name } } }"""
                ta_r = http_client.post("https://api.platform.opentargets.org/api/v4/graphql",
                                     json={"query": ta_q, "variables": {"id": efo}}, timeout=8)
                if ta_r and ta_r.ok:
                    d = (ta_r.json().get("data") or {}).get("disease") or {}
                    therapeutic_areas = [t.get("name", "") for t in (d.get("therapeuticAreas") or [])]
            developability = _dev.score(smiles, therapeutic_areas=therapeutic_areas)
        except Exception as e:
            logger.debug(f"developability error: {e}")

    return {
        "chembl_id":      chembl_id,
        "drug_name":      drug_name,
        "disease_name":   disease_name,
        "composite_score": round(composite_score, 4),
        "scores":         scores,
        "screens":        {
            "target":     0.25, "pathway":    0.20, "ppi":        0.20,
            "clinical":   0.15, "indication": 0.10, "regulatory": 0.10,
        },
        "disease_context": {
            "efo_id":      disease_info.get("disease_id", ""),
            "description": disease_info.get("description", "")[:400],
            "gene_count":  len(disease_genes),
            "top_genes":   disease_genes[:12],
            "pathway_count": len(disease_pathways),
            "top_pathways": [
                {"name": pw.get("pathway_name", ""), "id": pw.get("pathway_id", "")}
                for pw in disease_pathways[:5]
            ],
        },
        "drug_targets":       drug_genes[:15],
        "mechanistic":        mechanistic,
        "regulatory":         regulatory,
        "physchem":           physchem,
        "developability":     developability,
        "literature":         literature,
        "clinical_trials":    trials,
        "molecule_trials":      molecule_trials,
        "molecule_literature":  molecule_literature,
        "adverse_events":       adverse_events,
        "registry_links":       registry_links,
        "max_phase":          max_phase,
        "existing_indications": existing_ind,
    }
