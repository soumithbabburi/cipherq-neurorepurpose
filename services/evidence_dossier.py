"""
Evidence Dossier Service
Generates structured evidence reports for drug-disease pairs.
Supports 505(b)(2) regulatory documentation and clinical rationale.
"""

import logging
import re
from typing import Dict, List, Optional

import requests

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
        sr = requests.get(f"{NCBI_BASE}/esearch.fcgi",
            params={"db": "pubmed", "term": f"{drug}[tiab] AND {disease}[tiab]",
                    "retmax": n, "retmode": "json"}, timeout=8)
        ids = sr.json().get("esearchresult", {}).get("idlist", [])
        if ids:
            smr = requests.get(f"{NCBI_BASE}/esummary.fcgi",
                params={"db": "pubmed", "id": ",".join(ids), "retmode": "json"}, timeout=8)
            res = smr.json().get("result", {})
            for pid in ids:
                if pid in res:
                    papers.append({
                        "pmid":    pid,
                        "title":   res[pid].get("title", "")[:100],
                        "journal": res[pid].get("source", ""),
                        "year":    res[pid].get("pubdate", "")[:4],
                        "url":     f"https://pubmed.ncbi.nlm.nih.gov/{pid}/",
                    })
    except Exception:
        pass
    return papers


def _fetch_trials(drug: str, disease: str, n: int = 5) -> List[Dict]:
    trials = []
    try:
        r = requests.get(CT_BASE,
            params={"query.cond": disease, "query.term": drug,
                    "pageSize": n, "format": "json"}, timeout=10)
        if r.ok:
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
    except Exception:
        pass
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
                ta_r = requests.post("https://api.platform.opentargets.org/api/v4/graphql",
                                     json={"query": ta_q, "variables": {"id": efo}}, timeout=8)
                if ta_r.ok:
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
        "max_phase":          max_phase,
        "existing_indications": existing_ind,
    }
