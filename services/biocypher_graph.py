"""
Biomedical Knowledge Graph
Builds Disease → Gene → Protein → Pathway ← Compound evidence networks from:
  - DRKG (Drug Repurposing Knowledge Graph) — local cache when available
    Integrates DrugBank, Hetionet, GNBR, String, IntAct, DGIdb (5.8M triplets)
  - Open Targets Platform API (disease→gene evidence scores)
  - ChEMBL REST API (compound indications + mechanism/target data)
  - STRING DB (gene-gene protein interaction network)
"""
import logging
import re
import time
from typing import Dict, List, Optional, Tuple

import requests  # noqa: F401

from services import http_client

logger = logging.getLogger(__name__)

# ── DRKG cache (loaded once, kept in memory) ───────────────────────────────────
_drkg_cache: Optional[dict] = None
_drkg_loaded = False


def _get_drkg() -> dict:
    """Lazy-load DRKG cache."""
    global _drkg_cache, _drkg_loaded
    if not _drkg_loaded:
        try:
            from database.drkg_loader import load_cache
            _drkg_cache = load_cache() or None
        except Exception as e:
            logger.debug(f"DRKG cache load: {e}")
            _drkg_cache = None
        _drkg_loaded = True
    return _drkg_cache or {}


def drkg_available() -> bool:
    return bool(_get_drkg().get("_ts"))

OT_URL      = "https://api.platform.opentargets.org/api/v4/graphql"
CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"
STRING_URL  = "https://string-db.org/api/json"

# ── Node palette ───────────────────────────────────────────────────────────────
_COLORS = {
    "Disease":   "#e11d48",
    "Compound":  "#630ed4",
    "Gene":      "#059669",
    "Protein":   "#0284c7",
    "Pathway":   "#d97706",
}

# Gene symbol → human-readable protein name (major drug targets across all disease areas)
_PROTEIN_MAP: Dict[str, str] = {
    # CNS / Neurological targets
    "ACHE":    "Acetylcholinesterase",
    "BCHE":    "Butyrylcholinesterase",
    "CHRM1":   "Muscarinic M1 Receptor",
    "BACE1":   "β-Secretase 1",
    "BACE2":   "β-Secretase 2",
    "APP":     "Amyloid Precursor Protein",
    "PSEN1":   "Presenilin-1",
    "PSEN2":   "Presenilin-2",
    "MAPT":    "Microtubule-Assoc. Tau",
    "GSK3B":   "Glycogen Synthase Kinase 3β",
    "CDK5":    "Cyclin-Dependent Kinase 5",
    "APOE":    "Apolipoprotein E",
    "CLU":     "Clusterin",
    "SNCA":    "α-Synuclein",
    "LRRK2":   "LRRK2 Kinase",
    "PINK1":   "PINK1 Kinase",
    "PARK2":   "Parkin Ubiquitin Ligase",
    "HTT":     "Huntingtin Protein",
    "DRD1":    "Dopamine D1 Receptor",
    "DRD2":    "Dopamine D2 Receptor",
    "DRD3":    "Dopamine D3 Receptor",
    "DRD4":    "Dopamine D4 Receptor",
    "DRD5":    "Dopamine D5 Receptor",
    "HTR1A":   "5-HT1A Serotonin Receptor",
    "HTR2A":   "5-HT2A Serotonin Receptor",
    "HTR2C":   "5-HT2C Serotonin Receptor",
    "SLC6A4":  "Serotonin Transporter",
    "SLC6A3":  "Dopamine Transporter",
    "SLC18A2": "Vesicular Monoamine Transporter 2",
    "SLC6A2":  "Norepinephrine Transporter",
    "MAOB":    "Monoamine Oxidase B",
    "MAOA":    "Monoamine Oxidase A",
    "GRIN2A":  "NMDA Receptor NR2A",
    "GRIN2B":  "NMDA Receptor NR2B",
    "GRIA1":   "AMPA Receptor GluA1",
    "GABRA1":  "GABA-A Receptor α1",
    "GABRA2":  "GABA-A Receptor α2",
    "ADRA1A":  "Adrenergic α1A Receptor",
    "ADRA2A":  "Adrenergic α2A Receptor",
    "ADRB2":   "Adrenergic β2 Receptor",
    "HDAC1":   "Histone Deacetylase 1",
    "HDAC2":   "Histone Deacetylase 2",
    "SIRT1":   "Sirtuin-1 Deacetylase",
    "BECN1":   "Beclin-1 Autophagy Regulator",
    "TREM2":   "TREM2 Microglial Receptor",
    "BDNF":    "Brain-Derived Neurotrophic Factor",
    "NTRK2":   "TrkB Neurotrophin Receptor",
    "NGF":     "Nerve Growth Factor",
    # Oncology targets
    "EGFR":    "EGF Receptor",
    "HER2":    "HER2/ErbB2 Receptor",
    "KRAS":    "KRAS GTPase",
    "PIK3CA":  "PI3K p110α Catalytic Subunit",
    "PTEN":    "PTEN Phosphatase",
    "TP53":    "Tumor Protein p53",
    "BCL2":    "B-cell Lymphoma 2",
    "CDK4":    "Cyclin-Dependent Kinase 4",
    "CDK6":    "Cyclin-Dependent Kinase 6",
    "VEGFA":   "VEGF-A",
    "MAPK1":   "ERK2 MAP Kinase",
    "MAPK3":   "ERK1 MAP Kinase",
    "AKT1":    "AKT1 Kinase",
    "MTOR":    "mTOR Kinase",
    "HDAC1":   "Histone Deacetylase 1",
    # Cardiovascular / Metabolic targets
    "ACE":     "Angiotensin-Converting Enzyme",
    "PPARG":   "PPAR-γ Nuclear Receptor",
    "PPARA":   "PPAR-α Nuclear Receptor",
    "HMGCR":   "HMG-CoA Reductase",
    "DPP4":    "Dipeptidyl Peptidase 4",
    "SLC2A1":  "GLUT1 Glucose Transporter",
    "SLC2A2":  "GLUT2 Glucose Transporter",
    "SLC2A3":  "GLUT3 Glucose Transporter",
    "NOS1":    "Neuronal Nitric Oxide Synthase",
    "NOS3":    "Endothelial Nitric Oxide Synthase",
    "PRKAA1":  "AMPK α1 Subunit",
    "PRKAA2":  "AMPK α2 Subunit",
    # Autoimmune / Inflammatory targets
    "TNF":     "Tumor Necrosis Factor-α",
    "IL6":     "Interleukin-6",
    "IL1B":    "Interleukin-1β",
    "PTGS2":   "Cyclooxygenase-2 (COX-2)",
    "PTGS1":   "Cyclooxygenase-1 (COX-1)",
}

_NEURO_OVERVIEW = [
    ("Alzheimer disease",  "MONDO_0004975"),
    ("Parkinson disease",  "MONDO_0005180"),
    ("Multiple Sclerosis", "MONDO_0005301"),
    ("Huntington disease", "MONDO_0007739"),
]

# ── Static Reactome pathway knowledge base ─────────────────────────────────────
_PATHWAY_MAP: Dict[str, Tuple[str, str]] = {
    # Gene symbol → (pathway_name, reactome_stable_id)
    "ACHE":    ("Cholinergic Transmission",     "R-HSA-352230"),
    "BCHE":    ("Cholinergic Transmission",     "R-HSA-352230"),
    "CHRM1":   ("Cholinergic Transmission",     "R-HSA-352230"),
    "BACE1":   ("APP Cleavage",                 "R-HSA-9609736"),
    "BACE2":   ("APP Cleavage",                 "R-HSA-9609736"),
    "APP":     ("Amyloid Precursor Processing", "R-HSA-9609736"),
    "PSEN1":   ("Presenilin Signaling",         "R-HSA-1980143"),
    "PSEN2":   ("Presenilin Signaling",         "R-HSA-1980143"),
    "MAPT":    ("Tau Phosphorylation",          "R-HSA-9615710"),
    "GSK3B":   ("Tau Phosphorylation",          "R-HSA-9615710"),
    "CDK5":    ("Tau Phosphorylation",          "R-HSA-9615710"),
    "APOE":    ("Lipid Metabolism",             "R-HSA-556833"),
    "CLU":     ("Lipid Metabolism",             "R-HSA-556833"),
    "SNCA":    ("Alpha-Synuclein Signaling",    "R-HSA-983168"),
    "LRRK2":   ("LRRK2 Signaling",             "R-HSA-5637815"),
    "PINK1":   ("Mitophagy",                    "R-HSA-5205647"),
    "PARK2":   ("Mitophagy",                    "R-HSA-5205647"),
    "HTT":     ("Huntingtin Interactions",      "R-HSA-5633007"),
    "DRD1":    ("Dopamine Receptor Signaling",  "R-HSA-390651"),
    "DRD2":    ("Dopamine Receptor Signaling",  "R-HSA-390651"),
    "DRD3":    ("Dopamine Receptor Signaling",  "R-HSA-390651"),
    "DRD4":    ("Dopamine Receptor Signaling",  "R-HSA-390651"),
    "DRD5":    ("Dopamine Receptor Signaling",  "R-HSA-390651"),
    "HTR1A":   ("Serotonin Receptor Signaling", "R-HSA-390666"),
    "HTR2A":   ("Serotonin Receptor Signaling", "R-HSA-390666"),
    "HTR2C":   ("Serotonin Receptor Signaling", "R-HSA-390666"),
    "SLC6A4":  ("Serotonin Transport",          "R-HSA-442660"),
    "SLC6A3":  ("Dopamine Transport",           "R-HSA-442660"),
    "SLC18A2": ("Monoamine Transport",          "R-HSA-442660"),
    "ADRA1A":  ("Adrenergic Signaling",         "R-HSA-390697"),
    "ADRA2A":  ("Adrenergic Signaling",         "R-HSA-390697"),
    "ADRA2B":  ("Adrenergic Signaling",         "R-HSA-390697"),
    "ADRB2":   ("Adrenergic Signaling",         "R-HSA-390697"),
    "GABRA1":  ("GABA Receptor Signaling",      "R-HSA-977443"),
    "GABRA2":  ("GABA Receptor Signaling",      "R-HSA-977443"),
    "GRIN2A":  ("Glutamate Signaling",          "R-HSA-416993"),
    "GRIN2B":  ("Glutamate Signaling",          "R-HSA-416993"),
    "GRIA1":   ("AMPA Receptor Signaling",      "R-HSA-399719"),
    "PRKAA1":  ("AMPK-mTOR Signaling",          "R-HSA-165159"),
    "PRKAA2":  ("AMPK-mTOR Signaling",          "R-HSA-165159"),
    "MTOR":    ("AMPK-mTOR Signaling",          "R-HSA-165159"),
    "SIRT1":   ("Autophagy",                    "R-HSA-9612973"),
    "BECN1":   ("Autophagy",                    "R-HSA-9612973"),
    "TNF":     ("Neuroinflammation",            "R-HSA-5620971"),
    "IL6":     ("Neuroinflammation",            "R-HSA-5620971"),
    "IL1B":    ("Neuroinflammation",            "R-HSA-5620971"),
    "TREM2":   ("Microglial Signaling",         "R-HSA-5620971"),
    "NOS1":    ("Nitric Oxide Signaling",       "R-HSA-392154"),
    "NOS3":    ("Nitric Oxide Signaling",       "R-HSA-392154"),
    "MAPK1":   ("MAPK Signaling",               "R-HSA-5683057"),
    "MAPK3":   ("MAPK Signaling",               "R-HSA-5683057"),
    "AKT1":    ("PI3K-Akt Signaling",           "R-HSA-9612973"),
    "PIK3CA":  ("PI3K-Akt Signaling",           "R-HSA-9612973"),
    "PPARG":   ("PPAR Signaling",               "R-HSA-383280"),
    "PPARA":   ("PPAR Signaling",               "R-HSA-383280"),
    "HDAC1":   ("Histone Deacetylation",        "R-HSA-3214815"),
    "HDAC2":   ("Histone Deacetylation",        "R-HSA-3214815"),
    "SLC2A1":  ("Glucose Transport",            "R-HSA-189200"),
    "SLC2A3":  ("Glucose Transport",            "R-HSA-189200"),
    "VEGFA":   ("VEGF Signaling",               "R-HSA-194138"),
    "BDNF":    ("BDNF-TrkB Signaling",          "R-HSA-9013148"),
    "NTRK2":   ("BDNF-TrkB Signaling",          "R-HSA-9013148"),
    "NGF":     ("NGF Signaling",                "R-HSA-9013148"),
    "SLC6A2":  ("Norepinephrine Transport",     "R-HSA-442660"),
    "SLC2A2":  ("Glucose Transport",            "R-HSA-189200"),
}

# ── ChEMBL target cache (avoids re-fetching same targets across drugs) ─────────
_CHEMBL_TARGET_CACHE: Dict[str, dict] = {}


# ── Node / edge constructors ───────────────────────────────────────────────────

def _net_node(nid: str, name: str, kind: str, **extra) -> dict:
    color = _COLORS.get(kind, "#334155")
    size  = 72 if kind == "Disease" else (52 if kind == "Compound" else (36 if kind == "Protein" else 40))
    label = name[:22] + ("…" if len(name) > 22 else "")
    return {"data": {"id": nid, "label": label, "full_label": name,
                     "kind": kind, "color": color, "size": size, **extra}}


def _net_edge(src: str, tgt: str, label: str, weight: float = 1.0) -> dict:
    return {"data": {"source": src, "target": tgt, "label": label, "weight": weight}}


# Keep old names available for flask_app.py import
def _node(id_, label, kind, data=None):
    return _net_node(id_, label, kind, **(data or {}))


def _edge(source, target, label, weight=1.0):
    return _net_edge(source, target, label, weight)


# ── API helpers ────────────────────────────────────────────────────────────────

def _ot_disease_targets(efo_id: str, n: int = 25) -> List[Dict]:
    """Fetch top gene targets for a disease from Open Targets."""
    q = """
    query($id: String!, $n: Int!) {
      disease(efoId: $id) {
        name
        associatedTargets(page: {index: 0, size: $n}) {
          rows {
            target { approvedSymbol approvedName }
            score
          }
        }
      }
    }"""
    try:
        r = http_client.post(OT_URL,
            json={"query": q, "variables": {"id": efo_id, "n": n}},
            timeout=12, headers={"Content-Type": "application/json"})
        if r and r.ok:
            dis = r.json().get("data", {}).get("disease") or {}
            rows = dis.get("associatedTargets", {}).get("rows", [])
            return [
                {"gene_symbol": row["target"]["approvedSymbol"],
                 "gene_name":   row["target"]["approvedName"],
                 "score":       round(row["score"], 4)}
                for row in rows if row.get("score", 0) > 0.05
            ]
    except Exception as e:
        logger.debug(f"OT targets error: {e}")
    return []


_GENE_PW_CACHE: Dict[str, list] = {}

def _gene_pathways(gene: str, top: int = 2) -> List[Tuple[str, str]]:
    """Reactome pathways a target gene participates in → [(name, reactome_id)].
    Curated neuro map first (fast/relevant); else live from Open Targets."""
    g = (gene or "").upper().strip()
    if not g:
        return []
    if g in _PATHWAY_MAP:
        return [_PATHWAY_MAP[g]]
    if g in _GENE_PW_CACHE:
        return _GENE_PW_CACHE[g][:top]
    out: List[Tuple[str, str]] = []
    q = """query($s:String!){ search(queryString:$s, entityNames:["target"],
        page:{index:0,size:3}){ hits{ object{ __typename
        ... on Target{ approvedSymbol pathways{ pathwayId pathway } } } } } }"""
    try:
        r = http_client.post(OT_URL, json={"query": q, "variables": {"s": g}},
                             timeout=12, headers={"Content-Type": "application/json"})
        if r and r.ok:
            for h in (r.json().get("data", {}).get("search", {}) or {}).get("hits", []):
                obj = h.get("object") or {}
                if obj.get("__typename") == "Target" and (obj.get("approvedSymbol", "").upper() == g):
                    seen = set()
                    for pw in (obj.get("pathways") or []):
                        pid, pname = pw.get("pathwayId"), pw.get("pathway")
                        if pid and pname and pid not in seen:
                            seen.add(pid); out.append((pname, pid))
                    break
    except Exception as e:
        logger.debug(f"gene pathways error {g}: {e}")
    _GENE_PW_CACHE[g] = out
    return out[:top]


def _chembl_indications(disease: str, limit: int = 15) -> List[Dict]:
    """ChEMBL drug indications for a disease name."""
    seen, results = set(), []
    for param in ("mesh_heading__icontains", "efo_term__icontains"):
        if len(results) >= limit:
            break
        try:
            r = http_client.get(f"{CHEMBL_BASE}/drug_indication.json",
                params={param: disease, "limit": limit, "format": "json"}, timeout=10)
            if r and r.ok:
                for ind in r.json().get("drug_indications", []):
                    mid = ind.get("molecule_chembl_id", "")
                    if mid and mid not in seen:
                        seen.add(mid)
                        results.append({
                            "chembl_id":  mid,
                            "indication": ind.get("mesh_heading") or ind.get("efo_term") or disease,
                            "max_phase":  ind.get("max_phase_for_ind") or 0,
                        })
        except Exception:
            pass
    return results


def _compound_top_disease(chembl_id: str) -> str:
    """Top indication (disease name) for a compound, by highest clinical phase.
    Lets the repurposing knowledge graph build when a molecule is opened
    without a pre-selected disease (e.g. straight from a molecule search)."""
    if not chembl_id:
        return ""
    try:
        r = http_client.get(f"{CHEMBL_BASE}/drug_indication.json",
            params={"molecule_chembl_id": chembl_id, "limit": 50, "format": "json"}, timeout=10)
        if r and r.ok:
            inds = r.json().get("drug_indications", [])
            inds.sort(key=lambda d: d.get("max_phase_for_ind") or 0, reverse=True)
            for ind in inds:
                name = ind.get("mesh_heading") or ind.get("efo_term")
                if name:
                    return name
    except Exception:
        pass
    return ""


def build_compound_graph(drug_id: str, top_n: int = 14) -> Tuple[List[dict], dict]:
    """Molecule-centric repurposing knowledge graph — works with NO pre-selected
    disease. Built from the same reverse-repurposing engine that powers the
    candidate list, so graph and list always agree:

        Compound  --BINDS_TO-->  Target gene  --ASSOCIATED_WITH-->  Candidate disease

    A repurposing scientist starts from the molecule, not a disease — this is the
    graph they actually need.
    """
    try:
        from services.reverse_repurposing import screen_indications_for_drug
        r = screen_indications_for_drug(drug_id)
    except Exception as e:
        logger.debug(f"compound graph — screen error: {e}")
        return [], {}

    chembl  = r.get("chembl_id") or drug_id
    name    = r.get("drug") or drug_id
    targets = r.get("drug_targets") or []
    cands   = (r.get("candidates") or [])[:top_n]
    if not targets and not cands:
        return [], {}

    elements, legend, seen = [], {}, set()

    def add_node(el):
        nid = el["data"].get("id")
        if nid in seen:
            return
        seen.add(nid)
        elements.append(el)
        legend[el["data"]["kind"]] = el["data"].get("color")

    # Focal molecule at the centre
    add_node(_net_node("cmp_center", name, "Compound", chembl_id=chembl,
                       focal=True, size=64, source_db="ChEMBL"))

    # Mechanism of action per target (ChEMBL): gene -> action_type / mechanism text
    mech = {}
    try:
        for m in _chembl_mechanism_detail(chembl):
            mech[m["gene"].upper()] = m
    except Exception:
        pass

    # ── Targets as PROTEIN nodes: mechanism-of-action edge + the pathways they sit in ──
    gene_ids, gene_pw = {}, {}
    for g in targets:
        gU, gid = g.upper(), f"gene_{g}"
        gene_ids[gU] = gid
        m = mech.get(gU)
        action = _action_to_edge_label((m or {}).get("action_type", ""))
        add_node(_net_node(gid, g, "Protein", source_db="ChEMBL",
                           mechanism=(m or {}).get("mechanism", ""),
                           action_type=(m or {}).get("action_type", "")))
        elements.append(_net_edge("cmp_center", gid, action, 1.0))
        pw_ids = []
        for (pname, pid) in _gene_pathways(g, top=2):
            pwid = f"pw_{pid}"
            add_node(_net_node(pwid, pname, "Pathway", reactome_id=pid,
                               reactome_url=f"https://reactome.org/PathwayBrowser/#/{pid}",
                               source_db="Reactome", size=42))
            elements.append(_net_edge(gid, pwid, "PARTICIPATES_IN", 0.85))
            pw_ids.append(pwid)
        gene_pw[gU] = pw_ids

    # ── Candidate diseases via their target — evidence-bearing edges ──
    # Fuse the platform's validated signals INTO the graph so it reads as a decision
    # tool, not just a topology picture: each Disease node carries its Repurposing
    # Value Score (burden × unmet-need × market — "worth pursuing?") and each
    # drug→disease edge carries the validated DWPC plausibility P(treats). Fail-soft.
    try:
        from services.disease_value import value_for as _dvalue
    except Exception:
        _dvalue = None
    try:
        from services.repurposing_predictor import plausibility as _plaus
    except Exception:
        _plaus = None

    gene_dis = {}
    for c in cands:
        dis = c.get("disease")
        if not dis:
            continue
        did = f"dis_{c.get('efo_id') or dis}"
        sc  = float(c.get("score") or c.get("association_score") or 0.0)
        dv = _dvalue(dis, c.get("efo_id", "")) if _dvalue else None
        vscore = (dv or {}).get("value_score")
        vtier  = (dv or {}).get("tier")
        pl = _plaus(name, dis) if _plaus else None
        p_treats = (pl or {}).get("probability")
        # Node SIZE now reflects pharma value (fallback to association score); tier +
        # plausibility ride along for colour/tooltip in the front end.
        node_val = vscore if vscore is not None else min(sc, 1.0)
        add_node(_net_node(did, dis, "Disease", efo_id=c.get("efo_id"),
                           score=round(sc, 3), via_target=c.get("via_target"),
                           value_score=vscore, value_tier=vtier,
                           plausibility=p_treats,
                           therapeutic_areas=c.get("therapeutic_areas", []),
                           sources=c.get("sources", []),
                           trial_count=c.get("trial_count", 0),
                           max_trial_phase=c.get("max_trial_phase", 0),
                           source_db="Open Targets",
                           size=32 + node_val * 34))
        via = (c.get("via_target") or "").upper()
        if via and via in gene_ids:
            elements.append({"data": {
                "source": gene_ids[via], "target": did, "label": "ASSOCIATED_WITH",
                "weight": sc, "evidence_score": round(float(c.get("evidence_score", 0) or 0), 3),
                "plausibility": p_treats, "value_score": vscore,
                "sources": c.get("sources", []), "trial_count": c.get("trial_count", 0),
                "max_trial_phase": c.get("max_trial_phase", 0), "source_db": "Open Targets"}})
            gene_dis.setdefault(via, []).append((did, sc))
        else:
            e = _net_edge("cmp_center", did, "REPURPOSE_CANDIDATE", sc)
            e["data"]["plausibility"] = p_treats
            e["data"]["value_score"] = vscore
            elements.append(e)

    # ── Pathway → Disease (mechanistic bridge): the target's pathway is a candidate
    #    mechanism for the diseases that target associates with (top few, for clarity) ──
    for gU, pw_ids in gene_pw.items():
        top_dis = sorted(gene_dis.get(gU, []), key=lambda t: -t[1])[:5]
        for pwid in pw_ids:
            for (did, sc) in top_dis:
                elements.append(_net_edge(pwid, did, "IMPLICATED_IN", round(sc * 0.8, 3)))

    return elements, legend


def _chembl_molecules(chembl_ids: List[str]) -> Dict[str, Dict]:
    """Batch-fetch molecule names + phases."""
    results = {}
    for i in range(0, min(len(chembl_ids), 40), 20):
        chunk = chembl_ids[i:i + 20]
        try:
            r = http_client.get(f"{CHEMBL_BASE}/molecule.json",
                params={"molecule_chembl_id__in": ",".join(chunk), "limit": 20, "format": "json"},
                timeout=10)
            if r and r.ok:
                for mol in r.json().get("molecules", []):
                    mid = mol.get("molecule_chembl_id", "")
                    results[mid] = {
                        "name":      mol.get("pref_name") or "",
                        "max_phase": mol.get("max_phase") or 0,
                    }
        except Exception:
            pass
    return results


def _get_chembl_target(tid: str) -> Optional[dict]:
    """Cached fetch of a ChEMBL target record."""
    if tid not in _CHEMBL_TARGET_CACHE:
        try:
            r = http_client.get(f"{CHEMBL_BASE}/target/{tid}.json", timeout=6)
            _CHEMBL_TARGET_CACHE[tid] = r.json() if (r and r.ok) else {}
        except Exception:
            _CHEMBL_TARGET_CACHE[tid] = {}
    return _CHEMBL_TARGET_CACHE[tid] or None


def _action_to_edge_label(action: str) -> str:
    """Map ChEMBL action_type to a graph edge label."""
    a = (action or "").upper().replace(" ", "_")
    if a in {"INHIBITOR", "NEGATIVE_ALLOSTERIC_MODULATOR", "ANTAGONIST",
             "BLOCKER", "INVERSE_AGONIST", "CHANNEL_BLOCKER", "NEGATIVE_MODULATOR"}:
        return "INHIBITS"
    if a in {"ACTIVATOR", "AGONIST", "POSITIVE_ALLOSTERIC_MODULATOR",
             "PARTIAL_AGONIST", "FULL_AGONIST", "OPENER"}:
        return "ACTIVATES"
    return "BINDS_TO"


def _chembl_mechanism_detail(chembl_id: str) -> List[Dict]:
    """Return [{gene, action_type, mechanism, target_id}] from ChEMBL mechanism API."""
    results: List[Dict] = []
    seen: set = set()
    try:
        r = http_client.get(f"{CHEMBL_BASE}/mechanism.json",
            params={"molecule_chembl_id": chembl_id, "limit": 20, "format": "json"}, timeout=8)
        if not r or not r.ok:
            return []
        for mech in r.json().get("mechanisms", []):
            tid = mech.get("target_chembl_id", "")
            if not tid:
                continue
            action    = (mech.get("action_type") or "").replace("_", " ").title()
            mech_text = mech.get("mechanism_of_action") or ""
            tdata     = _get_chembl_target(tid)
            if not tdata:
                continue
            ttype = tdata.get("target_type", "")
            if ttype and ttype not in ("SINGLE PROTEIN", "PROTEIN COMPLEX",
                                       "PROTEIN FAMILY", "SELECTIVITY GROUP"):
                continue
            for comp in (tdata.get("target_components") or []):
                gene = ""
                for syn in (comp.get("target_component_synonyms") or []):
                    if syn.get("syn_type") == "GENE_SYMBOL":
                        gene = syn.get("component_synonym", "")
                        break
                if gene and gene not in seen:
                    seen.add(gene)
                    results.append({
                        "gene":        gene,
                        "action_type": action,
                        "mechanism":   mech_text[:120],
                        "target_id":   tid,
                    })
    except Exception:
        pass
    return results


def _chembl_drug_genes(chembl_id: str) -> List[str]:
    """Gene targets for a single compound via ChEMBL mechanism (simple list form)."""
    return [m["gene"] for m in _chembl_mechanism_detail(chembl_id)]


def _string_ppi(gene_symbols: List[str], species: int = 9606) -> Dict[str, List[str]]:
    """STRING protein-protein interaction adjacency."""
    if not gene_symbols:
        return {}
    try:
        r = http_client.get(f"{STRING_URL}/network",
            params={
                "identifiers":     "%0d".join(gene_symbols[:20]),
                "species":         species,
                "required_score":  700,
                "caller_identity": "neurorepurpose_platform",
            }, timeout=10)
        if r and r.ok:
            adj: Dict[str, List[str]] = {}
            for link in r.json():
                a = link.get("preferredName_A", "")
                b = link.get("preferredName_B", "")
                if a and b:
                    adj.setdefault(a, []).append(b)
                    adj.setdefault(b, []).append(a)
            return adj
    except Exception:
        pass
    return {}


# ── Legacy stub (kept for import compatibility) ────────────────────────────────

def build_graph(mesh_ids: List[str] = None, compound_id: int = None,
                max_compounds: int = 20, max_targets: int = 30) -> Tuple[List[dict], dict]:
    disease = (mesh_ids or [""])[0] if mesh_ids else ""
    return build_network_graph(disease=disease, limit=max_compounds + max_targets)


# ── Main graph builder ─────────────────────────────────────────────────────────

def build_network_graph(disease: str = "", limit: int = 80) -> Tuple[List[dict], dict]:
    """
    Build a Disease → Gene → Compound evidence network.
    Returns (cytoscape_elements, legend_dict).
    """
    elements:   List[dict] = []
    seen_nodes: set        = set()
    seen_edges: set        = set()

    def add_node(n: dict):
        nid = n["data"]["id"]
        if nid not in seen_nodes:
            seen_nodes.add(nid)
            elements.append(n)

    def add_edge(src: str, tgt: str, lbl: str, weight: float = 1.0):
        key = (src, tgt, lbl)
        if key not in seen_edges and src in seen_nodes and tgt in seen_nodes:
            seen_edges.add(key)
            elements.append(_net_edge(src, tgt, lbl, weight))

    if not disease:
        _build_overview_graph(add_node, add_edge)
    else:
        _build_indication_graph(disease, add_node, add_edge, elements, seen_nodes, limit)

    kinds_present = {n["data"]["kind"] for n in elements if "kind" in n.get("data", {})}
    legend = {k: v for k, v in _COLORS.items() if k in kinds_present}
    return elements, legend


def _build_overview_graph(add_node, add_edge):
    """
    Show neuro disease overview. Uses DRKG when cached (much richer),
    falls back to Open Targets live API.
    """
    drkg = _get_drkg()

    # ── DRKG path (pre-cached multi-source data) ──────────────────────────────
    if drkg.get("diseases"):
        dis_map    = drkg["diseases"]          # entity_id → {name, mesh_id}
        cmp_map    = drkg["compounds"]         # entity_id → {name, ...}
        gene_map   = drkg["genes"]             # entity_id → {symbol, ...}
        edges      = drkg.get("edges", {})
        treat_e    = edges.get("treat", [])
        dg_e       = edges.get("dis_gene", [])

        # Build per-disease compound/gene counts for picking top nodes
        dis_cmp_count: Dict[str, int] = {}
        dis_cmp_map:   Dict[str, List[str]] = {}
        for e in treat_e:
            t = e["t"]
            dis_cmp_count[t] = dis_cmp_count.get(t, 0) + 1
            dis_cmp_map.setdefault(t, []).append(e["h"])

        dis_gene_map: Dict[str, List[str]] = {}
        for e in dg_e:
            dis_gene_map.setdefault(e["h"], []).append(e["t"])

        # Top 4 diseases (by compound count)
        top_dis = sorted(dis_map.items(),
                         key=lambda x: dis_cmp_count.get(x[0], 0), reverse=True)[:4]

        for did_ent, dis_info in top_dis:
            mesh_id = dis_info.get("id_value", dis_info.get("entity_id", did_ent))
            did = f"dis_m{mesh_id}"
            add_node(_net_node(did, dis_info["name"], "Disease",
                               mesh_id=mesh_id))

            # Top 5 genes per disease
            for gent in (dis_gene_map.get(did_ent) or [])[:5]:
                gi = gene_map.get(gent, {})
                sym = gi.get("symbol", _eid_raw(gent))
                if not sym or sym.startswith("0"):
                    continue
                gid = f"gene_{sym}"
                add_node(_net_node(gid, sym, "Gene",
                                   gene_name=gi.get("gene_name", "")))
                add_edge(gid, did, "associated_with")

            # Top 5 compounds per disease (highest max_phase first)
            cmp_list = dis_cmp_map.get(did_ent, [])
            cmp_list_info = [
                (cent, cmp_map.get(cent, {}))
                for cent in cmp_list
                if cmp_map.get(cent, {}).get("name")
            ]
            cmp_list_info.sort(key=lambda x: x[1].get("max_phase", 0), reverse=True)
            for cent, ci in cmp_list_info[:5]:
                name = ci.get("name", "")
                if not name:
                    continue
                cid = f"cmp_{ci.get('chembl_id') or ci.get('drugbank_id') or cent}"
                add_node(_net_node(cid, name.title(), "Compound",
                                   chembl_id=ci.get("chembl_id", ""),
                                   phase=ci.get("max_phase", 0)))
                lbl = "treats" if ci.get("max_phase", 0) >= 4 else "repurpose_candidate"
                add_edge(cid, did, lbl)
        return

    # ── Fallback: Open Targets live API ──────────────────────────────────────
    try:
        from services.disease_ontology import resolve_disease
    except Exception:
        for dis_name, efo_id in _NEURO_OVERVIEW:
            add_node(_net_node(f"dis_{efo_id}", dis_name, "Disease", efo_id=efo_id))
        return

    for dis_name, efo_id in _NEURO_OVERVIEW:
        did = f"dis_{efo_id}"
        add_node(_net_node(did, dis_name, "Disease", efo_id=efo_id))
        try:
            info = resolve_disease(dis_name)
            for t in info.get("targets", [])[:5]:
                gsym = t["gene_symbol"]
                gid  = f"gene_{gsym}"
                add_node(_net_node(gid, gsym, "Gene",
                                   ot_score=round(t.get("score", 0), 3)))
                add_edge(gid, did, "associated_with", weight=t.get("score", 0.3))
        except Exception:
            pass


def _eid_raw(entity: str) -> str:
    return entity.split("::")[-1] if "::" in entity else entity


def _build_indication_graph(disease: str, add_node, add_edge, elements: list,
                             seen_nodes: set, limit: int):
    """
    Full disease-specific graph.
    Layer 1: Open Targets (disease→gene, authoritative)
    Layer 2: DRKG cache (additional compound-disease + compound-gene from 6 databases)
    Layer 3: ChEMBL live API fallback (when DRKG unavailable)
    Layer 4: STRING PPI (gene-gene interactions)
    """
    # ── 1. Resolve disease via Open Targets ───────────────────────────────────
    disease_info = {}
    try:
        from services.disease_ontology import resolve_disease
        disease_info = resolve_disease(disease)
    except Exception as e:
        logger.debug(f"Disease resolve: {e}")

    dis_name = disease_info.get("disease_name") or disease.title()
    efo_id   = disease_info.get("disease_id", "")
    desc     = (disease_info.get("description") or "")[:150]
    dis_id   = "dis_center"

    add_node(_net_node(dis_id, dis_name, "Disease",
                       efo_id=efo_id, description=desc))

    # ── 2. Gene nodes from Open Targets (primary evidence) ───────────────────
    gene_syms: List[str] = []
    gene_id_map: Dict[str, str] = {}   # UPPER(symbol) → node_id

    for t in disease_info.get("targets", [])[:22]:
        gsym  = t["gene_symbol"]
        gid   = f"gene_{gsym}"
        score = t.get("score", 0)
        add_node(_net_node(gid, gsym, "Gene",
                           gene_name=t.get("gene_name", ""),
                           ot_score=round(score, 3)))
        gene_id_map[gsym.upper()] = gid
        gene_syms.append(gsym)
        add_edge(gid, dis_id, "associated_with", weight=max(0.2, score))

    # ── 3a. DRKG compound candidates ─────────────────────────────────────────
    added_cmp_ids: List[str] = []   # ChEMBL IDs of added compound nodes
    drkg = _get_drkg()

    if drkg.get("diseases"):
        # Find DRKG disease entity matching this query
        dis_ent = _match_drkg_disease(disease, drkg["diseases"])

        if dis_ent:
            cmp_map  = drkg["compounds"]
            gene_map = drkg["genes"]
            edges    = drkg.get("edges", {})

            # Compound-disease treat edges for this disease
            treat_e = [e for e in edges.get("treat", []) if e["t"] == dis_ent]
            # Sort by max_phase descending
            treat_e.sort(key=lambda e: cmp_map.get(e["h"], {}).get("max_phase", 0),
                         reverse=True)

            drkg_cmp_added = 0
            for e in treat_e[:min(limit // 4, 20)]:
                cent = e["h"]
                ci   = cmp_map.get(cent, {})
                name = (ci.get("name") or "").strip()
                if not name:
                    continue
                phase = ci.get("max_phase", 0)
                try:
                    phase = int(float(phase or 0))
                except (TypeError, ValueError):
                    phase = 0
                chembl_id = ci.get("chembl_id", "")
                cid = f"cmp_{chembl_id or ci.get('drugbank_id') or cent}"
                add_node(_net_node(cid, name.title(), "Compound",
                                   chembl_id=chembl_id,
                                   phase=phase,
                                   db_source="DRKG"))
                lbl = "treats" if phase >= 4 else "repurpose_candidate"
                add_edge(cid, dis_id, lbl, weight=max(0.2, phase / 4.0))
                if chembl_id:
                    added_cmp_ids.append(chembl_id)
                drkg_cmp_added += 1

            # Additional disease-gene edges from DRKG (genes not in OT)
            dg_e = [e for e in edges.get("dis_gene", []) if e["h"] == dis_ent]
            for e in dg_e[:40]:
                gent = e["t"]
                gi   = gene_map.get(gent, {})
                sym  = gi.get("symbol", "")
                if not sym or sym.isdigit() or sym.upper() in gene_id_map:
                    continue
                gid = f"gene_{sym}"
                add_node(_net_node(gid, sym, "Gene",
                                   gene_name=gi.get("gene_name", ""),
                                   ot_score=0, db_source="DRKG"))
                gene_id_map[sym.upper()] = gid
                gene_syms.append(sym)
                add_edge(gid, dis_id, "associated_with", weight=0.25)

            # DRKG compound → gene target edges
            target_e = [e for e in edges.get("target", [])
                        if e["h"] in {e2["h"] for e2 in treat_e[:20]}]
            for e in target_e[:80]:
                gent = e["t"]
                gi   = gene_map.get(gent, {})
                sym  = gi.get("symbol", "")
                if not sym:
                    continue
                gid = f"gene_{sym}"
                cent = e["h"]
                ci   = cmp_map.get(cent, {})
                chembl_id = ci.get("chembl_id", "")
                cid = f"cmp_{chembl_id or ci.get('drugbank_id') or cent}"
                if cid not in seen_nodes:
                    continue
                if gid in seen_nodes:
                    add_edge(cid, gid, "targets", weight=0.8)
                elif sym.upper() in gene_id_map:
                    add_edge(cid, gene_id_map[sym.upper()], "targets", weight=0.8)

            logger.debug(f"DRKG: added {drkg_cmp_added} compounds for {dis_name}")

    # ── 3b. ChEMBL live API fallback (when DRKG has no data for this disease) ──
    if len(added_cmp_ids) < 5:
        max_cmp   = min(limit // 4, 14)
        inds      = _chembl_indications(disease, limit=max_cmp + 6)
        mol_ids   = [i["chembl_id"] for i in inds]
        mol_map   = _chembl_molecules(mol_ids[:30])
        ind_phase = {i["chembl_id"]: i["max_phase"] for i in inds}

        for ind in inds[:max_cmp]:
            mid  = ind["chembl_id"]
            mol  = mol_map.get(mid, {})
            name = mol.get("name", "").strip()
            if not name or f"cmp_{mid}" in seen_nodes:
                continue
            try:
                phase = int(float(mol.get("max_phase") or ind_phase.get(mid) or 0))
            except (TypeError, ValueError):
                phase = 0
            cid = f"cmp_{mid}"
            add_node(_net_node(cid, name, "Compound",
                               chembl_id=mid, phase=phase,
                               indication=ind.get("indication", disease),
                               db_source="ChEMBL"))
            lbl = "treats" if phase >= 4 else "repurpose_candidate"
            add_edge(cid, dis_id, lbl, weight=max(0.2, phase / 4.0))
            added_cmp_ids.append(mid)

        # ChEMBL mechanism → gene edges (top 8)
        for mid in added_cmp_ids[:8]:
            cid   = f"cmp_{mid}"
            for gsym in _chembl_drug_genes(mid):
                gid = f"gene_{gsym}"
                if gid in seen_nodes:
                    add_edge(cid, gid, "targets", weight=0.8)
                elif gsym.upper() in gene_id_map:
                    add_edge(cid, gene_id_map[gsym.upper()], "targets", weight=0.8)
                elif len(seen_nodes) < limit + 10:
                    add_node(_net_node(gid, gsym, "Gene", ot_score=0,
                                       db_source="ChEMBL", via_drug=True))
                    gene_id_map[gsym.upper()] = gid
                    add_edge(gid, dis_id, "associated_with", weight=0.15)
                    add_edge(cid, gid, "targets", weight=0.8)

    # ── 4. STRING PPI edges ───────────────────────────────────────────────────
    if gene_syms:
        try:
            from services.disease_ontology import get_ppi_network
            ppi = get_ppi_network(gene_syms[:15])
            for g1, neighbors in ppi.items():
                g1id = f"gene_{g1}"
                for g2 in neighbors:
                    g2id = f"gene_{g2}"
                    if g1id in seen_nodes and g2id in seen_nodes:
                        add_edge(g1id, g2id, "interacts", weight=0.4)
        except Exception as e:
            logger.debug(f"STRING PPI: {e}")


def _match_drkg_disease(disease_query: str, dis_map: dict) -> Optional[str]:
    """Find the DRKG disease entity ID that best matches a query string."""
    q = disease_query.lower()
    # Keyword map: query terms → MESH IDs
    _kw = {
        "alzheimer":          "D000544",
        "parkinson":          "D010300",
        "multiple sclerosis": "D009103",
        "sclerosis":          "D009103",
        " ms ":               "D009103",
        "als":                "D000690",
        "amyotrophic":        "D000690",
        "huntington":         "D006816",
        "epilep":             "D004827",
        "schizophr":          "D012559",
        "depress":            "D003866",
        "bipolar":            "D001714",
        "adhd":               "D001289",
        "attention deficit":  "D001289",
        "autis":              "D001321",
        "migraine":           "D008881",
    }
    for kw, mid in _kw.items():
        if kw in q:
            ent = f"Disease::MESH:{mid}"
            if ent in dis_map:
                return ent
    # Fuzzy fallback: check disease names in cache
    for ent, info in dis_map.items():
        name_words = set(info.get("name", "").lower().split())
        query_words = set(q.split())
        if len(name_words & query_words) >= 1:
            return ent
    return None


# ── Repurposing story graph ────────────────────────────────────────────────────

def build_repurpose_story_graph(disease: str, top_n: int = 3,
                                focal_compound: str = None) -> Tuple[List[dict], dict]:
    """
    Evidence-layered repurposing story graph with four biological layers.
    Every node connects back to Disease through at least one path:

    Layer C (Clinical):   Drug ──[TREATS / REPURPOSE_CANDIDATE]──► Disease
    Layer A (Molecular):  Gene ──[ASSOCIATED_WITH, OT score]──► Disease
                          Gene ──[ENCODES]──► Protein
                          Drug ──[BINDS_TO / INHIBITS / ACTIVATES]──► Protein
    Layer B (Pathway):    Protein ──[PARTICIPATES_IN]──► Pathway
                          Pathway ──[IMPLICATED_IN]──► Disease
    """
    elements:   List[dict] = []
    seen_nodes: set = set()
    seen_edges: set = set()
    gene_to_pw: Dict[str, list] = {}   # UPPER(gene) → [(pathway_name, reactome_id)] (dynamic, from OT/Reactome)
    pw_ids_added: set = set()
    PW_CAP = 12                         # max distinct pathway nodes in the graph

    def add_node(n: dict):
        nid = n["data"]["id"]
        if nid not in seen_nodes:
            seen_nodes.add(nid)
            elements.append(n)

    def add_edge(src: str, tgt: str, lbl: str, weight: float = 1.0, **extra):
        key = (src, tgt, lbl)
        if key not in seen_edges and src in seen_nodes and tgt in seen_nodes:
            seen_edges.add(key)
            elements.append({"data": {"source": src, "target": tgt,
                                      "label": lbl, "weight": weight, **extra}})

    def _attach_protein(gsym: str, gene_name: str = "", ot_score: float = 0) -> Optional[str]:
        """Add a Protein node + Gene→Protein ENCODES edge for ANY gene (dynamic).
        Protein name comes from Open Targets' approvedName, else a curated label."""
        gsym_upper = gsym.upper()
        # Prefer a clean curated name, else the OT gene/protein name, else a generic label
        prot_name = _PROTEIN_MAP.get(gsym_upper) or gene_name or f"{gsym} protein"
        gid = f"gene_{gsym}"
        prot_id = f"prot_{gsym}"
        if prot_id not in seen_nodes:
            add_node(_net_node(prot_id, prot_name, "Protein",
                               gene_symbol=gsym, gene_name=gene_name,
                               ot_score=ot_score, source_db="UniProt"))
        if gid in seen_nodes:
            add_edge(gid, prot_id, "ENCODES",
                     weight=0.95, source_db="UniProt",
                     evidence_score=0.95, evidence_type="Gene encodes protein")
        return prot_id

    def _attach_pathway(gsym_upper: str, pathway_src_id: str, dis_id: str,
                        max_per_gene: int = 2) -> None:
        """Add Pathway node(s) + PARTICIPATES_IN + IMPLICATED_IN edges for a gene,
        using the dynamic gene→pathway index (real Reactome pathways from Open Targets),
        falling back to the curated map only when OT has none."""
        pws = gene_to_pw.get(gsym_upper) or (
            [_PATHWAY_MAP[gsym_upper]] if gsym_upper in _PATHWAY_MAP else [])
        for pname, pid in pws[:max_per_gene]:
            if not pid:
                continue
            pw_id = f"pw_{pid.replace('-', '_')}"
            if pw_id not in pw_ids_added:
                if len(pw_ids_added) >= PW_CAP:
                    break
                add_node(_net_node(pw_id, pname, "Pathway",
                                   reactome_id=pid,
                                   reactome_url=f"https://reactome.org/PathwayBrowser/#/{pid}",
                                   source_db="Reactome", size=44))
                pw_ids_added.add(pw_id)
                add_edge(pw_id, dis_id, "IMPLICATED_IN",
                         weight=0.7, source_db="Reactome",
                         evidence_score=0.7,
                         evidence_type="Pathway implicated in disease",
                         pathway_id=pid)
            add_edge(pathway_src_id, pw_id, "PARTICIPATES_IN",
                     weight=0.85, source_db="Reactome", evidence_score=0.85,
                     evidence_type="Protein participates in pathway", pathway_id=pid)

    # ── 1. Run repurposing screen ──────────────────────────────────────────────
    try:
        from services.repurposing_engine import run_repurposing_screen
        result = run_repurposing_screen(disease, max_candidates=25)
    except Exception as e:
        logger.debug(f"Repurpose story graph — screen error: {e}")
        return [], {}

    candidates         = result.get("candidates", [])
    disease_genes_list = result.get("top_disease_genes", [])
    disease_genes      = {g.upper() for g in disease_genes_list}
    disease_info       = result.get("disease_info", {})

    # Dynamic gene→pathway index from Open Targets / Reactome (replaces the hardcoded
    # neuro-only _PATHWAY_MAP, so pathways populate for ANY disease/target).
    for pw in disease_info.get("pathways", []):
        pname = pw.get("pathway_name", "")
        pid   = pw.get("pathway_id", "")
        for g in pw.get("genes", []):
            gene_to_pw.setdefault(g.upper(), []).append((pname, pid))

    top_drugs = candidates[:top_n]

    # Ensure the focal compound is in the graph even when it is NOT among the
    # disease's own forward candidates. Otherwise, selecting a drug the platform
    # does not already rank for this disease (a genuine repurposing hypothesis)
    # would render a graph showing only OTHER drugs — never the one the user chose.
    if focal_compound and not any(c.get("chembl_id") == focal_compound for c in top_drugs):
        focal_cand = next((c for c in candidates if c.get("chembl_id") == focal_compound), None)
        if focal_cand is None:
            # Build the focal candidate from its own ChEMBL record so the graph is
            # ABOUT the selected drug; its mechanism edges (below) connect it to the
            # disease genes, or — honestly — show no bridge when there is no shared
            # target (a weak hypothesis reads as a drug with no molecular link).
            try:
                from services.reverse_repurposing import resolve_drug, canonical_pair_score
                info = resolve_drug(focal_compound)
                if info.get("chembl_id"):
                    try:
                        _sr = canonical_pair_score(info["chembl_id"], disease,
                                                   drug_genes=info.get("targets"))
                        _score = _sr.get("composite_score", 0)
                    except Exception:
                        _score = 0
                    focal_cand = {
                        "name":        info.get("name", ""),
                        "chembl_id":   info["chembl_id"],
                        "max_phase":   info.get("max_phase", 0),
                        "indications": "; ".join(k["name"] for k in info.get("known_indications", [])),
                        "targets":     ";".join(info.get("targets", [])),
                        "score":       _score,
                    }
            except Exception as e:
                logger.debug(f"focal compound build failed: {e}")
        if focal_cand:
            top_drugs = [focal_cand] + top_drugs[:top_n]

    if not top_drugs:
        return [], {}

    dis_name = disease_info.get("disease_name") or disease.title()
    efo_id   = disease_info.get("disease_id", "")
    dis_id   = "dis_center"

    # ── 2. Disease node ────────────────────────────────────────────────────────
    add_node(_net_node(dis_id, dis_name, "Disease",
                       efo_id=efo_id,
                       source_db="Open Targets",
                       description=(disease_info.get("description") or "")[:150]))

    # ── 3. Layer A/B: OT disease genes → Proteins → Pathways → Disease ───────
    ot_gene_scores: Dict[str, float] = {}
    ot_gene_names:  Dict[str, str]   = {}

    for t in disease_info.get("targets", [])[:8]:
        gsym       = t["gene_symbol"]
        gsym_upper = gsym.upper()
        ot_score   = round(t.get("score", 0), 3)
        ot_gene_scores[gsym_upper] = ot_score
        ot_gene_names[gsym_upper]  = t.get("gene_name", "")

        gid = f"gene_{gsym}"
        add_node(_net_node(gid, gsym, "Gene",
                           gene_name=t.get("gene_name", ""),
                           ot_score=ot_score,
                           bridge=False,
                           source_db="Open Targets"))

        # Gene → Disease (genetic association evidence)
        add_edge(gid, dis_id, "ASSOCIATED_WITH",
                 weight=max(0.2, ot_score),
                 source_db="Open Targets",
                 evidence_score=ot_score,
                 evidence_type="Genetic association",
                 source_url=(f"https://platform.opentargets.org/disease/{efo_id}/associations"
                             if efo_id else ""))

        # Gene → Protein (ENCODES)
        prot_id = _attach_protein(gsym, t.get("gene_name", ""), ot_score)

        # Protein (or Gene) → Pathway → Disease
        path_src = prot_id if prot_id else gid
        _attach_pathway(gsym_upper, path_src, dis_id)

    # ── 4. Layer A/C: Drug nodes + mechanism edges ────────────────────────────
    for rank, cand in enumerate(top_drugs):
        name = (cand.get("name") or "").strip()
        if not name:
            continue

        chembl_id    = cand.get("chembl_id", "")
        score        = round(float(cand.get("score") or cand.get("composite_score") or 0), 3)
        phase        = int(float(cand.get("max_phase") or 0))
        approved_for = (
            (cand.get("indications") or "").split(";")[0].strip()
            or (cand.get("mechanisms") or "").split(";")[0].strip()
        )

        cid = f"cmp_{chembl_id or name.lower().replace(' ', '_')}_{rank}"

        is_focal = bool(focal_compound and chembl_id and chembl_id == focal_compound)
        add_node(_net_node(cid, name.title(), "Compound",
                           chembl_id=chembl_id,
                           phase=phase,
                           score=score,
                           rank=rank + 1,
                           approved_for=approved_for,
                           source_db="ChEMBL",
                           focal=is_focal,
                           size=72 if is_focal else max(36, 60 - rank * 6)))

        # Layer C: Drug → Disease (clinical evidence)
        c_lbl = "TREATS" if phase >= 4 else "REPURPOSE_CANDIDATE"
        add_edge(cid, dis_id, c_lbl,
                 weight=max(0.3, score),
                 source_db="ChEMBL",
                 evidence_score=score,
                 max_phase=phase,
                 approved_for=approved_for,
                 source_id=chembl_id)

        # Layer A: ChEMBL mechanism → gene/protein targets
        mech_list = _chembl_mechanism_detail(chembl_id) if chembl_id else []

        if not mech_list:
            targets_raw = cand.get("targets", "") or cand.get("gene_targets", "") or ""
            for g in re.split(r"[;,/]", targets_raw):
                g = g.strip()
                if g:
                    mech_list.append({"gene": g, "action_type": "", "mechanism": "", "target_id": ""})

        drug_gene_map: Dict[str, Dict] = {}
        for m in mech_list:
            drug_gene_map[m["gene"].upper()] = m

        bridge = [k for k in drug_gene_map if k in disease_genes]
        others = [k for k in drug_gene_map if k not in disease_genes][:2]

        # Bridge genes — the molecular explanation for WHY this drug fits this disease
        for gsym_upper in bridge[:5]:
            m_info    = drug_gene_map[gsym_upper]
            gsym      = m_info["gene"]
            action    = m_info.get("action_type", "")
            mech_text = m_info.get("mechanism", "")
            gid       = f"gene_{gsym}"

            if gid not in seen_nodes:
                add_node(_net_node(gid, gsym, "Gene",
                                   bridge=True,
                                   ot_score=ot_gene_scores.get(gsym_upper, 0),
                                   gene_name=ot_gene_names.get(gsym_upper, ""),
                                   source_db="ChEMBL",
                                   note=f"Links {name.title()} to {dis_name}"))
                add_edge(gid, dis_id, "ASSOCIATED_WITH",
                         weight=max(0.3, ot_gene_scores.get(gsym_upper, 0.5)),
                         source_db="Open Targets",
                         evidence_score=ot_gene_scores.get(gsym_upper, 0.5),
                         evidence_type="Bridge gene")

                # Gene → Protein (ENCODES)
                prot_id = _attach_protein(gsym,
                                          ot_gene_names.get(gsym_upper, ""),
                                          ot_gene_scores.get(gsym_upper, 0))

                # Protein (or Gene) → Pathway → Disease
                path_src = prot_id if prot_id else gid
                _attach_pathway(gsym_upper, path_src, dis_id)
            else:
                # Mark existing Gene node as bridge
                for elem in elements:
                    if elem.get("data", {}).get("id") == gid:
                        elem["data"]["bridge"] = True
                        if not elem["data"].get("note"):
                            elem["data"]["note"] = f"Links {name.title()} to {dis_name}"
                        break
                # Ensure Protein exists for this bridge gene
                prot_id = _attach_protein(gsym,
                                          ot_gene_names.get(gsym_upper, ""),
                                          ot_gene_scores.get(gsym_upper, 0))

            # Drug → Protein (if protein exists) or Drug → Gene (fallback)
            prot_id_check = f"prot_{gsym}" if gsym_upper in _PROTEIN_MAP else None
            target_id = prot_id_check if (prot_id_check and prot_id_check in seen_nodes) else gid
            edge_lbl = _action_to_edge_label(action)
            add_edge(cid, target_id, edge_lbl,
                     weight=0.9,
                     source_db="ChEMBL",
                     evidence_score=0.92,
                     mechanism=action or "Direct interaction",
                     detail=mech_text[:100] if mech_text else "",
                     source_id=chembl_id)

        # Non-bridge targets — mechanism context
        for gsym_upper in others:
            m_info = drug_gene_map[gsym_upper]
            gsym   = m_info["gene"]
            action = m_info.get("action_type", "")
            gid    = f"gene_{gsym}"
            if gid not in seen_nodes:
                add_node(_net_node(gid, gsym, "Gene", bridge=False, source_db="ChEMBL"))
            prot_id = _attach_protein(gsym)
            target_id = f"prot_{gsym}" if (gsym.upper() in _PROTEIN_MAP and f"prot_{gsym}" in seen_nodes) else gid
            edge_lbl = _action_to_edge_label(action)
            add_edge(cid, target_id, edge_lbl,
                     weight=0.55,
                     source_db="ChEMBL",
                     evidence_score=0.55,
                     mechanism=action or "Direct interaction")

    kinds_present = {n["data"]["kind"] for n in elements if "kind" in n.get("data", {})}
    legend = {k: v for k, v in _COLORS.items() if k in kinds_present}
    return elements, legend


# ── Cytoscape stylesheet (legacy — kept for dash_app.py) ─────────────────────

CYTO_STYLESHEET = [
    {
        "selector": "node",
        "style": {
            "background-color": "data(color)",
            "border-width": 2,
            "border-color": "rgba(255,255,255,0.25)",
            "label": "data(label)",
            "color": "#ffffff",
            "font-size": "10px",
            "text-valign": "center",
            "text-halign": "center",
            "text-wrap": "wrap",
            "text-max-width": "72px",
            "width": "data(size)",
            "height": "data(size)",
        },
    },
    {
        "selector": "node[kind='Disease']",
        "style": {"font-size": "12px", "font-weight": "bold",
                  "border-width": 3, "border-color": "rgba(255,255,255,0.5)"},
    },
    {
        "selector": "node[kind='Compound']",
        "style": {"shape": "round-rectangle"},
    },
    {
        "selector": "node:selected",
        "style": {"border-width": 4, "border-color": "#ffffff"},
    },
    {
        "selector": "edge",
        "style": {
            "curve-style": "bezier",
            "line-color": "#cbd5e1",
            "target-arrow-color": "#94a3b8",
            "target-arrow-shape": "triangle",
            "arrow-scale": 1.0,
            "width": 1.5,
            "opacity": 0.6,
        },
    },
    {
        "selector": "edge[label='treats']",
        "style": {"line-color": "#e11d48", "target-arrow-color": "#e11d48",
                  "width": 2.5, "opacity": 0.85},
    },
    {
        "selector": "edge[label='repurpose_candidate']",
        "style": {"line-color": "#630ed4", "target-arrow-color": "#630ed4",
                  "width": 1.8, "opacity": 0.7, "line-style": "dashed"},
    },
    {
        "selector": "edge[label='targets']",
        "style": {"line-color": "#059669", "target-arrow-color": "#059669",
                  "width": 1.5, "opacity": 0.55},
    },
    {
        "selector": "edge[label='interacts']",
        "style": {"line-color": "#d97706", "target-arrow-color": "#d97706",
                  "width": 1.0, "opacity": 0.35, "line-style": "dotted"},
    },
]
