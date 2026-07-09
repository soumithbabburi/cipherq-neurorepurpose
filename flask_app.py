#!/usr/bin/env python3
"""RepurposeIQ Intelligence Platform — Flask Application"""

import logging, os, warnings, json, time
warnings.filterwarnings("ignore")
from typing import Dict, List, Optional
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd
import requests
from flask import Flask, render_template, request, jsonify, abort, send_from_directory

logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)

# ── Service imports ──────────────────────────────────────────────────────────
try:
    from services.neuro_db_service import (
        get_compound_activities, get_compound_indications, get_compound_targets,
        get_compounds_for_disease, get_stats, is_available as _db_available,
        search_compounds as db_search,
        find_repurposing_candidates,
        get_compound_by_chembl,
        get_data_footprint,
    )
    DB_OK = _db_available()
except Exception:
    DB_OK = False
    def get_compounds_for_disease(*a, **k): return []
    def db_search(*a, **k): return []
    def get_compound_targets(*a, **k): return []
    def get_compound_activities(*a, **k): return []
    def get_compound_indications(*a, **k): return []
    def get_stats(): return {}
    def find_repurposing_candidates(*a, **k): return []
    def get_compound_by_chembl(*a, **k): return None
    def get_data_footprint(): return {}

try:
    from services.disease_resolver import expand_mesh_ids, mesh_available, resolve_disease
    MESH_OK = mesh_available() if DB_OK else False
except Exception:
    MESH_OK = False
    def resolve_disease(*a, **k): return []
    def expand_mesh_ids(*a, **k): return []

try:
    from services.repurposing_scorer import (
        score_compound_list, score_compound as _score_compound,
        approved_chembls_for_disease)
    SCORER_OK = True
except Exception:
    SCORER_OK = False
    def score_compound_list(c, m):
        for x in c: x.setdefault("score", float(x.get("max_phase") or 0)/4); x.setdefault("score_breakdown", {})
        c.sort(key=lambda x: x["score"], reverse=True); return c
    def _score_compound(*a, **k): return {"overall": 0.0, "indication_score": 0.0, "target_score": 0.0, "activity_score": 0.0, "network_score": 0.0, "phase_bonus": 0.0}
    def approved_chembls_for_disease(*a, **k): return set()

try:
    from services.compound_validator import validate_and_deduplicate
except Exception:
    def validate_and_deduplicate(c, **k): return c

try:
    from services.pbpk_simulation import PBPKSimulator
    PBPK_OK = True
except Exception:
    PBPK_OK = False; PBPKSimulator = None

try:
    from services.docking_service import DockingService
    _dock_svc = DockingService()
    DOCK_OK = _dock_svc.available
    try:
        from local_diffdock import LocalDiffDock as _ldd
        LOCAL_DIFFDOCK_OK = _ldd().available
        LOCAL_DIFFDOCK_INSTRUCTIONS = _ldd.install_instructions()
    except Exception:
        LOCAL_DIFFDOCK_OK = False
        LOCAL_DIFFDOCK_INSTRUCTIONS = ""
except Exception:
    DOCK_OK = False; _dock_svc = None
    LOCAL_DIFFDOCK_OK = False; LOCAL_DIFFDOCK_INSTRUCTIONS = ""

try:
    from quantum_optimization_strategies import run_quantum_optimization
    QUANTUM_OK = True
except Exception:
    QUANTUM_OK = False

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors
    RDKIT_OK = True
except Exception:
    RDKIT_OK = False

try:
    import py3Dmol
    PY3DMOL_OK = True
except Exception:
    PY3DMOL_OK = False

try:
    from services.biocypher_graph import build_graph, build_network_graph
    GRAPH_OK = True
except Exception:
    GRAPH_OK = False
    def build_graph(*a, **k): return [], {}
    def build_network_graph(*a, **k): return [], {}

try:
    from services.repurposing_engine import run_repurposing_screen, score_compound_for_disease
    ENGINE_OK = True
except Exception:
    ENGINE_OK = False
    def run_repurposing_screen(*a, **k): return {"candidates": [], "error": "Engine unavailable"}

try:
    from services.evidence_dossier import generate_dossier
    DOSSIER_OK = True
except Exception:
    DOSSIER_OK = False
    def generate_dossier(*a, **k): return {"error": "Dossier service unavailable"}

try:
    from services.reverse_repurposing import screen_indications_for_drug, resolve_drug, canonical_pair_score
    REVERSE_OK = True
except Exception:
    REVERSE_OK = False
    def screen_indications_for_drug(*a, **k): return {"candidates": [], "error": "Reverse engine unavailable"}
    def resolve_drug(*a, **k): return {}
    def canonical_pair_score(*a, **k): return {}

try:
    from services import portfolio as _portfolio
    PORTFOLIO_OK = _portfolio.is_available()
except Exception:
    PORTFOLIO_OK = False
    _portfolio = None

try:
    from services.pathway_screen import (
        screen_pathway, suggest_pathways, resolve_pathway)
    PATHWAY_OK = True
except Exception:
    PATHWAY_OK = False
    def screen_pathway(*a, **k): return {"error": "Pathway engine unavailable", "candidates": []}
    def suggest_pathways(*a, **k): return []
    def resolve_pathway(*a, **k): return {}

# ── Flask app ────────────────────────────────────────────────────────────────
app = Flask(__name__)
app.secret_key = os.environ.get("SECRET_KEY", "neurorepurpose-cipherq-2026")
app.config['TEMPLATES_AUTO_RELOAD'] = True

# Shared HTTP session — connection pooling cuts per-request overhead
_http = requests.Session()
_http.headers.update({"User-Agent": "RepurposeIQ/1.0"})

# Simple TTL cache for external API results (avoids re-fetching on every request)
_api_cache: dict = {}
_CACHE_TTL = 3600  # 1 hour

def _cache_get(key: str):
    entry = _api_cache.get(key)
    if entry and time.time() < entry[1]:
        return entry[0]
    return None

def _cache_set(key: str, value, ttl: int = _CACHE_TTL):
    _api_cache[key] = (value, time.time() + ttl)


# ── Helpers ──────────────────────────────────────────────────────────────────
def _resolve_fetch(query: str):
    resolved = resolve_disease(query)
    if not resolved: return [], [], []
    ids = [r["mesh_id"] for r in resolved if r.get("mesh_id")]
    expanded = expand_mesh_ids(ids, include_children=True) or ids
    comps = get_compounds_for_disease(expanded, limit=80)
    comps = validate_and_deduplicate(comps, require_smiles=False)

    # Use repurposing engine to augment + re-score when DB returns < 8 results.
    # Feed the engine the RESOLVED canonical MeSH heading, not the raw query —
    # otherwise aliases/abbreviations ("PD", "alzheimers") never match ChEMBL and
    # silently return zero candidates even though resolution succeeded.
    canonical_name = resolved[0].get("heading") or query
    # Always run the engine: it augments the DB pool with NOVEL HetioNet-connected
    # candidates and applies the network signal, so every disease gets discovery
    # (not just confirmation) and consistent scoring. Fast + cached.
    if ENGINE_OK:
        try:
            screen = run_repurposing_screen(canonical_name, max_candidates=50, db_compounds=comps)
            engine_cands = screen.get("candidates", [])
            if engine_cands:
                # Merge: engine candidates include DB + ChEMBL + novelty, already scored
                comps = engine_cands
        except Exception as _e:
            logger.debug(f"Repurposing engine augment failed: {_e}")

    if SCORER_OK and comps and not all(str(mid).startswith("LOCAL:") for mid in expanded):
        if not all(c.get("score") for c in comps):
            comps = score_compound_list(comps, expanded)
    comps = [c for c in comps if _is_repurposable(c)]
    # Flag/boost molecules already in the company portfolio (505(b)(2)-ready)
    if PORTFOLIO_OK and comps:
        try:
            comps = _portfolio.annotate_candidates(comps)
        except Exception as _e:
            logger.debug(f"portfolio annotate failed: {_e}")
    # Area-aware developability for each compound vs the searched disease's route
    try:
        from services import developability as _devsvc
        if _devsvc.is_available() and comps:
            t_areas = _disease_areas(query)
            for c in comps:
                smi = c.get("smiles", "")
                if smi and not c.get("developability"):
                    dvr = _devsvc.score(smi, therapeutic_areas=t_areas)
                    if dvr.get("available"):
                        c["developability"]       = dvr
                        c["developability_score"]  = dvr.get("score")
    except Exception as _e:
        logger.debug(f"developability annotate failed: {_e}")
    # Flag candidates already APPROVED for this indication (the one true novelty
    # disqualifier). Local indication data is sparse and many DB rows lack a real
    # ChEMBL id, so this checks ChEMBL drug_indication phase directly — the same
    # authoritative source the dossier verdict uses. Bounded to the candidates a
    # user actually sees, run in parallel and cached.
    try:
        top = [c for c in comps[:18] if (c.get("chembl_id") or "").startswith("CHEMBL")]
        if top:
            with ThreadPoolExecutor(max_workers=6) as pool:
                flags = list(pool.map(
                    lambda c: _approved_for_disease(c.get("chembl_id"), canonical_name), top))
            for c, f in zip(top, flags):
                c["approved_here"] = f
        for c in comps:
            c.setdefault("approved_here", False)
    except Exception as _e:
        logger.debug(f"approved-here annotate failed: {_e}")
    return resolved, expanded, comps


def _approved_for_disease(chembl_id: str, disease_name: str) -> bool:
    """True only if the molecule is APPROVED (ChEMBL max_phase_for_ind >= 4) for
    a disease matching disease_name. The single genuine 505(b)(2) novelty
    disqualifier; prior lower-phase development is intentionally not counted."""
    if not chembl_id or not chembl_id.startswith("CHEMBL"):
        return False
    key = f"apprvhere:{chembl_id}:{(disease_name or '').lower()}"
    cached = _cache_get(key)
    if cached is not None:
        return cached
    val = False
    try:
        from services.regulatory_verdict import _matches
        r = _http.get("https://www.ebi.ac.uk/chembl/api/data/drug_indication.json",
                      params={"molecule_chembl_id": chembl_id, "limit": 100, "format": "json"},
                      timeout=6)
        if r.ok:
            for ind in r.json().get("drug_indications", []):
                label = ind.get("mesh_heading") or ind.get("efo_term") or ""
                try:
                    ph = float(ind.get("max_phase_for_ind") or 0)
                except (TypeError, ValueError):
                    ph = 0.0
                if ph >= 4 and label and _matches(disease_name, label):
                    val = True
                    break
    except Exception:
        pass
    _cache_set(key, val)
    return val


def _disease_areas(disease_name: str) -> List[str]:
    """Open Targets therapeutic areas for a disease (cached) — picks the developability profile."""
    key = f"areas:{disease_name.lower().strip()}"
    cached = _cache_get(key)
    if cached is not None:
        return cached
    areas: List[str] = []
    try:
        from services.disease_ontology import resolve_disease as ot_resolve
        efo = ot_resolve(disease_name).get("disease_id", "")
        if efo:
            q = "query($id:String!){ disease(efoId:$id){ therapeuticAreas { name } } }"
            r = _http.post("https://api.platform.opentargets.org/api/v4/graphql",
                           json={"query": q, "variables": {"id": efo}}, timeout=6)
            if r.ok:
                d = (r.json().get("data") or {}).get("disease") or {}
                areas = [t.get("name", "") for t in (d.get("therapeuticAreas") or [])]
    except Exception:
        pass
    _cache_set(key, areas)
    return areas


def _safe_float(v, default=0.0):
    try: return float(v)
    except: return default


_BIOLOGIC_SUFFIXES = (
    "mab",   # monoclonal antibodies: gantenerumab, lecanemab, rituximab
    "cept",  # receptor fusions: etanercept, abatacept
    "kin",   # cytokines: aldesleukin
    "tide",  # peptides/oligonucleotides
    "ase",   # enzyme biologics: alteplase, laronidase
    "gene",  # gene therapies
    "cel",   # cell therapies: sipuleucel
    "rix",   # vaccines
    "vax",   # vaccines
)

def _is_repurposable(c: dict) -> bool:
    """True if compound is a small molecule suitable for repurposing analysis."""
    mtype = (c.get("molecule_type") or "").lower()
    if mtype and "small molecule" not in mtype:
        return False
    name = (c.get("name") or "").lower().strip()
    if any(name.endswith(s) for s in _BIOLOGIC_SUFFIXES):
        return False
    if not c.get("smiles") and float(c.get("mw") or 0) > 1500:
        return False
    return True


def _get_compound(chembl_id: str):
    """Fetch compound by chembl_id or NR:id prefix. Returns dict or None."""
    if chembl_id in _compound_cache:
        c = _compound_cache[chembl_id]
        # If SMILES was missing when first cached, try online once
        if c and not c.get("smiles") and chembl_id not in _no_smiles_ids:
            smiles = _fetch_smiles_online(chembl_id, c.get("name", ""))
            if smiles:
                c["smiles"] = smiles
            else:
                _no_smiles_ids.add(chembl_id)
        return c
    c = get_compound_by_chembl(chembl_id) if DB_OK else None
    if not c:
        rows = db_search(chembl_id, limit=5)
        c = next((r for r in rows if r.get("chembl_id") == chembl_id), None)
    # Final fallback: ChEMBL REST API — molecule + indication fetched in parallel
    if not c and chembl_id and not chembl_id.startswith("NR:"):
        cached = _cache_get(f"chembl_mol:{chembl_id}")
        if cached:
            c = cached
        else:
            def _pf(v):
                try: return float(v) if v is not None else None
                except: return None
            def _pi(v):
                try: return int(v) if v is not None else None
                except: return None

            def _fetch_mol():
                try:
                    r = _http.get(
                        f"https://www.ebi.ac.uk/chembl/api/data/molecule/{chembl_id}.json",
                        timeout=5)
                    return r.json() if r.ok else None
                except Exception:
                    return None

            def _fetch_ind():
                try:
                    ir = _http.get(
                        "https://www.ebi.ac.uk/chembl/api/data/drug_indication.json",
                        params={"molecule_chembl_id": chembl_id, "limit": 5, "format": "json"},
                        timeout=5)
                    return ir.json() if ir.ok else None
                except Exception:
                    return None

            with ThreadPoolExecutor(max_workers=2) as pool:
                fut_mol = pool.submit(_fetch_mol)
                fut_ind = pool.submit(_fetch_ind)
                mol = fut_mol.result()
                ind_data = fut_ind.result()

            if mol:
                inds = []
                if ind_data:
                    inds = [i.get("mesh_heading", "") for i in
                            ind_data.get("drug_indications", []) if i.get("mesh_heading")]
                structs = mol.get("molecule_structures") or {}
                props = mol.get("molecule_properties") or {}
                c = {
                    "chembl_id":      chembl_id,
                    "name":           mol.get("pref_name") or chembl_id,
                    "max_phase":      mol.get("max_phase") or 0,
                    "smiles":         structs.get("canonical_smiles", ""),
                    "molecule_type":  mol.get("molecule_type", ""),
                    "mw":             _pf(props.get("full_mw")),
                    "alogp":          _pf(props.get("alogp")),
                    "psa":            _pf(props.get("psa")),
                    "hba":            _pi(props.get("hba")),
                    "hbd":            _pi(props.get("hbd")),
                    "rtb":            _pi(props.get("rtb")),
                    "ro5_violations": _pi(props.get("num_ro5_violations")),
                    "qed":            _pf(props.get("qed_weighted")),
                    "indications":    "; ".join(inds),
                    "targets":        "",
                }
                _cache_set(f"chembl_mol:{chembl_id}", c)
    if c:
        # Ensure SMILES is populated — fetch from ChEMBL/PubChem if the DB had NULL
        if not c.get("smiles"):
            smiles = _fetch_smiles_online(chembl_id, c.get("name", ""))
            if smiles:
                c["smiles"] = smiles
        _compound_cache[chembl_id] = c
    return c


_compound_cache: dict = {}
_tgt_gene_cache: dict = {}
_tgt_name_cache: dict = {}
_no_smiles_ids:  set  = set()   # IDs we've already tried online and found no SMILES


def _chembl_target_info(tid: str) -> tuple:
    """Return (gene_symbol, pref_name) for a ChEMBL target ID (both cached)."""
    if not tid:
        return ("", "")
    if tid in _tgt_gene_cache:
        return (_tgt_gene_cache[tid], _tgt_name_cache.get(tid, ""))
    gene = ""
    pref = ""
    try:
        r = _http.get(
            f"https://www.ebi.ac.uk/chembl/api/data/target/{tid}.json",
            timeout=4)
        if r.ok:
            data = r.json()
            pref = data.get("pref_name", "")
            for comp in data.get("target_components", []):
                for syn in (comp.get("target_component_synonyms") or []):
                    if syn.get("syn_type") == "GENE_SYMBOL":
                        gene = syn.get("component_synonym", "")
                        break
                if not gene:
                    gene = (comp.get("gene_names") or "").split(",")[0].strip()
                if gene:
                    break
    except Exception:
        pass
    _tgt_gene_cache[tid] = gene
    _tgt_name_cache[tid] = pref
    return (gene, pref)


def _chembl_gene_for_target(tid: str) -> str:
    """Wrapper kept for activities endpoint compatibility."""
    return _chembl_target_info(tid)[0]


def _fetch_chembl_mechanism_rows(chembl_id: str) -> list:
    """ChEMBL mechanisms → rows with gene_symbol + action_type (cached). Shared by
    the targets endpoint and the dossier's signature-connectivity computation."""
    if not chembl_id or chembl_id.startswith("NR:"):
        return []
    cached = _cache_get(f"chembl_mech:{chembl_id}")
    if cached is not None:
        return cached
    rows: list = []
    try:
        r = _http.get(
            "https://www.ebi.ac.uk/chembl/api/data/mechanism.json",
            params={"molecule_chembl_id": chembl_id, "limit": 50, "format": "json"},
            timeout=5)
        if r.ok:
            mechanisms = r.json().get("mechanisms", [])
            tids = [m.get("target_chembl_id", "") for m in mechanisms]
            with ThreadPoolExecutor(max_workers=min(8, max(len(tids), 1))) as pool:
                infos = list(pool.map(_chembl_target_info, tids))
            for m, (gene, pref) in zip(mechanisms, infos):
                rows.append({
                    "gene_symbol": gene,
                    "name":        pref or m.get("mechanism_of_action", ""),
                    "mechanism":   m.get("mechanism_of_action", ""),
                    "action_type": m.get("action_type", ""),
                    "confidence":  "High" if m.get("direct_interaction") else "Medium",
                    "target_id":   m.get("target_chembl_id", ""),
                    "max_phase":   m.get("max_phase", 0),
                })
            _cache_set(f"chembl_mech:{chembl_id}", rows)
    except Exception:
        pass
    return rows


def _fetch_smiles_online(chembl_id: str, name: str = "") -> str:
    """Fetch canonical SMILES from ChEMBL, then PubChem if ChEMBL misses."""
    # 1. ChEMBL molecule structures
    if chembl_id and not chembl_id.startswith("NR:"):
        cached = _cache_get(f"smiles_raw:{chembl_id}")
        if cached:
            return cached
        try:
            r = _http.get(
                f"https://www.ebi.ac.uk/chembl/api/data/molecule/{chembl_id}.json",
                timeout=6)
            if r.ok:
                structs = r.json().get("molecule_structures") or {}
                smiles = structs.get("canonical_smiles", "")
                if smiles:
                    _cache_set(f"smiles_raw:{chembl_id}", smiles)
                    return smiles
        except Exception:
            pass
    # 2. PubChem by compound name
    if name:
        try:
            import urllib.parse
            enc = urllib.parse.quote(name)
            r = _http.get(
                f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{enc}/property/IsomericSMILES/JSON",
                timeout=6)
            if r.ok:
                props = r.json().get("PropertyTable", {}).get("Properties", [])
                smiles = props[0].get("IsomericSMILES", "") if props else ""
                if smiles:
                    _cache_set(f"smiles_raw:{chembl_id}", smiles)
                    return smiles
        except Exception:
            pass
    return ""


def _clr_gradient(val, mx, palette=("#1e3a5f","#0e4e7a","#0a6e9e","#0ea5e9","#38bdf8","#7dd3fc")):
    if mx <= 0: return palette[0]
    n = max(0.0, min(1.0, val / mx))
    idx = min(int(n*(len(palette)-1)), len(palette)-2)
    t = n*(len(palette)-1) - idx
    c0, c1 = int(palette[idx][1:], 16), int(palette[idx+1][1:], 16)
    r = int((c0>>16)*(1-t)+(c1>>16)*t)
    g = int(((c0>>8)&0xFF)*(1-t)+((c1>>8)&0xFF)*t)
    b = int((c0&0xFF)*(1-t)+(c1&0xFF)*t)
    return f"#{r:02x}{g:02x}{b:02x}"


def human_body_svg(concs: dict) -> str:
    bg="#f8fafc"; ol="#cbd5e1"; tx="#0f172a"; sub="#475569"
    plasma=float(concs.get("plasma",0)); liver=float(concs.get("liver",0))
    brain=float(concs.get("brain",0)); kidney=float(concs.get("kidney",0))
    muscle=float(concs.get("muscle",0))
    mx=max(plasma,liver,brain,kidney,muscle,1.0)
    bc=_clr_gradient(brain,mx); lc=_clr_gradient(liver,mx)
    pc=_clr_gradient(plasma,mx); kc=_clr_gradient(kidney,mx)
    mc=_clr_gradient(muscle,mx)
    def f(v): return f"{v:.0f}" if v>=10 else f"{v:.1f}"
    return f"""<svg viewBox="0 0 300 560" xmlns="http://www.w3.org/2000/svg">
<rect width="300" height="560" fill="{bg}"/>
<ellipse cx="150" cy="55" rx="50" ry="48" fill="{ol}" opacity="0.2"/>
<ellipse cx="150" cy="52" rx="42" ry="40" fill="{bc}" opacity="0.9" stroke="{ol}" stroke-width="1.5"/>
<text x="150" y="47" text-anchor="middle" font-family="sans-serif" font-size="11" fill="{tx}" font-weight="700">Brain</text>
<text x="150" y="62" text-anchor="middle" font-family="sans-serif" font-size="9" fill="{sub}">{f(brain)} ng/mL</text>
<rect x="138" y="95" width="24" height="20" rx="6" fill="{ol}" opacity="0.35"/>
<path d="M 88 115 L 212 115 L 222 308 L 78 308 Z" fill="{ol}" opacity="0.07" stroke="{ol}" stroke-width="1.5"/>
<ellipse cx="120" cy="178" rx="26" ry="38" fill="#0f4c75" opacity="0.7" stroke="{ol}" stroke-width="1"/>
<text x="120" y="175" text-anchor="middle" font-family="sans-serif" font-size="8" fill="{tx}">L.Lung</text>
<ellipse cx="180" cy="178" rx="26" ry="38" fill="#0f4c75" opacity="0.7" stroke="{ol}" stroke-width="1"/>
<text x="180" y="175" text-anchor="middle" font-family="sans-serif" font-size="8" fill="{tx}">R.Lung</text>
<ellipse cx="150" cy="160" rx="20" ry="22" fill="{pc}" opacity="0.9" stroke="{ol}" stroke-width="1.5"/>
<text x="150" y="156" text-anchor="middle" font-family="sans-serif" font-size="9" fill="{tx}" font-weight="600">Heart</text>
<text x="150" y="169" text-anchor="middle" font-family="sans-serif" font-size="8" fill="{sub}">{f(plasma)}</text>
<ellipse cx="137" cy="237" rx="37" ry="26" fill="{lc}" opacity="0.9" stroke="{ol}" stroke-width="1.5"/>
<text x="137" y="234" text-anchor="middle" font-family="sans-serif" font-size="11" fill="{tx}" font-weight="700">Liver</text>
<text x="137" y="249" text-anchor="middle" font-family="sans-serif" font-size="8" fill="{sub}">{f(liver)} ng/mL</text>
<ellipse cx="110" cy="277" rx="18" ry="23" fill="{kc}" opacity="0.9" stroke="{ol}" stroke-width="1.5"/>
<text x="110" y="274" text-anchor="middle" font-family="sans-serif" font-size="7.5" fill="{tx}">L.Kidney</text>
<text x="190" y="274" text-anchor="middle" font-family="sans-serif" font-size="7.5" fill="{tx}">R.Kidney</text>
<ellipse cx="190" cy="277" rx="18" ry="23" fill="{kc}" opacity="0.9" stroke="{ol}" stroke-width="1.5"/>
<rect x="52" y="118" width="32" height="148" rx="16" fill="{mc}" opacity="0.75" stroke="{ol}" stroke-width="1.2"/>
<rect x="216" y="118" width="32" height="148" rx="16" fill="{mc}" opacity="0.75" stroke="{ol}" stroke-width="1.2"/>
<rect x="82" y="308" width="136" height="24" rx="8" fill="{ol}" opacity="0.12" stroke="{ol}" stroke-width="1"/>
<rect x="85" y="330" width="52" height="198" rx="24" fill="{mc}" opacity="0.75" stroke="{ol}" stroke-width="1.2"/>
<rect x="163" y="330" width="52" height="198" rx="24" fill="{mc}" opacity="0.75" stroke="{ol}" stroke-width="1.2"/>
</svg>"""


def generate_3d_html(smiles: str) -> str:
    if not RDKIT_OK:
        return "<body style='background:#f8fafc;color:#475569;font-family:sans-serif;padding:2rem'><p>Requires RDKit for 3D structure generation.</p></body>"
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None: raise ValueError("Invalid SMILES")
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        AllChem.MMFFOptimizeMolecule(mol)
        sdf = Chem.MolToMolBlock(mol)
        sdf_js = sdf.replace('\\', '\\\\').replace('`', '\\`').replace('${', '\\${')
        return f"""<!DOCTYPE html><html><head>
<script src="https://3dmol.org/build/3Dmol-min.js"></script>
<style>body{{margin:0;background:#f8fafc}}#v{{width:100%;height:460px;position:relative}}</style>
</head><body><div id="v"></div>
<script>
var viewer=$3Dmol.createViewer('v',{{backgroundColor:'#f8fafc'}});
viewer.addModel(`{sdf_js}`,'sdf');
viewer.setStyle({{}},{{stick:{{colorscheme:'cyanCarbon',radius:0.15}},sphere:{{scale:0.25}}}});
viewer.addSurface($3Dmol.VDW,{{opacity:0.25,colorscheme:'cyanCarbon'}});
viewer.zoomTo();viewer.render();
</script></body></html>"""
    except Exception as e:
        return f"<body style='background:#f8fafc;color:#475569;padding:2rem'><p>3D error: {e}</p></body>"


# ── Page routes ──────────────────────────────────────────────────────────────
@app.route("/")
def landing():
    return render_template("landing.html")


@app.route("/home")
def index():
    stats = get_stats() if DB_OK else {}
    return render_template("index.html", stats=stats, db_ok=DB_OK)


@app.route("/assets/<path:filename>")
def assets(filename):
    from pathlib import Path
    return send_from_directory(Path(__file__).parent / "assets", filename)


@app.route("/discover")
def discover():
    q = request.args.get("q", "").replace("+", " ").strip()
    results = []; disease_name = q; mesh_ids = []
    if q:
        resolved, expanded, comps = _resolve_fetch(q)
        results = comps
        disease_name = resolved[0]["heading"] if resolved else q
        mesh_ids = expanded
    return render_template("discover.html", results=results,
                           disease_name=disease_name, mesh_ids=json.dumps(mesh_ids),
                           query=q, db_ok=DB_OK)


@app.route("/analysis")
def analysis():
    chembl_id = request.args.get("id", "")
    disease = request.args.get("disease", "")
    dock_method = (
        "Local DiffDock" if LOCAL_DIFFDOCK_OK
        else "NVIDIA DiffDock" if DOCK_OK
        else "AutoDock Vina / RDKit estimate"
    )
    return render_template("analysis.html", chembl_id=chembl_id, disease=disease, db_ok=DB_OK,
                           dock_ok=DOCK_OK, pbpk_ok=PBPK_OK, quantum_ok=QUANTUM_OK,
                           rdkit_ok=RDKIT_OK, py3dmol_ok=PY3DMOL_OK,
                           local_diffdock_ok=LOCAL_DIFFDOCK_OK,
                           diffdock_instructions=LOCAL_DIFFDOCK_INSTRUCTIONS,
                           dock_method=dock_method)


@app.after_request
def _no_html_cache(resp):
    """Never let the browser serve a stale rendered page — inline page JS changes
    (e.g. the knowledge-graph logic) must take effect on the next navigation, not
    after a manual hard-refresh. Static assets are untouched."""
    if "text/html" in resp.headers.get("Content-Type", ""):
        resp.headers["Cache-Control"] = "no-store, max-age=0, must-revalidate"
        resp.headers["Pragma"] = "no-cache"
        resp.headers["Expires"] = "0"
    return resp


@app.route("/graph")
def graph():
    disease = request.args.get("disease", "")
    return render_template("graph.html", disease=disease, db_ok=DB_OK)


@app.route("/database")
def database():
    stats = get_stats() if DB_OK else {}
    footprint = _scrub_sources(get_data_footprint()) if DB_OK else {}
    return render_template("database.html", stats=stats, footprint=footprint, db_ok=DB_OK,
                           standalone=True)


# ── Source-name scrubber ──────────────────────────────────────────────────────
# Per product policy, the underlying *database* sources are not named in the UI.
# Generic, credibility-preserving substitutions applied to all dynamic text
# rendered on the trust pages (validation findings, data footprint). Regulatory
# standards (ALCOA+/GAMP 5) and methodology terms are intentionally kept.
import re as _re

_SOURCE_SUBS = [
    (_re.compile(r"IUPHAR\s*/\s*Guide to Pharmacology", _re.I), "an independent pharmacology authority"),
    (_re.compile(r"IUPHAR\s*/\s*GtoPdb", _re.I), "an independent pharmacology authority"),
    (_re.compile(r"Guide to Pharmacology", _re.I), "an independent pharmacology authority"),
    (_re.compile(r"\bIUPHAR\b", _re.I), "an independent authority"),
    (_re.compile(r"\bGtoPdb\b", _re.I), "the independent authority"),
    (_re.compile(r"\bChEMBL[ _]?33\b"), "the bioactivity database (v33)"),
    (_re.compile(r"\bchembl_33\b"), "the bioactivity database"),
    (_re.compile(r"(?<!p)\bChEMBL\b"), "the bioactivity database"),
    (_re.compile(r"\bHetionet[ _]?v?1\.0\b", _re.I), "the knowledge graph (v1.0)"),
    (_re.compile(r"\bhetionet_v1\.0\b", _re.I), "the knowledge graph"),
    (_re.compile(r"\bHetionet\b", _re.I), "the knowledge graph"),
    (_re.compile(r"\brepoDB\b", _re.I), "an external repurposing gold standard"),
    (_re.compile(r"\bOpen Targets\b"), "a target-association resource"),
    (_re.compile(r"\bMeSH\b"), "the medical disease vocabulary"),
    (_re.compile(r"\bReactome\b", _re.I), "curated pathway data"),
    (_re.compile(r"\bPathwayCommons\b", _re.I), "curated pathway data"),
    (_re.compile(r"\bDrugBank\b", _re.I), "a drug reference"),
    (_re.compile(r"\bBindingDB\b", _re.I), "an independent binding database"),
    (_re.compile(r"\bSTRING\b"), "a protein-interaction database"),
    (_re.compile(r"Brown\s*&(?:amp;)?\s*Patel[^.;)]*", _re.I), "an external gold standard"),
    (_re.compile(r"Himmelstein[^.;)]*", _re.I), "a published method"),
    (_re.compile(r"\bRephetio\b", _re.I), "a published metapath method"),
    (_re.compile(r"\bEMBL-EBI\b"), "the data provider"),
    # The validation surface speaks only for RepurposeIQ (no sibling-platform refs).
    (_re.compile(r"shared by\s*RepurposeIQ\s*\+\s*CompoundIQ\s*\(POZ\)", _re.I), "for RepurposeIQ"),
    (_re.compile(r"RepurposeIQ\s*\+\s*CompoundIQ\s*\(POZ\)", _re.I), "RepurposeIQ"),
    (_re.compile(r"CompoundIQ\s*\(POZ\)", _re.I), "RepurposeIQ"),
    (_re.compile(r"CompoundIQ\s*/\s*POZ", _re.I), "RepurposeIQ"),
    (_re.compile(r"\bCompoundIQ\b", _re.I), "RepurposeIQ"),
    (_re.compile(r"\s*\(POZ\)", _re.I), ""),
    (_re.compile(r"\bPOZ\b"), "RepurposeIQ"),
]


def _scrub_sources(obj):
    """Recursively replace database source names in any string within a dict/list."""
    if isinstance(obj, str):
        s = obj
        for rx, repl in _SOURCE_SUBS:
            s = rx.sub(repl, s)
        return s
    if isinstance(obj, dict):
        return {k: _scrub_sources(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_scrub_sources(v) for v in obj]
    return obj


def _load_json_artifact(filename: str) -> dict:
    """Load a read-only validation artifact from the validation/ folder."""
    from pathlib import Path
    f = Path(__file__).parent / "validation" / filename
    try:
        if f.exists():
            return json.loads(f.read_text(encoding="utf-8"))
    except Exception as e:
        logger.debug(f"validation artifact {filename} load failed: {e}")
    return {}


def _load_validation() -> dict:
    return _load_json_artifact("validation_results_bioactivity.json")


@app.route("/validation")
def validation_page():
    v = _load_validation()
    conc = _load_json_artifact("concordance_results.json")
    res = _load_json_artifact("resolution_results.json")
    mech = _load_json_artifact("mechanisms_results.json")
    pred = _load_json_artifact("predictions_results.json")
    cal = _load_json_artifact("calibration_results.json")
    kg = _load_json_artifact("kg_results.json")
    kgm = _load_json_artifact("kg_model_results.json")
    kge = _load_json_artifact("kg_ensemble_results.json")
    mpv = _load_json_artifact("metapath_results.json")
    footprint = get_data_footprint() if DB_OK else {}
    v, conc, res, mech, pred, cal, kg, kgm, kge, mpv, footprint = (
        _scrub_sources(x) for x in (v, conc, res, mech, pred, cal, kg, kgm, kge, mpv, footprint))
    return render_template("validation.html", v=v, conc=conc, res=res,
                           mech=mech, pred=pred, cal=cal, kg=kg, kgm=kgm, kge=kge, mpv=mpv,
                           footprint=footprint, db_ok=DB_OK, standalone=True)


@app.route("/api/validation")
def api_validation():
    return jsonify(_scrub_sources({"results": _load_validation(),
                    "concordance": _load_json_artifact("concordance_results.json"),
                    "resolution": _load_json_artifact("resolution_results.json"),
                    "mechanisms": _load_json_artifact("mechanisms_results.json"),
                    "predictions": _load_json_artifact("predictions_results.json"),
                    "calibration": _load_json_artifact("calibration_results.json"),
                    "knowledge_graph": _load_json_artifact("kg_results.json"),
                    "kg_model": _load_json_artifact("kg_model_results.json"),
                    "kg_ensemble": _load_json_artifact("kg_ensemble_results.json"),
                    "metapath": _load_json_artifact("metapath_results.json"),
                    "footprint": get_data_footprint() if DB_OK else {}}))


@app.route("/business-case")
def business_case():
    return render_template("business_case.html", db_ok=DB_OK)


@app.route("/docs")
def docs():
    return render_template("docs.html", db_ok=DB_OK)


@app.route("/settings")
def settings():
    return render_template("settings.html", db_ok=DB_OK)


@app.route("/api/system-status")
def api_system_status():
    """Live operational status for the Settings page — real flags, data
    freshness and active configuration (no live external pings; instant)."""
    import os as _os
    from datetime import datetime as _dt

    services = [
        {"name": "PostgreSQL database",        "ok": DB_OK,             "kind": "Data"},
        {"name": "MeSH disease resolver",      "ok": MESH_OK,           "kind": "Data"},
        {"name": "Repurposing scorer",         "ok": SCORER_OK,         "kind": "Engine"},
        {"name": "Discovery engine",           "ok": ENGINE_OK,         "kind": "Engine"},
        {"name": "Reverse (drug→indication)",  "ok": REVERSE_OK,        "kind": "Engine"},
        {"name": "Pathway-first screen",       "ok": PATHWAY_OK,        "kind": "Engine"},
        {"name": "Evidence dossier",           "ok": DOSSIER_OK,        "kind": "Engine"},
        {"name": "Knowledge graph",            "ok": GRAPH_OK,          "kind": "Engine"},
        {"name": "Portfolio matching",         "ok": PORTFOLIO_OK,      "kind": "Engine"},
        {"name": "PBPK simulation",            "ok": PBPK_OK,           "kind": "Compute"},
        {"name": "Docking (NVIDIA DiffDock)",  "ok": DOCK_OK,           "kind": "Compute"},
        {"name": "Docking (local DiffDock)",   "ok": LOCAL_DIFFDOCK_OK, "kind": "Compute"},
        {"name": "Quantum (GFN2-xTB)",         "ok": QUANTUM_OK,        "kind": "Compute"},
        {"name": "RDKit cheminformatics",      "ok": RDKIT_OK,          "kind": "Compute"},
        {"name": "3D rendering (py3Dmol)",     "ok": PY3DMOL_OK,        "kind": "Compute"},
    ]

    data_sources = [
        {"name": "ChEMBL",            "use": "Compounds, mechanisms, activities, indications"},
        {"name": "Open Targets",     "use": "Disease genes, target associations, pathways"},
        {"name": "ClinicalTrials.gov","use": "Trials per indication and per region (live)"},
        {"name": "PubMed / NCBI",    "use": "Literature evidence"},
        {"name": "FDA Orange Book",  "use": "US patents, exclusivity, generics"},
        {"name": "openFDA FAERS",    "use": "Post-market adverse-event signal"},
        {"name": "RCSB PDB",         "use": "Protein structures for docking"},
        {"name": "MeSH / NCBI",      "use": "Disease ontology and synonyms"},
    ]

    # Orange Book freshness
    ob = {"available": False}
    try:
        from services.orange_book import OB_DIR
        pat = OB_DIR / "patent.txt"
        if pat.exists():
            ts = _dt.fromtimestamp(pat.stat().st_mtime)
            ob = {"available": True, "updated": ts.strftime("%Y-%m-%d"),
                  "age_days": (_dt.now() - ts).days}
    except Exception:
        pass

    # Active scoring configuration
    weights = {}
    try:
        from services.repurposing_scorer import (
            W_INDICATION, W_TARGET, W_ACTIVITY, W_NETWORK)
        weights["discovery"] = {"Indication": W_INDICATION, "Target": W_TARGET,
                                "Activity": W_ACTIVITY, "Network": W_NETWORK}
    except Exception:
        pass
    try:
        from services.reverse_repurposing import EVIDENCE_WEIGHTS
        weights["reverse"] = dict(EVIDENCE_WEIGHTS)
    except Exception:
        pass

    # PoS calibrator status
    pos = {"analytic": True, "calibrator": False, "trustworthy": None, "cv_auc": None}
    try:
        import joblib
        from pathlib import Path as _P
        f = _P(__file__).parent / "data" / "pos_model.pkl"
        if f.exists():
            b = joblib.load(f)
            meta = b.get("meta", {}) if isinstance(b, dict) else {}
            pos.update({"calibrator": True, "trustworthy": meta.get("trustworthy"),
                        "cv_auc": meta.get("cv_auc")})
    except Exception:
        pass

    stats = {}
    try:
        if DB_OK:
            stats = get_stats() or {}
    except Exception:
        pass

    return jsonify(_scrub_sources({
        "services": services,
        "data_sources": data_sources,
        "orange_book": ob,
        "weights": weights,
        "pos_model": pos,
        "jurisdictions": ["US", "EU", "India"],
        "db_stats": stats,
        "cache_entries": len(_api_cache),
    }))


@app.route("/api/cache/clear", methods=["POST"])
def api_cache_clear():
    n = len(_api_cache)
    _api_cache.clear()
    return jsonify({"cleared": n})


# ── API routes ───────────────────────────────────────────────────────────────
@app.route("/api/stats")
def api_stats():
    return jsonify(get_stats() if DB_OK else {})


@app.route("/api/search", methods=["POST"])
def api_search():
    data = request.get_json() or {}
    q = data.get("query", "").strip()
    if not q: return jsonify({"error": "No query"}), 400
    resolved, expanded, comps = _resolve_fetch(q)
    disease_name = resolved[0]["heading"] if resolved else q
    return jsonify({"disease": disease_name, "mesh_ids": expanded, "compounds": comps})


@app.route("/api/compound/<chembl_id>")
def api_compound(chembl_id):
    try:
        c = _get_compound(chembl_id)
        if not c:
            return jsonify({"error": f"Compound {chembl_id} not found", "chembl_id": chembl_id}), 404
        smiles = c.get("smiles", "")
        if RDKIT_OK and smiles:
            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    c["qed"] = round(Descriptors.qed(mol), 3)
                    c["tpsa_rdkit"] = round(Descriptors.TPSA(mol), 1)
                    c["rot_bonds"] = Descriptors.NumRotatableBonds(mol)
                    c["arom_rings"] = Descriptors.NumAromaticRings(mol)
            except Exception:
                pass
        # Canonical repurposing score for the (drug, disease) pair — the SAME function
        # the Repurpose card and the Evidence Dossier use, so the number matches everywhere.
        disease = request.args.get("disease", "").strip()
        if disease and REVERSE_OK and not c.get("score"):
            try:
                ps = canonical_pair_score(
                    chembl_id=chembl_id, disease=disease,
                    max_phase=int(float(c.get("max_phase") or 0)),
                    indications=c.get("indications", "") or "",
                    drug_name=c.get("name", ""),
                )
                if ps:
                    c["score"] = ps["score"]
                    c["scores"] = ps["scores"]
                    c["score_breakdown"] = ps["scores"]
            except Exception as _e:
                logger.debug(f"api_compound canonical scoring failed: {_e}")
        return jsonify(c)
    except Exception as e:
        logger.error(f"api_compound {chembl_id}: {e}")
        return jsonify({"error": str(e), "chembl_id": chembl_id}), 500


@app.route("/api/compound/<chembl_id>/targets")
def api_targets(chembl_id):
    c = _get_compound(chembl_id)
    rows = get_compound_targets(c.get("id")) if c and c.get("id") else []
    if not rows:
        rows = _fetch_chembl_mechanism_rows(chembl_id)
    return jsonify(rows)


@app.route("/api/compound/<chembl_id>/activities")
def api_activities(chembl_id):
    c = _get_compound(chembl_id)
    rows = get_compound_activities(c.get("id")) if c and c.get("id") else []
    if not rows and chembl_id and not chembl_id.startswith("NR:"):
        cached = _cache_get(f"chembl_act:{chembl_id}")
        if cached is not None:
            return jsonify(cached)
        try:
            r = _http.get(
                "https://www.ebi.ac.uk/chembl/api/data/activity.json",
                params={"molecule_chembl_id": chembl_id,
                        "pchembl_value__isnull": "false",
                        "limit": 50, "format": "json"},
                timeout=5)
            if r.ok:
                seen: set = set()
                raw_acts = []
                tid_set: set = set()
                for a in r.json().get("activities", []):
                    tgt = a.get("target_pref_name", "")
                    key = (tgt, a.get("activity_type"))
                    if key in seen:
                        continue
                    seen.add(key)
                    pv = a.get("pchembl_value")
                    tid = a.get("target_chembl_id", "")
                    if tid:
                        tid_set.add(tid)
                    raw_acts.append({
                        "target_id":      tid,
                        "target_name":    tgt,
                        "pchembl_value":  float(pv) if pv else None,
                        "activity_type":  a.get("activity_type", ""),
                        "standard_value": a.get("standard_value"),
                        "standard_units": a.get("standard_units", ""),
                        "assay_type":     a.get("assay_type", ""),
                        "assay_chembl_id": a.get("assay_chembl_id", ""),
                    })
                # Resolve gene symbols for each unique target in parallel
                tid_list = list(tid_set)
                if tid_list:
                    with ThreadPoolExecutor(max_workers=min(8, len(tid_list))) as pool:
                        genes = list(pool.map(_chembl_gene_for_target, tid_list))
                    tid_gene = dict(zip(tid_list, genes))
                else:
                    tid_gene = {}
                for act in raw_acts:
                    tid = act.pop("target_id", "")
                    gene = tid_gene.get(tid, "") or act["target_name"][:20]
                    act["gene_symbol"] = gene
                    rows.append(act)
                rows.sort(key=lambda x: -(x.get("pchembl_value") or 0))
                _cache_set(f"chembl_act:{chembl_id}", rows)
        except Exception:
            pass
    return jsonify(rows)


@app.route("/api/compound/<chembl_id>/repurposing-candidates")
def api_repurposing_candidates(chembl_id):
    disease = request.args.get("disease", "")
    mesh_raw = request.args.get("mesh_ids", "")
    mesh_ids = [m for m in mesh_raw.split(",") if m] if mesh_raw else None
    # If disease provided but no mesh_ids, try to resolve
    if disease and not mesh_ids:
        try:
            resolved = resolve_disease(disease)
            if resolved:
                ids = [r["mesh_id"] for r in resolved if r.get("mesh_id")]
                mesh_ids = expand_mesh_ids(ids, include_children=True) or ids
        except Exception:
            pass
    candidates = find_repurposing_candidates(chembl_id, mesh_ids=mesh_ids, disease_name=disease, limit=20)
    candidates = [c for c in candidates if _is_repurposable(c)]
    return jsonify(candidates)


@app.route("/api/compound/<chembl_id>/indications")
def api_indications(chembl_id):
    try:
        c = _get_compound(chembl_id)
        rows = get_compound_indications(c.get("id")) if c and c.get("id") else []
        if not rows and chembl_id and not chembl_id.startswith("NR:"):
            cached = _cache_get(f"chembl_ind:{chembl_id}")
            if cached is not None:
                return jsonify(cached)

            # Fetch ChEMBL + Open Targets in parallel
            def _fetch_chembl_ind():
                try:
                    r = _http.get(
                        "https://www.ebi.ac.uk/chembl/api/data/drug_indication.json",
                        params={"molecule_chembl_id": chembl_id, "limit": 100, "format": "json"},
                        timeout=5)
                    return r.json() if r.ok else None
                except Exception:
                    return None

            def _fetch_ot_ind():
                if not chembl_id.startswith("CHEMBL"):
                    return None
                try:
                    ot_q = """query($id:String!){drug(chemblId:$id){indications{rows{disease{id name}maxPhaseForIndication}}}}"""
                    ot = _http.post(
                        "https://api.platform.opentargets.org/api/v4/graphql",
                        json={"query": ot_q, "variables": {"id": chembl_id}},
                        timeout=5)
                    return ot.json() if ot.ok else None
                except Exception:
                    return None

            with ThreadPoolExecutor(max_workers=2) as pool:
                fut_chembl = pool.submit(_fetch_chembl_ind)
                fut_ot     = pool.submit(_fetch_ot_ind)
                chembl_data = fut_chembl.result()
                ot_data     = fut_ot.result()

            seen_diseases: set = set()
            if chembl_data:
                for ind in chembl_data.get("drug_indications", []):
                    disease = ind.get("mesh_heading") or ind.get("efo_term") or ""
                    if not disease or disease in seen_diseases:
                        continue
                    seen_diseases.add(disease)
                    rows.append({
                        "disease":   disease,
                        "mesh_id":   ind.get("mesh_id", ""),
                        "efo_id":    ind.get("efo_id", ""),
                        "max_phase": _safe_float(ind.get("max_phase_for_ind"), 0),
                        "source":    "Curated",
                    })
            if ot_data:
                ind_rows = ((ot_data.get("data") or {}).get("drug", {}) or {})
                ind_rows = (ind_rows.get("indications") or {}).get("rows", [])
                for ind in ind_rows:
                    disease = (ind.get("disease") or {}).get("name", "")
                    if not disease or disease in seen_diseases:
                        continue
                    seen_diseases.add(disease)
                    rows.append({
                        "disease":   disease,
                        "mesh_id":   (ind.get("disease") or {}).get("id", ""),
                        "efo_id":    (ind.get("disease") or {}).get("id", ""),
                        "max_phase": _safe_float(ind.get("maxPhaseForIndication"), 0),
                        "source":    "Genetic",
                    })
            if rows:
                _cache_set(f"chembl_ind:{chembl_id}", rows)

        # Sanitize: ensure max_phase is always numeric, tree_numbers/entry_terms are lists
        for row in rows:
            row["max_phase"] = _safe_float(row.get("max_phase"), 0)
            for arr_key in ("tree_numbers", "entry_terms", "parent_ids"):
                if row.get(arr_key) is not None and not isinstance(row[arr_key], list):
                    row[arr_key] = list(row[arr_key])
        rows.sort(key=lambda x: -x.get("max_phase", 0))
        return jsonify(_scrub_sources(rows))
    except Exception as e:
        logger.error(f"api_indications error: {e}")
        return jsonify([])


@app.route("/api/compound/<chembl_id>/3d")
def api_3d(chembl_id):
    c = _get_compound(chembl_id)
    smiles = (c or {}).get("smiles", "")
    if not smiles:
        mtype = (c or {}).get("molecule_type", "")
        label = f"{mtype} — no small-molecule structure" if mtype and mtype.lower() not in ("small molecule", "") else "No 3D structure available"
        return f"<p style='color:#64748b;padding:1rem'>{label}.</p>", 200
    return generate_3d_html(smiles), 200


@app.route("/api/compound/<chembl_id>/sdf3d")
def api_sdf3d(chembl_id):
    """Return 3D SDF for a compound. Tries RDKit first, falls back to PubChem REST."""
    c = _get_compound(chembl_id)
    smiles = (c or {}).get("smiles", "") or request.args.get("smiles", "")
    name   = (c or {}).get("name", "")
    # RDKit path
    if RDKIT_OK and smiles:
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                mol = Chem.AddHs(mol)
                AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
                AllChem.MMFFOptimizeMolecule(mol)
                return jsonify({"sdf": Chem.MolToMolBlock(mol), "source": "rdkit"})
        except Exception:
            pass
    # PubChem fallback — try name first, then SMILES
    for kind, val in [("name", name), ("smiles", smiles)]:
        if not val:
            continue
        try:
            import urllib.parse
            r = requests.get(
                f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/{kind}/{urllib.parse.quote(str(val), safe='')}/SDF",
                params={"record_type": "3d"}, timeout=12)
            if r.ok and r.text.strip():
                return jsonify({"sdf": r.text, "source": "pubchem"})
        except Exception:
            pass
    return jsonify({"error": "3D structure unavailable"}), 404


@app.route("/api/compound/<chembl_id>/quantum")
def api_quantum(chembl_id):
    if not QUANTUM_OK: return jsonify({"error": "Quantum module not available"}), 503
    c = _get_compound(chembl_id)
    smiles = (c or {}).get("smiles", "")
    name   = (c or {}).get("name", "Drug")
    if not smiles:
        mtype = (c or {}).get("molecule_type", "")
        if mtype and mtype.lower() not in ("small molecule", ""):
            return jsonify({"error": f"{name} is a {mtype} — quantum optimization requires a small molecule."}), 400
        return jsonify({"error": "No SMILES available for this compound."}), 400
    try:
        res = run_quantum_optimization(name, smiles)
        return jsonify(res or {"success": False})
    except Exception as e:
        return jsonify({"error": str(e)}), 500


def _chembl_pchembl_affinity(chembl_id: str):
    """Best measured potency for a molecule → (ΔG kcal/mol, pChEMBL). For PBPK occupancy,
    a measured Ki/IC50 is far more meaningful than the conservative docking score."""
    if not chembl_id or chembl_id.startswith("NR:"):
        return None, None
    cached = _cache_get(f"pchembl:{chembl_id}")
    if cached is not None:
        return cached
    try:
        r = _http.get("https://www.ebi.ac.uk/chembl/api/data/activity.json",
                      params={"molecule_chembl_id": chembl_id, "pchembl_value__isnull": "false",
                              "limit": 100, "format": "json"}, timeout=8)
        if r.ok:
            vals = [float(a["pchembl_value"]) for a in r.json().get("activities", [])
                    if a.get("pchembl_value")]
            if vals:
                best = max(vals)                       # highest potency
                dG = round(-1.418 * best, 2)           # ΔG = RT·ln(Ki); RT(310K)·ln10 ≈ 1.418
                out = (dG, round(best, 2))
                _cache_set(f"pchembl:{chembl_id}", out)
                return out
    except Exception:
        pass
    return (None, None)


@app.route("/api/pbpk", methods=["POST"])
def api_pbpk():
    if not PBPK_OK: return jsonify({"error": "PBPK not available"}), 503
    data = request.get_json() or {}
    try:
        disease = (data.get("disease") or "").strip()
        sim = PBPKSimulator(disease)
        # Affinity for target-occupancy: explicit docking ΔG, else measured ChEMBL potency
        ba, ki_source = None, None
        aff = data.get("affinity")
        if aff not in (None, ""):
            ba = _safe_float(aff, None); ki_source = "docking ΔG"
        else:
            ba, pchembl = _chembl_pchembl_affinity(data.get("chembl_id", ""))
            if ba is not None:
                ki_source = f"measured potency (pChEMBL {pchembl})"
        res = sim.simulate_drug_exposure(
            drug_name    = data.get("name", "Drug"),
            dose_mg      = _safe_float(data.get("dose"),  100),
            route        = data.get("route", "oral"),
            duration_hours = _safe_float(data.get("hours"), 24),
            binding_affinity = ba,
            params={
                "mw":   _safe_float(data.get("mw"),   350),
                "logp": _safe_float(data.get("logp"),  2.5),
                "psa":  _safe_float(data.get("psa"),  80.0),
                "hba":  int(_safe_float(data.get("hba"),  4)),
                "hbd":  int(_safe_float(data.get("hbd"),  2)),
            },
        )
        if isinstance(res, dict):
            res["ki_source"] = ki_source
        if res.get("success"):
            pk = res.get("pk_metrics", {})
            pl = res["plasma_concentration_ng_ml"]
            li = res["liver_concentration_ng_ml"]
            br = res["brain_concentration_ng_ml"]
            tg = res["target_tissue_concentration_ng_ml"]
            peak = int(np.argmax(pl)) if pl else 0
            body_concs = {
                "plasma": pl[peak] if pl else 0,
                "liver":  li[peak] if li else 0,
                "brain":  br[peak] if br else 0,
                "kidney": (tg[peak] if tg else 0)*0.7,
                "muscle": (tg[peak] if tg else 0)*0.4,
            }
            res["body_svg"] = human_body_svg(body_concs)
            res["pk_metrics"] = pk
        return jsonify(res)
    except Exception as e:
        return jsonify({"error": str(e)}), 500


@app.route("/api/pbpk/routes", methods=["POST"])
def api_pbpk_routes():
    """Multi-route administration analysis: PBPK across every route + known/feasible
    flags + a recommended route for the target organ."""
    if not PBPK_OK: return jsonify({"error": "PBPK not available"}), 503
    data = request.get_json() or {}
    try:
        from services.pbpk_simulation import analyze_routes
        ba, _ = _chembl_pchembl_affinity(data.get("chembl_id", ""))
        res = analyze_routes(
            drug_name    = data.get("name", "Drug"),
            chembl_id    = data.get("chembl_id", ""),
            dose_mg      = _safe_float(data.get("dose"), 100),
            disease_name = (data.get("disease") or "").strip(),
            target_organ = data.get("target_organ", ""),
            binding_affinity = ba,
            params={
                "mw":   _safe_float(data.get("mw"),   350),
                "logp": _safe_float(data.get("logp"),  2.5),
                "psa":  _safe_float(data.get("psa"),  80.0),
                "hbd":  int(_safe_float(data.get("hbd"),  2)),
            },
        )
        return jsonify(res)
    except Exception as e:
        return jsonify({"error": str(e)}), 500


@app.route("/api/docking", methods=["POST"])
def api_docking():
    if not _dock_svc: return jsonify({"error": "Docking service unavailable"}), 503
    data = request.get_json() or {}
    try:
        res = _dock_svc.perform_docking(
            drug_name=data.get("name","Drug"),
            target_name=data.get("target","BACE1"),
            ligand_smiles=data.get("smiles",""),
        )
        return jsonify(res)
    except Exception as e:
        return jsonify({"error": str(e)}), 500


@app.route("/api/graph")
def api_graph():
    disease = request.args.get("disease", "").strip()
    limit   = min(int(request.args.get("limit", 80)), 200)
    try:
        elems, legend = build_network_graph(disease=disease, limit=limit)
        try:
            from services.biocypher_graph import drkg_available
            drkg_ok = drkg_available()
        except Exception:
            drkg_ok = False
        return jsonify({"elements": elems, "legend": legend, "drkg": drkg_ok})
    except Exception as e:
        return jsonify({"elements": [], "legend": {}, "error": str(e)})


@app.route("/api/repurpose-graph")
def api_repurpose_graph():
    """
    Story graph: top-3 repurposing candidates from our DB for a disease,
    connected via shared gene targets to explain WHY each drug is a candidate.
    """
    disease     = request.args.get("disease", "").strip()
    top_n       = min(int(request.args.get("n", 3)), 5)
    compound_id = request.args.get("compound", "").strip()
    # Molecule-first entry: a compound but no disease -> build a COMPOUND-CENTRIC
    # repurposing graph (molecule -> protein targets -> candidate diseases). Needs
    # no known indication — this is the graph a repurposing scientist actually wants.
    if compound_id and not disease:
        try:
            from services.biocypher_graph import build_compound_graph
            elems, legend = build_compound_graph(compound_id, top_n=14)
            if elems:
                return jsonify({"elements": elems, "legend": legend,
                                "compound": compound_id, "mode": "compound"})
        except Exception:
            logger.exception("compound graph error")
        return jsonify({"elements": [], "legend": {}, "error":
                        "No protein targets or candidate diseases found for this molecule."})
    if not disease:
        return jsonify({"elements": [], "legend": {}, "error": "disease param required"})
    try:
        from services.biocypher_graph import build_repurpose_story_graph
        elems, legend = build_repurpose_story_graph(
            disease, top_n=top_n, focal_compound=compound_id or None)
        return jsonify({"elements": elems, "legend": legend, "disease": disease})
    except Exception as e:
        logger.exception("repurpose-graph error")
        return jsonify({"elements": [], "legend": {}, "error": str(e)}), 500


@app.route("/novel-targets")
def novel_targets_page():
    """Novel-target discovery (P2) — inferred targets + drugs reachable via them."""
    disease = request.args.get("disease", "").strip()
    return render_template("novel_targets.html", disease=disease)


@app.route("/api/novel-targets")
def api_novel_targets():
    """Novel-target discovery (P2): infer targets NOT in Open Targets for a disease
    via PPI guilt-by-association, and the drugs reachable only through them."""
    disease = request.args.get("disease", "").strip()
    if not disease:
        return jsonify({"error": "disease param required", "novel_targets": [], "drugs": []})
    with_drugs = request.args.get("drugs", "1") not in ("0", "false", "no")
    try:
        if with_drugs:
            from services.novel_targets import drugs_via_novel_targets
            return jsonify(drugs_via_novel_targets(disease))
        from services.novel_targets import infer_novel_targets
        return jsonify(infer_novel_targets(disease))
    except Exception as e:
        logger.exception("novel-targets error")
        return jsonify({"error": str(e), "novel_targets": [], "drugs": []}), 500


@app.route("/api/resolve-entity")
def api_resolve_entity():
    """Classify a knowledge-graph search query as a MOLECULE or a DISEASE, so the
    graph search can be molecule-aware: a drug -> molecule graph (drug -> targets ->
    diseases); a disease -> the existing disease -> drug graph (unchanged)."""
    q = request.args.get("q", "").strip()
    if not q:
        return jsonify({"type": "none"})
    ql = q.lower()
    # Obvious disease phrasing -> disease (avoids fuzzy drug-name false positives)
    DISEASE_HINTS = ("disease", "syndrome", "cancer", "disorder", "itis", "aemia",
                     "emia", "pathy", "deficiency", "infection", "tumor", "tumour",
                     "carcinoma", "failure", "hypertension", "hypotension",
                     "sclerosis", "fibrosis", "neoplasm", "palsy", "dystrophy")
    if any(h in ql for h in DISEASE_HINTS):
        return jsonify({"type": "disease", "query": q})
    try:
        from services.reverse_repurposing import resolve_drug
        info = resolve_drug(q) or {}
        cid = info.get("chembl_id", "")
        targets = info.get("targets", []) or []
        if cid and targets:   # a real druggable molecule with protein targets
            return jsonify({"type": "drug", "chembl_id": cid,
                            "name": info.get("name", q), "n_targets": len(targets)})
    except Exception:
        logger.exception("resolve-entity error")
    return jsonify({"type": "disease", "query": q})


@app.route("/api/graph-narrative")
def api_graph_narrative():
    """
    Build a mechanistic narrative explaining how the focal compound might treat the disease.
    Uses compound data + graph edges to trace drug → target → pathway → disease.
    """
    compound_id = request.args.get("compound", "").strip()
    disease     = request.args.get("disease", "").strip()
    if not compound_id or not disease:
        return jsonify({"error": "compound and disease params required"}), 400

    c      = _get_compound(compound_id) or {}
    name   = c.get("name", compound_id)
    phase  = int(float(c.get("max_phase", 0) or 0))
    mechs  = c.get("mechanisms", "") or ""
    inds   = c.get("indications", "") or ""
    logp   = float(c.get("logp") or 2.0)
    psa    = float(c.get("psa")  or 80.0)
    mw     = float(c.get("mw")   or 350.0)
    score  = float(c.get("score") or 0.0)

    # Pull gene/protein targets and pathways from the knowledge graph
    graph_genes    = []
    graph_pathways = []
    try:
        from services.biocypher_graph import build_repurpose_story_graph
        elems, _ = build_repurpose_story_graph(disease, top_n=3, focal_compound=compound_id)
        node_map = {e["data"]["id"]: e["data"] for e in elems if "source" not in e["data"]}
        edges    = [e["data"] for e in elems if "source" in e["data"]]
        cpd_node = next((n for n in node_map.values() if n.get("chembl_id") == compound_id), {})
        cpd_id   = cpd_node.get("id", "")
        for edge in edges:
            if edge.get("source") == cpd_id or edge.get("target") == cpd_id:
                other_id = edge["target"] if edge["source"] == cpd_id else edge["source"]
                other    = node_map.get(other_id, {})
                kind     = other.get("kind", "")
                mech_lbl = edge.get("mechanism") or edge.get("label", "")
                if kind in ("Gene", "Protein"):
                    graph_genes.append((other.get("label", ""), other.get("gene_name", "") or other.get("full_label", ""), mech_lbl))
                elif kind == "Pathway":
                    graph_pathways.append(other.get("full_label") or other.get("label", ""))
    except Exception:
        pass

    # Fallback: fetch mechanisms + gene symbols directly from ChEMBL if graph gave nothing
    mech_actions = []
    gene_names   = []
    if not graph_genes:
        try:
            r = _http.get(
                "https://www.ebi.ac.uk/chembl/api/data/mechanism.json",
                params={"molecule_chembl_id": compound_id, "limit": 8, "format": "json"},
                timeout=5)
            if r.ok:
                mdata  = r.json().get("mechanisms", [])
                mech_actions = [m.get("mechanism_of_action", "") for m in mdata if m.get("mechanism_of_action")][:3]
                tids   = [m.get("target_chembl_id") for m in mdata if m.get("target_chembl_id")][:4]
                if tids:
                    with ThreadPoolExecutor(max_workers=4) as pool:
                        infos = list(pool.map(_chembl_target_info, tids))
                    gene_names = [(g, p) for g, p in infos if g][:4]
        except Exception:
            pass
    else:
        gene_names   = [(g[0], g[1]) for g in graph_genes[:4]]
        mech_actions = [g[2] for g in graph_genes if g[2] and g[2] not in ("BINDS_TO", "ASSOCIATED_WITH", "binds_to")][:3]
        if not mech_actions:
            mech_actions = [mechs.split(",")[0].strip()] if mechs else []

    # Human-readable labels
    phase_labels = {0: "preclinical compound", 1: "Phase 1 compound", 2: "Phase 2 compound",
                    3: "Phase 3 compound",     4: "FDA-approved drug"}
    phase_str  = phase_labels.get(phase, "investigational compound")
    gene_str   = ", ".join([g for g, _ in gene_names[:3]]) if gene_names else "relevant molecular targets"
    target_str = "; ".join([f"{g} ({p})" if p else g for g, p in gene_names[:3]]) if gene_names else gene_str
    mech_prim  = mech_actions[0] if mech_actions else (mechs.split(",")[0].strip() if mechs else "modulating disease-relevant pathways")
    ind_list   = [i.strip() for i in inds.split(",") if i.strip()][:3]
    approved_str = f"It is already approved for {', '.join(ind_list[:2])}. " if ind_list else ""

    # CNS penetration heuristic (simplified Lipinski / PSA rule)
    bbb_ok  = logp > 0.5 and psa < 90 and mw < 500
    bbb_txt = ""
    if bbb_ok:
        bbb_txt = (f"Its physicochemical profile (LogP {logp:.1f}, PSA {psa:.0f} Å², MW {mw:.0f} Da) "
                   "is consistent with blood–brain barrier penetration — important for neurological targets.")
    elif psa > 120:
        bbb_txt = (f"Its high polar surface area ({psa:.0f} Å²) may limit brain penetration; "
                   "formulation or prodrug strategies could be considered for CNS delivery.")

    pathway_txt = ""
    if graph_pathways:
        pathway_txt = (f"The knowledge graph links it to {', '.join(graph_pathways[:2])} — "
                       f"pathways directly implicated in {disease} pathogenesis. ")

    # Build hypothesis paragraph
    disease_short = disease.split("(")[0].strip()
    parts = [
        f"{name} is a {phase_str} that acts by {mech_prim}.",
        f"Its primary molecular targets are {target_str}.",
        approved_str,
        f"In {disease_short}, dysregulation of {gene_str.split(',')[0].strip() if gene_names else 'these molecular pathways'} "
        f"is a well-established driver of disease onset and progression.",
        pathway_txt,
        f"By blocking or modulating this axis, {name} has the potential to slow pathological signalling, "
        f"reduce neuroinflammation or cellular stress, and preserve neuronal function in {disease_short} patients.",
        bbb_txt,
    ]
    hypothesis = " ".join(p for p in parts if p.strip())

    # Evidence-chain bullets
    evidence_chain = []
    if gene_names:
        evidence_chain.append({"icon": "🧬", "text": f"Targets {', '.join([g for g, _ in gene_names[:2]])} — validated in {disease_short} literature"})
    if mech_actions:
        evidence_chain.append({"icon": "⚗️", "text": f"Mechanism: {mech_actions[0]}"})
    if phase >= 3:
        evidence_chain.append({"icon": "🏥", "text": f"Phase {phase} clinical data — safety and efficacy signals established"})
    elif phase >= 1:
        evidence_chain.append({"icon": "🏥", "text": f"Phase {phase} data available — early human pharmacology known"})
    if graph_pathways:
        evidence_chain.append({"icon": "🔗", "text": f"Graph connects drug via {graph_pathways[0]} — key disease pathway"})
    if ind_list:
        evidence_chain.append({"icon": "✅", "text": f"Approved in {ind_list[0]} — potential mechanism overlap with {disease_short}"})
    if score:
        evidence_chain.append({"icon": "📊", "text": f"Repurposing score {score*100:.0f}/100 — target overlap, clinical evidence & network analysis"})
    if bbb_txt:
        evidence_chain.append({"icon": "🧠", "text": bbb_txt})

    return jsonify({
        "drug_name":      name,
        "disease":        disease,
        "phase":          phase,
        "phase_label":    phase_str,
        "hypothesis":     hypothesis,
        "evidence_chain": evidence_chain,
        "targets":        [[g, p] for g, p in gene_names[:4]],
        "mechanisms":     mech_actions[:3],
        "caution":        ("This is a computational repurposing analysis based on pharmacological and network evidence. "
                           "Clinical validation is required before therapeutic application."),
    })


@app.route("/api/drkg-status")
def api_drkg_status():
    try:
        from database.drkg_loader import cache_available, CACHE_FILE
        avail = cache_available()
        size  = CACHE_FILE.stat().st_size // 1024 if avail else 0
        return jsonify({"available": avail, "size_kb": size})
    except Exception as e:
        return jsonify({"available": False, "error": str(e)})


@app.route("/api/drkg-import", methods=["POST"])
def api_drkg_import():
    """Trigger DRKG download + import in the background."""
    import threading
    def _run():
        try:
            from database.drkg_loader import build_cache
            from services import biocypher_graph as bg
            build_cache()
            # Reset in-memory cache so next graph request picks up new data
            bg._drkg_cache = None
            bg._drkg_loaded = False
            logger.info("DRKG import complete")
        except Exception as e:
            logger.error(f"DRKG import failed: {e}")
    t = threading.Thread(target=_run, daemon=True)
    t.start()
    return jsonify({"status": "started",
                    "message": "DRKG download started in background (~5–15 min). "
                               "Refresh the graph when complete."})


@app.route("/api/database/search")
def api_db_search():
    q = request.args.get("q", "").strip()
    if not q: return jsonify([])
    rows = db_search(q, limit=40)
    return jsonify(rows)


@app.route("/api/disease-targets")
def api_disease_targets():
    disease = request.args.get("disease", "").strip()
    if not disease:
        return jsonify([])
    try:
        from services.disease_ontology import resolve_disease as ot_resolve
        info = ot_resolve(disease)
        targets = info.get("targets", [])[:25]
        return jsonify([
            {
                "gene_symbol": t["gene_symbol"],
                "gene_name":   t.get("gene_name", ""),
                "score":       round(t.get("score", 0), 3),
            }
            for t in targets
        ])
    except Exception as e:
        logger.debug(f"disease-targets error: {e}")
        return jsonify([])


@app.route("/api/suggest-diseases")
def api_suggest_diseases():
    q = request.args.get("q", "").strip()
    try:
        from services.disease_resolver import suggest_diseases
        return jsonify(suggest_diseases(q, limit=8))
    except Exception:
        return jsonify([])


@app.route("/api/clinical")
def api_clinical():
    drug = request.args.get("drug", "")
    disease = request.args.get("disease", "")
    trials = []
    try:
        r = requests.get("https://clinicaltrials.gov/api/v2/studies",
            params={"query.cond": disease, "query.term": drug,
                    "pageSize": 8, "format": "json"}, timeout=10)
        if r.ok:
            for s in r.json().get("studies", []):
                pm = s.get("protocolSection", {})
                nct = pm.get("identificationModule", {}).get("nctId", "")
                trials.append({
                    "nct": nct,
                    "title": pm.get("identificationModule", {}).get("briefTitle","")[:80],
                    "status": pm.get("statusModule", {}).get("overallStatus",""),
                    "phase": ", ".join(pm.get("designModule",{}).get("phases",[])) or "N/A",
                    "url": f"https://clinicaltrials.gov/study/{nct}",
                })
    except Exception: pass
    return jsonify(trials)


@app.route("/api/literature")
def api_literature():
    drug = request.args.get("drug", "")
    disease = request.args.get("disease", "")
    papers = []
    try:
        base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
        sr = requests.get(f"{base}/esearch.fcgi",
            params={"db":"pubmed","term":f"{drug}[tiab] AND {disease}[tiab]",
                    "retmax":8,"retmode":"json"}, timeout=8)
        ids = sr.json().get("esearchresult",{}).get("idlist",[])
        if ids:
            smr = requests.get(f"{base}/esummary.fcgi",
                params={"db":"pubmed","id":",".join(ids),"retmode":"json"}, timeout=8)
            res = smr.json().get("result",{})
            papers = [{"pmid":p,"title":res[p].get("title","")[:80],
                       "journal":res[p].get("source",""),"year":res[p].get("pubdate","")[:4],
                       "url":f"https://pubmed.ncbi.nlm.nih.gov/{p}/"} for p in ids if p in res]
    except Exception: pass
    return jsonify(papers)


@app.route("/api/compound/<chembl_id>/properties")
def api_enhanced_properties(chembl_id):
    try:
      c = _get_compound(chembl_id)
    except Exception as e:
      logger.error(f"api_enhanced_properties _get_compound failed: {e}")
      return jsonify({"error": str(e)}), 500
    if not c: return jsonify({"error": f"Compound {chembl_id} not found"}), 404
    result = {k: c.get(k) for k in ("mw","alogp","psa","hba","hbd","ro5_violations","max_phase","qed","rtb","smiles","chembl_id","name")}
    smiles = c.get("smiles", "")
    if RDKIT_OK and smiles:
        try:
            from rdkit.Chem import rdMolDescriptors
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                mw   = result["mw"]   or round(Descriptors.MolWt(mol), 1)
                logp = result["alogp"] or round(Descriptors.MolLogP(mol), 2)
                psa  = result["psa"]  or round(Descriptors.TPSA(mol), 1)
                hbd  = result["hbd"]  if result["hbd"]  is not None else rdMolDescriptors.CalcNumHBD(mol)
                hba  = result["hba"]  if result["hba"]  is not None else rdMolDescriptors.CalcNumHBA(mol)
                rb   = result["rtb"]  if result["rtb"]  is not None else Descriptors.NumRotatableBonds(mol)
                qed  = round(Descriptors.qed(mol), 3)
                ap   = Descriptors.NumAromaticRings(mol)
                # BBB penetration (rule-based, 0–100)
                bbb = 0
                if psa < 60:   bbb += 30
                elif psa < 90: bbb += 15
                if mw < 400:   bbb += 25
                elif mw < 500: bbb += 10
                if 1 <= logp <= 3:     bbb += 25
                elif 0 <= logp <= 4.5: bbb += 10
                if hbd <= 2:   bbb += 20
                elif hbd <= 3: bbb += 10
                bbb_level = "High" if bbb >= 75 else "Moderate" if bbb >= 45 else "Low"
                # CNS-MPO (0–6 scale)
                cns = 0.0
                if logp <= 3:   cns += 1.0
                elif logp <= 5: cns += max(0, (5 - logp) / 2)
                logd = logp - 1.0
                if logd <= 2:   cns += 1.0
                elif logd <= 4: cns += max(0, (4 - logd) / 2)
                if mw <= 360:   cns += 1.0
                elif mw <= 500: cns += (500 - mw) / 140
                if psa <= 40:   cns += 1.0
                elif psa <= 90: cns += max(0, (90 - psa) / 50)
                if hbd == 0:   cns += 1.0
                elif hbd == 1: cns += 0.75
                elif hbd == 2: cns += 0.25
                cns += 0.75  # default pKa contribution
                # ESOL estimated aqueous solubility
                logs = 0.16 - 0.63 * logp - 0.0062 * mw + 0.066 * rb - 0.74 * ap
                sol = "High" if logs > -1 else "Good" if logs > -2 else "Moderate" if logs > -3 else "Low" if logs > -4 else "Poor"
                result.update({
                    "mw": mw, "alogp": logp, "psa": psa, "hba": hba, "hbd": hbd,
                    "rtb": rb, "qed": qed, "arom_rings": ap,
                    "bbb_score": bbb, "bbb_level": bbb_level,
                    "cns_mpo": round(cns, 2), "logs": round(logs, 2), "sol_class": sol,
                })
        except Exception as e:
            result["rdkit_error"] = str(e)
    return jsonify(result)


def _get_therapeutic_profile(disease: str) -> dict:
    """
    Map disease name to therapeutic area property profile.
    lipinski_rules: list of (label, key, op, threshold)
      key: "mw"|"logp"|"psa"|"hbd"|"hba"
      op:  "le" (≤ thr) | "ge" (≥ thr) | "range" (lo ≤ x ≤ hi, thr=(lo,hi))
    Evaluated in api_optimize with actual compound values.
    """
    d = (disease or "").lower()

    def matches(keywords):
        return any(kw in d for kw in keywords)

    if matches(["alzheimer", "parkinson", "epilep", "multiple sclerosis", "schizophreni",
                "depression", "bipolar", "amyotrophic", "huntington", "migraine",
                "dementia", "anxiety", "adhd", "autism", "psychosis", "neurolog",
                "neuropath", "cerebell", "meningit", "encephalit", "als ", "lewy body",
                "frontotemporal", "vascular dementia"]):
        return {
            "area": "CNS / Neurological", "name": "CNS MPO Rules",
            "mw_ideal": 340, "mw_max": 500,
            "logp_ideal": 2.0, "logp_lo": 1.0, "logp_hi": 3.0,
            "psa_ideal": 60,  "psa_max": 90,
            "hbd_ideal": 2,   "hbd_max": 3,
            "hba_ideal": 6,   "hba_max": 7,
            "rb_ideal": 5,    "rb_max": 8,
            "qed_ideal": 0.70,
            "accept_violations": False,
            "lipinski_rules": [
                ("MW ≤ 500 Da",  "mw",   "le",    500),
                ("PSA ≤ 90 Å²",  "psa",  "le",    90),
                ("LogP 1–3",     "logp", "range", (1.0, 3.0)),
                ("HBD ≤ 3",      "hbd",  "le",    3),
            ],
        }
    if matches(["cancer", "tumor", "tumour", "carcinoma", "leukemia", "leukaemia",
                "lymphoma", "melanoma", "sarcoma", "glioma", "glioblastoma",
                "myeloma", "oncol", "neoplasm", "metastasis", "pancreatic",
                "hepatocell", "adenocarcinoma", "renal cell", "breast cancer",
                "lung cancer", "colon cancer", "ovarian", "prostate cancer"]):
        return {
            "area": "Oncology", "name": "Beyond Ro5 (Oncology)",
            "mw_ideal": 700, "mw_max": 900,
            "logp_ideal": 3.0, "logp_lo": -1.0, "logp_hi": 6.0,
            "psa_ideal": 160, "psa_max": 250,
            "hbd_ideal": 5,   "hbd_max": 8,
            "hba_ideal": 12,  "hba_max": 15,
            "rb_ideal": 10,   "rb_max": 15,
            "qed_ideal": 0.30,
            "accept_violations": True,
            "lipinski_rules": [
                ("MW ≤ 900 Da",  "mw",   "le", 900),
                ("LogP ≤ 6",     "logp", "le", 6),
                ("HBD ≤ 8",      "hbd",  "le", 8),
                ("HBA ≤ 15",     "hba",  "le", 15),
            ],
        }
    if matches(["cardiac", "cardiovascular", "heart failure", "hypertension",
                "arrhythmia", "atrial fibrill", "coronary", "atherosclerosis",
                "myocardial", "stroke", "thrombosis", "angina"]):
        return {
            "area": "Cardiovascular", "name": "Lipinski + hERG Filter",
            "mw_ideal": 400, "mw_max": 500,
            "logp_ideal": 2.0, "logp_lo": 0.0, "logp_hi": 3.0,
            "psa_ideal": 80,  "psa_max": 120,
            "hbd_ideal": 3,   "hbd_max": 5,
            "hba_ideal": 8,   "hba_max": 10,
            "rb_ideal": 6,    "rb_max": 9,
            "qed_ideal": 0.60,
            "accept_violations": False,
            "lipinski_rules": [
                ("MW ≤ 500 Da",  "mw",   "le",    500),
                ("LogP 0–3",     "logp", "range", (0.0, 3.0)),
                ("HBD ≤ 5",      "hbd",  "le",    5),
                ("HBA ≤ 10",     "hba",  "le",    10),
            ],
        }
    if matches(["hiv", "aids", "tuberculosis", " tb ", "malaria", "antibiotic",
                "antiviral", "antifungal", "bacterial", "virus", "viral ", "fungal",
                "covid", "influenza", "hepatitis", "sepsis", "infection",
                "antiparasit", "leishmani", "trypanosoma"]):
        return {
            "area": "Infectious Disease", "name": "eNTRy / Lipinski Ro5",
            "mw_ideal": 450, "mw_max": 600,
            "logp_ideal": 1.0, "logp_lo": -1.0, "logp_hi": 3.0,
            "psa_ideal": 80,  "psa_max": 140,
            "hbd_ideal": 3,   "hbd_max": 5,
            "hba_ideal": 9,   "hba_max": 12,
            "rb_ideal": 6,    "rb_max": 10,
            "qed_ideal": 0.55,
            "accept_violations": False,
            "lipinski_rules": [
                ("MW ≤ 600 Da",  "mw",   "le", 600),
                ("LogP ≤ 3",     "logp", "le", 3.0),
                ("HBD ≤ 5",      "hbd",  "le", 5),
                ("HBA ≤ 12",     "hba",  "le", 12),
            ],
        }
    if matches(["diabetes", "diabetic", "obesity", "obese", "nash", "nafld",
                "dyslipidemia", "insulin", "metabolic syndrome", "thyroid",
                "hyperlipid", "hyperglycemi"]):
        return {
            "area": "Metabolic / Diabetes", "name": "Lipinski Ro5",
            "mw_ideal": 400, "mw_max": 500,
            "logp_ideal": 2.5, "logp_lo": 0.0, "logp_hi": 4.0,
            "psa_ideal": 80,  "psa_max": 120,
            "hbd_ideal": 3,   "hbd_max": 5,
            "hba_ideal": 8,   "hba_max": 10,
            "rb_ideal": 6,    "rb_max": 9,
            "qed_ideal": 0.60,
            "accept_violations": False,
            "lipinski_rules": [
                ("MW ≤ 500 Da",  "mw",   "le", 500),
                ("LogP ≤ 5",     "logp", "le", 5),
                ("HBD ≤ 5",      "hbd",  "le", 5),
                ("HBA ≤ 10",     "hba",  "le", 10),
            ],
        }
    if matches(["rheumatoid", "lupus", "crohn", "inflammatory bowel", "psoriasis",
                "autoimmune", "colitis", "sjogren", "scleroderma", "asthma",
                "atopic", "ankylosing", "ibd "]):
        return {
            "area": "Inflammatory / Autoimmune", "name": "Ro5 + GI Permeability",
            "mw_ideal": 450, "mw_max": 550,
            "logp_ideal": 2.5, "logp_lo": 0.0, "logp_hi": 4.0,
            "psa_ideal": 100, "psa_max": 140,
            "hbd_ideal": 4,   "hbd_max": 5,
            "hba_ideal": 9,   "hba_max": 11,
            "rb_ideal": 7,    "rb_max": 10,
            "qed_ideal": 0.60,
            "accept_violations": False,
            "lipinski_rules": [
                ("MW ≤ 550 Da",  "mw",   "le", 550),
                ("LogP ≤ 4",     "logp", "le", 4.0),
                ("HBD ≤ 5",      "hbd",  "le", 5),
                ("HBA ≤ 11",     "hba",  "le", 11),
            ],
        }
    if matches(["copd", "pulmonary fibrosis", "cystic fibrosis", "respiratory",
                "bronchitis", "pneumonia", "lung disease", "pulmonary hypertension"]):
        return {
            "area": "Respiratory", "name": "Lipinski Ro5",
            "mw_ideal": 400, "mw_max": 500,
            "logp_ideal": 2.5, "logp_lo": 0.0, "logp_hi": 5.0,
            "psa_ideal": 80,  "psa_max": 130,
            "hbd_ideal": 3,   "hbd_max": 5,
            "hba_ideal": 8,   "hba_max": 10,
            "rb_ideal": 6,    "rb_max": 9,
            "qed_ideal": 0.60,
            "accept_violations": False,
            "lipinski_rules": [
                ("MW ≤ 500 Da",  "mw",   "le", 500),
                ("LogP ≤ 5",     "logp", "le", 5),
                ("HBD ≤ 5",      "hbd",  "le", 5),
                ("HBA ≤ 10",     "hba",  "le", 10),
            ],
        }
    if matches(["macular", "glaucoma", "retinopathy", "ophthalm", "ocular",
                "age-related macular", "amd ", "diabetic retina"]):
        return {
            "area": "Ophthalmology", "name": "Ocular MPO",
            "mw_ideal": 320, "mw_max": 400,
            "logp_ideal": 1.5, "logp_lo": 0.0, "logp_hi": 3.0,
            "psa_ideal": 60,  "psa_max": 80,
            "hbd_ideal": 2,   "hbd_max": 4,
            "hba_ideal": 6,   "hba_max": 8,
            "rb_ideal": 5,    "rb_max": 8,
            "qed_ideal": 0.65,
            "accept_violations": False,
            "lipinski_rules": [
                ("MW ≤ 400 Da",  "mw",   "le",    400),
                ("LogP 0–3",     "logp", "range", (0.0, 3.0)),
                ("PSA ≤ 80 Å²",  "psa",  "le",    80),
                ("HBD ≤ 4",      "hbd",  "le",    4),
            ],
        }
    if matches(["neuropathic pain", "fibromyalgia", "analges", "opioid",
                "chronic pain", " pain "]):
        return {
            "area": "Pain / Analgesia", "name": "CNS-lite Rules",
            "mw_ideal": 380, "mw_max": 500,
            "logp_ideal": 2.0, "logp_lo": 0.0, "logp_hi": 4.0,
            "psa_ideal": 70,  "psa_max": 100,
            "hbd_ideal": 3,   "hbd_max": 4,
            "hba_ideal": 7,   "hba_max": 9,
            "rb_ideal": 6,    "rb_max": 9,
            "qed_ideal": 0.65,
            "accept_violations": False,
            "lipinski_rules": [
                ("MW ≤ 500 Da",  "mw",   "le",    500),
                ("PSA ≤ 100 Å²", "psa",  "le",    100),
                ("LogP 0–4",     "logp", "range", (0.0, 4.0)),
                ("HBD ≤ 4",      "hbd",  "le",    4),
            ],
        }
    if matches(["eczema", "acne", "atopic dermatitis", "rosacea", "dermatolog",
                "wound healing", "skin disorder"]):
        return {
            "area": "Dermatology", "name": "Skin Penetration Rules",
            "mw_ideal": 420, "mw_max": 550,
            "logp_ideal": 2.5, "logp_lo": 1.0, "logp_hi": 4.0,
            "psa_ideal": 80,  "psa_max": 120,
            "hbd_ideal": 3,   "hbd_max": 5,
            "hba_ideal": 9,   "hba_max": 11,
            "rb_ideal": 6,    "rb_max": 9,
            "qed_ideal": 0.60,
            "accept_violations": False,
            "lipinski_rules": [
                ("MW ≤ 550 Da",  "mw",   "le",    550),
                ("LogP 1–4",     "logp", "range", (1.0, 4.0)),
                ("HBD ≤ 5",      "hbd",  "le",    5),
                ("HBA ≤ 11",     "hba",  "le",    11),
            ],
        }
    if matches(["gerd", "ibs ", "irritable bowel", "helicobacter", "gastroparesis",
                "constipation", "ulcerative colitis", "gastrointestin"]):
        return {
            "area": "Gastrointestinal", "name": "GI-Topical Rules",
            "mw_ideal": 500, "mw_max": 700,
            "logp_ideal": 1.0, "logp_lo": -1.0, "logp_hi": 3.0,
            "psa_ideal": 160, "psa_max": 220,
            "hbd_ideal": 5,   "hbd_max": 7,
            "hba_ideal": 12,  "hba_max": 14,
            "rb_ideal": 8,    "rb_max": 12,
            "qed_ideal": 0.40,
            "accept_violations": True,
            "lipinski_rules": [
                ("MW ≤ 700 Da",  "mw",   "le", 700),
                ("LogP ≤ 3",     "logp", "le", 3.0),
                ("HBD ≤ 7",      "hbd",  "le", 7),
                ("HBA ≤ 14",     "hba",  "le", 14),
            ],
        }
    if matches(["rare disease", "orphan", "gaucher", "fabry", "pompe",
                "spinal muscular", " sma ", "hemophilia", "phenylketonuria"]):
        return {
            "area": "Rare / Orphan", "name": "Beyond Ro5 (Orphan)",
            "mw_ideal": 600, "mw_max": 900,
            "logp_ideal": 3.0, "logp_lo": -2.0, "logp_hi": 6.0,
            "psa_ideal": 150, "psa_max": 200,
            "hbd_ideal": 6,   "hbd_max": 8,
            "hba_ideal": 12,  "hba_max": 15,
            "rb_ideal": 10,   "rb_max": 15,
            "qed_ideal": 0.35,
            "accept_violations": True,
            "lipinski_rules": [
                ("MW ≤ 900 Da",  "mw",   "le", 900),
                ("LogP ≤ 6",     "logp", "le", 6),
                ("HBD ≤ 8",      "hbd",  "le", 8),
                ("HBA ≤ 15",     "hba",  "le", 15),
            ],
        }
    # Default: standard Lipinski Ro5
    return {
        "area": "General", "name": "Lipinski Ro5",
        "mw_ideal": 400, "mw_max": 500,
        "logp_ideal": 2.5, "logp_lo": 0.0, "logp_hi": 5.0,
        "psa_ideal": 80,  "psa_max": 120,
        "hbd_ideal": 3,   "hbd_max": 5,
        "hba_ideal": 8,   "hba_max": 10,
        "rb_ideal": 6,    "rb_max": 9,
        "qed_ideal": 0.60,
        "accept_violations": False,
        "lipinski_rules": [
            ("MW ≤ 500 Da",  "mw",   "le", 500),
            ("LogP ≤ 5",     "logp", "le", 5),
            ("HBD ≤ 5",      "hbd",  "le", 5),
            ("HBA ≤ 10",     "hba",  "le", 10),
        ],
    }


@app.route("/api/compound/<chembl_id>/optimize")
def api_optimize(chembl_id):
    if not RDKIT_OK: return jsonify({"error": "RDKit not available"}), 503
    c = _get_compound(chembl_id)
    if not c: return jsonify({"error": f"Compound {chembl_id} not found"}), 404
    smiles = c.get("smiles", "")
    if not smiles:
        mtype = c.get("molecule_type", "")
        if mtype and mtype.lower() not in ("small molecule", ""):
            return jsonify({"error": (
                f"{c.get('name', chembl_id)} is a {mtype} — SMILES-based optimization "
                "applies only to small molecules."
            )}), 400
        return jsonify({"error": "No SMILES available for this compound."}), 400
    try:
        from rdkit.Chem import rdMolDescriptors
        mol = Chem.MolFromSmiles(smiles)
        if mol is None: return jsonify({"error": "Invalid SMILES"}), 400
        mw   = round(Descriptors.MolWt(mol), 1)
        logp = round(Descriptors.MolLogP(mol), 2)
        psa  = round(Descriptors.TPSA(mol), 1)
        hbd  = rdMolDescriptors.CalcNumHBD(mol)
        hba  = rdMolDescriptors.CalcNumHBA(mol)
        rb   = Descriptors.NumRotatableBonds(mol)
        qed  = round(Descriptors.qed(mol), 3)
        ap   = Descriptors.NumAromaticRings(mol)

        disease = request.args.get("disease", "").strip()
        profile = _get_therapeutic_profile(disease)
        # Evaluate (label, key, op, threshold) rules with current compound values
        _vals = {"mw": mw, "logp": logp, "psa": psa, "hbd": hbd, "hba": hba}
        lipinski_rules = []
        for label, key, op, thr in profile["lipinski_rules"]:
            v = _vals.get(key, 0)
            if op == "le":
                passed = v <= thr
            elif op == "ge":
                passed = v >= thr
            elif op == "range":
                passed = thr[0] <= v <= thr[1]
            else:
                passed = False
            lipinski_rules.append({"rule": label, "pass": passed})
        lipinski_pass_all = all(r["pass"] for r in lipinski_rules)

        # Profile-specific suggestions
        suggestions = []
        psa_thr  = profile["psa_max"]
        mw_thr   = profile["mw_max"]
        logp_hi  = profile["logp_hi"]
        logp_lo  = profile["logp_lo"]
        hbd_thr  = profile["hbd_max"]

        if psa > psa_thr:
            action = ("Replace polar groups with non-polar bioisosteres; N-methylate amides"
                      if profile["area"] != "Gastrointestinal"
                      else "Reduce polar substituents; add lipophilic groups for GI absorption")
            suggestions.append({"property": "PSA", "current": f"{psa:.1f} Å²",
                                 "target": f"< {psa_thr} Å²",
                                 "impact": "BBB penetration" if "CNS" in profile["area"] else "Membrane permeability",
                                 "action": action})
        if mw > mw_thr:
            suggestions.append({"property": "Mol. Weight", "current": f"{mw:.1f} Da",
                                 "target": f"< {mw_thr} Da",
                                 "impact": "Oral bioavailability" if not profile["accept_violations"] else "Formulation & delivery",
                                 "action": "Remove non-essential substituents; truncate peripheral chains"})
        if logp > logp_hi:
            suggestions.append({"property": "LogP", "current": f"{logp:.2f}",
                                 "target": f"{logp_lo:.0f}–{logp_hi:.0f}",
                                 "impact": "Solubility & metabolic stability",
                                 "action": "Add polar groups (OH, NH₂); reduce aromatic ring count"})
        elif logp < logp_lo:
            suggestions.append({"property": "LogP", "current": f"{logp:.2f}",
                                 "target": f"{logp_lo:.0f}–{logp_hi:.0f}",
                                 "impact": "Membrane permeability",
                                 "action": "Add methyl/ethyl groups; replace amines with lipophilic bioisosteres"})
        if hbd > hbd_thr:
            suggestions.append({"property": "H-Bond Donors", "current": str(hbd),
                                 "target": f"≤ {hbd_thr}",
                                 "impact": "Membrane permeability & BBB" if "CNS" in profile["area"] else "Oral absorption",
                                 "action": "Cap -OH/-NH with esters or amides; N-methylate NH groups"})
        if rb > profile["rb_max"]:
            suggestions.append({"property": "Rotatable Bonds", "current": str(rb),
                                 "target": f"≤ {profile['rb_max']}",
                                 "impact": "Oral bioavailability (conformational entropy)",
                                 "action": "Cyclize flexible chains; add conformational constraints"})
        if ap > 3:
            suggestions.append({"property": "Aromatic Rings", "current": str(ap),
                                 "target": "≤ 3",
                                 "impact": "Metabolic stability & hERG risk",
                                 "action": "Replace aromatic rings with saturated bioisosteres"})
        if not suggestions:
            suggestions.append({"property": "Drug-likeness", "current": "All parameters within range",
                                 "target": "Maintain",
                                 "impact": f"Favorable profile for {profile['area']}",
                                 "action": "Focus on potency/selectivity optimization; physicochemical properties are already within target range"})

        # Radar — profile-aware labels and scoring
        area = profile["area"]
        if "Oncology" in area:
            radar_labels = ["Tumor Access", "Solubility", "Cell Uptake", "Drug-like", "Metabolic", "Selectivity"]
            radar_values = [
                min(100, max(0, 100 - max(0, psa - 80) * 0.5)),
                min(100, max(0, 100 - max(0, logp - 2) * 10)),
                min(100, max(0, 100 - max(0, mw - 400) * 0.12)),
                min(100, round(qed * 100)),
                min(100, max(0, 100 - rb * 4)),
                min(100, max(0, 80 - ap * 5)),
            ]
        elif "CNS" in area:
            radar_labels = ["BBB Score", "Solubility", "Oral BA", "Drug-like", "Metabolic", "CNS Safety"]
            radar_values = [
                min(100, max(0, 100 - max(0, psa - 20) * 1.5)),
                min(100, max(0, 100 - max(0, logp - 1) * 20)),
                min(100, max(0, 100 - max(0, mw - 200) * 0.35)),
                min(100, round(qed * 100)),
                min(100, max(0, 100 - rb * 6)),
                min(100, max(0, 90 - ap * 8)),
            ]
        elif "Cardiovascular" in area:
            radar_labels = ["Oral BA", "Solubility", "hERG Safety", "Drug-like", "Metabolic", "Selectivity"]
            radar_values = [
                min(100, max(0, 100 - max(0, mw - 200) * 0.25)),
                min(100, max(0, 100 - max(0, logp - 1) * 15)),
                min(100, max(0, 100 - max(0, logp - 1) * 12)),
                min(100, round(qed * 100)),
                min(100, max(0, 100 - rb * 5)),
                min(100, max(0, 80 - ap * 5)),
            ]
        else:
            radar_labels = ["BBB", "Solubility", "Oral BA", "Drug-like", "Metabolic", "Selectivity"]
            radar_values = [
                min(100, max(0, 100 - (psa - 30) * 1.2)),
                min(100, max(0, 100 - max(0, logp - 1) * 15)),
                min(100, max(0, 100 - max(0, mw - 200) * 0.25)),
                min(100, round(qed * 100)),
                min(100, max(0, 100 - rb * 5)),
                min(100, max(0, 80 - ap * 5)),
            ]
        radar = {"labels": radar_labels, "values": radar_values}

        targets = {
            "mw": profile["mw_ideal"], "logp": profile["logp_ideal"],
            "psa": profile["psa_ideal"], "hbd": profile["hbd_ideal"],
            "hba": profile["hba_ideal"], "rot_bonds": profile["rb_ideal"],
            "qed": profile["qed_ideal"],
        }
        profile_out = {
            "area": profile["area"], "name": profile["name"],
            "mw_max": profile["mw_max"], "logp_hi": profile["logp_hi"],
            "logp_lo": profile["logp_lo"], "psa_max": profile["psa_max"],
            "hbd_max": profile["hbd_max"], "hba_max": profile["hba_max"],
            "rb_max": profile["rb_max"],
            "accept_violations": profile["accept_violations"],
        }
        return jsonify({
            "name":   c.get("name", "Drug"),
            "smiles": smiles,
            "mw": mw, "logp": logp, "psa": psa,
            "hbd": hbd, "hba": hba, "rot_bonds": rb, "qed": qed, "arom_rings": ap,
            "lipinski": lipinski_rules,
            "lipinski_pass_all": lipinski_pass_all,
            "suggestions": suggestions,
            "radar": radar,
            "targets": targets,
            "profile": profile_out,
        })
    except Exception as e:
        return jsonify({"error": str(e)}), 500


@app.route("/api/optimize-smiles", methods=["POST"])
def api_optimize_smiles():
    """Compare original vs modified SMILES: property deltas, safety, repurposing narrative."""
    if not RDKIT_OK:
        return jsonify({"error": "RDKit not available"}), 503
    data = request.get_json() or {}
    orig_smiles  = data.get("original_smiles", "").strip()
    new_smiles   = data.get("smiles", "").strip()
    drug_name    = data.get("drug_name", "compound")
    disease      = data.get("disease", "")
    change_desc  = data.get("change_description", "")

    if not new_smiles:
        return jsonify({"error": "No SMILES provided"}), 400

    def _props(smi):
        from rdkit.Chem import rdMolDescriptors
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            raise ValueError(f"Invalid SMILES")
        mw   = round(Descriptors.MolWt(mol), 1)
        logp = round(Descriptors.MolLogP(mol), 2)
        psa  = round(Descriptors.TPSA(mol), 1)
        hbd  = rdMolDescriptors.CalcNumHBD(mol)
        hba  = rdMolDescriptors.CalcNumHBA(mol)
        rb   = Descriptors.NumRotatableBonds(mol)
        qed  = round(Descriptors.qed(mol), 3)
        bbb = 0
        if psa < 60:   bbb += 30
        elif psa < 90: bbb += 15
        if mw < 400:   bbb += 25
        elif mw < 500: bbb += 10
        if 1 <= logp <= 3:     bbb += 25
        elif 0 <= logp <= 4.5: bbb += 10
        if hbd <= 2:   bbb += 20
        elif hbd <= 3: bbb += 10
        bbb_level = "High" if bbb >= 75 else "Moderate" if bbb >= 45 else "Low"
        ro5 = sum([mw > 500, logp > 5, hbd > 5, hba > 10])
        return {"mw": mw, "logp": logp, "psa": psa, "hbd": hbd, "hba": hba,
                "rb": rb, "qed": qed, "bbb_score": bbb, "bbb_level": bbb_level,
                "ro5_violations": ro5}

    try:
        new_p = _props(new_smiles)
    except ValueError as e:
        return jsonify({"error": str(e)}), 400

    orig_p = {}
    if orig_smiles:
        try:
            orig_p = _props(orig_smiles)
        except ValueError:
            pass

    delta_keys = ["mw", "logp", "psa", "hbd", "hba", "rb", "qed", "bbb_score"]
    deltas = {k: round(new_p[k] - orig_p[k], 3) for k in delta_keys if k in orig_p and k in new_p}

    # Classify each delta as improvement, neutral, or regression
    improvements = []
    for prop, label, good_if_neg, unit in [
        ("psa",       "Polar Surface Area",  True,  " Å²"),
        ("mw",        "Molecular Weight",    True,  " Da"),
        ("hbd",       "H-Bond Donors",       True,  ""),
        ("hba",       "H-Bond Acceptors",    True,  ""),
        ("rb",        "Rotatable Bonds",     True,  ""),
        ("bbb_score", "BBB Score",           False, ""),
        ("qed",       "QED Drug-likeness",   False, ""),
        ("logp",      "LogP",                None,  ""),
    ]:
        if prop not in deltas:
            continue
        d = deltas[prop]
        if prop == "logp":
            new_val = new_p.get("logp", 0)
            if orig_p.get("logp") is not None:
                orig_good = 0 <= orig_p["logp"] <= 4
                new_good  = 0 <= new_val <= 4
                if new_good and not orig_good:
                    direction = "better"
                elif not new_good and orig_good:
                    direction = "worse"
                else:
                    direction = "neutral"
            else:
                direction = "neutral"
        elif prop == "bbb_score" or prop == "qed":
            direction = "better" if d > 2 else "worse" if d < -2 else "neutral"
        else:
            direction = "better" if d < 0 else "worse" if d > 0 else "neutral"
        if direction != "neutral":
            delta_str = f"{d:+.2f}" if isinstance(d, float) else f"{d:+d}"
            improvements.append({
                "property":  label,
                "delta":     delta_str,
                "unit":      unit,
                "direction": direction,
                "old":       orig_p.get(prop),
                "new":       new_p.get(prop),
            })

    safety = [
        {"check": "Lipinski Ro5",    "pass": new_p["ro5_violations"] == 0, "detail": f"{new_p['ro5_violations']} violation(s)"},
        {"check": "BBB Penetration", "pass": new_p["bbb_level"] != "Low",  "detail": new_p["bbb_level"]},
        {"check": "MW ≤ 500 Da",     "pass": new_p["mw"] <= 500,           "detail": f"{new_p['mw']} Da"},
        {"check": "LogP 0–5",        "pass": 0 <= new_p["logp"] <= 5,      "detail": str(new_p["logp"])},
        {"check": "HBD ≤ 5",         "pass": new_p["hbd"] <= 5,            "detail": str(new_p["hbd"])},
    ]

    # Narrative
    parts = []
    if change_desc:
        parts.append(f"By {change_desc.lower().rstrip('.')}")
    elif orig_smiles:
        parts.append("After this structural modification")
    improved_labels = [i["property"] for i in improvements if i["direction"] == "better"]
    worsened_labels = [i["property"] for i in improvements if i["direction"] == "worse"]
    if improved_labels:
        parts.append(f"key pharmacokinetic properties improved: {', '.join(improved_labels[:3])}")
    psa_d = deltas.get("psa", 0)
    bbb_d = deltas.get("bbb_score", 0)
    if psa_d < -5:
        parts.append(f"PSA reduced by {abs(psa_d):.1f} Å² — enhancing CNS membrane permeability")
    if bbb_d > 5:
        bbb_old = orig_p.get("bbb_level", "—")
        parts.append(f"BBB score improved by {int(bbb_d)} points ({bbb_old} → {new_p['bbb_level']} penetration)")
    if disease:
        dis_lo = disease.lower()
        if any(w in dis_lo for w in ("alzheimer", "parkinson", "huntington", "neuro", "depression", "schizo")):
            parts.append(f"supporting repurposing toward {disease} where CNS penetration is a critical requirement")
        else:
            parts.append(f"improving the pharmacokinetic profile for {disease} therapy")
    if worsened_labels:
        parts.append(f"Note: {worsened_labels[0]} increased — monitor this carefully in further optimization")
    narrative = ". ".join(parts).rstrip(".") + "." if parts else "Modified compound shows an altered physicochemical profile."

    # 3D SDF
    sdf = None
    try:
        mol = Chem.MolFromSmiles(new_smiles)
        if mol:
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
            AllChem.MMFFOptimizeMolecule(mol)
            sdf = Chem.MolToMolBlock(mol)
    except Exception:
        pass

    n_better = len([i for i in improvements if i["direction"] == "better"])
    n_worse  = len([i for i in improvements if i["direction"] == "worse"])
    verdict = "Improved" if n_better > n_worse + 1 else "Declined" if n_worse > n_better + 1 else "Mixed"

    return jsonify({
        "original_props": orig_p,
        "new_props":      new_p,
        "deltas":         deltas,
        "improvements":   improvements,
        "safety":         safety,
        "narrative":      narrative,
        "sdf":            sdf,
        "verdict":        verdict,
    })


def _patent_url(pat_id: str, prefix: str = "") -> str:
    """Return a working Google Patents URL for a given patent ID."""
    import re as _re
    pid = (pat_id or "").strip().replace(" ", "")
    if not pid:
        return ""
    # Already has a recognised country prefix → use as-is
    if _re.match(r'^(US|WO|EP|GB|JP|CA|AU|CN|KR|IN)\w', pid, _re.I):
        return f"https://patents.google.com/patent/{pid}/en"
    # Pure number → prepend caller-supplied prefix (e.g. "US") or default search
    if pid.isdigit() and prefix:
        return f"https://patents.google.com/patent/{prefix}{pid}/en"
    # Fallback: Google Patents search by ID
    return f"https://patents.google.com/patent/{pid}/en"


def _patent_expiry(grant_year: str) -> Optional[int]:
    try:
        return int(grant_year) + 20
    except (TypeError, ValueError):
        return None


def _patent_status(grant_year: str) -> str:
    """Estimate patent status. US utility patents expire ~20 yrs from filing (~17-20 from grant)."""
    try:
        y = int(grant_year)
        if y + 20 < 2026:
            return "expired"
        if y + 17 < 2026:
            return "likely_expired"
        return "active"
    except (TypeError, ValueError):
        return "unknown"


@app.route("/api/patents")
def api_patents():
    drug      = request.args.get("drug", "").strip()
    chembl_id = request.args.get("chembl_id", "").strip()
    patents   = []
    exclusivities = []
    seen_ids  = set()

    # 0. FDA Orange Book — authoritative patents + exclusivity with REAL expiry & status
    if drug:
        try:
            from services.orange_book import orange_book_protection
            ob = orange_book_protection(drug)
            for p in ob.get("patents", []):
                if p["id"] in seen_ids:
                    continue
                seen_ids.add(p["id"])
                ttl = (p.get("type", "") or "Patent").title()
                patents.append({
                    "id": p["id"],
                    "title": ttl + " patent" + (f" — {p['trade']}" if p.get("trade") else ""),
                    "assignee": p.get("applicant", ""), "year": (p.get("expiry_iso") or "")[:4],
                    "abstract": "", "url": p.get("url", ""), "source": "FDA Orange Book",
                    "status": p.get("status", "unknown"), "expiry": p.get("expiry"),
                    "expiry_iso": p.get("expiry_iso"), "type": p.get("type", ""),
                    "use_code": p.get("use_code", ""), "authoritative": True,
                })
            exclusivities = ob.get("exclusivities", [])
        except Exception as e:
            logger.debug(f"orange book error: {e}")

    # 1. ChEMBL mechanism → patent refs (most reliable source for drug-specific patents)
    if chembl_id and not chembl_id.startswith("NR:"):
        try:
            r = requests.get("https://www.ebi.ac.uk/chembl/api/data/mechanism.json",
                params={"molecule_chembl_id": chembl_id, "limit": 50, "format": "json"},
                timeout=10)
            if r.ok:
                for m in r.json().get("mechanisms", []):
                    for ref in (m.get("mechanism_refs") or []):
                        if ref.get("ref_type") == "Patent":
                            pid = (ref.get("ref_id") or "").strip()
                            if not pid or pid in seen_ids:
                                continue
                            seen_ids.add(pid)
                            patents.append({
                                "id":       pid,
                                "title":    (m.get("mechanism_of_action") or m.get("target_name") or f"Patent {pid}")[:150],
                                "assignee": "",
                                "year":     "",
                                "abstract": m.get("mechanism_of_action", ""),
                                "url":      _patent_url(pid),
                                "source":   "ChEMBL",
                                "status":   "unknown",
                                "expiry":   None,
                            })
        except Exception:
            pass

    # 2. ChEMBL drug indication patent refs
    if chembl_id and not chembl_id.startswith("NR:"):
        try:
            r = requests.get("https://www.ebi.ac.uk/chembl/api/data/drug_indication.json",
                params={"molecule_chembl_id": chembl_id, "limit": 50, "format": "json"},
                timeout=8)
            if r.ok:
                for ind in r.json().get("drug_indications", []):
                    for ref in (ind.get("indication_refs") or []):
                        if ref.get("ref_type") == "Patent":
                            pid = (ref.get("ref_id") or "").strip()
                            if not pid or pid in seen_ids:
                                continue
                            seen_ids.add(pid)
                            disease_name = ind.get("mesh_heading") or ind.get("efo_term") or ""
                            patents.append({
                                "id":       pid,
                                "title":    (f"Indication patent: {disease_name}" if disease_name else f"Patent {pid}")[:150],
                                "assignee": "",
                                "year":     "",
                                "abstract": f"Indication: {disease_name}",
                                "url":      _patent_url(pid),
                                "source":   "ChEMBL",
                                "status":   "unknown",
                                "expiry":   None,
                            })
        except Exception:
            pass

    # 3. ChEMBL activities from patent documents (SureChEMBL src_id=34)
    if chembl_id and not chembl_id.startswith("NR:"):
        try:
            r = requests.get(
                "https://www.ebi.ac.uk/chembl/api/data/activity.json",
                params={"molecule_chembl_id": chembl_id, "document_src_id": 34,
                        "limit": 20, "format": "json"},
                timeout=10)
            if r.ok:
                doc_ids = []
                for act in r.json().get("activities", []):
                    did = act.get("document_chembl_id", "")
                    if did and did not in doc_ids:
                        doc_ids.append(did)
                for did in doc_ids[:8]:
                    try:
                        dr = requests.get(
                            f"https://www.ebi.ac.uk/chembl/api/data/document/{did}.json",
                            timeout=5)
                        if dr.ok:
                            doc = dr.json()
                            # Only use explicit patent_id — never fall back to journal field
                            # (journal = "J. Med. Chem." for non-patent docs, produces broken URLs)
                            pat_id = (doc.get("patent_id") or "").strip()
                            if not pat_id:
                                continue
                            title  = (doc.get("title") or f"Patent {pat_id}")[:150]
                            year   = str(doc.get("year") or "")[:4]
                            if pat_id not in seen_ids:
                                seen_ids.add(pat_id)
                                patents.append({
                                    "id":       pat_id,
                                    "title":    title,
                                    "assignee": (doc.get("authors") or "")[:100],
                                    "year":     year,
                                    "abstract": "",
                                    "url":      _patent_url(pat_id),
                                    "source":   "SureChEMBL",
                                    "status":   _patent_status(year),
                                    "expiry":   _patent_expiry(year),
                                })
                    except Exception:
                        pass
        except Exception:
            pass

    # 4. PatentsView (USPTO) — uses v3 API (launched 2023)
    if drug:
        try:
            pv_r = requests.post(
                "https://api.patentsview.org/patents/query",
                timeout=12,
                headers={"Content-Type": "application/json"},
                json={
                    "q": {"_text_any": {"patent_title": drug}},
                    "f": ["patent_id", "patent_title", "patent_date",
                          "assignee_organization", "patent_abstract"],
                    "o": {"per_page": 20},
                    "s": [{"patent_date": "desc"}],
                })
            pv_list = []
            if pv_r.ok:
                body = pv_r.json()
                # v3 API uses "patents" key; v1 fallback also uses "patents"
                pv_list = body.get("patents") or []
            for p in pv_list:
                # v3 uses "patent_id", v1 used "patent_number"
                pid = (p.get("patent_id") or p.get("patent_number") or "").strip()
                if not pid or pid in seen_ids:
                    continue
                seen_ids.add(pid)
                # assignee_organization is a direct field in v3; in v1 it was nested
                org = (p.get("assignee_organization") or
                       "; ".join(a.get("assignee_organization","")
                                 for a in (p.get("assignees") or [])
                                 if a.get("assignee_organization")))[:100]
                year = (p.get("patent_date","") or "")[:4]
                patents.append({
                    "id":       pid,
                    "title":    (p.get("patent_title") or p.get("patent_title_text") or "")[:150],
                    "assignee": org,
                    "year":     year,
                    "abstract": (p.get("patent_abstract") or "")[:300],
                    "url":      _patent_url(pid, prefix="US"),
                    "source":   "USPTO",
                    "status":   _patent_status(year),
                    "expiry":   _patent_expiry(year),
                })
        except Exception:
            pass

    # Sort: authoritative (Orange Book) first, then active before expired, then year desc
    _order = {"active": 0, "likely_expired": 1, "expired": 2, "unknown": 3}
    patents.sort(key=lambda p: (0 if p.get("authoritative") else 1,
                                _order.get(p.get("status", "unknown"), 3),
                                -(int(p["year"]) if str(p.get("year", "")).isdigit() else 0)))
    return jsonify(_scrub_sources({
        "patents": patents, "exclusivities": exclusivities,
        "summary": {
            "total": len(patents),
            "active": sum(1 for p in patents if p.get("status") == "active"),
            "expired": sum(1 for p in patents if p.get("status") in ("expired", "likely_expired")),
            "orange_book": sum(1 for p in patents if p.get("authoritative")),
        },
    }))


@app.route("/api/repurposing-screen")
def api_repurposing_screen():
    disease = request.args.get("disease", "").strip()
    if not disease:
        return jsonify({"error": "disease parameter required"}), 400
    if not ENGINE_OK:
        return jsonify({"error": "Repurposing engine not available"}), 503
    try:
        resolved, expanded, db_comps = _resolve_fetch(disease)
        screen = run_repurposing_screen(disease, max_candidates=50, db_compounds=db_comps)
        screen.pop("_ts", None)
        if "candidates" in screen:
            screen["candidates"] = [c for c in screen["candidates"] if _is_repurposable(c)]
        return jsonify(screen)
    except Exception as e:
        return jsonify({"error": str(e)}), 500


@app.route("/api/compound/<chembl_id>/class-comparison")
def api_class_comparison(chembl_id):
    """Within-class ('why this drug') + out-of-class ('why this mechanism')
    differentiation for a (drug, disease) repurposing hypothesis."""
    disease = (request.args.get("disease") or "").strip()
    drug = (request.args.get("drug") or request.args.get("name") or "").strip()
    if not disease:
        return jsonify({"error": "disease parameter required"}), 400
    try:
        from services.class_comparison import compare_classes
        return jsonify(compare_classes(chembl_id, disease, drug_name=drug,
                                       smiles=(request.args.get("smiles") or "")))
    except Exception as e:
        logger.exception("class-comparison error")
        return jsonify({"error": str(e)}), 500


@app.route("/api/compound/<chembl_id>/dossier")
def api_dossier(chembl_id):
    disease = request.args.get("disease", "").strip()
    if not disease:
        return jsonify({"error": "disease parameter required"}), 400
    c = _get_compound(chembl_id)
    # Fallback: frontend passes name/smiles it already has in session
    if not c:
        fb_name   = request.args.get("name", "").strip()
        fb_smiles = request.args.get("smiles", "").strip()
        fb_phase  = _safe_float(request.args.get("max_phase", "0"), 0)
        if fb_name or fb_smiles:
            c = {"chembl_id": chembl_id, "name": fb_name or chembl_id,
                 "smiles": fb_smiles, "max_phase": int(fb_phase),
                 "indications": "", "targets": ""}
    if not c:
        return jsonify({"error": f"Compound {chembl_id} not found"}), 404

    # Only run repurposing screen when compound has a DB id (avoids slow failure for API compounds)
    screen = None
    if ENGINE_OK and c.get("id"):
        try:
            resolved, expanded, db_comps = _resolve_fetch(disease)
            screen = run_repurposing_screen(disease, max_candidates=50, db_compounds=db_comps)
        except Exception:
            pass

    dossier = generate_dossier(chembl_id, disease, compound=c, screen_result=screen)
    # Tag whether this molecule is in the company portfolio (505(b)(2) reference base)
    if PORTFOLIO_OK:
        try:
            entry = _portfolio.match(chembl_id=chembl_id, name=(c or {}).get("name", ""))
            if entry:
                dossier["portfolio"] = {
                    "in_portfolio": True,
                    "status": entry.get("status", ""),
                    "form":   entry.get("form", ""),
                    "note":   entry.get("note", ""),
                }
        except Exception:
            pass

    drug_name = dossier.get("drug_name") or c.get("name", "")
    phase = float(c.get("max_phase") or 0)

    # Asset Scoring — phase-progression Probability of Success (analytic survival
    # model; the composite score is the new-indication efficacy evidence).
    try:
        from services.pos_model import predict_progression
        dossier["pos"] = predict_progression(
            current_phase=phase,
            evidence_score=dossier.get("composite_score") or 0.0,
            is_repurposing=phase >= 4,
        )
    except Exception as _e:
        logger.debug(f"dossier PoS failed: {_e}")

    # Asset Acquisition Signal — FDA Orange Book (patents, exclusivity, generics).
    try:
        from services.acquisition_signal import acquisition_signal
        dossier["acquisition"] = acquisition_signal(drug_name)
    except Exception as _e:
        logger.debug(f"dossier acquisition failed: {_e}")

    # Signature-based connectivity — does the drug REVERSE the disease signature?
    # Disease genes from the dossier; drug target directions from ChEMBL mechanisms.
    try:
        from services.signature_engine import score_reversal
        disease_genes = (dossier.get("disease_context") or {}).get("top_genes") or []
        drug_targets = _fetch_chembl_mechanism_rows(chembl_id)
        if not drug_targets and c.get("id"):
            # Local DB mechanisms carry gene_symbol + action_type — no network needed.
            drug_targets = get_compound_targets(c.get("id"))
        if disease_genes and drug_targets:
            dossier["signature"] = score_reversal(disease_genes, drug_targets)
    except Exception as _e:
        logger.debug(f"dossier signature failed: {_e}")

    # Reconciled regulatory + novelty verdict (single source of truth) and the
    # consolidated US/EU/India landscape. known_indications carry max_phase so
    # novelty is phase-aware (approved-here vs merely studied/co-mentioned).
    known_inds = []
    try:
        if REVERSE_OK and drug_name:
            known_inds = (resolve_drug(drug_name) or {}).get("known_indications", []) or []
    except Exception as _e:
        logger.debug(f"resolve_drug for verdict failed: {_e}")
    try:
        from services.regulatory_verdict import assess as _assess
        dossier["verdict"] = _assess(drug_name, phase, known_inds, dossier.get("disease_name") or disease)
    except Exception as _e:
        logger.debug(f"dossier verdict failed: {_e}")
    try:
        from services.global_landscape import landscape as _landscape
        dossier["global_landscape"] = _landscape(drug_name, phase)
    except Exception as _e:
        logger.debug(f"dossier landscape failed: {_e}")

    return jsonify(dossier)


@app.route("/api/drug/<path:drug_id>/new-indications")
def api_new_indications(drug_id):
    """
    Reverse repurposing: given a molecule (name or ChEMBL ID), return ranked NEW
    candidate indications derived dynamically from its protein targets.
    Optional ?area=dermatology|ophthalmology filters by Open Targets therapeutic area.
    """
    if not REVERSE_OK:
        return jsonify({"error": "Reverse repurposing engine unavailable", "candidates": []}), 503
    area = request.args.get("area", "").strip() or None
    # exclude_oncology defaults to True; pass ?oncology=1 to include cancer indications
    exclude_oncology = request.args.get("oncology", "0").strip() not in ("1", "true", "yes")
    try:
        result = screen_indications_for_drug(drug_id, area_filter=area,
                                             exclude_oncology=exclude_oncology)
        # Mark the source molecule's portfolio status for the 505(b)(2) framing
        if PORTFOLIO_OK:
            entry = _portfolio.match(chembl_id=result.get("chembl_id", ""),
                                     name=result.get("drug", ""))
            result["portfolio"] = (
                {"in_portfolio": True, "status": entry.get("status", ""),
                 "form": entry.get("form", ""), "note": entry.get("note", "")}
                if entry else {"in_portfolio": False}
            )
        return jsonify(result)
    except Exception as e:
        logger.exception("new-indications error")
        return jsonify({"error": str(e), "candidates": []}), 500


@app.route("/repurpose")
def repurpose():
    """Drug → new indications page (reverse repurposing)."""
    drug = request.args.get("drug", "").strip()
    area = request.args.get("area", "").strip()
    return render_template("repurpose.html", drug=drug, area=area,
                           reverse_ok=REVERSE_OK, portfolio_ok=PORTFOLIO_OK, db_ok=DB_OK)


@app.route("/pathways")
def pathways():
    """Pathway-first repurposing page (pathway → drugs → new indications)."""
    pathway = request.args.get("pathway", "").strip()
    direction = request.args.get("direction", "either").strip()
    disease = request.args.get("disease", "").strip()
    return render_template("pathways.html", pathway=pathway, direction=direction,
                           disease=disease, pathway_ok=PATHWAY_OK, db_ok=DB_OK)


@app.route("/api/suggest-pathways")
def api_suggest_pathways():
    """Autocomplete for pathway names."""
    q = request.args.get("q", "").strip()
    if not q or not PATHWAY_OK:
        return jsonify([])
    try:
        return jsonify(suggest_pathways(q, limit=12))
    except Exception as e:
        logger.error(f"suggest-pathways: {e}")
        return jsonify([])


@app.route("/api/pathway-screen")
def api_pathway_screen():
    """Pathway-first repurposing screen.

    Query params:
      pathway   (required) — pathway name or id
      direction — suppress | activate | either   (desired effect on the pathway)
      disease   — optional; supplying it switches to Mode A (disease-anchored)
    """
    if not PATHWAY_OK:
        return jsonify({"error": "Pathway engine unavailable", "candidates": []}), 503
    pathway = request.args.get("pathway", "").strip()
    if not pathway:
        return jsonify({"error": "No pathway specified", "candidates": []}), 400
    direction = request.args.get("direction", "either").strip()
    disease = request.args.get("disease", "").strip() or None
    try:
        res = screen_pathway(pathway, direction=direction, disease=disease)
        return jsonify(res)
    except Exception as e:
        logger.error(f"pathway-screen '{pathway}': {e}")
        return jsonify({"error": str(e), "candidates": []}), 500


@app.route("/api/compound/<chembl_id>/developability")
def api_developability(chembl_id):
    """Area-aware developability for a compound vs the searched disease's route.
    Picks the profile (skin/ocular/CNS/oral) from the disease's therapeutic area
    instead of always scoring CNS/BBB."""
    disease = request.args.get("disease", "").strip()
    c = _get_compound(chembl_id)
    smiles = (c or {}).get("smiles", "") or request.args.get("smiles", "")
    if not smiles:
        return jsonify({"available": False})
    try:
        from services import developability as _dev
        areas = _disease_areas(disease) if disease else []
        return jsonify(_dev.score(smiles, area="", therapeutic_areas=areas))
    except Exception as e:
        logger.debug(f"developability endpoint error: {e}")
        return jsonify({"available": False, "error": str(e)})


@app.route("/api/compound/<chembl_id>/docking-poses", methods=["POST"])
def api_docking_poses(chembl_id):
    """Run docking and return SDF pose strings for 3Dmol.js rendering."""
    if not _dock_svc:
        return jsonify({"error": "Docking service unavailable"}), 503
    data = request.get_json() or {}
    c = _get_compound(chembl_id)
    if not c:
        return jsonify({"error": f"Compound {chembl_id} not found"}), 404
    try:
        target_name = data.get("target", "BACE1")
        method = data.get("method") if data.get("method") in ("diffdock", "boltz") else "fast"
        ligand_smiles = data.get("smiles") or c.get("smiles", "")
        # Cache identical (compound, target, ligand, method) runs — instant repeat
        cache_key = f"{chembl_id}|{target_name}|{method}|{hash(ligand_smiles)}"
        cached = _dock_result_cache.get(cache_key)
        if cached is not None:
            return jsonify(cached)
        # Reuse protein already fetched by /api/protein-pdb to avoid a 30-45s re-fetch
        cached_pdb = _protein_pdb_cache.get(target_name, "")
        res = _dock_svc.perform_docking(
            drug_name=data.get("name") or c.get("name", "Drug"),
            target_name=target_name,
            ligand_smiles=ligand_smiles,
            protein_pdb=cached_pdb or None,
            method=method,
        )
        # Ensure poses are SDF strings (already captured by local_diffdock)
        poses = res.get("poses", [])
        confs = res.get("confidence_scores", [])
        affs  = res.get("binding_affinities", [])
        # Include protein PDB in response so frontend can render protein+ligand complex
        protein_for_viewer = (
            res.get("protein_pdb")          # cleaned PDB from DiffDock run
            or cached_pdb                    # or the pre-fetched cached PDB
            or ""
        )
        payload = {
            "success":  res.get("success", False),
            "method":   res.get("docking_method", "Unknown"),
            "target":   res.get("target_name", ""),
            "drug":     res.get("drug_name", ""),
            "poses":    poses,
            "confidence_scores":  confs,
            "binding_affinities": affs,
            "protein_pdb": protein_for_viewer,
            "pockets":      res.get("pockets", []),
            "pose_pockets": res.get("pose_pockets", []),
            "pose_valid":   res.get("pose_valid", []),
            "pocket_source": res.get("pocket_source", ""),
            "structure_quality": res.get("structure_quality", "real"),
            "structure_warning": res.get("structure_warning", ""),
            "empirical_affinities": res.get("empirical_affinities", []),
            "boltz_affinity_pred": res.get("boltz_affinity_pred"),
            "boltz_binder_probability": res.get("boltz_binder_probability"),
            "consensus_note": res.get("consensus_note", ""),
            "note":     res.get("note", ""),
            "error":    res.get("error", ""),
        }
        if payload["success"]:
            _dock_result_cache[cache_key] = payload
        return jsonify(payload)
    except Exception as e:
        return jsonify({"error": str(e)}), 500


@app.route("/api/transform-smiles", methods=["POST"])
def api_transform_smiles():
    """Apply a named medicinal-chemistry transform to SMILES via RDKit SMARTS reactions."""
    data = request.get_json() or {}
    smiles = data.get("smiles", "").strip()
    action = data.get("action", "").lower()
    if not smiles:
        return jsonify({"error": "No SMILES"}), 400
    if not RDKIT_OK:
        return jsonify({"smiles": smiles, "applied": False, "note": "RDKit not available"})
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return jsonify({"error": "Invalid SMILES"}), 400

    # (keyword, SMARTS reaction, human description of what was done)
    TRANSFORMS = [
        ("n-methylat",              "[NH1;$(N-C=O):1]>>[N:1]C",  "N-methylated amide N-H → N-Me (reduces PSA, improves BBB)"),
        ("methylat amide",          "[NH1;$(N-C=O):1]>>[N:1]C",  "N-methylated amide N-H → N-Me (reduces PSA, improves BBB)"),
        ("cap amide",               "[NH1;$(N-C=O):1]>>[N:1]C",  "Capped amide N-H with methyl (reduces PSA)"),
        ("cap polar",               "[NH1;$(N-C=O):1]>>[N:1]C",  "Capped polar amide N-H (reduces PSA)"),
        ("non-polar bioisostere",   "[c:1][OH1]>>[c:1]F",        "Replaced phenol OH with F bioisostere (reduces PSA)"),
        ("replace polar",           "[c:1][OH1]>>[c:1]F",        "Replaced phenol OH with F bioisostere (reduces PSA)"),
        ("remove hydroxyl",         "[c:1][OH1]>>[c:1]F",        "Replaced phenol OH with F bioisostere (reduces PSA)"),
        ("remove oh",               "[c:1][OH1]>>[c:1]F",        "Replaced phenol OH with F bioisostere (reduces PSA)"),
        ("add methyl",              "[c;H1:1]>>[c:1]C",          "Added methyl group to aromatic ring (increases LogP)"),
        ("methyl group",            "[c;H1:1]>>[c:1]C",          "Added methyl group to aromatic ring"),
        ("remove non-essential",    "[c;H1:1]>>[c:1]",           "Simplified aromatic substitution"),
    ]

    rxn_smarts = None
    description = ""
    for keyword, smarts, desc in TRANSFORMS:
        if keyword in action:
            rxn_smarts = smarts
            description = desc
            break

    if not rxn_smarts:
        return jsonify({"smiles": smiles, "applied": False,
                        "note": "No automatic transform for this action — edit the Modified SMILES manually."})
    try:
        from rdkit.Chem import AllChem
        rxn = AllChem.ReactionFromSmarts(rxn_smarts)
        products = rxn.RunReactants((mol,))
        if not products:
            return jsonify({"smiles": smiles, "applied": False,
                            "note": "No matching groups found — edit Modified SMILES manually."})
        prod = products[0][0]
        Chem.SanitizeMol(prod)
        new_smi = Chem.MolToSmiles(prod)
        if new_smi == Chem.MolToSmiles(mol):
            return jsonify({"smiles": smiles, "applied": False,
                            "note": "Transform left molecule unchanged — edit Modified SMILES manually."})
        return jsonify({"smiles": new_smi, "applied": True, "description": description})
    except Exception as e:
        return jsonify({"smiles": smiles, "applied": False, "note": f"Transform error: {e}"})


_protein_pdb_cache: dict = {}
_dock_result_cache: dict = {}

@app.route("/api/protein-pdb")
def api_protein_pdb():
    """Fetch protein PDB for a target gene symbol (for 3D docking viewer)."""
    target = request.args.get("target", "").strip()
    if not target:
        return jsonify({"error": "target required"}), 400
    if target not in _protein_pdb_cache:
        pdb = ""
        try:
            from real_pdb_fetcher import RealPDBFetcher
            fetcher = RealPDBFetcher()
            pdb = fetcher.fetch_protein_structure(target) or ""
        except Exception as e:
            logger.debug(f"api_protein_pdb {target}: {e}")
        _protein_pdb_cache[target] = pdb
    pdb = _protein_pdb_cache[target]
    is_real = bool(pdb and "GENERIC PROTEIN" not in pdb and len(pdb) > 500)
    # Truncate very large PDBs to ATOM/HETATM records only (keep file compact for transfer)
    if is_real and len(pdb) > 300_000:
        lines = [ln for ln in pdb.splitlines() if ln.startswith(("ATOM","HETATM","TER","END"))]
        pdb = "\n".join(lines)
    return jsonify({"pdb": pdb if is_real else "", "real": is_real, "target": target})


_gene_info_cache: dict = {}

@app.route("/api/gene-info")
def api_gene_info():
    """
    Return gene-level evidence for the Drug→Gene→Disease chain shown in the Docking tab.
    Aggregates: OT association score, Reactome pathway, Open Targets function summary,
    and whether this compound is a known ChEMBL inhibitor of the gene.
    """
    gene    = request.args.get("gene", "").strip().upper()
    disease = request.args.get("disease", "").strip()
    chembl  = request.args.get("chembl_id", "").strip()
    if not gene:
        return jsonify({"error": "gene required"}), 400

    cache_key = f"{gene}|{disease}|{chembl}"
    if cache_key in _gene_info_cache:
        return jsonify(_gene_info_cache[cache_key])

    # ── Static pathway map (instant) ─────────────────────────────────────────
    from services.biocypher_graph import _PATHWAY_MAP
    pathway_name, pathway_id = _PATHWAY_MAP.get(gene, ("", ""))

    # ── Open Targets: disease→gene association score ──────────────────────────
    ot_score   = 0.0
    gene_name  = ""
    gene_func  = ""
    ot_disease = disease

    if disease:
        try:
            from services.disease_ontology import resolve_disease
            dis_info = resolve_disease(disease)
            ot_disease = dis_info.get("disease_name", disease)
            for t in dis_info.get("targets", []):
                if t.get("gene_symbol", "").upper() == gene:
                    ot_score  = round(t.get("score", 0), 4)
                    gene_name = t.get("gene_name", "")
                    break
        except Exception:
            pass

    # ── Open Targets: gene function summary ───────────────────────────────────
    if not gene_func:
        try:
            ot_q = """
            query($sym: String!) {
              targets(ensemblIds: [], queryString: $sym, page: {index: 0, size: 1}) {
                rows { id approvedName functionDescriptions }
              }
            }"""
            r = requests.post(
                "https://api.platform.opentargets.org/api/v4/graphql",
                json={"query": ot_q, "variables": {"sym": gene}},
                timeout=8, headers={"Content-Type": "application/json"})
            if r.ok:
                rows = r.json().get("data", {}).get("targets", {}).get("rows", [])
                if rows:
                    gene_name  = gene_name or rows[0].get("approvedName", "")
                    funcs      = rows[0].get("functionDescriptions") or []
                    gene_func  = (funcs[0] if funcs else "")[:300]
        except Exception:
            pass

    # ── ChEMBL: is this compound a known inhibitor / binder of this gene? ─────
    known_action = ""
    known_mech   = ""
    if chembl and not chembl.startswith("NR:"):
        try:
            r = requests.get("https://www.ebi.ac.uk/chembl/api/data/mechanism.json",
                params={"molecule_chembl_id": chembl, "limit": 50, "format": "json"},
                timeout=8)
            if r.ok:
                for m in r.json().get("mechanisms", []):
                    from services.biocypher_graph import _get_chembl_target
                    tdata = _get_chembl_target(m.get("target_chembl_id", ""))
                    if tdata:
                        for comp in (tdata.get("target_components") or []):
                            for syn in (comp.get("target_component_synonyms") or []):
                                if (syn.get("syn_type") == "GENE_SYMBOL"
                                        and syn.get("component_synonym","").upper() == gene):
                                    known_action = (m.get("action_type") or "").replace("_"," ").title()
                                    known_mech   = (m.get("mechanism_of_action") or "")[:200]
                                    break
                            if known_action:
                                break
                    if known_action:
                        break
        except Exception:
            pass

    result = {
        "gene":         gene,
        "gene_name":    gene_name,
        "gene_function": gene_func,
        "ot_score":     ot_score,
        "ot_disease":   ot_disease,
        "pathway_name": pathway_name,
        "pathway_id":   pathway_id,
        "reactome_url": f"https://reactome.org/PathwayBrowser/#/{pathway_id}" if pathway_id else "",
        "known_action": known_action,
        "known_mech":   known_mech,
        "is_known_target": bool(known_action),
    }
    _gene_info_cache[cache_key] = result
    return jsonify(result)


# ── Entry point ──────────────────────────────────────────────────────────────
if __name__ == "__main__":
    app.run(debug=False, host="127.0.0.1", port=5000, use_reloader=False, threaded=True)
