#!/usr/bin/env python3
"""
NeuroRepurpose Intelligence Platform
ChEMBL 33 · HetioNet · MeSH ontology · PBPK · Docking · Clinical Evidence
"""

import logging
import os
from pathlib import Path
from typing import Dict, List, Optional

import networkx as nx
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import requests
import streamlit as st

st.set_page_config(
    page_title="NeuroRepurpose",
    layout="wide",
    initial_sidebar_state="expanded",
)

logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)
REPO_ROOT = Path(__file__).parent

# ─── Core services ────────────────────────────────────────────────────────────

try:
    from services.neuro_db_service import (
        get_compounds_for_disease, search_compounds as db_search,
        get_compound_targets, get_compound_activities,
        get_compound_indications, get_hetionet_paths,
        get_stats, is_available as db_is_available,
        get_disease_top_compounds, get_available_diseases,
    )
    DB_AVAILABLE = db_is_available()
except Exception as _e:
    DB_AVAILABLE = False
    def get_compounds_for_disease(*a, **k): return []
    def db_search(*a, **k): return []
    def get_compound_targets(*a, **k): return []
    def get_compound_activities(*a, **k): return []
    def get_compound_indications(*a, **k): return []
    def get_hetionet_paths(*a, **k): return []
    def get_stats(): return {}
    def get_disease_top_compounds(*a, **k): return []
    def get_available_diseases(*a, **k): return []

try:
    from services.disease_resolver import (
        resolve_disease, expand_mesh_ids, suggest_diseases, mesh_available,
    )
    MESH_AVAILABLE = mesh_available() if DB_AVAILABLE else False
except Exception:
    MESH_AVAILABLE = False
    def resolve_disease(*a, **k): return []
    def expand_mesh_ids(*a, **k): return []
    def suggest_diseases(*a, **k): return []

try:
    from services.repurposing_scorer import score_compound_list
    SCORER_AVAILABLE = True
except Exception:
    SCORER_AVAILABLE = False
    def score_compound_list(c, m):
        for x in c: x.setdefault("score", float(x.get("max_phase") or 0)/4); x.setdefault("score_breakdown", {})
        c.sort(key=lambda x: x["score"], reverse=True); return c

try:
    from services.compound_validator import validate_and_deduplicate
except Exception:
    def validate_and_deduplicate(c, **k): return c

try:
    from services.pbpk_simulation import PBPKSimulator
    PBPK_AVAILABLE = True
except Exception:
    PBPK_AVAILABLE = False
    PBPKSimulator = None

try:
    from services.docking_service import DockingService
    _docking_svc = DockingService()
    DOCKING_AVAILABLE = _docking_svc.available
except Exception:
    DOCKING_AVAILABLE = False
    _docking_svc = None

try:
    from quantum_optimization_strategies import render_quantum_optimization_section
    QUANTUM_AVAILABLE = True
except Exception:
    QUANTUM_AVAILABLE = False

try:
    from drug_protein_3d_visualizer import DrugProtein3DVisualizer
    _viz = DrugProtein3DVisualizer()
    VIZ_3D_AVAILABLE = True
except Exception:
    VIZ_3D_AVAILABLE = False
    _viz = None


# ─── External data (Clinical Trials · PubMed) ─────────────────────────────────

@st.cache_data(ttl=3600, show_spinner=False)
def fetch_clinical_trials(drug_name: str, disease: str, n: int = 8) -> List[Dict]:
    try:
        r = requests.get(
            "https://clinicaltrials.gov/api/v2/studies",
            params={"query.cond": disease, "query.term": drug_name,
                    "pageSize": n, "format": "json"},
            timeout=12,
        )
        if r.status_code == 200:
            out = []
            for s in r.json().get("studies", []):
                pm = s.get("protocolSection", {})
                phases = pm.get("designModule", {}).get("phases", [])
                out.append({
                    "nct_id":  pm.get("identificationModule", {}).get("nctId", ""),
                    "title":   pm.get("identificationModule", {}).get("briefTitle", ""),
                    "status":  pm.get("statusModule", {}).get("overallStatus", ""),
                    "phase":   ", ".join(phases) if phases else "N/A",
                    "url":     f"https://clinicaltrials.gov/study/{pm.get('identificationModule',{}).get('nctId','')}",
                })
            return out
    except Exception:
        pass
    return []


@st.cache_data(ttl=3600, show_spinner=False)
def fetch_pubmed(drug_name: str, disease: str, n: int = 8) -> List[Dict]:
    try:
        base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
        q = f"{drug_name}[tiab] AND {disease}[tiab]"
        s = requests.get(f"{base}/esearch.fcgi",
                         params={"db": "pubmed", "term": q, "retmax": n, "retmode": "json"},
                         timeout=10)
        ids = s.json().get("esearchresult", {}).get("idlist", [])
        if not ids:
            return []
        sm = requests.get(f"{base}/esummary.fcgi",
                          params={"db": "pubmed", "id": ",".join(ids), "retmode": "json"},
                          timeout=10)
        res = sm.json().get("result", {})
        return [
            {"pmid": p, "title": res[p].get("title", ""), "journal": res[p].get("source", ""),
             "year": res[p].get("pubdate", "")[:4], "url": f"https://pubmed.ncbi.nlm.nih.gov/{p}/"}
            for p in ids if p in res
        ]
    except Exception:
        pass
    return []


# ─── Cached discovery pipeline ────────────────────────────────────────────────

@st.cache_data(ttl=1800, show_spinner=False)
def resolve_and_fetch(query: str):
    resolved = resolve_disease(query)
    if not resolved:
        return [], [], []
    mesh_ids = [r["mesh_id"] for r in resolved if r.get("mesh_id")]
    expanded = expand_mesh_ids(mesh_ids, include_children=True) or mesh_ids
    compounds = get_compounds_for_disease(expanded, limit=80)
    compounds = validate_and_deduplicate(compounds, require_smiles=False)
    return resolved, expanded, compounds


# ─── Theme ────────────────────────────────────────────────────────────────────

def apply_theme():
    st.markdown("""
    <style>
    /* ── Fonts ── */
    * { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', 'Inter', sans-serif !important; }
    code, .stCode * { font-family: 'JetBrains Mono', 'Fira Code', 'Cascadia Code', monospace !important; }

    /* ── Design tokens ── */
    :root {
        --bg:         #060b18;
        --card:       #0b1527;
        --raised:     #0f1e38;
        --surface:    #142040;
        --border:     #1c304f;
        --border-hi:  #2a4570;
        --cyan:       #0ea5e9;
        --cyan-dim:   rgba(14,165,233,0.12);
        --violet:     #8b5cf6;
        --violet-dim: rgba(139,92,246,0.12);
        --emerald:    #10b981;
        --amber:      #f59e0b;
        --rose:       #f43f5e;
        --text:       #eef2ff;
        --muted:      #6b7d9f;
        --dim:        #2a3d5c;
        --grad:       linear-gradient(135deg, #0ea5e9, #8b5cf6);
    }

    /* ── Global ── */
    .stApp, .main, [data-testid="stAppViewContainer"]  { background: var(--bg) !important; }
    [data-testid="stSidebar"]  { background: var(--card) !important; border-right: 1px solid var(--border); }
    .block-container           { padding: 1.5rem 2.5rem !important; max-width: 1500px; }
    section[data-testid="stSidebar"] > div { padding-top: 1rem; }

    /* ── Typography ── */
    h1,h2,h3,h4,h5,h6 { color: var(--text) !important; letter-spacing: -0.02em; }
    p, span, label, div, li { color: var(--text) !important; }
    .stMarkdown p { line-height: 1.65; color: var(--muted) !important; }
    small { color: var(--muted) !important; font-size: 0.78rem; }

    /* ── Inputs ── */
    .stTextInput input, .stTextArea textarea, .stNumberInput input {
        background: var(--raised) !important;
        border: 1px solid var(--border) !important;
        border-radius: 8px !important;
        color: var(--text) !important;
        font-size: 0.92rem !important;
        transition: border-color .2s, box-shadow .2s;
    }
    .stTextInput input:focus {
        border-color: var(--cyan) !important;
        box-shadow: 0 0 0 3px var(--cyan-dim) !important;
        outline: none !important;
    }
    .stSelectbox > div > div {
        background: var(--raised) !important;
        border-color: var(--border) !important;
        color: var(--text) !important;
        border-radius: 8px !important;
    }

    /* ── Buttons — primary ── */
    .stButton > button {
        background: var(--grad) !important;
        color: #fff !important;
        border: none !important;
        border-radius: 8px !important;
        font-weight: 600 !important;
        font-size: 0.85rem !important;
        padding: 0.5rem 1.4rem !important;
        letter-spacing: 0.01em;
        transition: opacity .18s, transform .12s, box-shadow .18s !important;
    }
    .stButton > button:hover {
        opacity: 0.85 !important;
        transform: translateY(-1px) !important;
        box-shadow: 0 4px 20px rgba(14,165,233,0.25) !important;
    }
    .stButton > button[kind="secondary"] {
        background: var(--raised) !important;
        border: 1px solid var(--border-hi) !important;
        color: var(--muted) !important;
    }
    .stButton > button[kind="secondary"]:hover {
        border-color: var(--cyan) !important;
        color: var(--cyan) !important;
        box-shadow: none !important;
        opacity: 1 !important;
    }

    /* ── Metrics ── */
    [data-testid="stMetric"] {
        background: var(--raised);
        border: 1px solid var(--border);
        border-radius: 12px;
        padding: 1rem 1.25rem;
    }
    [data-testid="stMetricValue"] { color: var(--cyan) !important; font-size: 1.6rem !important; font-weight: 700 !important; }
    [data-testid="stMetricLabel"] { color: var(--muted) !important; font-size: 0.78rem !important; text-transform: uppercase; letter-spacing: 0.06em; }

    /* ── Tabs ── */
    .stTabs [data-baseweb="tab-list"] {
        gap: 0;
        border-bottom: 1px solid var(--border) !important;
        background: transparent !important;
    }
    .stTabs [data-baseweb="tab"] {
        background: transparent !important;
        color: var(--muted) !important;
        border: none !important;
        border-radius: 0 !important;
        padding: 0.6rem 1.25rem !important;
        font-size: 0.85rem !important;
        font-weight: 500 !important;
    }
    .stTabs [aria-selected="true"] {
        color: var(--cyan) !important;
        border-bottom: 2px solid var(--cyan) !important;
        background: transparent !important;
    }

    /* ── DataFrame ── */
    [data-testid="stDataFrame"] { border-radius: 10px; overflow: hidden; border: 1px solid var(--border); }
    [data-testid="stDataFrame"] thead th { background: var(--raised) !important; color: var(--cyan) !important; font-size: 0.78rem !important; text-transform: uppercase; letter-spacing: 0.05em; }
    [data-testid="stDataFrame"] tbody td { font-size: 0.84rem !important; color: var(--text) !important; }
    [data-testid="stDataFrame"] tbody tr:nth-child(even) td { background: var(--raised) !important; }

    /* ── Progress / spinner ── */
    .stProgress > div > div { background: var(--grad) !important; border-radius: 4px; }
    .stSpinner > div { border-top-color: var(--cyan) !important; }
    [data-testid="stStatusWidget"] { display: none; }

    /* ── Expander ── */
    .streamlit-expanderHeader { background: var(--raised) !important; border-radius: 8px; color: var(--muted) !important; font-size: 0.85rem !important; }
    .streamlit-expanderContent { background: var(--card) !important; border: 1px solid var(--border); border-top: none; border-radius: 0 0 8px 8px; }

    hr { border-color: var(--border) !important; margin: 1.25rem 0; }

    /* ════════════════ DESIGN COMPONENTS ════════════════ */

    /* Hero */
    .nr-hero {
        background: linear-gradient(160deg, #0a1628 0%, #0f1e38 55%, #091424 100%);
        border: 1px solid var(--border);
        border-radius: 20px;
        padding: 2.75rem 3rem;
        margin-bottom: 2rem;
        position: relative;
        overflow: hidden;
    }
    .nr-hero::after {
        content: '';
        position: absolute; top: -100px; right: -80px;
        width: 320px; height: 320px; border-radius: 50%;
        background: radial-gradient(circle, rgba(14,165,233,0.06) 0%, transparent 65%);
        pointer-events: none;
    }
    .nr-hero::before {
        content: '';
        position: absolute; bottom: -60px; left: 30%;
        width: 220px; height: 220px; border-radius: 50%;
        background: radial-gradient(circle, rgba(139,92,246,0.05) 0%, transparent 65%);
        pointer-events: none;
    }
    .nr-hero-label {
        font-size: 0.68rem; font-weight: 600; letter-spacing: 0.18em;
        text-transform: uppercase; color: var(--cyan) !important;
        margin-bottom: 0.75rem;
    }
    .nr-hero-title {
        font-size: 2.2rem; font-weight: 800; line-height: 1.15;
        background: linear-gradient(90deg, #f0f4ff 30%, #0ea5e9 70%, #8b5cf6 100%);
        -webkit-background-clip: text; -webkit-text-fill-color: transparent;
        background-clip: text; margin-bottom: 0.75rem;
    }
    .nr-hero-sub { font-size: 0.95rem; color: var(--muted) !important; line-height: 1.6; max-width: 680px; }

    /* Standard card */
    .nr-card {
        background: var(--card);
        border: 1px solid var(--border);
        border-radius: 14px;
        padding: 1.25rem 1.5rem;
        margin-bottom: 0.75rem;
        transition: border-color .2s, box-shadow .2s;
        position: relative;
    }
    .nr-card:hover {
        border-color: var(--border-hi);
        box-shadow: 0 0 24px rgba(14,165,233,0.06);
    }

    /* Highlight card (top-3) */
    .nr-card-hi {
        background: linear-gradient(160deg, #0d1d38, #0f1e3a);
        border: 1px solid var(--border-hi);
        border-radius: 16px;
        padding: 1.5rem;
        margin-bottom: 0.75rem;
        position: relative;
        overflow: hidden;
        transition: border-color .2s, box-shadow .2s;
    }
    .nr-card-hi::before {
        content: '';
        position: absolute; top: 0; left: 0; right: 0; height: 1px;
        background: linear-gradient(90deg, transparent, var(--cyan), transparent);
    }
    .nr-card-hi:hover {
        border-color: rgba(14,165,233,0.5);
        box-shadow: 0 0 30px rgba(14,165,233,0.1);
    }

    /* Step / how-it-works card */
    .nr-step {
        background: var(--card);
        border: 1px solid var(--border);
        border-radius: 12px;
        padding: 1.25rem 1.25rem 1.25rem 1.5rem;
        border-left: 3px solid var(--cyan);
    }
    .nr-step-num { font-size: 0.7rem; font-weight: 700; color: var(--cyan) !important; letter-spacing: 0.1em; text-transform: uppercase; margin-bottom: 0.3rem; }
    .nr-step-title { font-size: 0.9rem; font-weight: 700; color: var(--text) !important; margin-bottom: 0.4rem; }
    .nr-step-body { font-size: 0.78rem; color: var(--muted) !important; line-height: 1.55; }

    /* Compound card */
    .cpd-card {
        background: var(--card);
        border: 1px solid var(--border);
        border-radius: 12px;
        padding: 1.2rem 1.5rem;
        margin-bottom: 0.6rem;
        transition: border-color .18s, box-shadow .18s;
    }
    .cpd-card:hover { border-color: var(--border-hi); box-shadow: 0 2px 20px rgba(14,165,233,0.05); }
    .cpd-name { font-size: 1rem; font-weight: 700; color: var(--text) !important; }
    .cpd-meta { font-size: 0.78rem; color: var(--muted) !important; margin-top: 2px; }

    /* Badges */
    .b-ph  { display:inline-block; background:var(--violet-dim); color:#a78bfa; border:1px solid rgba(139,92,246,.3); border-radius:5px; padding:2px 9px; font-size:0.72rem; font-weight:700; letter-spacing:.03em; }
    .b-hi  { display:inline-block; background:rgba(16,185,129,.12); color:#34d399; border:1px solid rgba(16,185,129,.3); border-radius:5px; padding:2px 9px; font-size:0.72rem; font-weight:700; }
    .b-mid { display:inline-block; background:rgba(245,158,11,.12); color:#fbbf24; border:1px solid rgba(245,158,11,.3); border-radius:5px; padding:2px 9px; font-size:0.72rem; font-weight:700; }
    .b-lo  { display:inline-block; background:rgba(244,63,94,.12);  color:#fb7185; border:1px solid rgba(244,63,94,.3);  border-radius:5px; padding:2px 9px; font-size:0.72rem; font-weight:700; }
    .b-tag { display:inline-block; background:var(--raised); color:var(--muted); border:1px solid var(--border); border-radius:5px; padding:2px 8px; font-size:0.7rem; font-weight:500; margin-right:4px; }

    /* Score bars */
    .sbar { display:flex; align-items:center; gap:8px; margin:3px 0; }
    .sbar-lbl { width:82px; font-size:0.7rem; color:var(--muted) !important; flex-shrink:0; text-align:right; }
    .sbar-track { flex:1; background:var(--dim); border-radius:3px; height:4px; overflow:hidden; }
    .sbar-fill  { height:4px; border-radius:3px; }
    .sbar-val   { width:30px; font-size:0.7rem; color:var(--text) !important; flex-shrink:0; font-variant-numeric:tabular-nums; }

    /* Section label */
    .sec-label {
        font-size: 0.65rem; font-weight: 700; letter-spacing: 0.15em;
        text-transform: uppercase; color: var(--muted) !important;
        margin-bottom: 0.75rem; padding-bottom: 0.5rem;
        border-bottom: 1px solid var(--border);
    }

    /* Status row */
    .st-row { display:flex; align-items:center; gap:6px; font-size:0.8rem; margin-bottom:4px; }
    .dot-on  { width:6px; height:6px; border-radius:50%; background:var(--emerald); flex-shrink:0; }
    .dot-off { width:6px; height:6px; border-radius:50%; background:var(--rose); flex-shrink:0; }
    .dot-na  { width:6px; height:6px; border-radius:50%; background:var(--dim); flex-shrink:0; }

    /* Evidence card (trials / papers) */
    .ev-card {
        background: var(--raised);
        border: 1px solid var(--border);
        border-radius: 10px;
        padding: 1rem 1.25rem;
        margin-bottom: 0.5rem;
    }
    .ev-title { font-size: 0.88rem; font-weight: 600; color: var(--text) !important; line-height: 1.4; margin-bottom: 0.3rem; }
    .ev-meta  { font-size: 0.75rem; color: var(--muted) !important; }
    .ev-link  { font-size: 0.75rem; color: var(--cyan) !important; text-decoration: none; }

    /* PBPK metric pill */
    .pk-pill {
        background: var(--raised);
        border: 1px solid var(--border);
        border-radius: 10px;
        padding: 0.75rem 1rem;
        text-align: center;
    }
    .pk-val  { font-size: 1.1rem; font-weight: 700; color: var(--cyan) !important; }
    .pk-lbl  { font-size: 0.68rem; color: var(--muted) !important; text-transform: uppercase; letter-spacing: 0.06em; margin-top: 2px; }

    /* Resolution pill */
    .res-pill {
        background: var(--raised);
        border: 1px solid var(--border-hi);
        border-radius: 10px;
        padding: 0.9rem 1.25rem;
        margin-bottom: 1.25rem;
        display: flex;
        align-items: center;
        gap: 1.25rem;
    }
    .res-heading { font-size: 1rem; font-weight: 700; color: var(--cyan) !important; }
    .res-meta    { font-size: 0.76rem; color: var(--muted) !important; margin-top: 1px; }

    /* Sidebar nav */
    .nav-active {
        background: var(--cyan-dim) !important;
        color: var(--cyan) !important;
        border-left: 3px solid var(--cyan);
        padding-left: calc(1rem - 3px);
    }
    </style>
    """, unsafe_allow_html=True)


# ─── Helpers ──────────────────────────────────────────────────────────────────

def _phase_lbl(phase) -> str:
    try: p = int(float(phase or 0))
    except Exception: p = 0
    return {4:"FDA Approved", 3:"Phase III", 2:"Phase II", 1:"Phase I"}.get(p, "Preclinical")


def _phase_badge(phase) -> str:
    return f"<span class='b-ph'>{_phase_lbl(phase)}</span>"


def _score_badge(s: float) -> str:
    cls = "b-hi" if s >= 0.60 else "b-mid" if s >= 0.35 else "b-lo"
    return f"<span class='{cls}'>{s:.0%}</span>"


def _bar_color(v: float) -> str:
    if v >= 0.60: return "#10b981"
    if v >= 0.35: return "#f59e0b"
    return "#f43f5e"


def _score_bars(bd: Dict, compact: bool = False) -> str:
    labels = [("Indication","indication_score"),("Target","target_score"),
              ("Activity","activity_score"),("Network","network_score")]
    rows = ""
    for lbl, key in labels:
        v = float(bd.get(key) or 0)
        pct = round(v * 100, 1)
        col = _bar_color(v)
        rows += (f"<div class='sbar'>"
                 f"<span class='sbar-lbl'>{lbl}</span>"
                 f"<div class='sbar-track'><div class='sbar-fill' style='width:{pct}%;background:{col};'></div></div>"
                 f"<span class='sbar-val'>{v:.0%}</span>"
                 f"</div>")
    return rows


def _status(label: str, ok: bool, na: bool = False):
    dot = "dot-na" if na else ("dot-on" if ok else "dot-off")
    st.markdown(f"<div class='st-row'><span class='{dot}'></span><span style='color:#6b7d9f;'>{label}</span></div>",
                unsafe_allow_html=True)


def _sec(title: str):
    st.markdown(f"<div class='sec-label'>{title}</div>", unsafe_allow_html=True)


# ─── Session ──────────────────────────────────────────────────────────────────

def _init():
    defs = {"page": "dashboard", "disease_query": "", "compounds": [],
            "selected_compound": None}
    for k, v in defs.items():
        if k not in st.session_state: st.session_state[k] = v


# ─── Sidebar ──────────────────────────────────────────────────────────────────

def _sidebar():
    with st.sidebar:
        st.markdown(
            "<div style='padding:1rem 0 1.5rem;'>"
            "<div style='font-size:1.35rem;font-weight:800;"
            "background:linear-gradient(90deg,#0ea5e9,#8b5cf6);"
            "-webkit-background-clip:text;-webkit-text-fill-color:transparent;"
            "background-clip:text;'>NeuroRepurpose</div>"
            "<div style='font-size:0.72rem;color:#6b7d9f;margin-top:2px;letter-spacing:.04em;'>"
            "Drug Repurposing Intelligence</div></div>",
            unsafe_allow_html=True,
        )

        _sec("Navigation")
        nav_items = [("dashboard","Dashboard"),("discover","Discover"),
                     ("network","Knowledge Network"),("database","Data Explorer")]
        for pid, label in nav_items:
            active = st.session_state.page == pid
            t = "primary" if active else "secondary"
            if st.button(label, key=f"nav_{pid}", use_container_width=True, type=t):
                st.session_state.page = pid
                st.rerun()

        st.markdown("---")
        _sec("Data Sources")
        _status("ChEMBL 33",    DB_AVAILABLE)
        _status("HetioNet",     DB_AVAILABLE)
        _status("MeSH Ontology",MESH_AVAILABLE)
        _status("Scoring Engine",SCORER_AVAILABLE)
        _status("PBPK Simulator",PBPK_AVAILABLE)
        _status("Molecular Docking",DOCKING_AVAILABLE)
        _status("3D Viewer",    VIZ_3D_AVAILABLE)

        if DB_AVAILABLE:
            st.markdown("---")
            _sec("Live Counts")
            try:
                stats = get_stats()
                for k, v in stats.items():
                    st.markdown(
                        f"<div style='display:flex;justify-content:space-between;"
                        f"font-size:0.8rem;margin-bottom:2px;'>"
                        f"<span style='color:#6b7d9f;'>{k.replace('_',' ').title()}</span>"
                        f"<span style='color:#f0f4ff;font-weight:600;'>{v:,}</span></div>",
                        unsafe_allow_html=True,
                    )
            except Exception:
                pass

        st.markdown("---")
        st.markdown("<div style='font-size:0.7rem;color:#3d5080;line-height:1.6;'>"
                    "ChEMBL 33 · HetioNet · MeSH NLM<br>RDKit · ClinicalTrials.gov · PubMed</div>",
                    unsafe_allow_html=True)


# ─── Dashboard ────────────────────────────────────────────────────────────────

def page_dashboard():
    st.markdown(
        "<div class='nr-hero'>"
        "<div class='nr-hero-label'>Drug Repurposing Platform</div>"
        "<div class='nr-hero-title'>NeuroRepurpose Intelligence</div>"
        "<div class='nr-hero-sub'>"
        "Identify new therapeutic uses for approved drugs using ChEMBL 33, HetioNet knowledge graphs, "
        "MeSH disease ontology, multi-signal evidence scoring, PBPK pharmacokinetics, and real clinical evidence. "
        "Supports any disease — not just neurology."
        "</div>"
        "</div>",
        unsafe_allow_html=True,
    )

    col_q, col_btn = st.columns([6, 1])
    with col_q:
        q = st.text_input("hero_q", placeholder="Search any disease — Alzheimer, MS, epilepsy, depression, ALS...",
                          label_visibility="collapsed")
    with col_btn:
        if st.button("Search", use_container_width=True) and q:
            st.session_state.disease_query = q
            st.session_state.page = "discover"
            st.rerun()

    st.markdown("---")

    stats = {}
    if DB_AVAILABLE:
        try: stats = get_stats()
        except Exception: pass

    c1,c2,c3,c4,c5 = st.columns(5)
    with c1: st.metric("Compounds",       f"{stats.get('compounds',0):,}")
    with c2: st.metric("Protein Targets", f"{stats.get('targets',0):,}")
    with c3: st.metric("Drug Indications",f"{stats.get('indications',0):,}")
    with c4: st.metric("KG Nodes",        f"{stats.get('hetionet_nodes',0):,}")
    with c5: st.metric("MeSH Diseases",   f"{stats.get('mesh_diseases',0):,}")

    if not DB_AVAILABLE:
        st.warning("Database unavailable. Run: `python database/importer.py`")
        return

    st.markdown("---")
    st.markdown("### Disease Areas")
    st.caption("All diseases loaded from the database. Click any to run drug discovery.")

    diseases = get_available_diseases(limit=20)
    if diseases:
        cols = st.columns(5)
        for i, d in enumerate(diseases[:20]):
            with cols[i % 5]:
                if st.button(d, key=f"ql_{i}", use_container_width=True, type="secondary"):
                    st.session_state.disease_query = d
                    st.session_state.page = "discover"
                    st.rerun()
    else:
        st.info("No diseases in database.")

    st.markdown("---")
    st.markdown("### Platform Pipeline")
    c1,c2,c3,c4 = st.columns(4)
    steps = [
        ("01","MeSH Resolution",
         "Free-text query matched against the full MeSH ontology — exact headings, entry terms, "
         "and abbreviation aliases. Expanded to parent and child disease concepts."),
        ("02","Evidence Scoring",
         "Four signals: clinical indication evidence (40%), target mechanistic overlap (30%), "
         "pChEMBL activity potency (20%), and HetioNet network path count (10%)."),
        ("03","PBPK Simulation",
         "One-compartment pharmacokinetic model predicts plasma, liver, brain, and kidney "
         "concentration-time curves from MW + LogP. Disease-specific tissue adjustments."),
        ("04","Clinical Evidence",
         "Real-time lookup of ClinicalTrials.gov studies and PubMed publications for "
         "each drug-disease pair. No cached or synthetic data."),
    ]
    for col,(num,title,body) in zip([c1,c2,c3,c4],steps):
        with col:
            st.markdown(
                f"<div class='nr-step'>"
                f"<div class='nr-step-num'>Step {num}</div>"
                f"<div class='nr-step-title'>{title}</div>"
                f"<div class='nr-step-body'>{body}</div>"
                f"</div>",
                unsafe_allow_html=True,
            )


# ─── Discover ─────────────────────────────────────────────────────────────────

def page_discover():
    st.markdown("## Discover")
    st.caption("Enter any disease name. MeSH ontology resolves abbreviations, synonyms, and related concepts.")

    col_inp, col_btn = st.columns([6, 1])
    with col_inp:
        disease_input = st.text_input(
            "discover_inp",
            value=st.session_state.disease_query,
            placeholder="e.g. Parkinson disease · MS · ALS · schizophrenia · type 2 diabetes...",
            label_visibility="collapsed",
        )
    with col_btn:
        run = st.button("Analyse", use_container_width=True)

    # Autocomplete suggestions
    if disease_input and len(disease_input) >= 2 and MESH_AVAILABLE:
        sugs = suggest_diseases(disease_input, limit=6)
        if sugs and disease_input.strip().lower() not in [s.lower() for s in sugs]:
            st.markdown("<div style='font-size:0.74rem;color:#6b7d9f;margin:4px 0 2px;'>Suggestions</div>",
                        unsafe_allow_html=True)
            scols = st.columns(min(len(sugs), 6))
            for i, sg in enumerate(sugs):
                with scols[i]:
                    if st.button(sg, key=f"sug_{i}", use_container_width=True, type="secondary"):
                        st.session_state.disease_query = sg
                        st.rerun()

    if run and disease_input.strip():
        st.session_state.disease_query = disease_input.strip()

    if not st.session_state.disease_query:
        st.info("Enter a disease name above to begin.")
        return

    with st.spinner("Resolving disease and retrieving compounds..."):
        try:
            resolved, expanded_ids, compounds = resolve_and_fetch(st.session_state.disease_query)
        except Exception as e:
            st.error(f"Query failed: {e}")
            return

    if not resolved:
        st.warning(f"No MeSH match found for '{st.session_state.disease_query}'.")
        if not MESH_AVAILABLE:
            st.info("MeSH table empty. Run: `python database/mesh_importer.py`")
        return

    primary = resolved[0]

    # Resolution banner
    st.markdown(
        f"<div class='res-pill'>"
        f"<div style='flex:1;'>"
        f"<div class='res-heading'>{primary.get('heading', st.session_state.disease_query)}</div>"
        f"<div class='res-meta'>"
        f"MeSH {primary.get('mesh_id','N/A')} &nbsp;·&nbsp; "
        f"{len(resolved)} record(s) matched &nbsp;·&nbsp; "
        f"{len(expanded_ids)} expanded IDs (parents + children)"
        f"</div></div>"
        f"</div>",
        unsafe_allow_html=True,
    )

    if not compounds:
        st.warning("No compounds found for this disease in ChEMBL.")
        return

    with st.spinner(f"Scoring {len(compounds)} candidates..."):
        try:
            scored = score_compound_list(list(compounds), expanded_ids)
        except Exception:
            scored = compounds
            for c in scored:
                c.setdefault("score", float(c.get("max_phase") or 0) / 4)
                c.setdefault("score_breakdown", {})

    # Top-3
    top3 = [c for c in scored if float(c.get("max_phase") or 0) >= 3][:3] or scored[:3]
    if top3:
        _sec("Top Candidates")
        tcols = st.columns(len(top3))
        for i, c in enumerate(top3):
            _top_card(c, tcols[i], i)

    st.markdown("---")

    # Filters
    fc1, fc2, fc3 = st.columns([1, 1, 2])
    with fc1: min_ph = st.selectbox("Min phase", [0,1,2,3,4], key="d_ph")
    with fc2: show_n = st.slider("Show", 5, min(80,len(scored)), min(30,len(scored)), key="d_n")
    with fc3: sort_by = st.selectbox("Sort", ["Score","Phase","Name"], key="d_sort")

    filt = [c for c in scored if float(c.get("max_phase") or 0) >= min_ph]
    if sort_by == "Name": filt.sort(key=lambda x: (x.get("name") or "").lower())
    elif sort_by == "Phase": filt.sort(key=lambda x: float(x.get("max_phase") or 0), reverse=True)
    filt = filt[:show_n]

    if len(filt) >= 4:
        with st.expander("Phase distribution", expanded=False):
            _phase_chart(filt)

    _sec(f"{len(filt)} Repurposing Candidates")
    for i, c in enumerate(filt):
        _compound_card(c, i)

    if filt:
        df_dl = pd.DataFrame([{
            "Name": c.get("name"), "ChEMBL": c.get("chembl_id"),
            "Phase": c.get("max_phase"),
            "Score": round(float(c.get("score") or 0), 4),
            "Indication": round(float((c.get("score_breakdown") or {}).get("indication_score") or 0), 4),
            "Target":     round(float((c.get("score_breakdown") or {}).get("target_score") or 0), 4),
            "Activity":   round(float((c.get("score_breakdown") or {}).get("activity_score") or 0), 4),
            "Network":    round(float((c.get("score_breakdown") or {}).get("network_score") or 0), 4),
            "Mechanisms": c.get("mechanisms",""), "SMILES": c.get("smiles",""),
        } for c in filt])
        st.download_button("Download CSV", data=df_dl.to_csv(index=False),
                           file_name=f"neurorepurpose_{primary.get('mesh_id','x')}.csv",
                           mime="text/csv")


def _top_card(c: Dict, col, idx: int):
    name = c.get("name") or "Unknown"
    score = float(c.get("score") or 0)
    bd = c.get("score_breakdown") or {}
    mechs = (c.get("mechanisms") or "—")[:90]
    with col:
        st.markdown(
            f"<div class='nr-card-hi'>"
            f"<div style='display:flex;justify-content:space-between;align-items:flex-start;margin-bottom:.6rem;'>"
            f"<span style='font-size:.95rem;font-weight:700;color:#eef2ff;'>{name}</span>"
            f"<span>{_score_badge(score)}</span>"
            f"</div>"
            f"<div style='margin-bottom:.5rem;'>{_phase_badge(c.get('max_phase'))}</div>"
            f"<div style='font-size:.75rem;color:#6b7d9f;margin-bottom:.7rem;line-height:1.4;'>{mechs}</div>"
            f"{_score_bars(bd)}"
            f"</div>",
            unsafe_allow_html=True,
        )
        if st.button("Full Analysis", key=f"ta_{idx}_{name[:5]}", use_container_width=True):
            st.session_state.selected_compound = c
            st.session_state.page = "analysis"
            st.rerun()


def _compound_card(c: Dict, idx: int):
    name    = c.get("name") or "Unknown"
    score   = float(c.get("score") or 0)
    phase   = c.get("max_phase", 0)
    mechs   = (c.get("mechanisms") or "—")[:140]
    targets = (c.get("targets") or "—")[:120]
    bd      = c.get("score_breakdown") or {}

    col_main, col_act = st.columns([7, 1])
    with col_main:
        st.markdown(
            f"<div class='cpd-card'>"
            f"<div style='display:flex;justify-content:space-between;align-items:flex-start;margin-bottom:.45rem;'>"
            f"  <span class='cpd-name'>{name} &nbsp; {_phase_badge(phase)}</span>"
            f"  <span>{_score_badge(score)}</span>"
            f"</div>"
            f"<div class='cpd-meta'><b style='color:#c4cde8;'>Mechanism:</b> {mechs}</div>"
            f"<div class='cpd-meta' style='margin-bottom:.55rem;'><b style='color:#c4cde8;'>Targets:</b> {targets}</div>"
            f"{_score_bars(bd)}"
            f"</div>",
            unsafe_allow_html=True,
        )
    with col_act:
        st.markdown("<div style='height:1.2rem;'></div>", unsafe_allow_html=True)
        if st.button("Analyse", key=f"ca_{idx}_{name[:6]}", use_container_width=True):
            st.session_state.selected_compound = c
            st.session_state.page = "analysis"
            st.rerun()


def _phase_chart(compounds: List[Dict]):
    pm = {4:"FDA Approved",3:"Phase III",2:"Phase II",1:"Phase I",0:"Preclinical"}
    cnt: Dict = {}
    for c in compounds:
        lbl = pm.get(int(float(c.get("max_phase") or 0)), "Preclinical")
        cnt[lbl] = cnt.get(lbl, 0) + 1
    df = pd.DataFrame(list(cnt.items()), columns=["Phase","Count"])
    fig = px.pie(df, values="Count", names="Phase", hole=.52,
                 color_discrete_sequence=["#0ea5e9","#8b5cf6","#10b981","#f59e0b","#f43f5e"])
    fig.update_layout(paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
                      font_color="#f0f4ff", margin=dict(t=10,b=10,l=10,r=10),
                      legend=dict(font=dict(color="#6b7d9f")))
    st.plotly_chart(fig, use_container_width=True)


# ─── Compound Analysis ────────────────────────────────────────────────────────

def page_analysis():
    c = st.session_state.selected_compound
    if not c:
        st.info("Select a compound from the Discover page.")
        if st.button("Go to Discover"):
            st.session_state.page = "discover"
            st.rerun()
        return

    name      = c.get("name") or "Unknown"
    smiles    = c.get("smiles") or ""
    chembl_id = c.get("chembl_id") or ""
    cid       = c.get("id")
    score     = float(c.get("score") or 0)
    bd        = c.get("score_breakdown") or {}
    phase     = c.get("max_phase", 0)

    # Header
    col_h, col_meta = st.columns([4, 1])
    with col_h:
        st.markdown(f"## {name}")
        if chembl_id:
            st.caption(f"ChEMBL ID: {chembl_id}")
        mech = c.get("mechanisms") or "—"
        st.markdown(f"<div style='font-size:.88rem;color:#6b7d9f;margin-top:.25rem;'>Mechanism: "
                    f"<span style='color:#c4cde8;'>{mech[:260]}</span></div>", unsafe_allow_html=True)
    with col_meta:
        st.markdown(
            f"<div style='text-align:right;padding-top:.5rem;'>"
            f"{_phase_badge(phase)}"
            f"<br><br>{_score_badge(score)} overall"
            f"</div>",
            unsafe_allow_html=True,
        )
        if bd:
            st.markdown(f"<div style='margin-top:.75rem;'>{_score_bars(bd)}</div>", unsafe_allow_html=True)

    st.markdown("---")

    # Build tab list
    tab_labels = ["Properties","Targets","Activities","Indications",
                  "Clinical Evidence","PBPK Simulation","Docking","Network"]
    if QUANTUM_AVAILABLE: tab_labels.append("Quantum")
    if VIZ_3D_AVAILABLE:  tab_labels.append("3D Structure")
    tabs = st.tabs(tab_labels)
    ti = 0

    # ── Properties ─────────────────────────────────────────────────────────
    with tabs[ti]:
        _tab_properties(c, smiles)
    ti += 1

    # ── Targets ────────────────────────────────────────────────────────────
    with tabs[ti]:
        tgts = get_compound_targets(int(cid)) if (DB_AVAILABLE and cid) else []
        if tgts:
            df_t = pd.DataFrame(tgts)
            keep = [col for col in ["name","gene_symbol","target_type","mechanism","action_type","confidence","organism"] if col in df_t.columns]
            st.dataframe(df_t[keep], use_container_width=True)
            if "confidence" in df_t.columns and len(df_t) >= 2:
                top_t = df_t.dropna(subset=["confidence"]).sort_values("confidence").tail(12)
                nc = "gene_symbol" if "gene_symbol" in top_t.columns else "name"
                fig = px.bar(top_t, x="confidence", y=nc, orientation="h",
                             color="confidence",
                             color_continuous_scale=[[0,"#f43f5e"],[.5,"#f59e0b"],[1,"#10b981"]],
                             title="Target confidence")
                fig.update_layout(paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
                                  font_color="#f0f4ff", margin=dict(t=30,b=10))
                st.plotly_chart(fig, use_container_width=True)
        else:
            st.info("No target data in database for this compound.")
    ti += 1

    # ── Activities ─────────────────────────────────────────────────────────
    with tabs[ti]:
        acts = get_compound_activities(int(cid)) if (DB_AVAILABLE and cid) else []
        if acts:
            df_a = pd.DataFrame(acts)
            st.dataframe(df_a, use_container_width=True)
            if "pchembl_value" in df_a.columns:
                dp = df_a.dropna(subset=["pchembl_value"]).sort_values("pchembl_value", ascending=False).head(15)
                if not dp.empty:
                    nc = "gene_symbol" if "gene_symbol" in dp.columns else "target_name"
                    fig = px.bar(dp.sort_values("pchembl_value"), x="pchembl_value", y=nc,
                                 orientation="h",
                                 color="pchembl_value",
                                 color_continuous_scale=[[0,"#8b5cf6"],[.5,"#f59e0b"],[1,"#10b981"]],
                                 title="pChEMBL values (higher = more potent)")
                    fig.add_vline(x=6, line_dash="dash", line_color="#f59e0b",
                                  annotation=dict(text="1 µM", font_color="#f59e0b"))
                    fig.add_vline(x=8, line_dash="dash", line_color="#10b981",
                                  annotation=dict(text="10 nM", font_color="#10b981"))
                    fig.update_layout(paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
                                      font_color="#f0f4ff")
                    st.plotly_chart(fig, use_container_width=True)
        else:
            st.info("No quantitative activity data available.")
    ti += 1

    # ── Indications ────────────────────────────────────────────────────────
    with tabs[ti]:
        inds = get_compound_indications(int(cid)) if (DB_AVAILABLE and cid) else []
        if inds:
            df_i = pd.DataFrame(inds)
            keep = [col for col in ["disease","mesh_id","max_phase","source"] if col in df_i.columns]
            st.dataframe(df_i[keep], use_container_width=True)
        elif c.get("indications"):
            for i_str in str(c["indications"]).split(";"):
                st.markdown(f"- {i_str.strip()}")
        else:
            st.info("No indication data available.")
    ti += 1

    # ── Clinical Evidence ──────────────────────────────────────────────────
    with tabs[ti]:
        _tab_clinical_evidence(name)
    ti += 1

    # ── PBPK Simulation ────────────────────────────────────────────────────
    with tabs[ti]:
        _tab_pbpk(c, smiles)
    ti += 1

    # ── Docking ────────────────────────────────────────────────────────────
    with tabs[ti]:
        _tab_docking(name, smiles, cid)
    ti += 1

    # ── Network ────────────────────────────────────────────────────────────
    with tabs[ti]:
        _tab_network(name)
    ti += 1

    # ── Quantum ────────────────────────────────────────────────────────────
    if QUANTUM_AVAILABLE:
        with tabs[ti]:
            render_quantum_optimization_section(drug_name=name, smiles=smiles, drug_data=c)
        ti += 1

    # ── 3D Structure ───────────────────────────────────────────────────────
    if VIZ_3D_AVAILABLE:
        with tabs[ti]:
            if not smiles:
                st.warning(f"No SMILES available for {name}.")
            elif _viz:
                tgts2 = get_compound_targets(int(cid)) if (DB_AVAILABLE and cid) else []
                tgt_name = tgts2[0]["name"] if tgts2 else ""
                pdb = None
                if tgt_name:
                    try:
                        from real_pdb_fetcher import RealPDBFetcher
                        pdb = RealPDBFetcher().fetch_pdb(tgt_name)
                    except Exception: pass
                _viz.render_visualization(name, tgt_name or "", smiles, pdb)
        ti += 1


def _tab_properties(c: Dict, smiles: str):
    mw = c.get("mw"); logp = c.get("alogp"); psa = c.get("psa")
    hba = c.get("hba"); hbd = c.get("hbd"); qed = None

    if smiles and not all([mw, logp, psa]):
        try:
            from rdkit import Chem
            from rdkit.Chem import Descriptors, rdMolDescriptors, QED as rdQED
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                mw   = mw   or round(Descriptors.ExactMolWt(mol), 1)
                logp = logp or round(Descriptors.MolLogP(mol), 2)
                psa  = psa  or round(Descriptors.TPSA(mol), 1)
                hba  = hba  or rdMolDescriptors.CalcNumHBA(mol)
                hbd  = hbd  or rdMolDescriptors.CalcNumHBD(mol)
                qed  = round(rdQED.qed(mol), 3)
        except Exception: pass

    col1,col2,col3,col4,col5 = st.columns(5)
    with col1: st.metric("Mol Weight (Da)", f"{mw:.1f}" if mw else "N/A")
    with col2: st.metric("LogP",            f"{logp:.2f}" if logp is not None else "N/A")
    with col3: st.metric("TPSA (sq A)",     f"{psa:.1f}" if psa else "N/A")
    with col4: st.metric("QED",             f"{qed:.3f}" if qed is not None else "N/A")
    with col5: st.metric("HBA / HBD",       f"{hba} / {hbd}" if hba is not None else "N/A")

    viols = sum([bool(mw and mw>500), bool(logp is not None and logp>5),
                 bool(hbd is not None and hbd>5), bool(hba is not None and hba>10)])
    col = "#10b981" if viols == 0 else "#f59e0b" if viols == 1 else "#f43f5e"
    txt = "Passes Lipinski Rule of 5" if viols == 0 else f"Lipinski: {viols} violation(s)"
    st.markdown(f"<div style='font-size:.85rem;color:{col};margin:.75rem 0;'>{txt}</div>",
                unsafe_allow_html=True)

    ro5 = [("Molecular weight", f"{mw:.1f} Da" if mw else "N/A", bool(mw and mw>500)),
           ("LogP",             f"{logp:.2f}" if logp is not None else "N/A", bool(logp is not None and logp>5)),
           ("H-bond donors",    str(hbd) if hbd is not None else "N/A", bool(hbd is not None and hbd>5)),
           ("H-bond acceptors", str(hba) if hba is not None else "N/A", bool(hba is not None and hba>10))]
    for prop, val, fail in ro5:
        col2 = "#f43f5e" if fail else "#10b981"
        mk   = "FAIL" if fail else "PASS"
        st.markdown(
            f"<div style='display:flex;justify-content:space-between;background:#0f1e38;"
            f"border-radius:7px;padding:5px 14px;margin-bottom:3px;font-size:.82rem;'>"
            f"<span style='color:#6b7d9f;'>{prop}</span>"
            f"<span style='color:#eef2ff;'>{val}</span>"
            f"<span style='color:{col2};font-weight:700;font-size:.7rem;'>{mk}</span>"
            f"</div>",
            unsafe_allow_html=True,
        )

    if smiles:
        st.markdown("<div style='margin-top:1rem;'></div>", unsafe_allow_html=True)
        st.code(smiles, language=None)


def _tab_clinical_evidence(drug_name: str):
    disease_query = st.session_state.get("disease_query", "")
    if not disease_query:
        st.info("Run a disease search first to load clinical evidence.")
        return

    col_t, col_p = st.columns(2)

    with col_t:
        _sec("Clinical Trials  (ClinicalTrials.gov)")
        with st.spinner("Fetching trials..."):
            trials = fetch_clinical_trials(drug_name, disease_query, n=8)

        if trials:
            for t in trials:
                status_col = {"RECRUITING":"#10b981","COMPLETED":"#0ea5e9",
                              "ACTIVE, NOT RECRUITING":"#f59e0b"}.get(t.get("status",""), "#6b7d9f")
                st.markdown(
                    f"<div class='ev-card'>"
                    f"<div class='ev-title'>{t.get('title','—')}</div>"
                    f"<div class='ev-meta'>"
                    f"<span style='color:{status_col};font-weight:600;'>{t.get('status','—')}</span>"
                    f" &nbsp;·&nbsp; {t.get('phase','N/A')}"
                    f" &nbsp;·&nbsp; <a class='ev-link' href='{t.get('url','#')}' target='_blank'>"
                    f"{t.get('nct_id','')}</a>"
                    f"</div></div>",
                    unsafe_allow_html=True,
                )
        else:
            st.info(f"No trials found for {drug_name} in {disease_query}.")

    with col_p:
        _sec("Publications  (PubMed)")
        with st.spinner("Fetching papers..."):
            papers = fetch_pubmed(drug_name, disease_query, n=8)

        if papers:
            for p in papers:
                st.markdown(
                    f"<div class='ev-card'>"
                    f"<div class='ev-title'>{p.get('title','—')}</div>"
                    f"<div class='ev-meta'>"
                    f"{p.get('journal','—')} &nbsp;·&nbsp; {p.get('year','')}"
                    f" &nbsp;·&nbsp; <a class='ev-link' href='{p.get('url','#')}' target='_blank'>"
                    f"PMID {p.get('pmid','')}</a>"
                    f"</div></div>",
                    unsafe_allow_html=True,
                )
        else:
            st.info(f"No PubMed papers found for {drug_name} + {disease_query}.")


def _tab_pbpk(c: Dict, smiles: str):
    _sec("PBPK Pharmacokinetic Simulation")

    if not PBPK_AVAILABLE:
        st.info("PBPK module unavailable. Check services/pbpk_simulation.py.")
        return

    mw   = float(c.get("mw") or 0)
    logp = float(c.get("alogp") or 0)
    name = c.get("name") or "Unknown"

    # Compute from RDKit if needed
    if (not mw or not logp) and smiles:
        try:
            from rdkit import Chem
            from rdkit.Chem import Descriptors
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                mw   = mw   or round(Descriptors.ExactMolWt(mol), 1)
                logp = logp or round(Descriptors.MolLogP(mol), 2)
        except Exception: pass

    if not mw:
        st.info("Molecular weight not available. Cannot run PBPK simulation.")
        return

    # Controls
    disease_query = st.session_state.get("disease_query", "Unknown")
    col1, col2, col3, col4 = st.columns(4)
    with col1: dose    = st.number_input("Dose (mg)", 1.0, 2000.0, 100.0, step=10.0)
    with col2: route   = st.selectbox("Route", ["oral","IV"])
    with col3: dur     = st.number_input("Duration (h)", 4.0, 168.0, 24.0, step=4.0)
    with col4: binding = st.number_input("Binding affinity (kcal/mol)", -20.0, 0.0, -8.0, step=0.5)

    if st.button("Run Simulation", key="run_pbpk"):
        with st.spinner("Running compartmental PK model..."):
            try:
                sim = PBPKSimulator(disease_name=disease_query)
                result = sim.simulate_drug_exposure(
                    drug_name=name, molecular_weight=mw, logp=logp,
                    dose_mg=dose, route=route, duration_hours=dur,
                    binding_affinity=binding,
                )
                st.session_state["pbpk_result"] = result
            except Exception as e:
                st.error(f"Simulation failed: {e}")

    result = st.session_state.get("pbpk_result")
    if result and result.get("success"):
        pk = result.get("pk_metrics", {})
        adme = result.get("adme_parameters", {})

        # PK metrics row
        pkc = st.columns(4)
        metrics = [
            ("Cmax", f"{pk.get('cmax_ng_ml',0):.1f}", "ng/mL"),
            ("Tmax", f"{pk.get('tmax_hours',0):.2f}", "hours"),
            ("AUC",  f"{pk.get('auc_ng_h_ml',0):.0f}", "ng·h/mL"),
            ("t½",   f"{pk.get('t_half_hours',0):.2f}", "hours"),
        ]
        for col, (lbl, val, unit) in zip(pkc, metrics):
            with col:
                st.markdown(
                    f"<div class='pk-pill'><div class='pk-val'>{val}</div>"
                    f"<div class='pk-lbl'>{lbl} ({unit})</div></div>",
                    unsafe_allow_html=True,
                )

        st.markdown("<div style='height:.5rem;'></div>", unsafe_allow_html=True)

        # Concentration-time chart
        times = result.get("time_hours", [])
        fig = go.Figure()
        compartments = [
            ("plasma_concentration_ng_ml", "Plasma", "#0ea5e9"),
            ("liver_concentration_ng_ml",  "Liver",  "#f59e0b"),
            ("brain_concentration_ng_ml",  "Brain",  "#8b5cf6"),
        ]
        for key, label, color in compartments:
            vals = result.get(key, [])
            if vals:
                fig.add_trace(go.Scatter(x=times, y=vals, name=label,
                                         line=dict(color=color, width=2)))
        fig.update_layout(
            title=f"Concentration-Time Profile — {name} ({route}, {dose:.0f} mg)",
            xaxis_title="Time (hours)", yaxis_title="Concentration (ng/mL)",
            paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="#0b1527",
            font_color="#f0f4ff",
            legend=dict(bgcolor="rgba(0,0,0,0)", font_color="#6b7d9f"),
            xaxis=dict(gridcolor="#1c304f"), yaxis=dict(gridcolor="#1c304f"),
        )
        st.plotly_chart(fig, use_container_width=True)

        # ADME + Safety
        col_a, col_s = st.columns(2)
        with col_a:
            _sec("ADME Parameters")
            adme_items = [("Bioavailability (F)",f"{adme.get('f',0):.0%}"),
                          ("Vol. Distribution (Vd)",f"{adme.get('vd',0):.2f} L/kg"),
                          ("Clearance (CL)",f"{adme.get('cl',0):.3f} L/h/kg"),
                          ("Fraction Unbound",f"{adme.get('fu',0):.0%}"),
                          ("Brain Kp",f"{adme.get('kp_brain',0):.2f}"),
                          ("Liver Kp",f"{adme.get('kp_liver',0):.2f}")]
            for lbl, val in adme_items:
                st.markdown(
                    f"<div style='display:flex;justify-content:space-between;"
                    f"background:#0f1e38;border-radius:7px;padding:5px 14px;"
                    f"margin-bottom:3px;font-size:.82rem;'>"
                    f"<span style='color:#6b7d9f;'>{lbl}</span>"
                    f"<span style='color:#eef2ff;font-weight:600;'>{val}</span></div>",
                    unsafe_allow_html=True,
                )
        with col_s:
            safety = result.get("safety_assessment", {})
            _sec("Safety Assessment")
            margin = safety.get("safety_margin","Unknown")
            margin_col = "#10b981" if margin == "Good" else "#f59e0b"
            st.markdown(
                f"<div style='font-size:.9rem;font-weight:700;color:{margin_col};margin-bottom:.5rem;'>"
                f"Safety Margin: {margin}</div>",
                unsafe_allow_html=True,
            )
            for w in safety.get("warnings", []):
                st.markdown(f"<div style='font-size:.82rem;color:#6b7d9f;margin-bottom:3px;'>- {w}</div>",
                            unsafe_allow_html=True)
            tw = safety.get("therapeutic_window","")
            if tw:
                st.markdown(f"<div style='font-size:.82rem;color:#6b7d9f;margin-top:.5rem;'>"
                            f"Therapeutic window: {tw}</div>", unsafe_allow_html=True)


def _tab_docking(name: str, smiles: str, cid):
    _sec("Molecular Docking")

    if not DOCKING_AVAILABLE:
        st.info(
            "Molecular docking requires an NVIDIA BioNeMo API key or AutoDock Vina installation. "
            "Set the `NVIDIA_API_KEY` environment variable to enable DiffDock."
        )
        return

    if not smiles:
        st.warning(f"No SMILES structure for {name}. Cannot run docking.")
        return

    tgts = get_compound_targets(int(cid)) if (DB_AVAILABLE and cid) else []
    target_names = [t["name"] for t in tgts if t.get("name")] or []

    col1, col2 = st.columns(2)
    with col1:
        if target_names:
            selected_target = st.selectbox("Select target", target_names)
        else:
            selected_target = st.text_input("Target protein name", placeholder="e.g. AChE, BACE1, GSK3B")
    with col2:
        n_poses = st.slider("Number of poses", 1, 20, 10)

    if st.button("Run Docking", key="run_dock") and selected_target:
        with st.spinner(f"Running docking: {name} -> {selected_target}..."):
            try:
                result = _docking_svc.perform_docking(
                    drug_name=name, target_name=selected_target,
                    ligand_smiles=smiles,
                )
                st.session_state["dock_result"] = result
            except Exception as e:
                st.error(f"Docking failed: {e}")

    dr = st.session_state.get("dock_result")
    if dr:
        if not dr.get("success"):
            st.warning(f"Docking returned no result: {dr.get('error','Unknown error')}")
            return

        affs = dr.get("binding_affinities", [])
        poses = dr.get("poses", [])

        if affs:
            c1, c2, c3 = st.columns(3)
            with c1: st.metric("Best Affinity",   f"{min(affs):.2f} kcal/mol")
            with c2: st.metric("Mean Affinity",   f"{sum(affs)/len(affs):.2f} kcal/mol")
            with c3: st.metric("Poses Generated", len(poses))

            fig = go.Figure(go.Bar(
                x=list(range(1, len(affs)+1)), y=affs,
                marker_color=["#0ea5e9" if a == min(affs) else "#8b5cf6" for a in affs],
            ))
            fig.update_layout(
                title="Docking Pose Binding Affinities",
                xaxis_title="Pose", yaxis_title="Affinity (kcal/mol)",
                paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="#0b1527",
                font_color="#f0f4ff",
                xaxis=dict(gridcolor="#1c304f"), yaxis=dict(gridcolor="#1c304f"),
            )
            st.plotly_chart(fig, use_container_width=True)

            conf = dr.get("confidence_scores", [])
            if conf:
                df_dock = pd.DataFrame({
                    "Pose": list(range(1, len(conf)+1)),
                    "Binding Affinity (kcal/mol)": affs[:len(conf)],
                    "Confidence": conf,
                })
                st.dataframe(df_dock, use_container_width=True)


def _tab_network(compound_name: str):
    query = st.session_state.get("disease_query", "")
    mesh_ids = []
    if query:
        try:
            _, mesh_ids, _ = resolve_and_fetch(query)
        except Exception:
            pass

    paths = []
    if DB_AVAILABLE and mesh_ids:
        with st.spinner("Loading HetioNet paths..."):
            try: paths = get_hetionet_paths(compound_name, mesh_ids)
            except Exception as e: st.warning(f"Network query failed: {e}")

    if not paths:
        st.info("No HetioNet paths found for this compound. "
                "The compound may not be in the HetioNet node set.")
        return

    G = nx.DiGraph()
    for p in paths[:60]:
        s, t = p.get("source_name",""), p.get("target_name","")
        if s and t:
            G.add_node(s, kind=p.get("source_kind",""))
            G.add_node(t, kind=p.get("target_kind",""))
            G.add_edge(s, t, me=p.get("metaedge",""))

    if not G.nodes: st.info("No network paths."); return

    pos = nx.spring_layout(G, seed=42, k=1.2)
    kc  = {"Compound":"#0ea5e9","Gene":"#10b981","Disease":"#f43f5e",
           "Anatomy":"#f59e0b","Pathway":"#8b5cf6","Biological Process":"#a78bfa"}

    nx_, ny_, nt_, nc_ = [], [], [], []
    for node, attr in G.nodes(data=True):
        x, y = pos[node]
        nx_.append(x); ny_.append(y)
        nt_.append(node[:30])
        nc_.append(kc.get(attr.get("kind",""), "#6b7d9f"))

    ex_, ey_ = [], []
    for u, v in G.edges():
        x0, y0 = pos[u]; x1, y1 = pos[v]
        ex_ += [x0,x1,None]; ey_ += [y0,y1,None]

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=ex_, y=ey_, mode="lines",
                              line=dict(color="#1c304f", width=1), hoverinfo="none"))
    fig.add_trace(go.Scatter(x=nx_, y=ny_, mode="markers+text",
                              marker=dict(size=13, color=nc_,
                                          line=dict(width=1, color="#060b18")),
                              text=nt_, textposition="top center",
                              textfont=dict(size=8, color="#6b7d9f"),
                              hoverinfo="text"))
    fig.update_layout(showlegend=False, hovermode="closest", height=500,
                      paper_bgcolor="#060b18", plot_bgcolor="#060b18",
                      font_color="#f0f4ff",
                      xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                      yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                      margin=dict(l=10,r=10,t=10,b=10))
    st.plotly_chart(fig, use_container_width=True)

    leg_cols = st.columns(len(kc))
    for col, (kind, color) in zip(leg_cols, kc.items()):
        col.markdown(f"<span style='color:{color};font-size:.75rem;'>&#9679; {kind}</span>",
                     unsafe_allow_html=True)
    st.caption(f"{G.number_of_nodes()} nodes · {G.number_of_edges()} edges")


# ─── Knowledge Network page ───────────────────────────────────────────────────

def page_network():
    st.markdown("## Knowledge Network")
    st.caption("Drug-gene-disease connections from HetioNet knowledge graph.")

    col1, col2 = st.columns([5, 1])
    with col1:
        center = st.text_input("net_q", value=st.session_state.disease_query or "",
                               placeholder="Enter disease or drug name...",
                               label_visibility="collapsed")
    with col2:
        build = st.button("Build", use_container_width=True)

    if not center:
        st.info("Enter a disease or drug name to build the network.")
        return

    if not DB_AVAILABLE:
        st.warning("Database unavailable.")
        return

    with st.spinner("Building network..."):
        resolved = resolve_disease(center) if center else []
        mesh_ids = [r["mesh_id"] for r in resolved if r.get("mesh_id")]
        expanded = expand_mesh_ids(mesh_ids) if mesh_ids else []
        compounds = (get_compounds_for_disease(expanded or mesh_ids, limit=20)
                     if (expanded or mesh_ids) else db_search(center, limit=20))
        compounds = validate_and_deduplicate(compounds, require_smiles=False)[:20]

    if not compounds:
        st.warning("No compounds found."); return

    G = nx.Graph()
    cl = resolved[0]["heading"] if resolved else center
    G.add_node(cl, kind="Disease" if resolved else "Query")

    for c in compounds:
        n = c.get("name") or "";
        if not n: continue
        G.add_node(n, kind="Compound")
        G.add_edge(cl, n)
        for tgt in (c.get("targets") or "").split(";")[:2]:
            tgt = tgt.strip()
            if tgt and len(tgt) > 2:
                G.add_node(tgt, kind="Gene")
                G.add_edge(n, tgt)

    pos = nx.spring_layout(G, seed=42, k=1.4)
    kc  = {"Disease":"#f43f5e","Query":"#0ea5e9","Compound":"#8b5cf6","Gene":"#10b981"}

    nx_, ny_, nt_, nc_, ns_ = [], [], [], [], []
    for node, attr in G.nodes(data=True):
        x, y = pos[node]
        nx_.append(x); ny_.append(y); nt_.append(node[:28])
        nc_.append(kc.get(attr.get("kind",""), "#6b7d9f"))
        ns_.append(22 if attr.get("kind") in ("Disease","Query") else 14 if attr.get("kind")=="Compound" else 9)

    ex_, ey_ = [], []
    for u, v in G.edges():
        x0,y0=pos[u]; x1,y1=pos[v]
        ex_+=[x0,x1,None]; ey_+=[y0,y1,None]

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=ex_, y=ey_, mode="lines",
                              line=dict(color="#1c304f", width=1), hoverinfo="none"))
    fig.add_trace(go.Scatter(x=nx_, y=ny_, mode="markers+text",
                              marker=dict(size=ns_, color=nc_,
                                          line=dict(width=1, color="#060b18")),
                              text=nt_, textposition="top center",
                              textfont=dict(size=8, color="#6b7d9f"),
                              hoverinfo="text"))
    fig.update_layout(showlegend=False, hovermode="closest", height=580,
                      paper_bgcolor="#060b18", plot_bgcolor="#060b18",
                      font_color="#f0f4ff",
                      xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                      yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                      margin=dict(l=10,r=10,t=10,b=10))
    st.plotly_chart(fig, use_container_width=True)

    leg_cols = st.columns(4)
    for col, (kind, color) in zip(leg_cols, kc.items()):
        col.markdown(f"<span style='color:{color};font-size:.78rem;'>&#9679; {kind}</span>",
                     unsafe_allow_html=True)
    st.caption(f"{G.number_of_nodes()} nodes · {G.number_of_edges()} edges")


# ─── Data Explorer ────────────────────────────────────────────────────────────

def page_database():
    st.markdown("## Data Explorer")

    if not DB_AVAILABLE:
        st.warning("Database unavailable. Run: `python database/importer.py`")
        return

    tabs = st.tabs(["Compounds","Targets","Diseases","Statistics"])

    with tabs[0]:
        col1, col2 = st.columns([4, 1])
        with col1:
            q = st.text_input("cmp_q", placeholder="Name, mechanism, target...",
                              label_visibility="collapsed")
        with col2:
            lim = st.number_input("Rows", 10, 500, 50, key="cmp_lim")
        compounds = db_search(q or "", limit=int(lim)) if q else []
        if not compounds and not q:
            try:
                from database.schema import get_connection
                with get_connection() as conn:
                    cur = conn.cursor()
                    cur.execute("SELECT c.chembl_id,c.name,c.smiles,c.max_phase,cp.mw,cp.alogp "
                                "FROM compounds c LEFT JOIN compound_properties cp ON cp.compound_id=c.id "
                                "ORDER BY c.max_phase DESC NULLS LAST LIMIT %s", (lim,))
                    compounds = [{"chembl_id":r[0],"name":r[1],"smiles":r[2],
                                  "max_phase":r[3],"mw":r[4],"alogp":r[5]} for r in cur.fetchall()]
            except Exception as e:
                st.error(f"Query error: {e}")
        if compounds:
            df = pd.DataFrame(compounds)
            keep = [c for c in ["chembl_id","name","max_phase","mw","alogp","psa","mechanisms"] if c in df.columns]
            st.dataframe(df[keep], use_container_width=True, height=430)
            st.download_button("Download CSV", data=df.to_csv(index=False),
                               file_name="compounds.csv", mime="text/csv")
        else:
            st.info("No results.")

    with tabs[1]:
        try:
            from database.schema import get_connection
            with get_connection() as conn:
                cur = conn.cursor()
                cur.execute("SELECT chembl_tid,name,target_type,gene_symbol,organism "
                            "FROM targets ORDER BY name LIMIT 500")
                rows = cur.fetchall()
            if rows:
                df_t = pd.DataFrame(rows, columns=["ChEMBL TID","Name","Type","Gene","Organism"])
                st.dataframe(df_t, use_container_width=True, height=430)
            else:
                st.info("No targets in database.")
        except Exception as e:
            st.error(f"Query error: {e}")

    with tabs[2]:
        try:
            from database.schema import get_connection
            with get_connection() as conn:
                cur = conn.cursor()
                cur.execute("SELECT mesh_id,heading,array_length(tree_numbers,1) AS trees,"
                            "array_length(entry_terms,1) AS synonyms "
                            "FROM mesh_diseases ORDER BY heading LIMIT 500")
                rows = cur.fetchall()
            if rows:
                df_d = pd.DataFrame(rows, columns=["MeSH ID","Heading","Tree Entries","Synonyms"])
                st.dataframe(df_d, use_container_width=True, height=430)
            else:
                st.info("MeSH table empty. Run: `python database/mesh_importer.py`")
        except Exception as e:
            st.error(f"Query error: {e}")

    with tabs[3]:
        try:
            stats = get_stats()
            c1,c2,c3 = st.columns(3)
            for i,(k,v) in enumerate(stats.items()):
                with [c1,c2,c3][i%3]:
                    st.metric(k.replace("_"," ").title(), f"{v:,}")
            from database.schema import get_connection
            with get_connection() as conn:
                cur = conn.cursor()
                cur.execute("SELECT kind,COUNT(*) AS n FROM hetionet_nodes GROUP BY kind ORDER BY n DESC")
                kr = cur.fetchall()
            if kr:
                st.markdown("---")
                st.markdown("### HetioNet Node Types")
                df_k = pd.DataFrame(kr, columns=["Kind","Count"])
                fig = px.bar(df_k, x="Kind", y="Count",
                             color="Count",
                             color_continuous_scale=[[0,"#8b5cf6"],[1,"#0ea5e9"]])
                fig.update_layout(paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="#0b1527",
                                  font_color="#f0f4ff",
                                  xaxis=dict(gridcolor="#1c304f"), yaxis=dict(gridcolor="#1c304f"))
                st.plotly_chart(fig, use_container_width=True)
        except Exception as e:
            st.error(f"Stats error: {e}")


# ─── Router ───────────────────────────────────────────────────────────────────

def main():
    _init()
    apply_theme()
    _sidebar()
    page = st.session_state.page
    if   page == "dashboard": page_dashboard()
    elif page == "discover":  page_discover()
    elif page == "analysis":  page_analysis()
    elif page == "network":   page_network()
    elif page == "database":  page_database()
    else:                     page_dashboard()


if __name__ == "__main__":
    main()
