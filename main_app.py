#!/usr/bin/env python3
"""
NeuroRepurpose Intelligence Platform
ChEMBL 33 · Knowledge Graph · MeSH · PBPK · Docking · Clinical Evidence
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

st.set_page_config(page_title="NeuroRepurpose", layout="wide",
                   initial_sidebar_state="expanded")

logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)
REPO_ROOT = Path(__file__).parent

# ─── Service imports ─────────────────────────────────────────────────────────

try:
    from services.neuro_db_service import (
        get_compounds_for_disease, search_compounds as db_search,
        get_compound_targets, get_compound_activities,
        get_compound_indications, get_hetionet_paths,
        get_stats, is_available as db_is_available,
        get_available_diseases,
    )
    DB_AVAILABLE = db_is_available()
except Exception:
    DB_AVAILABLE = False
    def get_compounds_for_disease(*a, **k): return []
    def db_search(*a, **k): return []
    def get_compound_targets(*a, **k): return []
    def get_compound_activities(*a, **k): return []
    def get_compound_indications(*a, **k): return []
    def get_hetionet_paths(*a, **k): return []
    def get_stats(): return {}
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
        for x in c:
            x.setdefault("score", float(x.get("max_phase") or 0) / 4)
            x.setdefault("score_breakdown", {})
        c.sort(key=lambda x: x["score"], reverse=True)
        return c

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
    from quantum_optimization_strategies import run_quantum_optimization
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


# ─── External data ────────────────────────────────────────────────────────────

@st.cache_data(ttl=3600, show_spinner=False)
def fetch_trials(drug: str, disease: str, n: int = 8) -> List[Dict]:
    try:
        r = requests.get(
            "https://clinicaltrials.gov/api/v2/studies",
            params={"query.cond": disease, "query.term": drug,
                    "pageSize": n, "format": "json"},
            timeout=12,
        )
        if r.status_code == 200:
            out = []
            for s in r.json().get("studies", []):
                pm = s.get("protocolSection", {})
                phases = pm.get("designModule", {}).get("phases", [])
                out.append({
                    "nct_id": pm.get("identificationModule", {}).get("nctId", ""),
                    "title":  pm.get("identificationModule", {}).get("briefTitle", ""),
                    "status": pm.get("statusModule", {}).get("overallStatus", ""),
                    "phase":  ", ".join(phases) or "N/A",
                    "url":    f"https://clinicaltrials.gov/study/{pm.get('identificationModule',{}).get('nctId','')}",
                })
            return out
    except Exception:
        pass
    return []


@st.cache_data(ttl=3600, show_spinner=False)
def fetch_papers(drug: str, disease: str, n: int = 8) -> List[Dict]:
    try:
        base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
        sr = requests.get(f"{base}/esearch.fcgi",
            params={"db": "pubmed",
                    "term": f"{drug}[tiab] AND {disease}[tiab]",
                    "retmax": n, "retmode": "json"}, timeout=10)
        ids = sr.json().get("esearchresult", {}).get("idlist", [])
        if not ids:
            return []
        smr = requests.get(f"{base}/esummary.fcgi",
            params={"db": "pubmed", "id": ",".join(ids), "retmode": "json"}, timeout=10)
        res = smr.json().get("result", {})
        return [
            {"pmid": p, "title": res[p].get("title", ""), "journal": res[p].get("source", ""),
             "year": res[p].get("pubdate", "")[:4],
             "url": f"https://pubmed.ncbi.nlm.nih.gov/{p}/"}
            for p in ids if p in res
        ]
    except Exception:
        pass
    return []


# ─── Cached pipeline ──────────────────────────────────────────────────────────

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


# ─── Theme CSS ────────────────────────────────────────────────────────────────

_DARK_VARS = """
    --bg:          #04080f;
    --card:        #080f1e;
    --raised:      #0d1729;
    --surface:     #122038;
    --border:      #1a2e4a;
    --border-hi:   #243d62;
    --cyan:        #38bdf8;
    --cyan-glow:   rgba(56,189,248,0.12);
    --violet:      #a78bfa;
    --violet-glow: rgba(167,139,250,0.12);
    --emerald:     #34d399;
    --amber:       #fbbf24;
    --rose:        #fb7185;
    --text:        #f0f4ff;
    --text-2:      #c5cfe8;
    --muted:       #7c8daa;
    --dim:         #1a2d48;
    --invert:      0;
"""

_LIGHT_VARS = """
    --bg:          #f4f7fb;
    --card:        #ffffff;
    --raised:      #eef2f8;
    --surface:     #e4eaf4;
    --border:      #d0daea;
    --border-hi:   #b8c8de;
    --cyan:        #0284c7;
    --cyan-glow:   rgba(2,132,199,0.10);
    --violet:      #6d28d9;
    --violet-glow: rgba(109,40,217,0.10);
    --emerald:     #059669;
    --amber:       #b45309;
    --rose:        #be123c;
    --text:        #0c1629;
    --text-2:      #253655;
    --muted:       #4e6080;
    --dim:         #dde6f0;
    --invert:      1;
"""


def apply_theme():
    light = st.session_state.get("light_mode", False)
    vars_ = _LIGHT_VARS if light else _DARK_VARS
    plot_bg = "rgba(240,244,251,0.6)" if light else "#080f1e"
    grid_c  = "#d0daea" if light else "#1a2e4a"
    sb_bg   = "#ffffff" if light else "#06101f"

    st.markdown(f"""
    <style>
    :root {{ {vars_} }}

    * {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', 'Inter', Helvetica, Arial, sans-serif !important; }}
    code, pre, .stCode * {{ font-family: 'JetBrains Mono', 'Fira Code', 'Cascadia Code', monospace !important; }}

    /* ── Layout ── */
    .stApp, .main, [data-testid="stAppViewContainer"] {{ background: var(--bg) !important; }}
    [data-testid="stSidebar"] {{ background: {sb_bg} !important; border-right: 1px solid var(--border); }}
    .block-container {{ padding: 1.5rem 2.5rem !important; max-width: 1540px; }}

    /* ── Base text ── */
    h1,h2,h3,h4,h5,h6 {{ color: var(--text) !important; letter-spacing: -0.025em; font-weight: 700; }}
    p, span, div, li   {{ color: var(--text) !important; }}
    .stMarkdown p       {{ color: var(--muted) !important; line-height: 1.7; }}
    small               {{ color: var(--muted) !important; font-size: 0.77rem; }}

    /* ── Inputs ── */
    .stTextInput input, .stTextArea textarea, .stNumberInput input {{
        background: var(--raised) !important;
        border: 1px solid var(--border) !important;
        border-radius: 9px !important;
        color: var(--text) !important;
        font-size: 0.9rem !important;
        padding: 0.55rem 0.9rem !important;
        transition: border-color .18s, box-shadow .18s;
    }}
    .stTextInput input:focus {{
        border-color: var(--cyan) !important;
        box-shadow: 0 0 0 3px var(--cyan-glow) !important;
        outline: none !important;
    }}
    .stSelectbox > div > div {{
        background: var(--raised) !important;
        border: 1px solid var(--border) !important;
        border-radius: 9px !important;
        color: var(--text) !important;
    }}
    .stSlider [data-testid="stSlider"] {{ padding: 0; }}

    /* ── Buttons ── */
    .stButton > button {{
        background: linear-gradient(135deg, var(--cyan), var(--violet)) !important;
        color: #fff !important;
        border: none !important;
        border-radius: 9px !important;
        font-weight: 600 !important;
        font-size: 0.85rem !important;
        padding: 0.5rem 1.4rem !important;
        letter-spacing: 0.015em;
        transition: opacity .15s, transform .1s, box-shadow .15s !important;
        box-shadow: 0 1px 3px rgba(0,0,0,.15);
    }}
    .stButton > button:hover {{
        opacity: 0.84 !important;
        transform: translateY(-1px) !important;
        box-shadow: 0 6px 20px var(--cyan-glow) !important;
    }}
    .stButton > button[kind="secondary"] {{
        background: var(--raised) !important;
        border: 1px solid var(--border-hi) !important;
        color: var(--muted) !important;
        box-shadow: none !important;
    }}
    .stButton > button[kind="secondary"]:hover {{
        border-color: var(--cyan) !important;
        color: var(--cyan) !important;
        box-shadow: none !important;
        opacity: 1 !important;
        transform: none !important;
    }}

    /* ── Metrics ── */
    [data-testid="stMetric"] {{
        background: var(--card);
        border: 1px solid var(--border);
        border-radius: 12px;
        padding: 1rem 1.2rem;
    }}
    [data-testid="stMetricValue"] {{
        color: var(--cyan) !important;
        font-size: 1.55rem !important;
        font-weight: 700 !important;
    }}
    [data-testid="stMetricLabel"] {{
        color: var(--muted) !important;
        font-size: 0.72rem !important;
        text-transform: uppercase;
        letter-spacing: 0.07em;
    }}

    /* ── Tabs ── */
    .stTabs [data-baseweb="tab-list"] {{
        gap: 0;
        border-bottom: 1px solid var(--border) !important;
        background: transparent !important;
        padding: 0;
    }}
    .stTabs [data-baseweb="tab"] {{
        background: transparent !important;
        color: var(--muted) !important;
        border: none !important;
        border-radius: 0 !important;
        padding: 0.55rem 1.1rem !important;
        font-size: 0.83rem !important;
        font-weight: 500 !important;
        margin: 0;
    }}
    .stTabs [aria-selected="true"] {{
        color: var(--cyan) !important;
        border-bottom: 2px solid var(--cyan) !important;
        background: transparent !important;
        font-weight: 600 !important;
    }}

    /* ── DataFrame ── */
    [data-testid="stDataFrame"] {{
        border-radius: 10px;
        overflow: hidden;
        border: 1px solid var(--border);
    }}
    [data-testid="stDataFrame"] thead th {{
        background: var(--raised) !important;
        color: var(--cyan) !important;
        font-size: 0.75rem !important;
        text-transform: uppercase;
        letter-spacing: 0.06em;
    }}
    [data-testid="stDataFrame"] tbody td {{
        font-size: 0.83rem !important;
        color: var(--text) !important;
    }}
    [data-testid="stDataFrame"] tbody tr:nth-child(even) td {{
        background: var(--raised) !important;
    }}

    /* ── Expander ── */
    .streamlit-expanderHeader {{
        background: var(--raised) !important;
        border: 1px solid var(--border);
        border-radius: 9px;
        color: var(--muted) !important;
        font-size: 0.84rem !important;
    }}
    .streamlit-expanderContent {{
        background: var(--card) !important;
        border: 1px solid var(--border);
        border-top: none;
        border-radius: 0 0 9px 9px;
    }}

    hr {{ border-color: var(--border) !important; margin: 1.1rem 0; }}

    /* ════════════ COMPONENTS ════════════ */

    /* Hero */
    .nr-hero {{
        background: linear-gradient(160deg, var(--raised) 0%, var(--surface) 55%, var(--raised) 100%);
        border: 1px solid var(--border);
        border-radius: 20px;
        padding: 2.75rem 3rem;
        margin-bottom: 2rem;
        position: relative;
        overflow: hidden;
    }}
    .nr-hero::after {{
        content: '';
        position: absolute; top: -120px; right: -80px;
        width: 340px; height: 340px; border-radius: 50%;
        background: radial-gradient(circle, var(--cyan-glow) 0%, transparent 65%);
        pointer-events: none;
    }}
    .nr-hero::before {{
        content: '';
        position: absolute; bottom: -80px; left: 25%;
        width: 240px; height: 240px; border-radius: 50%;
        background: radial-gradient(circle, var(--violet-glow) 0%, transparent 65%);
        pointer-events: none;
    }}
    .nr-hero-chip {{
        display: inline-block;
        background: var(--cyan-glow);
        color: var(--cyan) !important;
        border: 1px solid rgba(56,189,248,.25);
        border-radius: 20px;
        padding: 3px 12px;
        font-size: 0.72rem;
        font-weight: 600;
        letter-spacing: 0.12em;
        text-transform: uppercase;
        margin-bottom: 1rem;
    }}
    .nr-hero-title {{
        font-size: 2.4rem;
        font-weight: 800;
        line-height: 1.12;
        background: linear-gradient(90deg, var(--text) 20%, var(--cyan) 60%, var(--violet) 100%);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        background-clip: text;
        margin-bottom: 0.8rem;
    }}
    .nr-hero-sub {{
        font-size: 0.95rem;
        color: var(--muted) !important;
        line-height: 1.65;
        max-width: 660px;
    }}

    /* Search wrapper */
    .nr-search-wrap {{
        background: var(--card);
        border: 1px solid var(--border-hi);
        border-radius: 14px;
        padding: 1.25rem 1.5rem;
        margin-bottom: 1.5rem;
        box-shadow: 0 2px 12px rgba(0,0,0,.06);
    }}

    /* Resolution pill */
    .nr-resolve {{
        display: flex;
        align-items: center;
        gap: 1rem;
        background: var(--card);
        border: 1px solid var(--border-hi);
        border-left: 3px solid var(--cyan);
        border-radius: 10px;
        padding: 0.85rem 1.25rem;
        margin-bottom: 1.25rem;
    }}
    .nr-resolve-heading {{ font-size: 0.95rem; font-weight: 700; color: var(--cyan) !important; }}
    .nr-resolve-meta    {{ font-size: 0.74rem; color: var(--muted) !important; margin-top: 1px; }}

    /* Compound card */
    .cpd {{
        display: flex;
        gap: 1.25rem;
        background: var(--card);
        border: 1px solid var(--border);
        border-radius: 12px;
        padding: 1.1rem 1.4rem;
        margin-bottom: 0.55rem;
        transition: border-color .18s, box-shadow .18s;
        position: relative;
    }}
    .cpd:hover {{ border-color: var(--border-hi); box-shadow: 0 2px 16px rgba(0,0,0,.06); }}
    .cpd-left  {{ flex: 1; min-width: 0; }}
    .cpd-right {{ flex-shrink: 0; width: 160px; display: flex; flex-direction: column; gap: 5px; justify-content: center; }}
    .cpd-name  {{ font-size: 0.95rem; font-weight: 700; color: var(--text) !important; }}
    .cpd-mech  {{ font-size: 0.77rem; color: var(--muted) !important; margin-top: 2px; line-height: 1.4; }}

    /* Top-candidate card */
    .top-cpd {{
        background: var(--card);
        border: 1px solid var(--border-hi);
        border-radius: 14px;
        padding: 1.4rem;
        margin-bottom: 0.75rem;
        position: relative;
        overflow: hidden;
        transition: border-color .18s, box-shadow .18s;
    }}
    .top-cpd::before {{
        content: '';
        position: absolute; top: 0; left: 0; right: 0; height: 2px;
        background: linear-gradient(90deg, var(--cyan), var(--violet));
    }}
    .top-cpd:hover {{
        border-color: var(--cyan);
        box-shadow: 0 4px 24px var(--cyan-glow);
    }}

    /* Score bars */
    .sbar {{ display:flex; align-items:center; gap:7px; margin:2px 0; }}
    .sbar-l {{ width:76px; font-size:0.69rem; color:var(--muted) !important; flex-shrink:0; text-align:right; }}
    .sbar-t {{ flex:1; background:var(--dim); border-radius:3px; height:4px; overflow:hidden; }}
    .sbar-f {{ height:4px; border-radius:3px; transition:width .3s; }}
    .sbar-v {{ width:28px; font-size:0.69rem; color:var(--text-2) !important; flex-shrink:0; font-variant-numeric:tabular-nums; }}

    /* Badges */
    .b-ph  {{ display:inline-flex;align-items:center;background:var(--violet-glow);color:var(--violet);border:1px solid rgba(167,139,250,.3);border-radius:5px;padding:1px 8px;font-size:0.7rem;font-weight:700;letter-spacing:.03em; }}
    .b-hi  {{ display:inline-flex;align-items:center;background:rgba(52,211,153,.11);color:var(--emerald);border:1px solid rgba(52,211,153,.3);border-radius:5px;padding:1px 8px;font-size:0.7rem;font-weight:700; }}
    .b-mid {{ display:inline-flex;align-items:center;background:rgba(251,191,36,.11);color:var(--amber);border:1px solid rgba(251,191,36,.3);border-radius:5px;padding:1px 8px;font-size:0.7rem;font-weight:700; }}
    .b-lo  {{ display:inline-flex;align-items:center;background:rgba(251,113,133,.11);color:var(--rose);border:1px solid rgba(251,113,133,.3);border-radius:5px;padding:1px 8px;font-size:0.7rem;font-weight:700; }}

    /* Section label */
    .sec {{ font-size:0.64rem;font-weight:700;letter-spacing:.14em;text-transform:uppercase;color:var(--muted) !important;padding-bottom:0.45rem;border-bottom:1px solid var(--border);margin-bottom:0.85rem; }}

    /* Status dot */
    .srow {{ display:flex;align-items:center;gap:6px;font-size:0.79rem;margin-bottom:4px; }}
    .dot-on  {{ width:6px;height:6px;border-radius:50%;background:var(--emerald);flex-shrink:0; }}
    .dot-off {{ width:6px;height:6px;border-radius:50%;background:var(--rose);flex-shrink:0; }}

    /* Evidence card */
    .ev {{
        background: var(--raised);
        border: 1px solid var(--border);
        border-radius: 9px;
        padding: 0.9rem 1.1rem;
        margin-bottom: 0.45rem;
    }}
    .ev-title {{ font-size:0.85rem;font-weight:600;color:var(--text) !important;line-height:1.4;margin-bottom:3px; }}
    .ev-meta  {{ font-size:0.73rem;color:var(--muted) !important; }}
    .ev-link  {{ font-size:0.73rem;color:var(--cyan) !important; }}

    /* PBPK metric pill */
    .pk-pill {{
        background: var(--raised);
        border: 1px solid var(--border);
        border-radius: 10px;
        padding: 0.8rem 0.5rem;
        text-align: center;
    }}
    .pk-val {{ font-size:1.05rem;font-weight:700;color:var(--cyan) !important; }}
    .pk-lbl {{ font-size:0.64rem;color:var(--muted) !important;text-transform:uppercase;letter-spacing:.06em;margin-top:2px; }}

    /* Info / step card */
    .step-card {{
        background: var(--card);
        border: 1px solid var(--border);
        border-left: 3px solid var(--cyan);
        border-radius: 10px;
        padding: 1rem 1.25rem;
    }}
    .step-num   {{ font-size:0.65rem;font-weight:700;color:var(--cyan) !important;letter-spacing:.1em;text-transform:uppercase;margin-bottom:.2rem; }}
    .step-title {{ font-size:0.88rem;font-weight:700;color:var(--text) !important;margin-bottom:.3rem; }}
    .step-body  {{ font-size:0.77rem;color:var(--muted) !important;line-height:1.55; }}

    /* Row KV pair */
    .kv {{
        display: flex;
        justify-content: space-between;
        background: var(--raised);
        border-radius: 7px;
        padding: 5px 13px;
        margin-bottom: 3px;
        font-size: 0.82rem;
    }}
    .kv-k {{ color: var(--muted) !important; }}
    .kv-v {{ color: var(--text) !important; font-weight: 600; }}
    .kv-fail {{ color: var(--rose) !important; font-weight: 700; font-size: 0.7rem; }}
    .kv-pass {{ color: var(--emerald) !important; font-weight: 700; font-size: 0.7rem; }}

    /* Toggle row */
    .theme-row {{
        display: flex;
        align-items: center;
        justify-content: space-between;
        padding: 0.5rem 0;
        font-size: 0.82rem;
        color: var(--muted) !important;
    }}
    </style>
    """, unsafe_allow_html=True)


# ─── Helpers ──────────────────────────────────────────────────────────────────

def _phase_lbl(phase) -> str:
    try: p = int(float(phase or 0))
    except Exception: p = 0
    return {4:"FDA Approved",3:"Phase III",2:"Phase II",1:"Phase I"}.get(p,"Preclinical")

def _phase_b(phase) -> str:
    return f"<span class='b-ph'>{_phase_lbl(phase)}</span>"

def _score_b(s: float) -> str:
    cls = "b-hi" if s >= .60 else "b-mid" if s >= .35 else "b-lo"
    return f"<span class='{cls}'>{s:.0%}</span>"

def _bar_col(v: float) -> str:
    return "var(--emerald)" if v >= .60 else "var(--amber)" if v >= .35 else "var(--rose)"

def _score_bars(bd: Dict) -> str:
    labels = [("Indication","indication_score"),("Target","target_score"),
              ("Activity","activity_score"),("Network","network_score")]
    out = ""
    for lbl, key in labels:
        v = float(bd.get(key) or 0)
        out += (f"<div class='sbar'>"
                f"<span class='sbar-l'>{lbl}</span>"
                f"<div class='sbar-t'><div class='sbar-f' style='width:{v*100:.1f}%;background:{_bar_col(v)};'></div></div>"
                f"<span class='sbar-v'>{v:.0%}</span>"
                f"</div>")
    return out

def _status(label: str, ok: bool):
    d = "dot-on" if ok else "dot-off"
    st.markdown(f"<div class='srow'><span class='{d}'></span><span style='color:var(--muted);'>{label}</span></div>",
                unsafe_allow_html=True)

def _sec(t: str):
    st.markdown(f"<div class='sec'>{t}</div>", unsafe_allow_html=True)

def _kv(k: str, v: str, fail: bool = None):
    mark = ""
    if fail is True:  mark = "<span class='kv-fail'>FAIL</span>"
    if fail is False: mark = "<span class='kv-pass'>PASS</span>"
    st.markdown(f"<div class='kv'><span class='kv-k'>{k}</span><span class='kv-v'>{v}</span>{mark}</div>",
                unsafe_allow_html=True)

def _plot_style(fig, light: bool):
    bg   = "rgba(240,244,251,0.6)" if light else "#0d1729"
    grid = "#d0daea" if light else "#1a2e4a"
    fc   = "#0c1629" if light else "#f0f4ff"
    fig.update_layout(paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor=bg,
                      font_color=fc, legend=dict(font=dict(color=fc)),
                      xaxis=dict(gridcolor=grid), yaxis=dict(gridcolor=grid))
    return fig


# ─── Session init ─────────────────────────────────────────────────────────────

def _init():
    defs = {"page":"dashboard","disease_query":"","selected_compound":None,
            "light_mode":False}
    for k, v in defs.items():
        if k not in st.session_state: st.session_state[k] = v


# ─── Sidebar ──────────────────────────────────────────────────────────────────

def _sidebar():
    light = st.session_state.get("light_mode", False)
    with st.sidebar:
        st.markdown(
            "<div style='padding:1.1rem 0 1.5rem;'>"
            "<div style='font-size:1.3rem;font-weight:800;"
            "background:linear-gradient(90deg,var(--cyan),var(--violet));"
            "-webkit-background-clip:text;-webkit-text-fill-color:transparent;"
            "background-clip:text;'>NeuroRepurpose</div>"
            "<div style='font-size:0.7rem;color:var(--muted);margin-top:2px;letter-spacing:.05em;'>"
            "Drug Repurposing Intelligence</div></div>",
            unsafe_allow_html=True,
        )

        # Theme toggle
        st.markdown("<div style='font-size:0.7rem;color:var(--muted);text-transform:uppercase;"
                    "letter-spacing:.1em;margin-bottom:4px;'>Appearance</div>", unsafe_allow_html=True)
        col_dk, col_lt = st.columns(2)
        with col_dk:
            if st.button("Dark", key="thm_dk", use_container_width=True,
                         type="secondary" if light else "primary"):
                st.session_state.light_mode = False
                st.rerun()
        with col_lt:
            if st.button("Light", key="thm_lt", use_container_width=True,
                         type="primary" if light else "secondary"):
                st.session_state.light_mode = True
                st.rerun()

        st.markdown("<div style='height:.5rem;'></div>", unsafe_allow_html=True)
        st.markdown("---")
        _sec("Navigation")

        nav = [("dashboard","Dashboard"),("discover","Discover"),
               ("graph","Knowledge Graph"),("database","Data Explorer")]
        for pid, lbl in nav:
            active = st.session_state.page == pid
            t = "primary" if active else "secondary"
            if st.button(lbl, key=f"nav_{pid}", use_container_width=True, type=t):
                st.session_state.page = pid
                st.rerun()

        st.markdown("---")
        _sec("Data Sources")
        _status("ChEMBL 33",         DB_AVAILABLE)
        _status("MeSH Ontology",     MESH_AVAILABLE)
        _status("Scoring Engine",    SCORER_AVAILABLE)
        _status("PBPK Simulator",    PBPK_AVAILABLE)
        _status("Molecular Docking", DOCKING_AVAILABLE)
        _status("Quantum / RDKit",   QUANTUM_AVAILABLE)
        _status("3D Viewer",         VIZ_3D_AVAILABLE)

        if DB_AVAILABLE:
            st.markdown("---")
            _sec("Live Counts")
            try:
                stats = get_stats()
                for k, v in stats.items():
                    st.markdown(
                        f"<div style='display:flex;justify-content:space-between;"
                        f"font-size:0.78rem;margin-bottom:2px;'>"
                        f"<span style='color:var(--muted);'>{k.replace('_',' ').title()}</span>"
                        f"<span style='color:var(--text);font-weight:600;'>{v:,}</span></div>",
                        unsafe_allow_html=True,
                    )
            except Exception: pass

        st.markdown("---")
        st.markdown("<div style='font-size:0.68rem;color:var(--muted);line-height:1.7;'>"
                    "ChEMBL 33 · MeSH NLM · RDKit<br>"
                    "ClinicalTrials.gov · PubMed</div>", unsafe_allow_html=True)


# ─── Dashboard ────────────────────────────────────────────────────────────────

def page_dashboard():
    light = st.session_state.get("light_mode", False)

    st.markdown(
        "<div class='nr-hero'>"
        "<div class='nr-hero-chip'>Drug Repurposing Platform</div>"
        "<div class='nr-hero-title'>NeuroRepurpose Intelligence</div>"
        "<div class='nr-hero-sub'>"
        "Identify new therapeutic applications for approved compounds using ChEMBL 33, "
        "disease knowledge graphs, MeSH ontology, multi-signal evidence scoring, "
        "PBPK pharmacokinetics, and real clinical trial data."
        "</div>"
        "</div>",
        unsafe_allow_html=True,
    )

    col_q, col_btn = st.columns([6, 1])
    with col_q:
        q = st.text_input("hero_q",
            placeholder="Search any disease — Alzheimer, MS, epilepsy, depression, ALS...",
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
    with c4: st.metric("Graph Nodes",     f"{stats.get('hetionet_nodes',0):,}")
    with c5: st.metric("MeSH Diseases",   f"{stats.get('mesh_diseases',0):,}")

    if not DB_AVAILABLE:
        st.warning("Database unavailable. Run: `python database/importer.py`")
        return

    st.markdown("---")
    st.markdown("### Disease Areas")
    st.caption("Loaded live from the database. Click any to run drug discovery.")
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
    st.markdown("### Analysis Pipeline")
    c1,c2,c3,c4 = st.columns(4)
    steps = [
        ("01","MeSH Disease Resolution",
         "Free-text matched against the full MeSH ontology — headings, entry terms, "
         "abbreviation aliases (AD, MS, ALS, ADHD...). Expanded to parent and child concepts."),
        ("02","Evidence-Based Scoring",
         "Four signals ranked per compound: clinical indication evidence (40%), "
         "target mechanistic overlap (30%), pChEMBL activity potency (20%), "
         "and knowledge graph path count (10%)."),
        ("03","PBPK Pharmacokinetics",
         "One-compartment PK model predicts plasma, liver, brain, and kidney "
         "concentration-time curves from MW and LogP. Cmax, AUC, t½, safety margin."),
        ("04","Clinical Evidence",
         "Live ClinicalTrials.gov and PubMed lookup per drug-disease pair. "
         "Real NCT IDs, enrollment status, phase, and peer-reviewed publication links."),
    ]
    for col, (num, title, body) in zip([c1,c2,c3,c4], steps):
        with col:
            st.markdown(
                f"<div class='step-card'>"
                f"<div class='step-num'>Step {num}</div>"
                f"<div class='step-title'>{title}</div>"
                f"<div class='step-body'>{body}</div>"
                f"</div>",
                unsafe_allow_html=True,
            )


# ─── Discover ─────────────────────────────────────────────────────────────────

def page_discover():
    light = st.session_state.get("light_mode", False)
    st.markdown("## Discover")
    st.caption("Enter any disease name. The platform resolves abbreviations, synonyms, and related MeSH concepts automatically.")

    col_inp, col_btn = st.columns([6, 1])
    with col_inp:
        inp = st.text_input("disc_inp",
            value=st.session_state.disease_query,
            placeholder="e.g. Parkinson disease · ALS · MS · schizophrenia · type 2 diabetes...",
            label_visibility="collapsed")
    with col_btn:
        run = st.button("Analyse", use_container_width=True)

    if inp and len(inp) >= 2 and MESH_AVAILABLE:
        sugs = suggest_diseases(inp, limit=5)
        if sugs and inp.strip().lower() not in [s.lower() for s in sugs]:
            st.markdown("<div style='font-size:0.73rem;color:var(--muted);margin:3px 0 2px;'>Suggestions</div>",
                        unsafe_allow_html=True)
            scols = st.columns(min(len(sugs), 5))
            for i, sg in enumerate(sugs):
                with scols[i]:
                    if st.button(sg, key=f"sug_{i}", use_container_width=True, type="secondary"):
                        st.session_state.disease_query = sg
                        st.rerun()

    if run and inp.strip():
        st.session_state.disease_query = inp.strip()
    if not st.session_state.disease_query:
        st.info("Enter a disease name above to begin.")
        return

    with st.spinner("Resolving disease and retrieving compounds..."):
        try:
            resolved, expanded_ids, compounds = resolve_and_fetch(st.session_state.disease_query)
        except Exception as e:
            st.error(f"Query failed: {e}"); return

    if not resolved:
        st.warning(f"No MeSH match for '{st.session_state.disease_query}'.")
        if not MESH_AVAILABLE:
            st.info("MeSH table empty — run: `python database/mesh_importer.py`")
        return

    primary = resolved[0]
    st.markdown(
        f"<div class='nr-resolve'>"
        f"<div><div class='nr-resolve-heading'>{primary.get('heading', st.session_state.disease_query)}</div>"
        f"<div class='nr-resolve-meta'>"
        f"MeSH {primary.get('mesh_id','N/A')} &nbsp;·&nbsp; "
        f"{len(resolved)} record(s) &nbsp;·&nbsp; "
        f"{len(expanded_ids)} expanded IDs (parents + children)"
        f"</div></div></div>",
        unsafe_allow_html=True,
    )

    if not compounds:
        st.warning("No compounds found for this disease.")
        return

    with st.spinner(f"Scoring {len(compounds)} candidates across 4 evidence signals..."):
        try:
            scored = score_compound_list(list(compounds), expanded_ids)
        except Exception:
            scored = compounds
            for c in scored:
                c.setdefault("score", float(c.get("max_phase") or 0) / 4)
                c.setdefault("score_breakdown", {})

    # Top 3
    top3 = [c for c in scored if float(c.get("max_phase") or 0) >= 3][:3] or scored[:3]
    if top3:
        _sec("Top Candidates")
        tcols = st.columns(len(top3))
        for i, c in enumerate(top3):
            _top_card(c, tcols[i], i)

    st.markdown("---")

    fc1, fc2, fc3 = st.columns([1, 1, 2])
    with fc1: min_ph = st.selectbox("Min Phase", [0,1,2,3,4], key="d_ph")
    with fc2: show_n = st.slider("Show", 5, min(80,len(scored)), min(30,len(scored)), key="d_n")
    with fc3: sort_by = st.selectbox("Sort", ["Score","Phase","Name"], key="d_sort")

    filt = [c for c in scored if float(c.get("max_phase") or 0) >= min_ph]
    if sort_by == "Name": filt.sort(key=lambda x: (x.get("name") or "").lower())
    elif sort_by == "Phase": filt.sort(key=lambda x: float(x.get("max_phase") or 0), reverse=True)
    filt = filt[:show_n]

    if len(filt) >= 4:
        with st.expander("Phase distribution", expanded=False):
            _phase_donut(filt, light)

    _sec(f"{len(filt)} Repurposing Candidates — ranked by evidence score")
    for i, c in enumerate(filt):
        _compound_row(c, i)

    if filt:
        df_dl = pd.DataFrame([{
            "Name": c.get("name"), "ChEMBL": c.get("chembl_id"),
            "Phase": c.get("max_phase"),
            "Score": round(float(c.get("score") or 0), 4),
            "Indication Score": round(float((c.get("score_breakdown") or {}).get("indication_score") or 0), 4),
            "Target Score":     round(float((c.get("score_breakdown") or {}).get("target_score") or 0), 4),
            "Activity Score":   round(float((c.get("score_breakdown") or {}).get("activity_score") or 0), 4),
            "Network Score":    round(float((c.get("score_breakdown") or {}).get("network_score") or 0), 4),
            "Mechanisms": c.get("mechanisms",""), "SMILES": c.get("smiles",""),
        } for c in filt])
        st.download_button("Download CSV", data=df_dl.to_csv(index=False),
                           file_name=f"repurposing_{primary.get('mesh_id','x')}.csv",
                           mime="text/csv")


def _top_card(c: Dict, col, idx: int):
    name  = c.get("name") or "Unknown"
    score = float(c.get("score") or 0)
    bd    = c.get("score_breakdown") or {}
    mechs = (c.get("mechanisms") or "—")[:90]
    with col:
        st.markdown(
            f"<div class='top-cpd'>"
            f"<div style='display:flex;justify-content:space-between;align-items:flex-start;"
            f"margin-bottom:.55rem;'>"
            f"<span style='font-size:.93rem;font-weight:700;color:var(--text);'>{name}</span>"
            f"<span>{_score_b(score)}</span>"
            f"</div>"
            f"<div style='margin-bottom:.4rem;'>{_phase_b(c.get('max_phase'))}</div>"
            f"<div style='font-size:.74rem;color:var(--muted);line-height:1.4;margin-bottom:.6rem;'>{mechs}</div>"
            f"{_score_bars(bd)}"
            f"</div>",
            unsafe_allow_html=True,
        )
        if st.button("Full Analysis", key=f"ta_{idx}_{name[:5]}", use_container_width=True):
            st.session_state.selected_compound = c
            st.session_state.page = "analysis"
            st.rerun()


def _compound_row(c: Dict, idx: int):
    name    = c.get("name") or "Unknown"
    score   = float(c.get("score") or 0)
    phase   = c.get("max_phase", 0)
    mechs   = (c.get("mechanisms") or "—")[:140]
    targets = (c.get("targets") or "—")[:110]
    bd      = c.get("score_breakdown") or {}

    col_main, col_act = st.columns([7, 1])
    with col_main:
        st.markdown(
            f"<div class='cpd'>"
            f"<div class='cpd-left'>"
            f"<div style='display:flex;align-items:center;gap:.5rem;margin-bottom:.3rem;'>"
            f"<span class='cpd-name'>{name}</span>{_phase_b(phase)}"
            f"</div>"
            f"<div class='cpd-mech'><b style='color:var(--text-2);'>Mechanism:</b> {mechs}</div>"
            f"<div class='cpd-mech'><b style='color:var(--text-2);'>Targets:</b> {targets}</div>"
            f"</div>"
            f"<div class='cpd-right'>"
            f"<div style='text-align:right;margin-bottom:.35rem;'>{_score_b(score)}</div>"
            f"{_score_bars(bd)}"
            f"</div>"
            f"</div>",
            unsafe_allow_html=True,
        )
    with col_act:
        st.markdown("<div style='height:1.15rem;'></div>", unsafe_allow_html=True)
        if st.button("Analyse", key=f"ca_{idx}_{name[:6]}", use_container_width=True):
            st.session_state.selected_compound = c
            st.session_state.page = "analysis"
            st.rerun()


def _phase_donut(compounds: List[Dict], light: bool):
    pm = {4:"FDA Approved",3:"Phase III",2:"Phase II",1:"Phase I",0:"Preclinical"}
    cnt: Dict = {}
    for c in compounds:
        lbl = pm.get(int(float(c.get("max_phase") or 0)), "Preclinical")
        cnt[lbl] = cnt.get(lbl, 0) + 1
    df = pd.DataFrame(list(cnt.items()), columns=["Phase","Count"])
    fig = px.pie(df, values="Count", names="Phase", hole=.52,
                 color_discrete_sequence=["#38bdf8","#a78bfa","#34d399","#fbbf24","#fb7185"])
    fig = _plot_style(fig, light)
    fig.update_layout(margin=dict(t=10,b=10,l=10,r=10),
                      legend=dict(font=dict(size=11)))
    st.plotly_chart(fig, use_container_width=True)


# ─── Compound Analysis ────────────────────────────────────────────────────────

def page_analysis():
    light = st.session_state.get("light_mode", False)
    c = st.session_state.selected_compound
    if not c:
        st.info("Select a compound from Discover.")
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
    # Per-compound state keys so switching compound resets results
    pbpk_key  = f"pbpk_{cid}"
    dock_key  = f"dock_{cid}"

    col_h, col_meta = st.columns([4, 1])
    with col_h:
        st.markdown(f"## {name}")
        if chembl_id:
            st.caption(f"ChEMBL ID: {chembl_id}")
        mech = c.get("mechanisms") or "—"
        st.markdown(
            f"<div style='font-size:.86rem;color:var(--muted);margin-top:.2rem;'>"
            f"Mechanism of action: <span style='color:var(--text-2);'>{mech[:280]}</span></div>",
            unsafe_allow_html=True,
        )
    with col_meta:
        st.markdown(
            f"<div style='text-align:right;padding-top:.4rem;'>"
            f"{_phase_b(phase)}"
            f"<br><br>{_score_b(score)} overall score"
            f"</div>",
            unsafe_allow_html=True,
        )
        if bd:
            st.markdown(f"<div style='margin-top:.6rem;'>{_score_bars(bd)}</div>",
                        unsafe_allow_html=True)

    st.markdown("---")

    # Build tabs
    tab_labels = ["Properties","Bioactivity","Indications","Clinical Evidence",
                  "PBPK Simulation","Docking"]
    if QUANTUM_AVAILABLE: tab_labels.append("Optimization")
    if VIZ_3D_AVAILABLE:  tab_labels.append("3D Structure")
    tabs = st.tabs(tab_labels)
    ti = 0

    # ── Properties ──────────────────────────────────────────────────────────
    with tabs[ti]:
        _properties_tab(c, smiles, light)
    ti += 1

    # ── Bioactivity (Targets + Activities) ───────────────────────────────────
    with tabs[ti]:
        _bioactivity_tab(cid, light)
    ti += 1

    # ── Indications ─────────────────────────────────────────────────────────
    with tabs[ti]:
        inds = get_compound_indications(int(cid)) if (DB_AVAILABLE and cid) else []
        if inds:
            df_i = pd.DataFrame(inds)
            keep = [col for col in ["disease","mesh_id","max_phase","source"] if col in df_i.columns]
            st.dataframe(df_i[keep], use_container_width=True)
        elif c.get("indications"):
            for ind in str(c["indications"]).split(";"):
                st.markdown(f"- {ind.strip()}")
        else:
            st.info("No indication data in database.")
    ti += 1

    # ── Clinical Evidence ────────────────────────────────────────────────────
    with tabs[ti]:
        _clinical_tab(name)
    ti += 1

    # ── PBPK ────────────────────────────────────────────────────────────────
    with tabs[ti]:
        _pbpk_tab(c, smiles, light, pbpk_key)
    ti += 1

    # ── Docking ─────────────────────────────────────────────────────────────
    with tabs[ti]:
        _docking_tab(name, smiles, cid, light, dock_key)
    ti += 1

    # ── Quantum Optimization ─────────────────────────────────────────────────
    if QUANTUM_AVAILABLE:
        with tabs[ti]:
            _quantum_tab(name, smiles, light)
        ti += 1

    # ── 3D Structure ─────────────────────────────────────────────────────────
    if VIZ_3D_AVAILABLE:
        with tabs[ti]:
            _3d_tab(name, smiles, cid)
        ti += 1


# ── Tab: Properties ───────────────────────────────────────────────────────────

def _properties_tab(c: Dict, smiles: str, light: bool):
    mw = c.get("mw"); logp = c.get("alogp"); psa = c.get("psa")
    hba = c.get("hba"); hbd = c.get("hbd"); rtb = c.get("rtb")
    ro5 = c.get("ro5_violations"); qed = None

    if smiles:
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
                rtb  = rtb  or rdMolDescriptors.CalcNumRotatableBonds(mol)
                qed  = round(rdQED.qed(mol), 3)
        except Exception: pass

    c1,c2,c3,c4,c5,c6 = st.columns(6)
    def _m(col, lbl, val): col.metric(lbl, val)
    _m(c1,"MW (Da)",      f"{mw:.1f}" if mw else "N/A")
    _m(c2,"LogP",         f"{logp:.2f}" if logp is not None else "N/A")
    _m(c3,"TPSA (sq A)",  f"{psa:.1f}" if psa else "N/A")
    _m(c4,"QED",          f"{qed:.3f}" if qed is not None else "N/A")
    _m(c5,"HBA / HBD",    f"{hba} / {hbd}" if hba is not None else "N/A")
    _m(c6,"Rotatable Bonds", str(rtb) if rtb is not None else "N/A")

    st.markdown("<div style='height:.5rem;'></div>", unsafe_allow_html=True)

    viols = sum([bool(mw and mw>500), bool(logp is not None and logp>5),
                 bool(hbd is not None and hbd>5), bool(hba is not None and hba>10)])
    col = "var(--emerald)" if viols==0 else "var(--amber)" if viols==1 else "var(--rose)"
    txt  = "Passes Lipinski Rule of 5 — good oral bioavailability potential" \
           if viols==0 else f"Lipinski: {viols} violation(s)"
    st.markdown(f"<div style='font-size:.84rem;color:{col};margin-bottom:.85rem;'>{txt}</div>",
                unsafe_allow_html=True)

    col_l, col_r = st.columns(2)
    with col_l:
        _sec("Lipinski Rule of 5")
        _kv("Molecular weight", f"{mw:.1f} Da" if mw else "N/A", bool(mw and mw>500))
        _kv("LogP (lipophilicity)", f"{logp:.2f}" if logp is not None else "N/A",
            bool(logp is not None and logp>5))
        _kv("H-bond donors", str(hbd) if hbd is not None else "N/A",
            bool(hbd is not None and hbd>5))
        _kv("H-bond acceptors", str(hba) if hba is not None else "N/A",
            bool(hba is not None and hba>10))
    with col_r:
        _sec("Drug-likeness")
        _kv("QED score",           f"{qed:.3f}" if qed else "N/A")
        _kv("TPSA",                f"{psa:.1f} sq A" if psa else "N/A")
        _kv("Rotatable bonds",     str(rtb) if rtb is not None else "N/A")
        _kv("Ro5 violations (DB)", str(int(ro5)) if ro5 is not None else "N/A")

    if smiles:
        st.markdown("<div style='margin-top:1rem;'></div>", unsafe_allow_html=True)
        _sec("SMILES")
        st.code(smiles, language=None)


# ── Tab: Bioactivity ─────────────────────────────────────────────────────────

def _bioactivity_tab(cid, light: bool):
    tgts = get_compound_targets(int(cid)) if (DB_AVAILABLE and cid) else []
    acts = get_compound_activities(int(cid)) if (DB_AVAILABLE and cid) else []

    col_t, col_a = st.columns(2)

    with col_t:
        _sec(f"Protein Targets ({len(tgts)})")
        if tgts:
            df_t = pd.DataFrame(tgts)
            keep = [col for col in ["gene_symbol","name","target_type","mechanism",
                                     "action_type","confidence"] if col in df_t.columns]
            st.dataframe(df_t[keep], use_container_width=True, height=260)

            if "confidence" in df_t.columns and len(df_t) >= 2:
                top = df_t.dropna(subset=["confidence"]).sort_values("confidence").tail(10)
                nc  = "gene_symbol" if "gene_symbol" in top.columns else "name"
                fig = px.bar(top, x="confidence", y=nc, orientation="h",
                             color="confidence",
                             color_continuous_scale=[[0,"#fb7185"],[.5,"#fbbf24"],[1,"#34d399"]])
                fig.update_layout(title="Target confidence", margin=dict(t=30,b=0))
                fig = _plot_style(fig, light)
                st.plotly_chart(fig, use_container_width=True)
        else:
            st.info("No targets found in database for this compound.")

    with col_a:
        _sec(f"pChEMBL Activities ({len(acts)})")
        if acts:
            df_a = pd.DataFrame(acts)
            keep = [col for col in ["gene_symbol","target_name","activity_type",
                                     "pchembl_value","standard_value","standard_units"]
                    if col in df_a.columns]
            st.dataframe(df_a[keep], use_container_width=True, height=260)

            if "pchembl_value" in df_a.columns:
                dp = df_a.dropna(subset=["pchembl_value"]) \
                         .sort_values("pchembl_value", ascending=False).head(12)
                if not dp.empty:
                    nc = "gene_symbol" if "gene_symbol" in dp.columns else "target_name"
                    fig = px.bar(dp.sort_values("pchembl_value"), x="pchembl_value", y=nc,
                                 orientation="h",
                                 color="pchembl_value",
                                 color_continuous_scale=[[0,"#a78bfa"],[.5,"#fbbf24"],[1,"#34d399"]])
                    fig.add_vline(x=6, line_dash="dash", line_color="#fbbf24",
                                  annotation=dict(text="1 µM",font_color="#fbbf24"))
                    fig.add_vline(x=8, line_dash="dash", line_color="#34d399",
                                  annotation=dict(text="10 nM",font_color="#34d399"))
                    fig.update_layout(title="pChEMBL (higher = more potent)",
                                      margin=dict(t=30,b=0))
                    fig = _plot_style(fig, light)
                    st.plotly_chart(fig, use_container_width=True)
        else:
            st.info("No quantitative activity data available.")


# ── Tab: Clinical Evidence ────────────────────────────────────────────────────

def _clinical_tab(drug_name: str):
    dq = st.session_state.get("disease_query", "")
    if not dq:
        st.info("Run a disease search first to load clinical evidence.")
        return

    col_t, col_p = st.columns(2)

    with col_t:
        _sec("Clinical Trials  (ClinicalTrials.gov)")
        with st.spinner("Fetching..."):
            trials = fetch_trials(drug_name, dq, n=8)
        if trials:
            for t in trials:
                sc = {"RECRUITING":"var(--emerald)","COMPLETED":"var(--cyan)",
                      "ACTIVE, NOT RECRUITING":"var(--amber)"}.get(t.get("status",""), "var(--muted)")
                nct = t.get("nct_id","")
                url = t.get("url","#")
                st.markdown(
                    f"<div class='ev'>"
                    f"<div class='ev-title'>{t.get('title','—')}</div>"
                    f"<div class='ev-meta'>"
                    f"<span style='color:{sc};font-weight:600;'>{t.get('status','—')}</span>"
                    f" &nbsp;·&nbsp; {t.get('phase','N/A')}"
                    f" &nbsp;·&nbsp; <a class='ev-link' href='{url}' target='_blank'>{nct}</a>"
                    f"</div></div>",
                    unsafe_allow_html=True,
                )
        else:
            st.info(f"No trials found for {drug_name} in {dq}.")

    with col_p:
        _sec("Publications  (PubMed / NCBI)")
        with st.spinner("Fetching..."):
            papers = fetch_papers(drug_name, dq, n=8)
        if papers:
            for p in papers:
                st.markdown(
                    f"<div class='ev'>"
                    f"<div class='ev-title'>{p.get('title','—')}</div>"
                    f"<div class='ev-meta'>"
                    f"{p.get('journal','—')} &nbsp;·&nbsp; {p.get('year','')}"
                    f" &nbsp;·&nbsp; <a class='ev-link' href='{p.get('url','#')}' "
                    f"target='_blank'>PMID {p.get('pmid','')}</a>"
                    f"</div></div>",
                    unsafe_allow_html=True,
                )
        else:
            st.info(f"No PubMed papers found for {drug_name} + {dq}.")


# ── Tab: PBPK ─────────────────────────────────────────────────────────────────

def _pbpk_tab(c: Dict, smiles: str, light: bool, state_key: str):
    _sec("PBPK Pharmacokinetic Simulation")
    if not PBPK_AVAILABLE:
        st.info("PBPK unavailable — check `services/pbpk_simulation.py`.")
        return

    name = c.get("name") or "Unknown"
    mw   = float(c.get("mw") or 0)
    logp = float(c.get("alogp") or 0)

    if smiles and not (mw and logp):
        try:
            from rdkit import Chem
            from rdkit.Chem import Descriptors
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                mw   = mw   or round(Descriptors.ExactMolWt(mol), 1)
                logp = logp or round(Descriptors.MolLogP(mol), 2)
        except Exception: pass

    if not mw:
        st.info("Molecular weight unavailable — cannot run PBPK simulation.")
        return

    st.markdown(
        f"<div style='font-size:.82rem;color:var(--muted);margin-bottom:.75rem;'>"
        f"Using MW = {mw:.1f} Da, LogP = {logp:.2f}. "
        f"Adjust dose and route below then click Run.</div>",
        unsafe_allow_html=True,
    )

    col1,col2,col3,col4 = st.columns(4)
    with col1: dose    = st.number_input("Dose (mg)",       1.0, 2000.0, 100.0, 10.0, key=f"pbpk_dose_{state_key}")
    with col2: route   = st.selectbox("Route",              ["oral","IV"],               key=f"pbpk_route_{state_key}")
    with col3: dur     = st.number_input("Duration (h)",    4.0,  168.0,  24.0,  4.0,  key=f"pbpk_dur_{state_key}")
    with col4: binding = st.number_input("Binding (kcal/mol)", -20.0, 0.0, -8.0, 0.5,  key=f"pbpk_ba_{state_key}")

    dq = st.session_state.get("disease_query", "Unknown")

    if st.button("Run Simulation", key=f"run_pbpk_{state_key}"):
        with st.spinner("Running compartmental PK model..."):
            try:
                sim = PBPKSimulator(disease_name=dq)
                result = sim.simulate_drug_exposure(
                    drug_name=name, molecular_weight=mw, logp=logp,
                    dose_mg=dose, route=route, duration_hours=dur,
                    binding_affinity=binding,
                )
                st.session_state[state_key] = result
            except Exception as e:
                st.error(f"Simulation failed: {e}")

    result = st.session_state.get(state_key)
    if not result or not result.get("success"):
        return

    pk   = result.get("pk_metrics", {})
    adme = result.get("adme_parameters", {})

    # PK metrics row
    pkc = st.columns(4)
    for col, (lbl, val, unit) in zip(pkc, [
        ("Cmax",  f"{pk.get('cmax_ng_ml',0):.1f}",   "ng/mL"),
        ("Tmax",  f"{pk.get('tmax_hours',0):.2f}",    "h"),
        ("AUC",   f"{pk.get('auc_ng_h_ml',0):.0f}",  "ng·h/mL"),
        ("t½",    f"{pk.get('t_half_hours',0):.2f}",  "h"),
    ]):
        with col:
            st.markdown(f"<div class='pk-pill'><div class='pk-val'>{val}</div>"
                        f"<div class='pk-lbl'>{lbl} ({unit})</div></div>",
                        unsafe_allow_html=True)

    st.markdown("<div style='height:.6rem;'></div>", unsafe_allow_html=True)

    # Concentration-time chart
    times = result.get("time_hours", [])
    fig = go.Figure()
    for key, label, color in [
        ("plasma_concentration_ng_ml","Plasma","#38bdf8"),
        ("liver_concentration_ng_ml", "Liver", "#fbbf24"),
        ("brain_concentration_ng_ml", "Brain", "#a78bfa"),
    ]:
        vals = result.get(key, [])
        if vals:
            fig.add_trace(go.Scatter(x=times, y=vals, name=label,
                                     line=dict(color=color, width=2.2)))
    fig.update_layout(
        title=f"Concentration-time — {name}  ({route}, {dose:.0f} mg)",
        xaxis_title="Time (h)", yaxis_title="Concentration (ng/mL)",
    )
    fig = _plot_style(fig, light)
    st.plotly_chart(fig, use_container_width=True)

    col_a, col_s = st.columns(2)
    with col_a:
        _sec("ADME Parameters")
        for lbl, val in [
            ("Bioavailability (F)",    f"{adme.get('f',0):.0%}"),
            ("Vol. distribution (Vd)", f"{adme.get('vd',0):.2f} L/kg"),
            ("Clearance (CL)",         f"{adme.get('cl',0):.3f} L/h/kg"),
            ("Fraction unbound (fu)",  f"{adme.get('fu',0):.0%}"),
            ("Brain partition (Kp)",   f"{adme.get('kp_brain',0):.2f}"),
            ("Liver partition (Kp)",   f"{adme.get('kp_liver',0):.2f}"),
        ]:
            _kv(lbl, val)

    with col_s:
        safety = result.get("safety_assessment", {})
        _sec("Safety Assessment")
        margin = safety.get("safety_margin","Unknown")
        mc = "var(--emerald)" if margin=="Good" else "var(--amber)"
        st.markdown(f"<div style='font-size:.88rem;font-weight:700;color:{mc};"
                    f"margin-bottom:.45rem;'>Safety Margin: {margin}</div>",
                    unsafe_allow_html=True)
        for w in safety.get("warnings", []):
            st.markdown(f"<div style='font-size:.8rem;color:var(--muted);margin-bottom:3px;'>"
                        f"- {w}</div>", unsafe_allow_html=True)
        tw = safety.get("therapeutic_window","")
        if tw:
            st.markdown(f"<div style='font-size:.8rem;color:var(--muted);margin-top:.5rem;'>"
                        f"Therapeutic window: {tw}</div>", unsafe_allow_html=True)


# ── Tab: Docking ──────────────────────────────────────────────────────────────

def _docking_tab(name: str, smiles: str, cid, light: bool, state_key: str):
    _sec("Molecular Docking")
    if not DOCKING_AVAILABLE:
        st.info(
            "Molecular docking requires an NVIDIA BioNeMo API key (DiffDock) "
            "or AutoDock Vina installation. "
            "Set `NVIDIA_API_KEY` environment variable to enable."
        )
        return
    if not smiles:
        st.warning(f"No SMILES structure available for {name}.")
        return

    tgts = get_compound_targets(int(cid)) if (DB_AVAILABLE and cid) else []
    tgt_names = [t["name"] for t in tgts if t.get("name")]

    col1, col2 = st.columns(2)
    with col1:
        if tgt_names:
            sel_tgt = st.selectbox("Target protein", tgt_names, key=f"dock_tgt_{state_key}")
        else:
            sel_tgt = st.text_input("Target protein name",
                placeholder="e.g. BACE1, AChE, GSK3B", key=f"dock_tgt_{state_key}")
    with col2:
        n_poses = st.slider("Poses to generate", 1, 20, 10, key=f"dock_poses_{state_key}")

    if st.button("Run Docking", key=f"run_dock_{state_key}") and sel_tgt:
        with st.spinner(f"Running docking: {name} → {sel_tgt}..."):
            try:
                result = _docking_svc.perform_docking(
                    drug_name=name, target_name=sel_tgt, ligand_smiles=smiles)
                st.session_state[state_key] = result
            except Exception as e:
                st.error(f"Docking failed: {e}")

    dr = st.session_state.get(state_key)
    if not dr: return
    if not dr.get("success"):
        st.warning(f"Docking result: {dr.get('error','No result returned')}")
        return

    affs  = dr.get("binding_affinities", [])
    poses = dr.get("poses", [])
    confs = dr.get("confidence_scores", [])

    if affs:
        c1,c2,c3 = st.columns(3)
        with c1: st.metric("Best Affinity",  f"{min(affs):.2f} kcal/mol")
        with c2: st.metric("Mean Affinity",  f"{sum(affs)/len(affs):.2f} kcal/mol")
        with c3: st.metric("Poses",          len(poses))

        fig = go.Figure(go.Bar(
            x=list(range(1, len(affs)+1)), y=affs,
            marker_color=["#38bdf8" if a==min(affs) else "#a78bfa" for a in affs],
            name="Affinity",
        ))
        fig.update_layout(title="Pose binding affinities (kcal/mol)",
                          xaxis_title="Pose", yaxis_title="Affinity (kcal/mol)")
        fig = _plot_style(fig, light)
        st.plotly_chart(fig, use_container_width=True)

        if confs:
            df_d = pd.DataFrame({
                "Pose":             list(range(1, len(confs)+1)),
                "Affinity (kcal/mol)": affs[:len(confs)],
                "Confidence":       confs,
            })
            st.dataframe(df_d, use_container_width=True)


# ── Tab: Quantum Optimization ─────────────────────────────────────────────────

def _quantum_tab(name: str, smiles: str, light: bool):
    _sec("Molecular Optimization")
    if not smiles:
        st.info("No SMILES structure available.")
        return

    with st.spinner("Running RDKit force-field optimization and electronic property calculation..."):
        try:
            results = run_quantum_optimization(name, smiles)
        except Exception as e:
            st.error(f"Optimization failed: {e}")
            return

    opt   = results.get("optimization", {})
    props = results.get("quantum_properties", {})

    col_o, col_e = st.columns(2)
    with col_o:
        _sec("3D Conformer Optimization")
        if opt.get("optimized"):
            ff = opt.get("ff_type","UFF")
            st.markdown(f"<div style='font-size:.84rem;color:var(--emerald);"
                        f"margin-bottom:.5rem;'>Optimized with {ff} force field</div>",
                        unsafe_allow_html=True)
            _kv("Energy before optimization", f"{opt.get('energy_before_kcal','—')} kcal/mol")
            _kv("Energy after optimization",  f"{opt.get('energy_after_kcal','—')} kcal/mol")
            _kv("Energy reduction (delta)",   f"{opt.get('energy_delta_kcal','—')} kcal/mol")
            _kv("Conformers generated",       str(opt.get("conformers","—")))
            _kv("Total atoms (with H)",       str(opt.get("num_atoms","—")))
        else:
            st.warning(f"Optimization limited: {opt.get('reason','unknown')}")

    with col_e:
        _sec("Electronic Properties (RDKit Estimates)")
        if props:
            _kv("HOMO-LUMO gap (proxy)",  f"{props.get('homo_lumo_gap_eV',0):.2f} eV")
            _kv("Dipole moment (proxy)",  f"{props.get('dipole_moment_proxy_D',0):.2f} D")
            _kv("Polarisability",         f"{props.get('polarisability_A3',0):.1f} Å³")
            _kv("Max partial charge",     f"{props.get('max_partial_charge',0):.4f}")
            _kv("Min partial charge",     f"{props.get('min_partial_charge',0):.4f}")
            _kv("QED (drug-likeness)",    f"{props.get('qed',0):.4f}")
        else:
            st.info("Electronic properties unavailable.")

    if props:
        with st.expander("Full property data"):
            st.json(props)


# ── Tab: 3D Structure ─────────────────────────────────────────────────────────

def _3d_tab(name: str, smiles: str, cid):
    _sec("3D Molecular Structure")
    if not smiles:
        st.warning(f"No SMILES available for {name}.")
        return
    if not (_viz and VIZ_3D_AVAILABLE):
        st.info("3D viewer requires py3Dmol and stmol.")
        st.code(smiles, language=None)
        return

    tgts     = get_compound_targets(int(cid)) if (DB_AVAILABLE and cid) else []
    tgt_name = tgts[0]["name"] if tgts else ""
    pdb      = None
    if tgt_name:
        try:
            from real_pdb_fetcher import RealPDBFetcher
            pdb = RealPDBFetcher().fetch_pdb(tgt_name)
        except Exception: pass

    _viz.render_visualization(name, tgt_name or "", smiles, pdb)


# ─── Knowledge Graph page ─────────────────────────────────────────────────────

def page_graph():
    light = st.session_state.get("light_mode", False)
    st.markdown("## Knowledge Graph")
    st.caption("Explore drug-gene-disease connections from the biomedical knowledge graph.")

    col1, col2 = st.columns([5, 1])
    with col1:
        center = st.text_input("graph_q",
            value=st.session_state.disease_query or "",
            placeholder="Enter disease or drug name...",
            label_visibility="collapsed")
    with col2:
        build = st.button("Build", use_container_width=True)

    if not center:
        st.info("Enter a disease or drug name to build the graph.")
        return
    if not DB_AVAILABLE:
        st.warning("Database unavailable.")
        return

    with st.spinner("Building knowledge graph..."):
        resolved = resolve_disease(center) if center else []
        mesh_ids = [r["mesh_id"] for r in resolved if r.get("mesh_id")]
        expanded = expand_mesh_ids(mesh_ids) if mesh_ids else []
        compounds = (get_compounds_for_disease(expanded or mesh_ids, limit=20)
                     if (expanded or mesh_ids) else db_search(center, limit=20))
        compounds = validate_and_deduplicate(compounds, require_smiles=False)[:20]

    if not compounds:
        st.warning("No compounds found.")
        return

    G = nx.Graph()
    cl = resolved[0]["heading"] if resolved else center
    G.add_node(cl, kind="Disease" if resolved else "Query")

    for c in compounds:
        n = c.get("name") or ""
        if not n: continue
        G.add_node(n, kind="Compound")
        G.add_edge(cl, n)
        for tgt in (c.get("targets") or "").split(";")[:2]:
            tgt = tgt.strip()
            if tgt and len(tgt) > 2:
                G.add_node(tgt, kind="Gene")
                G.add_edge(n, tgt)

    pos = nx.spring_layout(G, seed=42, k=1.4)
    kc  = {"Disease":"#fb7185","Query":"#38bdf8","Compound":"#a78bfa","Gene":"#34d399"}
    ns_map = {"Disease":24,"Query":22,"Compound":13,"Gene":9}

    nx_, ny_, nt_, nc_, ns_ = [], [], [], [], []
    for node, attr in G.nodes(data=True):
        x, y = pos[node]
        nx_.append(x); ny_.append(y); nt_.append(node[:28])
        nc_.append(kc.get(attr.get("kind",""), "#7c8daa"))
        ns_.append(ns_map.get(attr.get("kind",""), 11))

    ex_, ey_ = [], []
    for u, v in G.edges():
        x0,y0=pos[u]; x1,y1=pos[v]
        ex_+=[x0,x1,None]; ey_+=[y0,y1,None]

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=ex_, y=ey_, mode="lines",
        line=dict(color="#1a2e4a" if not light else "#d0daea", width=1),
        hoverinfo="none"))
    fig.add_trace(go.Scatter(x=nx_, y=ny_, mode="markers+text",
        marker=dict(size=ns_, color=nc_, line=dict(width=1, color="#04080f" if not light else "#f4f7fb")),
        text=nt_, textposition="top center",
        textfont=dict(size=8, color="#7c8daa"),
        hoverinfo="text"))
    fig.update_layout(
        showlegend=False, hovermode="closest", height=600,
        margin=dict(l=10,r=10,t=10,b=10),
    )
    fig = _plot_style(fig, light)
    fig.update_layout(xaxis=dict(showgrid=False,zeroline=False,showticklabels=False),
                      yaxis=dict(showgrid=False,zeroline=False,showticklabels=False))
    st.plotly_chart(fig, use_container_width=True)

    leg = st.columns(4)
    for col, (kind, color) in zip(leg, kc.items()):
        col.markdown(f"<span style='color:{color};font-size:.78rem;'>&#9679; {kind}</span>",
                     unsafe_allow_html=True)
    st.caption(f"{G.number_of_nodes()} nodes · {G.number_of_edges()} edges")


# ─── Data Explorer ────────────────────────────────────────────────────────────

def page_database():
    light = st.session_state.get("light_mode", False)
    st.markdown("## Data Explorer")
    if not DB_AVAILABLE:
        st.warning("Database unavailable."); return

    tabs = st.tabs(["Compounds","Targets","Diseases","Statistics"])

    with tabs[0]:
        col1, col2 = st.columns([4,1])
        with col1: q = st.text_input("cq", placeholder="Name, mechanism, target...", label_visibility="collapsed")
        with col2: lim = st.number_input("Rows", 10, 500, 50, key="clim")
        compounds = db_search(q, limit=int(lim)) if q else []
        if not compounds and not q:
            try:
                from database.schema import get_connection
                with get_connection() as conn:
                    cur = conn.cursor()
                    cur.execute("SELECT c.chembl_id,c.name,c.smiles,c.max_phase,cp.mw,cp.alogp "
                                "FROM compounds c LEFT JOIN compound_properties cp "
                                "ON cp.compound_id=c.id ORDER BY c.max_phase DESC NULLS LAST LIMIT %s",(lim,))
                    compounds = [{"chembl_id":r[0],"name":r[1],"smiles":r[2],
                                  "max_phase":r[3],"mw":r[4],"alogp":r[5]} for r in cur.fetchall()]
            except Exception as e: st.error(f"Error: {e}")
        if compounds:
            df = pd.DataFrame(compounds)
            keep = [c for c in ["chembl_id","name","max_phase","mw","alogp","psa","mechanisms"] if c in df.columns]
            st.dataframe(df[keep], use_container_width=True, height=430)
            st.download_button("Download CSV", data=df.to_csv(index=False),
                               file_name="compounds.csv", mime="text/csv")
        else: st.info("No results.")

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
            else: st.info("No targets.")
        except Exception as e: st.error(f"Error: {e}")

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
            else: st.info("MeSH table empty — run: `python database/mesh_importer.py`")
        except Exception as e: st.error(f"Error: {e}")

    with tabs[3]:
        try:
            stats = get_stats()
            c1,c2,c3 = st.columns(3)
            for i,(k,v) in enumerate(stats.items()):
                with [c1,c2,c3][i%3]: st.metric(k.replace("_"," ").title(), f"{v:,}")
            from database.schema import get_connection
            with get_connection() as conn:
                cur = conn.cursor()
                cur.execute("SELECT kind,COUNT(*) n FROM hetionet_nodes GROUP BY kind ORDER BY n DESC")
                kr = cur.fetchall()
            if kr:
                st.markdown("---")
                _sec("Graph Node Types")
                df_k = pd.DataFrame(kr, columns=["Type","Count"])
                fig = px.bar(df_k, x="Type", y="Count", color="Count",
                             color_continuous_scale=[[0,"#a78bfa"],[1,"#38bdf8"]])
                fig = _plot_style(fig, light)
                st.plotly_chart(fig, use_container_width=True)
        except Exception as e: st.error(f"Stats error: {e}")


# ─── Router ───────────────────────────────────────────────────────────────────

def main():
    _init()
    apply_theme()
    _sidebar()
    page = st.session_state.page
    if   page == "dashboard": page_dashboard()
    elif page == "discover":  page_discover()
    elif page == "analysis":  page_analysis()
    elif page == "graph":     page_graph()
    elif page == "database":  page_database()
    else:                     page_dashboard()


if __name__ == "__main__":
    main()
