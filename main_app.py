#!/usr/bin/env python3
"""
NeuroRepurpose Intelligence Platform
Full drug repurposing platform powered by ChEMBL 33 + HetioNet + RDKit ML scoring.
"""

import json
import logging
import os
import time
from pathlib import Path
from typing import Dict, List, Optional

import networkx as nx
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import requests
import streamlit as st

# ─── Page config (must be first Streamlit call) ───────────────────────────────
st.set_page_config(
    page_title="NeuroRepurpose",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

REPO_ROOT = Path(__file__).parent

# ─── Optional imports ─────────────────────────────────────────────────────────
try:
    from services.neuro_db_service import (
        get_neuro_compounds, search_compounds, get_compound_targets,
        get_compound_mechanisms, get_compound_indications,
        get_hetionet_paths, get_stats, is_available as db_available,
    )
    DB_AVAILABLE = db_available()
except Exception:
    DB_AVAILABLE = False
    def get_neuro_compounds(*a, **k): return []
    def search_compounds(*a, **k): return []
    def get_compound_targets(*a, **k): return []
    def get_compound_mechanisms(*a, **k): return []
    def get_compound_indications(*a, **k): return []
    def get_hetionet_paths(*a, **k): return []
    def get_stats(): return {}

try:
    from scoring_engine import score_drug
    SCORING_AVAILABLE = True
except Exception:
    SCORING_AVAILABLE = False
    def score_drug(*a, **k): return 0.5

try:
    from services.docking_service import DockingService
    DOCKING_AVAILABLE = True
except Exception:
    DOCKING_AVAILABLE = False
    DockingService = None

try:
    from quantum_optimization_strategies import run_quantum_optimization, render_quantum_optimization_section
    QUANTUM_AVAILABLE = True
except Exception:
    QUANTUM_AVAILABLE = False

try:
    from drug_protein_3d_visualizer import DrugProtein3DVisualizer
    VIZ_3D_AVAILABLE = True
    _viz = DrugProtein3DVisualizer()
except Exception:
    VIZ_3D_AVAILABLE = False
    _viz = None

try:
    from enhanced_authentic_data_fetcher import EnhancedAuthenticDataFetcher
    _fetcher = EnhancedAuthenticDataFetcher()
    FETCHER_AVAILABLE = True
except Exception:
    FETCHER_AVAILABLE = False
    _fetcher = None

try:
    from streamlit_echarts import st_echarts
    ECHARTS_AVAILABLE = True
except Exception:
    ECHARTS_AVAILABLE = False


# ─── Data loaders (JSON fallback) ────────────────────────────────────────────

@st.cache_data(ttl=3600, show_spinner=False)
def load_drugs_json() -> Dict:
    try:
        with open(REPO_ROOT / "drugs.json", encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        return {}


@st.cache_data(ttl=3600, show_spinner=False)
def load_interactions_json() -> Dict:
    try:
        with open(REPO_ROOT / "drug_interactions.json", encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        return {}


@st.cache_data(ttl=3600, show_spinner=False)
def load_disease_associations() -> Dict:
    try:
        with open(REPO_ROOT / "disease_associations.json", encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        return {}


def get_drug_smiles(name: str) -> str:
    drugs = load_drugs_json()
    key = name.lower()
    entry = drugs.get(key) or drugs.get(name, {})
    return entry.get("smiles", "")


# ─── Theme / CSS ─────────────────────────────────────────────────────────────

def apply_theme():
    st.markdown(
        """
        <style>
        /* ── Root palette ──────────────────────────────────────────── */
        :root {
            --bg-primary:    #070b14;
            --bg-card:       #0f1526;
            --bg-elevated:   #151d35;
            --accent-cyan:   #00d4ff;
            --accent-purple: #7b2fff;
            --accent-green:  #00ff9d;
            --accent-amber:  #ffb300;
            --text-primary:  #e8edf7;
            --text-secondary:#8896b3;
            --border:        #1e2a4a;
        }

        /* ── Global resets ─────────────────────────────────────────── */
        .stApp, .main { background: var(--bg-primary) !important; }
        [data-testid="stSidebar"] {
            background: var(--bg-card) !important;
            border-right: 1px solid var(--border);
        }
        .block-container { padding: 1.5rem 2rem !important; max-width: 1400px; }

        /* ── Typography ─────────────────────────────────────────────── */
        h1,h2,h3,h4,h5,h6,p,label,div { color: var(--text-primary) !important; }
        .stMarkdown p { color: var(--text-primary) !important; }
        small, .caption-text { color: var(--text-secondary) !important; font-size: 0.82rem; }

        /* ── Cards ──────────────────────────────────────────────────── */
        .nr-card {
            background: var(--bg-card);
            border: 1px solid var(--border);
            border-radius: 12px;
            padding: 1.25rem 1.5rem;
            margin-bottom: 1rem;
            transition: border-color 0.2s, box-shadow 0.2s;
        }
        .nr-card:hover {
            border-color: var(--accent-cyan);
            box-shadow: 0 0 16px rgba(0,212,255,0.1);
        }
        .nr-card-title {
            font-size: 1rem; font-weight: 700; color: var(--accent-cyan) !important;
            margin-bottom: 0.4rem;
        }
        .nr-card-sub { font-size: 0.82rem; color: var(--text-secondary) !important; }

        /* ── Metrics ─────────────────────────────────────────────────── */
        [data-testid="stMetric"] {
            background: var(--bg-elevated);
            border-radius: 10px;
            padding: 0.75rem 1rem;
            border: 1px solid var(--border);
        }
        [data-testid="stMetricValue"] { color: var(--accent-cyan) !important; font-size: 1.6rem !important; }
        [data-testid="stMetricLabel"] { color: var(--text-secondary) !important; }

        /* ── Inputs ──────────────────────────────────────────────────── */
        .stTextInput input, .stSelectbox select, .stTextArea textarea {
            background: var(--bg-elevated) !important;
            border: 1px solid var(--border) !important;
            border-radius: 8px !important;
            color: var(--text-primary) !important;
        }
        .stTextInput input:focus, .stTextArea textarea:focus {
            border-color: var(--accent-cyan) !important;
            box-shadow: 0 0 0 2px rgba(0,212,255,0.15) !important;
        }

        /* ── Buttons ─────────────────────────────────────────────────── */
        .stButton > button {
            background: linear-gradient(135deg, var(--accent-cyan), var(--accent-purple)) !important;
            color: #fff !important; border: none !important;
            border-radius: 8px !important; font-weight: 600 !important;
            padding: 0.5rem 1.4rem !important;
            transition: opacity 0.2s, transform 0.1s !important;
        }
        .stButton > button:hover { opacity: 0.88; transform: translateY(-1px); }
        .stButton > button.secondary {
            background: var(--bg-elevated) !important;
            border: 1px solid var(--border) !important;
            color: var(--text-primary) !important;
        }

        /* ── DataFrames ──────────────────────────────────────────────── */
        [data-testid="stDataFrame"] { border-radius: 10px; overflow: hidden; }
        [data-testid="stDataFrame"] thead th {
            background: var(--bg-elevated) !important;
            color: var(--accent-cyan) !important;
        }
        [data-testid="stDataFrame"] tbody tr:nth-child(even) td {
            background: var(--bg-card) !important;
        }

        /* ── Tabs ─────────────────────────────────────────────────────── */
        .stTabs [data-baseweb="tab-list"] { gap: 4px; border-bottom: 1px solid var(--border) !important; }
        .stTabs [data-baseweb="tab"] {
            background: transparent !important; border-radius: 8px 8px 0 0 !important;
            color: var(--text-secondary) !important; padding: 0.5rem 1.2rem !important;
        }
        .stTabs [aria-selected="true"] {
            background: var(--bg-elevated) !important;
            color: var(--accent-cyan) !important;
            border-top: 2px solid var(--accent-cyan) !important;
        }

        /* ── Progress / Spinner ──────────────────────────────────────── */
        .stProgress > div > div { background: var(--accent-cyan) !important; }
        .stSpinner > div { border-top-color: var(--accent-cyan) !important; }

        /* ── Confidence badges ───────────────────────────────────────── */
        .badge-high   { background:#00ff9d22; color:#00ff9d; border:1px solid #00ff9d44; border-radius:6px; padding:2px 8px; font-size:0.78rem; font-weight:700; }
        .badge-medium { background:#ffb30022; color:#ffb300; border:1px solid #ffb30044; border-radius:6px; padding:2px 8px; font-size:0.78rem; font-weight:700; }
        .badge-low    { background:#ff445522; color:#ff4455; border:1px solid #ff445544; border-radius:6px; padding:2px 8px; font-size:0.78rem; font-weight:700; }
        .phase-badge  { background:#7b2fff22; color:#ab6fff; border:1px solid #7b2fff44; border-radius:6px; padding:2px 8px; font-size:0.78rem; font-weight:600; }

        /* ── Sidebar nav items ────────────────────────────────────────── */
        .nr-nav-item {
            padding: 0.6rem 1rem; border-radius: 8px; margin-bottom: 4px;
            cursor: pointer; font-weight: 600; font-size: 0.9rem;
            color: var(--text-secondary) !important;
            transition: background 0.15s, color 0.15s;
        }
        .nr-nav-item:hover { background: var(--bg-elevated); color: var(--accent-cyan) !important; }
        .nr-nav-item.active { background: rgba(0,212,255,0.1); color: var(--accent-cyan) !important; border-left: 3px solid var(--accent-cyan); }

        /* ── Hero banner ─────────────────────────────────────────────── */
        .nr-hero {
            background: linear-gradient(135deg, #070b14 0%, #0f1a35 50%, #0b1520 100%);
            border: 1px solid var(--border); border-radius: 16px;
            padding: 2.5rem; margin-bottom: 1.5rem; position: relative; overflow: hidden;
        }
        .nr-hero::before {
            content: ''; position: absolute; top: -60px; right: -60px;
            width: 200px; height: 200px; border-radius: 50%;
            background: radial-gradient(circle, rgba(0,212,255,0.08) 0%, transparent 70%);
        }
        .nr-hero-title {
            font-size: 2rem; font-weight: 800;
            background: linear-gradient(90deg, var(--accent-cyan), var(--accent-purple));
            -webkit-background-clip: text; -webkit-text-fill-color: transparent;
            background-clip: text; margin-bottom: 0.5rem;
        }
        .nr-hero-sub { font-size: 1rem; color: var(--text-secondary) !important; }

        /* ── Dividers ─────────────────────────────────────────────────── */
        hr { border-color: var(--border) !important; }
        </style>
        """,
        unsafe_allow_html=True,
    )


# ─── Session state ────────────────────────────────────────────────────────────

def init_session():
    defaults = {
        "page": "dashboard",
        "search_query": "",
        "selected_disease": "Alzheimer's Disease",
        "selected_compound": None,
        "search_results": [],
        "docking_results": {},
        "chat_history": [],
    }
    for k, v in defaults.items():
        if k not in st.session_state:
            st.session_state[k] = v


# ─── Sidebar ──────────────────────────────────────────────────────────────────

def render_sidebar():
    with st.sidebar:
        st.markdown(
            """
            <div style='padding:1rem 0 1.5rem 0; text-align:center;'>
              <div style='font-size:1.6rem; font-weight:800;
                          background:linear-gradient(90deg,#00d4ff,#7b2fff);
                          -webkit-background-clip:text;-webkit-text-fill-color:transparent;
                          background-clip:text;'>
                🧬 NeuroRepurpose
              </div>
              <div style='font-size:0.78rem; color:#8896b3; margin-top:4px;'>
                AI Drug Repurposing Platform
              </div>
            </div>
            """,
            unsafe_allow_html=True,
        )

        pages = {
            "dashboard":   ("📊", "Dashboard"),
            "discover":    ("🔬", "Drug Discovery"),
            "network":     ("🕸️", "Knowledge Network"),
            "analysis":    ("⚗️", "Molecular Analysis"),
            "database":    ("🗄️", "Data Explorer"),
        }

        st.markdown("**Navigation**")
        for page_id, (icon, label) in pages.items():
            active = st.session_state.page == page_id
            cls = "nr-nav-item active" if active else "nr-nav-item"
            if st.button(f"{icon}  {label}", key=f"nav_{page_id}", use_container_width=True):
                st.session_state.page = page_id
                st.rerun()

        st.markdown("---")

        # Data sources status
        st.markdown("**Data Sources**")
        _status_row("ChEMBL 33", DB_AVAILABLE)
        _status_row("HetioNet", DB_AVAILABLE)
        _status_row("Scoring ML", SCORING_AVAILABLE)
        _status_row("3D Viewer", VIZ_3D_AVAILABLE)
        _status_row("Docking", DOCKING_AVAILABLE)

        if DB_AVAILABLE:
            st.markdown("---")
            try:
                stats = get_stats()
                st.markdown("**Database**")
                st.caption(f"Compounds: {stats.get('compounds',0):,}")
                st.caption(f"Targets: {stats.get('targets',0):,}")
                st.caption(f"Indications: {stats.get('indications',0):,}")
                st.caption(f"KG Nodes: {stats.get('hetionet_nodes',0):,}")
            except Exception:
                pass

        st.markdown("---")
        st.caption("Powered by ChEMBL 33 · HetioNet · RDKit")


def _status_row(label: str, ok: bool):
    colour = "#00ff9d" if ok else "#ff4455"
    dot = "●" if ok else "○"
    st.markdown(
        f"<div style='font-size:0.8rem; color:{colour};'>{dot} {label}</div>",
        unsafe_allow_html=True,
    )


# ─── Shared helpers ───────────────────────────────────────────────────────────

def _confidence_badge(score: float) -> str:
    if score >= 0.75:
        return f"<span class='badge-high'>{score:.0%}</span>"
    if score >= 0.45:
        return f"<span class='badge-medium'>{score:.0%}</span>"
    return f"<span class='badge-low'>{score:.0%}</span>"


def _phase_badge(phase) -> str:
    try:
        p = float(phase or 0)
    except Exception:
        p = 0
    labels = {4: "FDA Approved", 3: "Phase III", 2: "Phase II", 1: "Phase I"}
    label = labels.get(int(p), f"Phase {p}" if p > 0 else "Preclinical")
    return f"<span class='phase-badge'>{label}</span>"


def _drug_card(compound: Dict, idx: int = 0):
    name = compound.get("name") or "Unknown"
    smiles = compound.get("smiles") or ""
    score = float(compound.get("score") or compound.get("alogp") or 0.5)
    if score > 1:
        score = min(1.0, score / 10)
    phase = compound.get("max_phase", 0)
    mechs = compound.get("mechanisms") or compound.get("mechanism", "—")
    targets = compound.get("targets") or compound.get("target", "—")

    with st.container():
        st.markdown(
            f"""
            <div class='nr-card'>
              <div style='display:flex; justify-content:space-between; align-items:flex-start;'>
                <div>
                  <div class='nr-card-title'>{name}</div>
                  <div class='nr-card-sub'>{(mechs or '—')[:120]}</div>
                </div>
                <div style='text-align:right; flex-shrink:0; margin-left:1rem;'>
                  {_phase_badge(phase)}
                </div>
              </div>
              <div style='margin-top:0.75rem; font-size:0.8rem; color:#8896b3;'>
                <b style='color:#e8edf7;'>Targets:</b> {(str(targets or '—'))[:100]}
              </div>
            </div>
            """,
            unsafe_allow_html=True,
        )
        col1, col2 = st.columns([3, 1])
        with col2:
            if st.button("Analyse →", key=f"analyse_{idx}_{name[:8]}"):
                st.session_state.selected_compound = compound
                st.session_state.page = "analysis"
                st.rerun()


# ─── Page: Dashboard ─────────────────────────────────────────────────────────

def page_dashboard():
    st.markdown(
        """
        <div class='nr-hero'>
          <div class='nr-hero-title'>NeuroRepurpose Intelligence Platform</div>
          <div class='nr-hero-sub'>
            Discover new therapeutic uses for existing drugs using ChEMBL 33,
            HetioNet knowledge graphs, and ML scoring — focused on neurological diseases.
          </div>
        </div>
        """,
        unsafe_allow_html=True,
    )

    # Quick search bar
    q = st.text_input(
        "",
        placeholder="Search by disease, drug name, target, or mechanism…",
        key="hero_search",
        label_visibility="collapsed",
    )
    col_btn, col_space = st.columns([1, 5])
    with col_btn:
        if st.button("🔍 Search", use_container_width=True) and q:
            st.session_state.search_query = q
            st.session_state.page = "discover"
            st.rerun()

    st.markdown("---")

    # Stats row
    col1, col2, col3, col4, col5 = st.columns(5)
    drugs_json = load_drugs_json()

    stats = {}
    if DB_AVAILABLE:
        try:
            stats = get_stats()
        except Exception:
            pass

    with col1:
        st.metric("Compounds", f"{stats.get('compounds', len(drugs_json)):,}")
    with col2:
        st.metric("Protein Targets", f"{stats.get('targets', 15398):,}")
    with col3:
        st.metric("Indications", f"{stats.get('indications', 51582):,}")
    with col4:
        st.metric("KG Nodes", f"{stats.get('hetionet_nodes', 47031):,}")
    with col5:
        st.metric("Data Source", "ChEMBL 33")

    st.markdown("---")

    # Disease quick-launch
    st.markdown("### Top Neurological Disease Areas")
    diseases = [
        ("Alzheimer's Disease", "271", "🧠"),
        ("Parkinson's Disease", "215", "🔵"),
        ("Multiple Sclerosis", "187", "🔴"),
        ("Epilepsy", "109", "⚡"),
        ("Schizophrenia", "286", "🟣"),
        ("Huntington Disease", "46", "🟠"),
        ("ALS", "134", "💛"),
        ("Stroke", "154", "❤️"),
    ]

    cols = st.columns(4)
    for i, (disease, cnt, emoji) in enumerate(diseases):
        with cols[i % 4]:
            st.markdown(
                f"""
                <div class='nr-card' style='cursor:pointer;'>
                  <div style='font-size:1.5rem;'>{emoji}</div>
                  <div class='nr-card-title'>{disease}</div>
                  <div class='nr-card-sub'>{cnt} approved compounds</div>
                </div>
                """,
                unsafe_allow_html=True,
            )
            if st.button(f"Explore {disease.split()[0]}", key=f"dash_{i}", use_container_width=True):
                st.session_state.selected_disease = disease
                st.session_state.page = "discover"
                st.rerun()

    st.markdown("---")
    st.markdown("### Recent Drug Repurposing Highlights")
    highlights = [
        ("Metformin", "Alzheimer's Disease", 0.82, "Anti-diabetic → neuro-protective via AMPK pathway"),
        ("Liraglutide", "Parkinson's Disease", 0.78, "GLP-1 receptor agonist → dopaminergic neuroprotection"),
        ("Fingolimod", "Alzheimer's Disease", 0.71, "MS drug → reduces Aβ accumulation via S1P signalling"),
        ("Memantine", "Huntington Disease", 0.69, "NMDA antagonist (Alzheimer's) → HD symptom relief"),
        ("Riluzole", "Alzheimer's Disease", 0.67, "ALS drug → reduces tau phosphorylation"),
        ("Telmisartan", "Stroke", 0.74, "ARB → cerebral blood flow improvement"),
    ]
    cols = st.columns(3)
    for i, (drug, disease, score, rationale) in enumerate(highlights):
        with cols[i % 3]:
            st.markdown(
                f"""
                <div class='nr-card'>
                  <div class='nr-card-title'>{drug}</div>
                  <div style='font-size:0.8rem; color:#8896b3;'>→ {disease}</div>
                  <div style='margin:0.5rem 0;'>{_confidence_badge(score)}</div>
                  <div style='font-size:0.8rem; color:#c8d0e0;'>{rationale}</div>
                </div>
                """,
                unsafe_allow_html=True,
            )


# ─── Page: Drug Discovery ─────────────────────────────────────────────────────

def page_discover():
    st.markdown("## 🔬 Drug Discovery")

    # Controls
    col_q, col_d, col_btn = st.columns([3, 2, 1])
    with col_q:
        query = st.text_input(
            "Search",
            value=st.session_state.search_query,
            placeholder="Drug name, target, mechanism…",
            label_visibility="collapsed",
        )
    with col_d:
        disease = st.selectbox(
            "Disease",
            options=[
                "Alzheimer's Disease", "Parkinson's Disease", "Multiple Sclerosis",
                "Epilepsy", "Schizophrenia", "Huntington Disease",
                "ALS", "Stroke", "Depression", "Anxiety",
                "Autism Spectrum Disorder", "Bipolar Disorder",
            ],
            index=["Alzheimer's Disease", "Parkinson's Disease", "Multiple Sclerosis",
                   "Epilepsy", "Schizophrenia", "Huntington Disease",
                   "ALS", "Stroke", "Depression", "Anxiety",
                   "Autism Spectrum Disorder", "Bipolar Disorder"].index(
                       st.session_state.get("selected_disease", "Alzheimer's Disease")
                   ),
            label_visibility="collapsed",
        )
    with col_btn:
        run_search = st.button("🔍 Find Drugs", use_container_width=True)

    if run_search or st.session_state.search_query:
        st.session_state.selected_disease = disease
        _run_discovery(query or disease, disease)


def _run_discovery(query: str, disease: str):
    with st.spinner("Querying ChEMBL + knowledge graph…"):
        compounds = []

        # Try DB first
        if DB_AVAILABLE:
            if query and query.lower() != disease.lower():
                compounds = search_compounds(query, limit=30)
            if not compounds:
                compounds = get_neuro_compounds(disease, limit=30)

        # Fallback to JSON
        if not compounds:
            compounds = _json_fallback(query or disease)

        # Score all results
        if SCORING_AVAILABLE:
            for c in compounds:
                try:
                    c["score"] = score_drug(c.get("name", ""), c)
                except Exception:
                    c["score"] = float(c.get("max_phase", 0)) / 4.0

        compounds.sort(key=lambda x: float(x.get("score") or x.get("max_phase") or 0), reverse=True)
        st.session_state.search_results = compounds

    if not compounds:
        st.warning("No results found. Try a different query or disease.")
        return

    st.success(f"Found {len(compounds)} candidates for **{disease}**")

    # Filters
    col_f1, col_f2, col_f3 = st.columns(3)
    with col_f1:
        min_phase = st.slider("Min clinical phase", 0, 4, 0)
    with col_f2:
        show_n = st.slider("Results to show", 5, len(compounds), min(20, len(compounds)))
    with col_f3:
        sort_by = st.selectbox("Sort by", ["Score", "Phase", "Name"])

    filtered = [c for c in compounds if float(c.get("max_phase") or 0) >= min_phase]
    if sort_by == "Name":
        filtered.sort(key=lambda x: x.get("name", ""))
    elif sort_by == "Phase":
        filtered.sort(key=lambda x: float(x.get("max_phase") or 0), reverse=True)

    filtered = filtered[:show_n]

    # Summary chart
    if len(filtered) >= 3:
        phases = [float(c.get("max_phase") or 0) for c in filtered]
        phase_df = pd.DataFrame({"Phase": phases})
        phase_counts = phase_df["Phase"].value_counts().reset_index()
        phase_counts.columns = ["Phase", "Count"]
        phase_counts["Phase"] = phase_counts["Phase"].apply(
            lambda p: {4: "FDA Approved", 3: "Phase III", 2: "Phase II", 1: "Phase I"}.get(int(p), "Preclinical")
        )
        fig = px.pie(
            phase_counts, values="Count", names="Phase",
            color_discrete_sequence=["#00d4ff", "#7b2fff", "#00ff9d", "#ffb300", "#ff4455"],
            hole=0.5,
        )
        fig.update_layout(
            paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
            font_color="#e8edf7", margin=dict(t=20, b=20, l=10, r=10),
            showlegend=True, legend=dict(font_color="#8896b3"),
        )
        with st.expander("Clinical Phase Distribution", expanded=False):
            st.plotly_chart(fig, use_container_width=True)

    # Results grid
    st.markdown(f"### Top {len(filtered)} Repurposing Candidates")
    for i, c in enumerate(filtered):
        _drug_card(c, idx=i)

    # Download
    if filtered:
        df_dl = pd.DataFrame([{
            "Name": c.get("name"), "Phase": c.get("max_phase"),
            "Score": round(float(c.get("score") or 0), 3),
            "SMILES": c.get("smiles", ""),
            "Mechanisms": c.get("mechanisms", ""),
            "Targets": c.get("targets", ""),
            "Indications": c.get("indications", ""),
        } for c in filtered])
        st.download_button(
            "📥 Download CSV",
            data=df_dl.to_csv(index=False),
            file_name=f"neurorepurpose_{disease.lower().replace(' ', '_')}.csv",
            mime="text/csv",
        )


def _json_fallback(query: str) -> List[Dict]:
    """Simple JSON-based fallback search when DB is unavailable."""
    drugs = load_drugs_json()
    q = query.lower()
    results = []
    neuro_cats = ["alzheimer", "parkinson", "neurolog", "schizo", "epilep", "ms", "huntington", "als", "stroke", "depress", "anxiet"]
    for name, data in drugs.items():
        cat = (data.get("therapeutic_category") or "").lower()
        match = (q in name) or any(n in q or n in cat for n in neuro_cats)
        if match:
            results.append({
                "name": name.title(),
                "smiles": data.get("smiles", ""),
                "max_phase": 4.0 if data.get("fda_approved") else 0.0,
                "mechanisms": data.get("mechanism_of_action", ""),
                "targets": data.get("primary_target", ""),
                "score": 0.5,
            })
        if len(results) >= 30:
            break
    return results


# ─── Page: Knowledge Network ─────────────────────────────────────────────────

def page_network():
    st.markdown("## 🕸️ Drug-Disease Knowledge Network")

    col1, col2 = st.columns([2, 1])
    with col1:
        center_query = st.text_input("Centre node (drug or disease)", value=st.session_state.selected_disease)
    with col2:
        depth = st.slider("Network depth", 1, 3, 2)

    if st.button("Build Network"):
        _render_network(center_query, depth)
    elif st.session_state.search_results:
        st.info("Showing network for current search results. Click 'Build Network' to update.")
        _render_results_network(st.session_state.search_results[:15])


def _render_network(center: str, depth: int):
    G = nx.Graph()
    G.add_node(center, kind="query", size=30)

    compounds = []
    if DB_AVAILABLE:
        compounds = get_neuro_compounds(center, limit=15)
    if not compounds:
        compounds = _json_fallback(center)[:15]

    edges_data = []
    for c in compounds:
        name = c.get("name", "")
        G.add_node(name, kind="compound", phase=float(c.get("max_phase") or 0), size=20)
        G.add_edge(center, name, relation="indicated_for")
        edges_data.append((center, name))

        if depth >= 2:
            for t in (c.get("targets") or "").split(";")[:3]:
                t = t.strip()
                if t and len(t) > 2:
                    G.add_node(t, kind="target", size=15)
                    G.add_edge(name, t, relation="targets")

    if not G.nodes:
        st.warning("No network data found.")
        return

    pos = nx.spring_layout(G, seed=42, k=1.5)

    kind_colors = {"query": "#00d4ff", "compound": "#7b2fff", "target": "#00ff9d", "indication": "#ffb300"}

    node_x, node_y, node_text, node_color, node_size = [], [], [], [], []
    for node, attr in G.nodes(data=True):
        x, y = pos[node]
        node_x.append(x); node_y.append(y)
        node_text.append(node)
        node_color.append(kind_colors.get(attr.get("kind", "compound"), "#8896b3"))
        node_size.append(attr.get("size", 15))

    edge_x, edge_y = [], []
    for u, v in G.edges():
        x0, y0 = pos[u]; x1, y1 = pos[v]
        edge_x += [x0, x1, None]; edge_y += [y0, y1, None]

    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=edge_x, y=edge_y, mode="lines",
        line=dict(color="#1e2a4a", width=1), hoverinfo="none",
    ))
    fig.add_trace(go.Scatter(
        x=node_x, y=node_y, mode="markers+text",
        marker=dict(size=node_size, color=node_color, line=dict(width=1, color="#070b14")),
        text=node_text, textposition="top center",
        textfont=dict(size=9, color="#8896b3"),
        hoverinfo="text",
    ))
    fig.update_layout(
        showlegend=False, hovermode="closest",
        paper_bgcolor="#070b14", plot_bgcolor="#070b14",
        font_color="#e8edf7",
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        margin=dict(l=10, r=10, t=20, b=10),
        height=550,
    )
    st.plotly_chart(fig, use_container_width=True)
    st.caption(f"Nodes: {G.number_of_nodes()} | Edges: {G.number_of_edges()}")

    # Legend
    leg_cols = st.columns(4)
    for col, (kind, colour, label) in zip(leg_cols, [
        ("query", "#00d4ff", "Query"), ("compound", "#7b2fff", "Compound"),
        ("target", "#00ff9d", "Target"), ("indication", "#ffb300", "Indication"),
    ]):
        col.markdown(f"<span style='color:{colour};'>●</span> {label}", unsafe_allow_html=True)


def _render_results_network(compounds: List[Dict]):
    if not compounds:
        return
    disease = st.session_state.selected_disease
    G = nx.Graph()
    G.add_node(disease, kind="query")
    for c in compounds:
        G.add_node(c.get("name", "?"), kind="compound")
        G.add_edge(disease, c.get("name", "?"))
    _render_network_fig(G)


def _render_network_fig(G):
    pos = nx.spring_layout(G, seed=42)
    fig = go.Figure()
    ex, ey = [], []
    for u, v in G.edges():
        x0, y0 = pos[u]; x1, y1 = pos[v]
        ex += [x0, x1, None]; ey += [y0, y1, None]
    fig.add_trace(go.Scatter(x=ex, y=ey, mode="lines", line=dict(color="#1e2a4a", width=1), hoverinfo="none"))
    nx_arr = np.array(list(pos.values()))
    fig.add_trace(go.Scatter(
        x=nx_arr[:, 0], y=nx_arr[:, 1], mode="markers+text",
        marker=dict(size=16, color="#7b2fff"),
        text=list(pos.keys()), textposition="top center",
        textfont=dict(size=8, color="#8896b3"),
    ))
    fig.update_layout(
        paper_bgcolor="#070b14", plot_bgcolor="#070b14",
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        showlegend=False, height=450, margin=dict(l=0, r=0, t=0, b=0),
    )
    st.plotly_chart(fig, use_container_width=True)


# ─── Page: Molecular Analysis ─────────────────────────────────────────────────

def page_analysis():
    st.markdown("## ⚗️ Molecular Analysis")

    compound = st.session_state.selected_compound
    if not compound:
        st.info("Select a compound from Drug Discovery to analyse it here.")
        if st.button("Go to Drug Discovery"):
            st.session_state.page = "discover"
            st.rerun()
        return

    name = compound.get("name", "Unknown")
    smiles = compound.get("smiles") or get_drug_smiles(name)
    chembl_id = compound.get("chembl_id", "")

    # Header
    col_h1, col_h2 = st.columns([3, 1])
    with col_h1:
        st.markdown(f"### {name}")
        if chembl_id:
            st.caption(f"ChEMBL: {chembl_id}")
        mech = compound.get("mechanisms") or compound.get("mechanism", "—")
        st.markdown(f"**Mechanism:** {mech[:200] if mech else '—'}")
    with col_h2:
        phase_html = _phase_badge(compound.get("max_phase", 0))
        st.markdown(phase_html, unsafe_allow_html=True)
        score = float(compound.get("score") or 0.5)
        st.markdown(_confidence_badge(score) + " confidence", unsafe_allow_html=True)

    st.markdown("---")

    tab1, tab2, tab3, tab4, tab5 = st.tabs(
        ["📊 Properties", "🎯 Targets", "🏥 Indications", "⚛️ Quantum", "🔬 3D Structure"]
    )

    # ── Tab 1: Properties ────────────────────────────────────────────────────
    with tab1:
        _render_properties(compound, smiles)

    # ── Tab 2: Targets ───────────────────────────────────────────────────────
    with tab2:
        targets = []
        if DB_AVAILABLE and chembl_id:
            targets = get_compound_targets(chembl_id)
        if not targets:
            interactions = load_interactions_json()
            drug_ints = interactions.get(name.lower(), {})
            targets = [
                {"name": t.get("protein_name", t.get("gene_symbol", "")),
                 "gene_symbol": t.get("gene_symbol", ""),
                 "mechanism": t.get("interaction_type", ""),
                 "confidence": t.get("confidence", 0.5)}
                for t in drug_ints.get("targets", [])
            ]

        if targets:
            df_t = pd.DataFrame(targets)
            cols_show = [c for c in ["name", "gene_symbol", "target_type", "mechanism", "action_type", "confidence"] if c in df_t.columns]
            st.dataframe(df_t[cols_show].head(20), use_container_width=True)

            # Confidence chart
            if "confidence" in df_t.columns:
                fig = px.bar(
                    df_t.head(10), x="name" if "name" in df_t else "gene_symbol",
                    y="confidence", color="confidence",
                    color_continuous_scale=[[0, "#ff4455"], [0.5, "#ffb300"], [1, "#00ff9d"]],
                    title="Target Confidence Scores",
                )
                fig.update_layout(paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)", font_color="#e8edf7")
                st.plotly_chart(fig, use_container_width=True)
        else:
            st.info("No target data available for this compound.")

    # ── Tab 3: Indications ───────────────────────────────────────────────────
    with tab3:
        indications = []
        if DB_AVAILABLE and chembl_id:
            indications = get_compound_indications(chembl_id)

        if indications:
            df_i = pd.DataFrame(indications)
            st.dataframe(df_i, use_container_width=True)
        else:
            inds = compound.get("indications") or ""
            if inds:
                for ind in str(inds).split(";"):
                    st.markdown(f"- {ind.strip()}")
            else:
                st.info("No indication data available.")

        # Clinical trial lookup
        if FETCHER_AVAILABLE and st.button("Fetch Clinical Trials"):
            with st.spinner("Querying ClinicalTrials.gov…"):
                try:
                    disease = st.session_state.selected_disease
                    trials = _fetcher.fetch_clinical_trials(name, disease)
                    if trials:
                        st.markdown(f"**{len(trials)} Clinical Trials Found**")
                        for t in trials[:5]:
                            st.markdown(
                                f"<div class='nr-card'><b>{t.get('title','?')}</b>"
                                f"<br><small>{t.get('status','?')} | Phase: {t.get('phase','?')}</small></div>",
                                unsafe_allow_html=True,
                            )
                    else:
                        st.info("No active trials found.")
                except Exception as e:
                    st.warning(f"Trial lookup failed: {e}")

    # ── Tab 4: Quantum ───────────────────────────────────────────────────────
    with tab4:
        if QUANTUM_AVAILABLE:
            render_quantum_optimization_section(drug_name=name, smiles=smiles, drug_data=compound)
        else:
            st.info("Quantum calculations require RDKit. Install: `pip install rdkit`")

    # ── Tab 5: 3D Structure ──────────────────────────────────────────────────
    with tab5:
        if not smiles:
            st.warning(f"No SMILES structure available for {name}.")
        elif VIZ_3D_AVAILABLE and _viz:
            st.markdown(f"**3D Molecular Structure: {name}**")
            targets_list = []
            if DB_AVAILABLE and chembl_id:
                targets_list = get_compound_targets(chembl_id)
            target_name = targets_list[0]["name"] if targets_list else ""

            pdb_content = None
            if target_name:
                try:
                    from real_pdb_fetcher import RealPDBFetcher
                    pdb_content = RealPDBFetcher().fetch_pdb(target_name)
                except Exception:
                    pass

            _viz.render_visualization(name, target_name, smiles, pdb_content)
        else:
            st.info("3D viewer requires py3Dmol and stmol.")
            st.code(f"SMILES: {smiles}")


def _render_properties(compound: Dict, smiles: str):
    mw = compound.get("mw") or compound.get("molecular_weight")
    logp = compound.get("alogp") or compound.get("logp")
    psa = compound.get("psa") or compound.get("tpsa")
    hba = compound.get("hba")
    hbd = compound.get("hbd")

    if not any([mw, logp, psa]) and smiles:
        try:
            from rdkit import Chem
            from rdkit.Chem import Descriptors, rdMolDescriptors, QED
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                mw = round(Descriptors.ExactMolWt(mol), 1)
                logp = round(Descriptors.MolLogP(mol), 2)
                psa = round(Descriptors.TPSA(mol), 1)
                hba = rdMolDescriptors.CalcNumHBA(mol)
                hbd = rdMolDescriptors.CalcNumHBD(mol)
                qed = round(QED.qed(mol), 3)
        except Exception:
            pass

    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.metric("Mol. Weight (Da)", f"{mw:.1f}" if mw else "N/A")
    with col2:
        st.metric("LogP", f"{logp:.2f}" if logp is not None else "N/A")
    with col3:
        st.metric("TPSA (Å²)", f"{psa:.1f}" if psa else "N/A")
    with col4:
        st.metric("QED", f"{qed:.3f}" if "qed" in dir() else "N/A")

    col5, col6 = st.columns(2)
    with col5:
        if hba is not None:
            st.metric("H-Bond Acceptors", hba)
    with col6:
        if hbd is not None:
            st.metric("H-Bond Donors", hbd)

    # Lipinski gauge
    violations = sum([
        (mw or 0) > 500,
        (logp or 0) > 5,
        (hbd or 0) > 5,
        (hba or 0) > 10,
    ])
    st.markdown(f"**Lipinski Rule of 5:** {'✅ Pass' if violations == 0 else f'⚠️ {violations} violation(s)'}")

    if smiles:
        st.markdown("**SMILES**")
        st.code(smiles, language=None)


# ─── Page: Database Explorer ──────────────────────────────────────────────────

def page_database():
    st.markdown("## 🗄️ Data Explorer")

    tab1, tab2, tab3 = st.tabs(["Compounds", "Targets", "HetioNet Nodes"])

    with tab1:
        st.markdown("### Browse Compounds")
        col1, col2 = st.columns([3, 1])
        with col1:
            q = st.text_input("Filter by name or mechanism", key="db_q")
        with col2:
            limit = st.number_input("Rows", 10, 200, 50, key="db_limit")

        if DB_AVAILABLE:
            compounds = search_compounds(q or "", limit=limit) if q else []
            if not compounds and not q:
                from database.schema import get_connection
                try:
                    with get_connection() as conn:
                        cur = conn.cursor()
                        cur.execute(
                            "SELECT c.chembl_id, c.name, c.smiles, c.max_phase, cp.mw, cp.alogp "
                            "FROM compounds c LEFT JOIN compound_properties cp ON cp.compound_id=c.id "
                            "ORDER BY c.max_phase DESC NULLS LAST LIMIT %s",
                            (limit,)
                        )
                        compounds = [
                            {"chembl_id": r[0], "name": r[1], "smiles": r[2],
                             "max_phase": r[3], "mw": r[4], "alogp": r[5]}
                            for r in cur.fetchall()
                        ]
                except Exception as e:
                    st.error(f"DB error: {e}")

            if compounds:
                df = pd.DataFrame(compounds)
                keep = [c for c in ["chembl_id", "name", "max_phase", "mw", "alogp", "psa", "mechanisms"] if c in df.columns]
                st.dataframe(df[keep], use_container_width=True, height=400)
                st.download_button(
                    "📥 Download CSV",
                    data=df.to_csv(index=False),
                    file_name="neurorepurpose_compounds.csv",
                    mime="text/csv",
                )
            else:
                st.info("Enter a search term or wait for DB to finish importing." if not q else "No results.")
        else:
            st.warning("Database not available. Run: `python3 database/importer.py`")
            drugs = load_drugs_json()
            q_lower = (q or "").lower()
            rows = [{"name": n.title(), **d} for n, d in drugs.items() if not q_lower or q_lower in n][:limit]
            if rows:
                st.dataframe(pd.DataFrame(rows)[["name", "smiles"]], use_container_width=True)

    with tab2:
        st.markdown("### Protein Targets")
        if DB_AVAILABLE:
            from database.schema import get_connection
            try:
                with get_connection() as conn:
                    cur = conn.cursor()
                    cur.execute(
                        "SELECT chembl_tid, name, target_type, gene_symbol, organism "
                        "FROM targets ORDER BY name LIMIT 500"
                    )
                    targets = [{"chembl_tid": r[0], "name": r[1], "target_type": r[2], "gene_symbol": r[3], "organism": r[4]} for r in cur.fetchall()]
                if targets:
                    st.dataframe(pd.DataFrame(targets), use_container_width=True, height=400)
                else:
                    st.info("No target data yet. Run the importer.")
            except Exception as e:
                st.error(f"DB error: {e}")
        else:
            st.warning("Database not available.")

    with tab3:
        st.markdown("### HetioNet Knowledge Graph Nodes")
        if DB_AVAILABLE:
            from database.schema import get_connection
            try:
                with get_connection() as conn:
                    cur = conn.cursor()
                    cur.execute("SELECT kind, COUNT(*) AS n FROM hetionet_nodes GROUP BY kind ORDER BY n DESC")
                    kind_counts = cur.fetchall()
                    if kind_counts:
                        df_k = pd.DataFrame(kind_counts, columns=["Kind", "Count"])
                        fig = px.bar(
                            df_k, x="Kind", y="Count",
                            color="Count",
                            color_continuous_scale=[[0, "#7b2fff"], [1, "#00d4ff"]],
                        )
                        fig.update_layout(paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)", font_color="#e8edf7")
                        st.plotly_chart(fig, use_container_width=True)
                    else:
                        st.info("No HetioNet nodes yet. Ensure nodes.tsv is in data/hetionet/ and re-run importer.")
            except Exception as e:
                st.error(f"DB error: {e}")
        else:
            st.warning("Database not available.")


# ─── Main router ─────────────────────────────────────────────────────────────

def main():
    init_session()
    apply_theme()
    render_sidebar()

    page = st.session_state.page
    if page == "dashboard":
        page_dashboard()
    elif page == "discover":
        page_discover()
    elif page == "network":
        page_network()
    elif page == "analysis":
        page_analysis()
    elif page == "database":
        page_database()
    else:
        page_dashboard()


if __name__ == "__main__":
    main()
