#!/usr/bin/env python3
"""
NeuroRepurpose Intelligence Platform
Drug repurposing powered by ChEMBL 33 + HetioNet + MeSH ontology + RDKit scoring.
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
import streamlit as st

st.set_page_config(
    page_title="NeuroRepurpose",
    page_icon=None,
    layout="wide",
    initial_sidebar_state="expanded",
)

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

REPO_ROOT = Path(__file__).parent

# ─── Service imports ──────────────────────────────────────────────────────────

try:
    from services.neuro_db_service import (
        get_compounds_for_disease,
        search_compounds as db_search_compounds,
        get_compound_targets,
        get_compound_activities,
        get_compound_indications,
        get_hetionet_paths,
        get_stats,
        is_available as db_is_available,
        get_disease_top_compounds,
        get_available_diseases,
        get_compound_by_chembl,
    )
    DB_AVAILABLE = db_is_available()
except Exception as _e:
    logger.warning(f"neuro_db_service unavailable: {_e}")
    DB_AVAILABLE = False
    def get_compounds_for_disease(*a, **k): return []
    def db_search_compounds(*a, **k): return []
    def get_compound_targets(*a, **k): return []
    def get_compound_activities(*a, **k): return []
    def get_compound_indications(*a, **k): return []
    def get_hetionet_paths(*a, **k): return []
    def get_stats(): return {}
    def get_disease_top_compounds(*a, **k): return []
    def get_available_diseases(*a, **k): return []
    def get_compound_by_chembl(*a, **k): return None

try:
    from services.disease_resolver import (
        resolve_disease,
        expand_mesh_ids,
        suggest_diseases,
        mesh_available,
    )
    MESH_AVAILABLE = mesh_available() if DB_AVAILABLE else False
except Exception as _e:
    logger.warning(f"disease_resolver unavailable: {_e}")
    MESH_AVAILABLE = False
    def resolve_disease(*a, **k): return []
    def expand_mesh_ids(*a, **k): return []
    def suggest_diseases(*a, **k): return []

try:
    from services.repurposing_scorer import score_compound_list
    SCORER_AVAILABLE = True
except Exception as _e:
    logger.warning(f"repurposing_scorer unavailable: {_e}")
    SCORER_AVAILABLE = False
    def score_compound_list(compounds, mesh_ids):
        for c in compounds:
            c["score"] = float(c.get("max_phase") or 0) / 4.0
            c["score_breakdown"] = {}
        compounds.sort(key=lambda x: x["score"], reverse=True)
        return compounds

try:
    from services.compound_validator import validate_and_deduplicate
    VALIDATOR_AVAILABLE = True
except Exception as _e:
    logger.warning(f"compound_validator unavailable: {_e}")
    VALIDATOR_AVAILABLE = False
    def validate_and_deduplicate(compounds, **k): return compounds

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


# ─── CSS / theme ──────────────────────────────────────────────────────────────

def apply_theme():
    st.markdown(
        """
        <style>
        :root {
            --bg-primary:    #070b14;
            --bg-card:       #0f1526;
            --bg-elevated:   #151d35;
            --accent-cyan:   #00d4ff;
            --accent-purple: #7b2fff;
            --accent-green:  #00ff9d;
            --accent-amber:  #ffb300;
            --accent-red:    #ff4455;
            --text-primary:  #e8edf7;
            --text-secondary:#8896b3;
            --border:        #1e2a4a;
        }

        .stApp, .main                   { background: var(--bg-primary) !important; }
        [data-testid="stSidebar"]       { background: var(--bg-card) !important; border-right: 1px solid var(--border); }
        .block-container                { padding: 1.5rem 2rem !important; max-width: 1440px; }

        h1,h2,h3,h4,h5,h6              { color: var(--text-primary) !important; }
        p, label, div                  { color: var(--text-primary) !important; }
        .stMarkdown p                  { color: var(--text-primary) !important; }

        /* Cards */
        .nr-card {
            background: var(--bg-card);
            border: 1px solid var(--border);
            border-radius: 12px;
            padding: 1.25rem 1.5rem;
            margin-bottom: 0.75rem;
            transition: border-color 0.18s, box-shadow 0.18s;
        }
        .nr-card:hover { border-color: rgba(0,212,255,0.5); box-shadow: 0 0 18px rgba(0,212,255,0.08); }
        .nr-card-top {
            background: linear-gradient(135deg, #0f1526, #13213a);
            border: 1px solid rgba(0,212,255,0.25);
            border-radius: 14px;
            padding: 1.5rem;
            margin-bottom: 0.75rem;
        }

        /* Metrics */
        [data-testid="stMetric"]            { background: var(--bg-elevated); border-radius: 10px; padding: 0.75rem 1rem; border: 1px solid var(--border); }
        [data-testid="stMetricValue"]       { color: var(--accent-cyan) !important; font-size: 1.55rem !important; }
        [data-testid="stMetricLabel"]       { color: var(--text-secondary) !important; }

        /* Inputs */
        .stTextInput input, .stTextArea textarea, .stSelectbox select {
            background: var(--bg-elevated) !important;
            border: 1px solid var(--border) !important;
            border-radius: 8px !important;
            color: var(--text-primary) !important;
        }
        .stTextInput input:focus { border-color: var(--accent-cyan) !important; box-shadow: 0 0 0 2px rgba(0,212,255,0.12) !important; }

        /* Buttons */
        .stButton > button {
            background: linear-gradient(135deg, #00d4ff, #7b2fff) !important;
            color: #fff !important; border: none !important;
            border-radius: 8px !important; font-weight: 600 !important;
            padding: 0.5rem 1.4rem !important;
            transition: opacity 0.2s, transform 0.1s !important;
        }
        .stButton > button:hover { opacity: 0.85; transform: translateY(-1px); }
        .stButton > button[kind="secondary"] {
            background: var(--bg-elevated) !important;
            border: 1px solid var(--border) !important;
            color: var(--text-primary) !important;
        }

        /* DataFrame */
        [data-testid="stDataFrame"] { border-radius: 10px; overflow: hidden; }
        [data-testid="stDataFrame"] thead th { background: var(--bg-elevated) !important; color: var(--accent-cyan) !important; }
        [data-testid="stDataFrame"] tbody tr:nth-child(even) td { background: var(--bg-card) !important; }

        /* Tabs */
        .stTabs [data-baseweb="tab-list"] { gap: 4px; border-bottom: 1px solid var(--border) !important; }
        .stTabs [data-baseweb="tab"]      { background: transparent !important; border-radius: 8px 8px 0 0 !important; color: var(--text-secondary) !important; padding: 0.5rem 1.2rem !important; }
        .stTabs [aria-selected="true"]    { background: var(--bg-elevated) !important; color: var(--accent-cyan) !important; border-top: 2px solid var(--accent-cyan) !important; }

        /* Progress */
        .stProgress > div > div { background: var(--accent-cyan) !important; }
        .stSpinner > div        { border-top-color: var(--accent-cyan) !important; }

        /* Badges */
        .badge-high   { background:rgba(0,255,157,0.12); color:#00ff9d; border:1px solid rgba(0,255,157,0.3); border-radius:5px; padding:2px 8px; font-size:0.75rem; font-weight:700; }
        .badge-medium { background:rgba(255,179,0,0.12); color:#ffb300; border:1px solid rgba(255,179,0,0.3); border-radius:5px; padding:2px 8px; font-size:0.75rem; font-weight:700; }
        .badge-low    { background:rgba(255,68,85,0.12);  color:#ff4455; border:1px solid rgba(255,68,85,0.3);  border-radius:5px; padding:2px 8px; font-size:0.75rem; font-weight:700; }
        .phase-badge  { background:rgba(123,47,255,0.15); color:#ab6fff; border:1px solid rgba(123,47,255,0.3); border-radius:5px; padding:2px 8px; font-size:0.75rem; font-weight:600; }

        /* Hero */
        .nr-hero {
            background: linear-gradient(135deg, #070b14 0%, #0f1a35 60%, #0b1520 100%);
            border: 1px solid var(--border); border-radius: 16px;
            padding: 2.5rem 3rem; margin-bottom: 1.5rem; position: relative; overflow: hidden;
        }
        .nr-hero::before {
            content: ''; position: absolute; top: -80px; right: -80px;
            width: 260px; height: 260px; border-radius: 50%;
            background: radial-gradient(circle, rgba(0,212,255,0.07) 0%, transparent 70%);
            pointer-events: none;
        }
        .nr-hero-title {
            font-size: 2rem; font-weight: 800;
            background: linear-gradient(90deg, #00d4ff, #7b2fff);
            -webkit-background-clip: text; -webkit-text-fill-color: transparent;
            background-clip: text; margin-bottom: 0.5rem; line-height: 1.2;
        }
        .nr-hero-sub { font-size: 1rem; color: var(--text-secondary) !important; }

        /* Score bars */
        .score-bar-wrap   { display:flex; align-items:center; gap:8px; margin:3px 0; }
        .score-bar-label  { width:90px; font-size:0.73rem; color:var(--text-secondary); flex-shrink:0; }
        .score-bar-track  { flex:1; background:#1e2a4a; border-radius:4px; height:5px; overflow:hidden; }
        .score-bar-fill   { height:5px; border-radius:4px; transition:width 0.4s; }
        .score-bar-value  { width:32px; font-size:0.73rem; text-align:right; color:var(--text-primary); flex-shrink:0; }

        /* Status dot */
        .status-row { display:flex; align-items:center; gap:6px; font-size:0.8rem; margin-bottom:3px; }
        .dot-ok  { width:7px; height:7px; border-radius:50%; background:#00ff9d; flex-shrink:0; }
        .dot-err { width:7px; height:7px; border-radius:50%; background:#ff4455; flex-shrink:0; }

        hr { border-color: var(--border) !important; }
        </style>
        """,
        unsafe_allow_html=True,
    )


# ─── Shared helpers ───────────────────────────────────────────────────────────

def _phase_label(phase) -> str:
    try:
        p = int(float(phase or 0))
    except Exception:
        p = 0
    return {4: "FDA Approved", 3: "Phase III", 2: "Phase II", 1: "Phase I"}.get(p, "Preclinical")


def _phase_badge(phase) -> str:
    return f"<span class='phase-badge'>{_phase_label(phase)}</span>"


def _confidence_badge(score: float) -> str:
    if score >= 0.65:
        return f"<span class='badge-high'>{score:.0%}</span>"
    if score >= 0.35:
        return f"<span class='badge-medium'>{score:.0%}</span>"
    return f"<span class='badge-low'>{score:.0%}</span>"


def _bar_color(v: float) -> str:
    if v >= 0.65:
        return "#00ff9d"
    if v >= 0.35:
        return "#ffb300"
    return "#ff4455"


def _score_bars_html(breakdown: Dict) -> str:
    labels = [
        ("Indication",  "indication_score"),
        ("Target",      "target_score"),
        ("Activity",    "activity_score"),
        ("Network",     "network_score"),
    ]
    html = ""
    for label, key in labels:
        v = float(breakdown.get(key) or 0)
        pct = round(v * 100, 1)
        color = _bar_color(v)
        html += (
            f"<div class='score-bar-wrap'>"
            f"<span class='score-bar-label'>{label}</span>"
            f"<div class='score-bar-track'><div class='score-bar-fill' style='width:{pct}%;background:{color};'></div></div>"
            f"<span class='score-bar-value'>{v:.0%}</span>"
            f"</div>"
        )
    return html


def _status_row(label: str, ok: bool):
    dot_cls = "dot-ok" if ok else "dot-err"
    color = "#8896b3" if ok else "#8896b3"
    st.markdown(
        f"<div class='status-row'><span class='{dot_cls}'></span><span style='color:#8896b3;'>{label}</span></div>",
        unsafe_allow_html=True,
    )


# ─── Session state ────────────────────────────────────────────────────────────

def _init_session():
    defaults = {
        "page":              "dashboard",
        "disease_query":     "",
        "compounds":         [],
        "selected_compound": None,
    }
    for k, v in defaults.items():
        if k not in st.session_state:
            st.session_state[k] = v


# ─── Sidebar ──────────────────────────────────────────────────────────────────

def render_sidebar():
    with st.sidebar:
        st.markdown(
            """
            <div style='padding:1.2rem 0 1.6rem 0;'>
              <div style='font-size:1.45rem; font-weight:800;
                          background:linear-gradient(90deg,#00d4ff,#7b2fff);
                          -webkit-background-clip:text;-webkit-text-fill-color:transparent;
                          background-clip:text;'>
                NeuroRepurpose
              </div>
              <div style='font-size:0.75rem; color:#8896b3; margin-top:3px; letter-spacing:0.03em;'>
                Drug Repurposing Intelligence
              </div>
            </div>
            """,
            unsafe_allow_html=True,
        )

        pages = [
            ("dashboard",  "Dashboard"),
            ("discover",   "Discover"),
            ("network",    "Knowledge Network"),
            ("database",   "Data Explorer"),
        ]

        st.markdown("<div style='font-size:0.7rem; color:#8896b3; letter-spacing:0.08em; text-transform:uppercase; margin-bottom:6px;'>Navigation</div>", unsafe_allow_html=True)
        for page_id, label in pages:
            active = st.session_state.page == page_id
            style = "background:rgba(0,212,255,0.1); color:#00d4ff; border-left:3px solid #00d4ff;" if active else "color:#8896b3;"
            if st.button(
                label,
                key=f"nav_{page_id}",
                use_container_width=True,
                type="primary" if active else "secondary",
            ):
                st.session_state.page = page_id
                st.rerun()

        st.markdown("---")
        st.markdown("<div style='font-size:0.7rem; color:#8896b3; letter-spacing:0.08em; text-transform:uppercase; margin-bottom:6px;'>Data Sources</div>", unsafe_allow_html=True)
        _status_row("ChEMBL 33", DB_AVAILABLE)
        _status_row("HetioNet", DB_AVAILABLE)
        _status_row("MeSH Ontology", MESH_AVAILABLE)
        _status_row("Scoring Engine", SCORER_AVAILABLE)
        _status_row("3D Viewer", VIZ_3D_AVAILABLE)

        if DB_AVAILABLE:
            st.markdown("---")
            try:
                stats = get_stats()
                st.markdown("<div style='font-size:0.7rem; color:#8896b3; letter-spacing:0.08em; text-transform:uppercase; margin-bottom:6px;'>Database</div>", unsafe_allow_html=True)
                items = [
                    ("Compounds",   stats.get("compounds", 0)),
                    ("Targets",     stats.get("targets", 0)),
                    ("Indications", stats.get("indications", 0)),
                    ("KG Nodes",    stats.get("hetionet_nodes", 0)),
                    ("Diseases",    stats.get("mesh_diseases", 0)),
                ]
                for label, val in items:
                    st.markdown(
                        f"<div style='font-size:0.8rem; display:flex; justify-content:space-between; color:#8896b3; margin-bottom:2px;'>"
                        f"<span>{label}</span><span style='color:#e8edf7;'>{val:,}</span></div>",
                        unsafe_allow_html=True,
                    )
            except Exception:
                pass

        st.markdown("---")
        st.markdown("<div style='font-size:0.72rem; color:#8896b3;'>ChEMBL 33 · HetioNet · RDKit · MeSH</div>", unsafe_allow_html=True)


# ─── Page: Dashboard ──────────────────────────────────────────────────────────

def page_dashboard():
    st.markdown(
        """
        <div class='nr-hero'>
          <div class='nr-hero-title'>NeuroRepurpose Intelligence Platform</div>
          <div class='nr-hero-sub'>
            Identify new therapeutic indications for existing drugs using ChEMBL 33,
            HetioNet knowledge graphs, MeSH disease ontology, and a multi-signal scoring engine.
            Supports any disease — not just neurology.
          </div>
        </div>
        """,
        unsafe_allow_html=True,
    )

    # Quick search
    col_q, col_btn = st.columns([5, 1])
    with col_q:
        q = st.text_input(
            "disease_search_hero",
            placeholder="Enter disease name — e.g. Alzheimer disease, multiple sclerosis, ALS, depression...",
            label_visibility="collapsed",
            key="hero_q",
        )
    with col_btn:
        if st.button("Search", use_container_width=True, key="hero_btn") and q:
            st.session_state.disease_query = q
            st.session_state.page = "discover"
            st.rerun()

    st.markdown("---")

    # Stats row
    if DB_AVAILABLE:
        try:
            stats = get_stats()
        except Exception:
            stats = {}
    else:
        stats = {}

    col1, col2, col3, col4, col5 = st.columns(5)
    with col1: st.metric("Compounds",      f"{stats.get('compounds', 0):,}")
    with col2: st.metric("Protein Targets", f"{stats.get('targets', 0):,}")
    with col3: st.metric("Drug Indications", f"{stats.get('indications', 0):,}")
    with col4: st.metric("KG Nodes",        f"{stats.get('hetionet_nodes', 0):,}")
    with col5: st.metric("MeSH Diseases",   f"{stats.get('mesh_diseases', 0):,}")

    if not DB_AVAILABLE:
        st.warning("Database unavailable. Run: `python database/importer.py`")
        return

    st.markdown("---")
    st.markdown("### Disease Areas — Quick Launch")
    st.caption("Loaded from database. Click any disease to run drug discovery.")

    diseases = get_available_diseases(limit=24)
    if diseases:
        cols = st.columns(4)
        for i, disease in enumerate(diseases[:16]):
            with cols[i % 4]:
                if st.button(disease, key=f"ql_{i}", use_container_width=True):
                    st.session_state.disease_query = disease
                    st.session_state.page = "discover"
                    st.rerun()
    else:
        st.info("No diseases found in database. Ensure the importer has run.")

    st.markdown("---")
    st.markdown("### How It Works")
    col_a, col_b, col_c, col_d = st.columns(4)
    steps = [
        ("1. Disease Resolution",
         "Your query is matched against the full MeSH disease ontology — exact heading, "
         "entry terms, and synonyms. The canonical MeSH ID is then expanded to include "
         "parent and child disease concepts."),
        ("2. Compound Retrieval",
         "All compounds in ChEMBL 33 with an indication matching any of the expanded "
         "MeSH IDs are retrieved. Duplicate salt forms and hydrates are collapsed "
         "using InChIKey deduplication."),
        ("3. Multi-Signal Scoring",
         "Each compound is scored on four signals: indication evidence (40%), "
         "target mechanistic relevance (30%), pChEMBL activity strength (20%), "
         "and HetioNet network path count (10%)."),
        ("4. Ranked Results",
         "Results are ranked by overall score. Each card shows the full evidence "
         "breakdown so you can evaluate the confidence and rationale for each "
         "repurposing hypothesis."),
    ]
    for col, (title, body) in zip([col_a, col_b, col_c, col_d], steps):
        with col:
            st.markdown(
                f"<div class='nr-card'>"
                f"<div style='font-size:0.85rem; font-weight:700; color:#00d4ff; margin-bottom:0.5rem;'>{title}</div>"
                f"<div style='font-size:0.8rem; color:#8896b3; line-height:1.5;'>{body}</div>"
                f"</div>",
                unsafe_allow_html=True,
            )


# ─── Page: Discover ───────────────────────────────────────────────────────────

@st.cache_data(ttl=1800, show_spinner=False)
def _resolve_and_fetch(disease_query: str):
    """
    Resolve disease query to MeSH IDs, expand, fetch and validate compounds.
    Cached for 30 minutes keyed on the query string.
    Returns (resolved_records, expanded_mesh_ids, compounds)
    """
    resolved = resolve_disease(disease_query)
    if not resolved:
        return [], [], []

    mesh_ids = [r["mesh_id"] for r in resolved if r.get("mesh_id")]
    if not mesh_ids:
        return resolved, [], []

    expanded = expand_mesh_ids(mesh_ids, include_children=True)
    if not expanded:
        expanded = mesh_ids

    compounds = get_compounds_for_disease(expanded, limit=80)
    compounds = validate_and_deduplicate(compounds, require_smiles=False)
    return resolved, expanded, compounds


def page_discover():
    st.markdown("## Drug Discovery")
    st.caption("Enter any disease name. The platform resolves it against MeSH and retrieves ranked repurposing candidates.")

    # Disease input row
    col_inp, col_btn = st.columns([5, 1])
    with col_inp:
        disease_input = st.text_input(
            "disease_input",
            value=st.session_state.disease_query,
            placeholder="e.g. Parkinson disease, MS, epilepsy, schizophrenia, type 2 diabetes...",
            label_visibility="collapsed",
            key="discover_disease_input",
        )
    with col_btn:
        run = st.button("Analyse", use_container_width=True, key="discover_run_btn")

    # Autocomplete suggestions
    if disease_input and len(disease_input) >= 2 and MESH_AVAILABLE:
        suggestions = suggest_diseases(disease_input, limit=6)
        if suggestions and disease_input.strip().lower() not in [s.lower() for s in suggestions]:
            st.markdown(
                "<div style='font-size:0.78rem; color:#8896b3; margin-bottom:4px;'>Suggestions:</div>",
                unsafe_allow_html=True,
            )
            sug_cols = st.columns(min(len(suggestions), 6))
            for i, sug in enumerate(suggestions):
                with sug_cols[i]:
                    if st.button(sug, key=f"sug_{i}", use_container_width=True):
                        st.session_state.disease_query = sug
                        st.rerun()

    if run and disease_input.strip():
        st.session_state.disease_query = disease_input.strip()

    if not st.session_state.disease_query:
        st.info("Enter a disease name above to begin.")
        return

    # Load results
    with st.spinner("Resolving disease and retrieving compounds..."):
        try:
            resolved_records, expanded_ids, compounds = _resolve_and_fetch(st.session_state.disease_query)
        except Exception as e:
            st.error(f"Query failed: {e}")
            return

    if not resolved_records:
        st.warning(
            f"No MeSH match found for '{st.session_state.disease_query}'. "
            "Try a different spelling or use the suggestions above."
        )
        if not MESH_AVAILABLE:
            st.info("MeSH table is empty. Run: `python database/mesh_importer.py`")
        return

    # Resolution info
    primary = resolved_records[0]
    st.markdown(
        f"<div style='background:#0f1526; border:1px solid #1e2a4a; border-radius:10px; padding:1rem 1.25rem; margin-bottom:1rem;'>"
        f"<div style='font-size:0.72rem; color:#8896b3; text-transform:uppercase; letter-spacing:0.06em; margin-bottom:4px;'>Resolved disease</div>"
        f"<div style='font-size:1rem; font-weight:700; color:#00d4ff;'>{primary.get('heading', st.session_state.disease_query)}</div>"
        f"<div style='font-size:0.78rem; color:#8896b3; margin-top:3px;'>"
        f"MeSH ID: {primary.get('mesh_id', 'N/A')} &nbsp;|&nbsp; "
        f"{len(resolved_records)} record(s) matched &nbsp;|&nbsp; "
        f"{len(expanded_ids)} expanded IDs (parents + children)"
        f"</div>"
        f"</div>",
        unsafe_allow_html=True,
    )

    if not compounds:
        st.warning("No compounds found for this disease in the database.")
        st.caption("The disease was resolved but no ChEMBL indications match. Try a broader term.")
        return

    # Score compounds
    with st.spinner(f"Scoring {len(compounds)} compounds across 4 evidence signals..."):
        try:
            scored = score_compound_list(list(compounds), expanded_ids)
        except Exception as e:
            logger.warning(f"Scoring failed: {e}")
            scored = compounds
            for c in scored:
                c.setdefault("score", float(c.get("max_phase") or 0) / 4.0)
                c.setdefault("score_breakdown", {})

    st.markdown("---")

    # Top-3 highlight
    top3 = [c for c in scored if float(c.get("max_phase") or 0) >= 3][:3]
    if not top3:
        top3 = scored[:3]

    if top3:
        st.markdown("### Top Candidates")
        t_cols = st.columns(len(top3))
        for i, c in enumerate(top3):
            _top_compound_card(c, t_cols[i], idx=i)

    st.markdown("---")

    # Filters
    col_f1, col_f2, col_f3, col_f4 = st.columns([1, 1, 1, 2])
    with col_f1:
        min_phase = st.selectbox("Min phase", [0, 1, 2, 3, 4], index=0, key="disc_phase")
    with col_f2:
        show_n = st.slider("Show results", 5, min(80, len(scored)), min(30, len(scored)), key="disc_n")
    with col_f3:
        sort_by = st.selectbox("Sort by", ["Score", "Phase", "Name"], key="disc_sort")

    filtered = [c for c in scored if float(c.get("max_phase") or 0) >= min_phase]
    if sort_by == "Name":
        filtered.sort(key=lambda x: (x.get("name") or "").lower())
    elif sort_by == "Phase":
        filtered.sort(key=lambda x: float(x.get("max_phase") or 0), reverse=True)

    filtered = filtered[:show_n]

    # Phase distribution chart
    if len(filtered) >= 4:
        with st.expander("Phase distribution", expanded=False):
            _render_phase_chart(filtered)

    # Results list
    st.markdown(f"### {len(filtered)} Repurposing Candidates")
    for i, c in enumerate(filtered):
        _compound_list_card(c, idx=i)

    # Download
    if filtered:
        df_dl = pd.DataFrame([{
            "Name":        c.get("name"),
            "ChEMBL ID":   c.get("chembl_id"),
            "Phase":       c.get("max_phase"),
            "Score":       round(float(c.get("score") or 0), 4),
            "Indication":  round(float((c.get("score_breakdown") or {}).get("indication_score") or 0), 4),
            "Target":      round(float((c.get("score_breakdown") or {}).get("target_score") or 0), 4),
            "Activity":    round(float((c.get("score_breakdown") or {}).get("activity_score") or 0), 4),
            "Network":     round(float((c.get("score_breakdown") or {}).get("network_score") or 0), 4),
            "Mechanisms":  c.get("mechanisms", ""),
            "Targets":     c.get("targets", ""),
            "SMILES":      c.get("smiles", ""),
        } for c in filtered])
        st.download_button(
            "Download CSV",
            data=df_dl.to_csv(index=False),
            file_name=f"neurorepurpose_{primary.get('mesh_id','unknown')}.csv",
            mime="text/csv",
        )


def _top_compound_card(c: Dict, col, idx: int):
    name     = c.get("name") or "Unknown"
    score    = float(c.get("score") or 0)
    phase    = c.get("max_phase", 0)
    mechs    = (c.get("mechanisms") or "—")[:100]
    breakdown = c.get("score_breakdown") or {}

    with col:
        st.markdown(
            f"<div class='nr-card-top'>"
            f"<div style='display:flex; justify-content:space-between; align-items:flex-start; margin-bottom:0.6rem;'>"
            f"  <div style='font-size:0.95rem; font-weight:700; color:#e8edf7;'>{name}</div>"
            f"  <div>{_confidence_badge(score)}</div>"
            f"</div>"
            f"<div style='margin-bottom:0.5rem;'>{_phase_badge(phase)}</div>"
            f"<div style='font-size:0.75rem; color:#8896b3; margin-bottom:0.7rem;'>{mechs}</div>"
            f"{_score_bars_html(breakdown)}"
            f"</div>",
            unsafe_allow_html=True,
        )
        if st.button("Full analysis", key=f"top_analyse_{idx}_{name[:6]}", use_container_width=True):
            st.session_state.selected_compound = c
            st.session_state.page = "analysis"
            st.rerun()


def _compound_list_card(c: Dict, idx: int):
    name     = c.get("name") or "Unknown"
    score    = float(c.get("score") or 0)
    phase    = c.get("max_phase", 0)
    mechs    = (c.get("mechanisms") or "—")[:140]
    targets  = (c.get("targets") or "—")[:120]
    breakdown = c.get("score_breakdown") or {}

    with st.container():
        col_main, col_action = st.columns([6, 1])
        with col_main:
            st.markdown(
                f"<div class='nr-card'>"
                f"<div style='display:flex; justify-content:space-between; align-items:flex-start; margin-bottom:0.5rem;'>"
                f"  <div>"
                f"    <span style='font-size:0.95rem; font-weight:700; color:#e8edf7;'>{name}</span>"
                f"    &nbsp; {_phase_badge(phase)}"
                f"  </div>"
                f"  <div>{_confidence_badge(score)}</div>"
                f"</div>"
                f"<div style='font-size:0.78rem; color:#8896b3; margin-bottom:0.4rem;'>"
                f"  <b style='color:#c8d0e0;'>Mechanism:</b> {mechs}"
                f"</div>"
                f"<div style='font-size:0.78rem; color:#8896b3; margin-bottom:0.6rem;'>"
                f"  <b style='color:#c8d0e0;'>Targets:</b> {targets}"
                f"</div>"
                f"{_score_bars_html(breakdown)}"
                f"</div>",
                unsafe_allow_html=True,
            )
        with col_action:
            st.markdown("<div style='height:1.2rem;'></div>", unsafe_allow_html=True)
            if st.button("Analyse", key=f"analyse_{idx}_{name[:8]}", use_container_width=True):
                st.session_state.selected_compound = c
                st.session_state.page = "analysis"
                st.rerun()


def _render_phase_chart(compounds: List[Dict]):
    phase_map = {4: "FDA Approved", 3: "Phase III", 2: "Phase II", 1: "Phase I", 0: "Preclinical"}
    counts: Dict[str, int] = {}
    for c in compounds:
        label = phase_map.get(int(float(c.get("max_phase") or 0)), "Preclinical")
        counts[label] = counts.get(label, 0) + 1
    df_p = pd.DataFrame(list(counts.items()), columns=["Phase", "Count"])
    fig = px.pie(
        df_p, values="Count", names="Phase", hole=0.5,
        color_discrete_sequence=["#00d4ff", "#7b2fff", "#00ff9d", "#ffb300", "#ff4455"],
    )
    fig.update_layout(
        paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
        font_color="#e8edf7", margin=dict(t=10, b=10, l=10, r=10),
        legend=dict(font=dict(color="#8896b3")),
    )
    st.plotly_chart(fig, use_container_width=True)


# ─── Page: Compound Analysis ──────────────────────────────────────────────────

def page_analysis():
    compound = st.session_state.selected_compound
    if not compound:
        st.info("Select a compound from the Discover page.")
        if st.button("Go to Discover"):
            st.session_state.page = "discover"
            st.rerun()
        return

    name      = compound.get("name") or "Unknown"
    smiles    = compound.get("smiles") or ""
    chembl_id = compound.get("chembl_id") or ""
    cid       = compound.get("id")
    score     = float(compound.get("score") or 0)
    breakdown = compound.get("score_breakdown") or {}
    mechs     = compound.get("mechanisms") or "—"

    # Header
    col_h, col_meta = st.columns([3, 1])
    with col_h:
        st.markdown(f"## {name}")
        if chembl_id:
            st.caption(f"ChEMBL ID: {chembl_id}")
        st.markdown(f"**Mechanism of action:** {mechs[:300] if mechs else '—'}")
    with col_meta:
        st.markdown(
            f"<div style='text-align:right;'>"
            f"{_phase_badge(compound.get('max_phase', 0))}"
            f"<br><br>"
            f"{_confidence_badge(score)} overall score"
            f"</div>",
            unsafe_allow_html=True,
        )
        if breakdown:
            st.markdown(_score_bars_html(breakdown), unsafe_allow_html=True)

    st.markdown("---")

    tab_labels = ["Properties", "Targets", "Activities", "Indications", "Network"]
    if QUANTUM_AVAILABLE:
        tab_labels.append("Quantum")
    if VIZ_3D_AVAILABLE:
        tab_labels.append("3D Structure")

    tabs = st.tabs(tab_labels)
    tab_idx = 0

    # ── Properties ───────────────────────────────────────────────────────────
    with tabs[tab_idx]:
        _render_properties(compound, smiles)
    tab_idx += 1

    # ── Targets ──────────────────────────────────────────────────────────────
    with tabs[tab_idx]:
        targets = []
        if DB_AVAILABLE and cid:
            targets = get_compound_targets(int(cid))
        if targets:
            df_t = pd.DataFrame(targets)
            keep = [col for col in ["name", "gene_symbol", "target_type", "mechanism", "action_type", "confidence", "organism"] if col in df_t.columns]
            st.dataframe(df_t[keep], use_container_width=True)
            if "confidence" in df_t.columns and len(df_t) >= 2:
                top_t = df_t.dropna(subset=["confidence"]).head(12)
                name_col = "gene_symbol" if "gene_symbol" in top_t.columns else "name"
                fig = px.bar(
                    top_t.sort_values("confidence"), x="confidence", y=name_col,
                    orientation="h", color="confidence",
                    color_continuous_scale=[[0, "#ff4455"], [0.5, "#ffb300"], [1, "#00ff9d"]],
                    title="Target confidence",
                )
                fig.update_layout(paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)", font_color="#e8edf7", margin=dict(t=30, b=10))
                st.plotly_chart(fig, use_container_width=True)
        else:
            st.info("No target data available.")
    tab_idx += 1

    # ── Activities ────────────────────────────────────────────────────────────
    with tabs[tab_idx]:
        activities = []
        if DB_AVAILABLE and cid:
            activities = get_compound_activities(int(cid))
        if activities:
            df_a = pd.DataFrame(activities)
            st.dataframe(df_a, use_container_width=True)
            if "pchembl_value" in df_a.columns:
                df_plot = df_a.dropna(subset=["pchembl_value"]).head(15).sort_values("pchembl_value", ascending=False)
                if not df_plot.empty:
                    name_col = "gene_symbol" if "gene_symbol" in df_plot.columns else "target_name"
                    fig = px.bar(
                        df_plot, x="pchembl_value", y=name_col,
                        orientation="h",
                        color="pchembl_value",
                        color_continuous_scale=[[0, "#7b2fff"], [0.5, "#ffb300"], [1, "#00ff9d"]],
                        title="pChEMBL values (higher = more potent)",
                    )
                    fig.add_vline(x=6, line_dash="dash", line_color="#ffb300", annotation_text="IC50 1uM")
                    fig.add_vline(x=8, line_dash="dash", line_color="#00ff9d", annotation_text="IC50 10nM")
                    fig.update_layout(paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)", font_color="#e8edf7")
                    st.plotly_chart(fig, use_container_width=True)
        else:
            st.info("No quantitative activity data available for this compound.")
    tab_idx += 1

    # ── Indications ───────────────────────────────────────────────────────────
    with tabs[tab_idx]:
        indications = []
        if DB_AVAILABLE and cid:
            indications = get_compound_indications(int(cid))
        if indications:
            df_i = pd.DataFrame(indications)
            keep = [col for col in ["disease", "mesh_id", "max_phase", "source"] if col in df_i.columns]
            st.dataframe(df_i[keep], use_container_width=True)
        elif compound.get("indications"):
            for ind in str(compound["indications"]).split(";"):
                st.markdown(f"- {ind.strip()}")
        else:
            st.info("No indication data available.")
    tab_idx += 1

    # ── Network ───────────────────────────────────────────────────────────────
    with tabs[tab_idx]:
        _render_compound_network(name, compound)
    tab_idx += 1

    # ── Quantum ───────────────────────────────────────────────────────────────
    if QUANTUM_AVAILABLE:
        with tabs[tab_idx]:
            render_quantum_optimization_section(drug_name=name, smiles=smiles, drug_data=compound)
        tab_idx += 1

    # ── 3D Structure ──────────────────────────────────────────────────────────
    if VIZ_3D_AVAILABLE:
        with tabs[tab_idx]:
            if not smiles:
                st.warning(f"No SMILES available for {name}.")
            elif _viz:
                targets_list = get_compound_targets(int(cid)) if (DB_AVAILABLE and cid) else []
                target_name = targets_list[0]["name"] if targets_list else ""
                pdb_content = None
                if target_name:
                    try:
                        from real_pdb_fetcher import RealPDBFetcher
                        pdb_content = RealPDBFetcher().fetch_pdb(target_name)
                    except Exception:
                        pass
                _viz.render_visualization(name, target_name or "", smiles, pdb_content)
        tab_idx += 1


def _render_properties(compound: Dict, smiles: str):
    mw   = compound.get("mw")
    logp = compound.get("alogp")
    psa  = compound.get("psa")
    hba  = compound.get("hba")
    hbd  = compound.get("hbd")
    qed  = None

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
        except Exception:
            pass

    col1, col2, col3, col4, col5 = st.columns(5)
    with col1: st.metric("Mol Weight (Da)", f"{mw:.1f}" if mw else "N/A")
    with col2: st.metric("LogP", f"{logp:.2f}" if logp is not None else "N/A")
    with col3: st.metric("TPSA (sq A)", f"{psa:.1f}" if psa else "N/A")
    with col4: st.metric("QED", f"{qed:.3f}" if qed is not None else "N/A")
    with col5: st.metric("HBA / HBD", f"{hba} / {hbd}" if hba is not None else "N/A")

    violations = sum([
        bool(mw and mw > 500),
        bool(logp is not None and logp > 5),
        bool(hbd is not None and hbd > 5),
        bool(hba is not None and hba > 10),
    ])
    ro5_text = "Passes Lipinski Rule of 5" if violations == 0 else f"Lipinski violations: {violations}"
    ro5_color = "#00ff9d" if violations == 0 else "#ffb300" if violations == 1 else "#ff4455"
    st.markdown(f"<div style='margin:0.75rem 0; font-size:0.85rem; color:{ro5_color};'>{ro5_text}</div>", unsafe_allow_html=True)

    if smiles:
        st.markdown("**SMILES**")
        st.code(smiles, language=None)

    ro5_items = [
        ("Molecular weight", f"{mw:.1f} Da" if mw else "N/A", bool(mw and mw > 500)),
        ("LogP", f"{logp:.2f}" if logp is not None else "N/A", bool(logp is not None and logp > 5)),
        ("H-bond donors", str(hbd) if hbd is not None else "N/A", bool(hbd is not None and hbd > 5)),
        ("H-bond acceptors", str(hba) if hba is not None else "N/A", bool(hba is not None and hba > 10)),
    ]
    st.markdown("**Lipinski Rule of 5 Breakdown**")
    for prop, val, violated in ro5_items:
        color = "#ff4455" if violated else "#00ff9d"
        mark  = "FAIL" if violated else "PASS"
        st.markdown(
            f"<div style='font-size:0.82rem; display:flex; justify-content:space-between; "
            f"background:#0f1526; border-radius:6px; padding:4px 12px; margin-bottom:3px;'>"
            f"<span style='color:#8896b3;'>{prop}</span>"
            f"<span style='color:#e8edf7;'>{val}</span>"
            f"<span style='color:{color}; font-weight:700; font-size:0.72rem;'>{mark}</span>"
            f"</div>",
            unsafe_allow_html=True,
        )


def _render_compound_network(compound_name: str, compound: Dict):
    query = st.session_state.get("disease_query", "")
    mesh_ids = []
    if query:
        try:
            _, mesh_ids, _ = _resolve_and_fetch(query)
        except Exception:
            pass

    paths = []
    if DB_AVAILABLE and mesh_ids:
        with st.spinner("Loading HetioNet paths..."):
            try:
                paths = get_hetionet_paths(compound_name, mesh_ids)
            except Exception as e:
                st.warning(f"Network query failed: {e}")

    if not paths:
        st.info("No HetioNet paths found for this compound. The compound may not be in the HetioNet node set.")
        return

    G = nx.DiGraph()
    for p in paths[:60]:
        src = p.get("source_name", "")
        tgt = p.get("target_name", "")
        edge = p.get("metaedge", "")
        src_kind = p.get("source_kind", "")
        tgt_kind = p.get("target_kind", "")
        if src and tgt:
            G.add_node(src, kind=src_kind)
            G.add_node(tgt, kind=tgt_kind)
            G.add_edge(src, tgt, metaedge=edge)

    if not G.nodes:
        st.info("No network paths to render.")
        return

    pos = nx.spring_layout(G, seed=42, k=1.2)
    kind_color = {
        "Compound": "#00d4ff", "Gene": "#00ff9d", "Disease": "#ff4455",
        "Anatomy": "#ffb300", "Pathway": "#7b2fff", "Biological Process": "#ab6fff",
    }
    default_color = "#8896b3"

    node_x, node_y, node_text, node_color = [], [], [], []
    for node, attr in G.nodes(data=True):
        x, y = pos[node]
        node_x.append(x); node_y.append(y)
        node_text.append(node[:30])
        node_color.append(kind_color.get(attr.get("kind", ""), default_color))

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
        marker=dict(size=14, color=node_color, line=dict(width=1, color="#070b14")),
        text=node_text, textposition="top center",
        textfont=dict(size=8, color="#8896b3"),
        hoverinfo="text",
    ))
    fig.update_layout(
        showlegend=False, hovermode="closest", height=500,
        paper_bgcolor="#070b14", plot_bgcolor="#070b14", font_color="#e8edf7",
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        margin=dict(l=10, r=10, t=10, b=10),
    )
    st.plotly_chart(fig, use_container_width=True)

    # Legend
    leg_cols = st.columns(len(kind_color))
    for col, (kind, color) in zip(leg_cols, kind_color.items()):
        col.markdown(f"<span style='color:{color}; font-size:0.78rem;'>&#9679; {kind}</span>", unsafe_allow_html=True)

    st.caption(f"HetioNet: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges shown")


# ─── Page: Knowledge Network ──────────────────────────────────────────────────

def page_network():
    st.markdown("## Knowledge Network")
    st.caption("Explore drug-gene-disease connections in the HetioNet knowledge graph.")

    col_inp, col_btn = st.columns([4, 1])
    with col_inp:
        center = st.text_input(
            "network_center",
            value=st.session_state.disease_query or "",
            placeholder="Enter disease or drug name...",
            label_visibility="collapsed",
            key="net_center",
        )
    with col_btn:
        build = st.button("Build Network", use_container_width=True)

    if not build and not center:
        st.info("Enter a disease or drug name to build the network.")
        return

    if not DB_AVAILABLE:
        st.warning("Database unavailable.")
        return

    with st.spinner("Building network..."):
        resolved = resolve_disease(center) if center else []
        mesh_ids = [r["mesh_id"] for r in resolved if r.get("mesh_id")]
        expanded = expand_mesh_ids(mesh_ids) if mesh_ids else []

        compounds = get_compounds_for_disease(expanded or mesh_ids, limit=20) if expanded or mesh_ids else db_search_compounds(center, limit=20)
        compounds = validate_and_deduplicate(compounds, require_smiles=False)[:20]

    if not compounds:
        st.warning("No compounds found. Try a different term.")
        return

    G = nx.Graph()
    center_label = resolved[0]["heading"] if resolved else center
    G.add_node(center_label, kind="Disease" if resolved else "Query")

    for c in compounds:
        name = c.get("name") or ""
        if not name:
            continue
        G.add_node(name, kind="Compound", phase=float(c.get("max_phase") or 0))
        G.add_edge(center_label, name)
        for tgt in (c.get("targets") or "").split(";")[:2]:
            tgt = tgt.strip()
            if tgt and len(tgt) > 2:
                G.add_node(tgt, kind="Gene")
                G.add_edge(name, tgt)

    pos = nx.spring_layout(G, seed=42, k=1.4)
    kind_color = {"Disease": "#ff4455", "Query": "#00d4ff", "Compound": "#7b2fff", "Gene": "#00ff9d"}

    node_x, node_y, node_text, node_color, node_size = [], [], [], [], []
    for node, attr in G.nodes(data=True):
        x, y = pos[node]
        node_x.append(x); node_y.append(y)
        node_text.append(node[:28])
        node_color.append(kind_color.get(attr.get("kind", ""), "#8896b3"))
        node_size.append(24 if attr.get("kind") in ("Disease", "Query") else 14 if attr.get("kind") == "Compound" else 10)

    edge_x, edge_y = [], []
    for u, v in G.edges():
        x0, y0 = pos[u]; x1, y1 = pos[v]
        edge_x += [x0, x1, None]; edge_y += [y0, y1, None]

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=edge_x, y=edge_y, mode="lines", line=dict(color="#1e2a4a", width=1), hoverinfo="none"))
    fig.add_trace(go.Scatter(
        x=node_x, y=node_y, mode="markers+text",
        marker=dict(size=node_size, color=node_color, line=dict(width=1, color="#070b14")),
        text=node_text, textposition="top center",
        textfont=dict(size=8, color="#8896b3"),
        hoverinfo="text",
    ))
    fig.update_layout(
        showlegend=False, hovermode="closest", height=560,
        paper_bgcolor="#070b14", plot_bgcolor="#070b14", font_color="#e8edf7",
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        margin=dict(l=10, r=10, t=10, b=10),
    )
    st.plotly_chart(fig, use_container_width=True)

    leg_cols = st.columns(4)
    for col, (kind, color) in zip(leg_cols, kind_color.items()):
        col.markdown(f"<span style='color:{color}; font-size:0.8rem;'>&#9679; {kind}</span>", unsafe_allow_html=True)

    st.caption(f"{G.number_of_nodes()} nodes, {G.number_of_edges()} edges")


# ─── Page: Data Explorer ──────────────────────────────────────────────────────

def page_database():
    st.markdown("## Data Explorer")

    if not DB_AVAILABLE:
        st.warning("Database unavailable. Run: `python database/importer.py`")
        return

    tab_compounds, tab_targets, tab_diseases, tab_stats = st.tabs(
        ["Compounds", "Targets", "Diseases", "Database Stats"]
    )

    with tab_compounds:
        st.markdown("### Browse Compounds")
        col_q, col_lim = st.columns([3, 1])
        with col_q:
            q = st.text_input("Filter", placeholder="Name, mechanism, or target...", key="db_q", label_visibility="collapsed")
        with col_lim:
            limit = st.number_input("Rows", 10, 500, 50, key="db_lim")

        if q:
            compounds = db_search_compounds(q, limit=int(limit))
        else:
            compounds = db_search_compounds("", limit=int(limit)) or []
            if not compounds:
                try:
                    from database.schema import get_connection
                    with get_connection() as conn:
                        cur = conn.cursor()
                        cur.execute(
                            "SELECT c.chembl_id, c.name, c.smiles, c.max_phase, cp.mw, cp.alogp "
                            "FROM compounds c LEFT JOIN compound_properties cp ON cp.compound_id=c.id "
                            "ORDER BY c.max_phase DESC NULLS LAST LIMIT %s", (limit,)
                        )
                        compounds = [
                            {"chembl_id": r[0], "name": r[1], "smiles": r[2], "max_phase": r[3], "mw": r[4], "alogp": r[5]}
                            for r in cur.fetchall()
                        ]
                except Exception as e:
                    st.error(f"Query error: {e}")

        if compounds:
            df = pd.DataFrame(compounds)
            keep = [c for c in ["chembl_id", "name", "max_phase", "mw", "alogp", "psa", "mechanisms", "targets"] if c in df.columns]
            st.dataframe(df[keep], use_container_width=True, height=420)
            st.download_button(
                "Download CSV",
                data=df.to_csv(index=False),
                file_name="neurorepurpose_compounds.csv",
                mime="text/csv",
            )
        else:
            st.info("No results. Try a search term or check the importer.")

    with tab_targets:
        st.markdown("### Protein Targets")
        try:
            from database.schema import get_connection
            with get_connection() as conn:
                cur = conn.cursor()
                cur.execute(
                    "SELECT chembl_tid, name, target_type, gene_symbol, organism "
                    "FROM targets ORDER BY name LIMIT 500"
                )
                rows = cur.fetchall()
            if rows:
                df_t = pd.DataFrame(rows, columns=["ChEMBL TID", "Name", "Type", "Gene", "Organism"])
                st.dataframe(df_t, use_container_width=True, height=420)
            else:
                st.info("No targets in database.")
        except Exception as e:
            st.error(f"Query error: {e}")

    with tab_diseases:
        st.markdown("### MeSH Disease Ontology")
        try:
            from database.schema import get_connection
            with get_connection() as conn:
                cur = conn.cursor()
                cur.execute(
                    "SELECT mesh_id, heading, array_length(tree_numbers,1) AS tree_count, "
                    "array_length(entry_terms,1) AS synonym_count "
                    "FROM mesh_diseases ORDER BY heading LIMIT 500"
                )
                rows = cur.fetchall()
            if rows:
                df_d = pd.DataFrame(rows, columns=["MeSH ID", "Heading", "Tree Entries", "Synonyms"])
                st.dataframe(df_d, use_container_width=True, height=420)
            else:
                st.info("MeSH table is empty. Run: `python database/mesh_importer.py`")
        except Exception as e:
            st.error(f"Query error: {e}")

    with tab_stats:
        st.markdown("### Database Statistics")
        try:
            stats = get_stats()
            col1, col2, col3 = st.columns(3)
            items = list(stats.items())
            for i, (table, count) in enumerate(items):
                with [col1, col2, col3][i % 3]:
                    st.metric(table.replace("_", " ").title(), f"{count:,}")

            # HetioNet node type breakdown
            from database.schema import get_connection
            with get_connection() as conn:
                cur = conn.cursor()
                cur.execute("SELECT kind, COUNT(*) AS n FROM hetionet_nodes GROUP BY kind ORDER BY n DESC")
                kind_rows = cur.fetchall()
            if kind_rows:
                st.markdown("### HetioNet Node Types")
                df_k = pd.DataFrame(kind_rows, columns=["Kind", "Count"])
                fig = px.bar(
                    df_k, x="Kind", y="Count",
                    color="Count",
                    color_continuous_scale=[[0, "#7b2fff"], [1, "#00d4ff"]],
                )
                fig.update_layout(paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)", font_color="#e8edf7")
                st.plotly_chart(fig, use_container_width=True)
        except Exception as e:
            st.error(f"Stats error: {e}")


# ─── Main router ──────────────────────────────────────────────────────────────

def main():
    _init_session()
    apply_theme()
    render_sidebar()

    page = st.session_state.page
    if page == "dashboard":
        page_dashboard()
    elif page == "discover":
        page_discover()
    elif page == "analysis":
        page_analysis()
    elif page == "network":
        page_network()
    elif page == "database":
        page_database()
    else:
        page_dashboard()


if __name__ == "__main__":
    main()
