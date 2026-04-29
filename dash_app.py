#!/usr/bin/env python3
"""
NeuroRepurpose Intelligence Platform — Dash Application
Replaces Streamlit; no framework branding, full design control.
"""

import json
import logging
import os
import warnings
warnings.filterwarnings("ignore", message=".*use_container_width.*")
warnings.filterwarnings("ignore", category=DeprecationWarning, module="streamlit")
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import requests

import dash
from dash import ALL, MATCH, Input, Output, State, callback, ctx, dcc, html, no_update
import dash_cytoscape as cyto

cyto.load_extra_layouts()

logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)

# ── Service imports ────────────────────────────────────────────────────────────
try:
    from services.neuro_db_service import (
        get_compound_activities, get_compound_indications, get_compound_targets,
        get_compounds_for_disease, get_stats, is_available as db_available,
        search_compounds as db_search,
    )
    DB_OK = db_available()
except Exception:
    DB_OK = False
    def get_compounds_for_disease(*a, **k): return []
    def db_search(*a, **k): return []
    def get_compound_targets(*a, **k): return []
    def get_compound_activities(*a, **k): return []
    def get_compound_indications(*a, **k): return []
    def get_stats(): return {}

try:
    from services.disease_resolver import (
        expand_mesh_ids, mesh_available, resolve_disease, suggest_diseases,
    )
    MESH_OK = mesh_available() if DB_OK else False
except Exception:
    MESH_OK = False
    def resolve_disease(*a, **k): return []
    def expand_mesh_ids(*a, **k): return []
    def suggest_diseases(*a, **k): return []

try:
    from services.repurposing_scorer import score_compound_list
    SCORER_OK = True
except Exception:
    SCORER_OK = False
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
    PBPK_OK = True
except Exception:
    PBPK_OK = False
    PBPKSimulator = None

try:
    from services.docking_service import DockingService
    _dock_svc = DockingService()
    DOCK_OK = _dock_svc.available
except Exception:
    DOCK_OK = False
    _dock_svc = None

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
    from services.biocypher_graph import build_graph, CYTO_STYLESHEET
    GRAPH_OK = True
except Exception:
    GRAPH_OK = False
    def build_graph(*a, **k): return [], {}
    CYTO_STYLESHEET = []


# ── Helpers ────────────────────────────────────────────────────────────────────

def _clr(val: float, lo: float, hi: float, palette: tuple) -> str:
    """Map value in [lo,hi] to hex colour from a gradient palette (dark → bright)."""
    if hi <= lo:
        return palette[0]
    n = max(0.0, min(1.0, (val - lo) / (hi - lo)))
    idx = min(int(n * (len(palette) - 1)), len(palette) - 2)
    t = n * (len(palette) - 1) - idx
    c0, c1 = int(palette[idx][1:], 16), int(palette[idx + 1][1:], 16)
    r = int((c0 >> 16) * (1 - t) + (c1 >> 16) * t)
    g = int(((c0 >> 8) & 0xFF) * (1 - t) + ((c1 >> 8) & 0xFF) * t)
    b = int((c0 & 0xFF) * (1 - t) + (c1 & 0xFF) * t)
    return f"#{r:02x}{g:02x}{b:02x}"


def _score_color(s: float) -> str:
    if s >= 0.65: return "#10b981"
    if s >= 0.40: return "#0ea5e9"
    if s >= 0.20: return "#f59e0b"
    return "#64748b"


def _phase_class(p) -> str:
    p = int(p or 0)
    return f"phase-{min(p, 4)}"


def _phase_label(p) -> str:
    p = int(p or 0)
    labels = {4: "Phase 4 / Approved", 3: "Phase 3", 2: "Phase 2", 1: "Phase 1"}
    return labels.get(p, "Preclinical")


def _plotly_theme(dark: bool) -> dict:
    bg = "#020817" if dark else "#ffffff"
    grid = "#1e3a5f" if dark else "#e2e8f0"
    text = "#94a3b8" if dark else "#64748b"
    return dict(
        paper_bgcolor=bg, plot_bgcolor=bg,
        font_color=text, font_size=11,
        xaxis=dict(gridcolor=grid, linecolor=grid, zerolinecolor=grid),
        yaxis=dict(gridcolor=grid, linecolor=grid, zerolinecolor=grid),
        margin=dict(l=40, r=20, t=30, b=40),
    )


def _empty_fig(msg: str, dark: bool) -> go.Figure:
    fig = go.Figure()
    fig.add_annotation(text=msg, xref="paper", yref="paper", x=0.5, y=0.5,
                       showarrow=False, font=dict(size=13, color="#64748b"))
    fig.update_layout(**_plotly_theme(dark), height=250)
    return fig


def _resolve_and_fetch(query: str):
    resolved = resolve_disease(query)
    if not resolved:
        return [], [], []
    mesh_ids = [r["mesh_id"] for r in resolved if r.get("mesh_id")]
    expanded = expand_mesh_ids(mesh_ids, include_children=True) or mesh_ids
    compounds = get_compounds_for_disease(expanded, limit=80)
    compounds = validate_and_deduplicate(compounds, require_smiles=False)
    if SCORER_OK:
        compounds = score_compound_list(compounds, expanded)
    return resolved, expanded, compounds


def _fetch_trials(drug: str, disease: str, n: int = 8) -> List[Dict]:
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
                nct = pm.get("identificationModule", {}).get("nctId", "")
                out.append({
                    "nct_id": nct,
                    "title": pm.get("identificationModule", {}).get("briefTitle", "")[:90],
                    "status": pm.get("statusModule", {}).get("overallStatus", ""),
                    "phase": ", ".join(phases) or "N/A",
                    "url": f"https://clinicaltrials.gov/study/{nct}",
                })
            return out
    except Exception:
        pass
    return []


def _fetch_papers(drug: str, disease: str, n: int = 8) -> List[Dict]:
    try:
        base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
        sr = requests.get(f"{base}/esearch.fcgi",
            params={"db": "pubmed", "term": f"{drug}[tiab] AND {disease}[tiab]",
                    "retmax": n, "retmode": "json"}, timeout=10)
        ids = sr.json().get("esearchresult", {}).get("idlist", [])
        if not ids:
            return []
        smr = requests.get(f"{base}/esummary.fcgi",
            params={"db": "pubmed", "id": ",".join(ids), "retmode": "json"}, timeout=10)
        res = smr.json().get("result", {})
        return [{"pmid": p, "title": res[p].get("title", "")[:90],
                 "journal": res[p].get("source", ""), "year": res[p].get("pubdate", "")[:4],
                 "url": f"https://pubmed.ncbi.nlm.nih.gov/{p}/"} for p in ids if p in res]
    except Exception:
        pass
    return []


def _generate_3d_html(smiles: str, dark: bool = True) -> str:
    """Return standalone HTML string for py3Dmol viewer (embedded in iframe)."""
    if not RDKIT_OK or not PY3DMOL_OK:
        return "<body style='background:#020817;color:#64748b;font-family:sans-serif;padding:2rem'><p>3D viewer requires RDKit + py3Dmol.</p></body>"
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return "<body style='background:#020817;color:#64748b;font-family:sans-serif;padding:2rem'><p>Invalid SMILES.</p></body>"
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        AllChem.MMFFOptimizeMolecule(mol)
        molblock = Chem.MolToMolBlock(mol)
        bg = "#020817" if dark else "#f8fafc"
        view = py3Dmol.view(width=760, height=440, linked=False)
        view.addModel(molblock, "sdf")
        view.setStyle({}, {"stick": {"colorscheme": "cyanCarbon", "radius": 0.14}})
        view.addSurface(py3Dmol.VDW, {"opacity": 0.28, "colorscheme": "cyanCarbon"})
        view.setBackgroundColor(bg)
        view.zoomTo()
        return view._make_html()
    except Exception as e:
        return f"<body style='background:#020817;color:#64748b;padding:2rem'><p>3D error: {e}</p></body>"


def _human_body_svg(concs: dict, dark: bool = True) -> str:
    """SVG anatomical body diagram with concentration-based coloring."""
    bg = "#020817" if dark else "#f8fafc"
    outline = "#1e3a5f" if dark else "#cbd5e1"
    txt = "#f0f9ff" if dark else "#0f172a"
    sub = "#94a3b8" if dark else "#64748b"

    plasma  = float(concs.get("plasma", 0))
    liver   = float(concs.get("liver", 0))
    brain   = float(concs.get("brain", 0))
    kidney  = float(concs.get("kidney", 0))
    muscle  = float(concs.get("muscle", 0))

    mx = max(plasma, liver, brain, kidney, muscle, 1.0)

    _pal = ("#1e3a5f", "#0e4e7a", "#0a6e9e", "#0ea5e9", "#38bdf8", "#7dd3fc")

    def c(val): return _clr(val, 0, mx, _pal)

    bc = c(brain);  lc = c(liver);  pc = c(plasma)
    kc = c(kidney); mc = c(muscle)

    def fmt(v): return f"{v:.0f}" if v >= 10 else f"{v:.1f}"

    return f"""<svg viewBox="0 0 300 560" xmlns="http://www.w3.org/2000/svg">
<rect width="300" height="560" fill="{bg}"/>

<!-- Head silhouette -->
<ellipse cx="150" cy="55" rx="50" ry="48" fill="{outline}" opacity="0.25"/>

<!-- Brain compartment -->
<ellipse cx="150" cy="52" rx="42" ry="40" fill="{bc}" opacity="0.9" stroke="{outline}" stroke-width="1.5"/>
<text x="150" y="47" text-anchor="middle" font-family="sans-serif" font-size="11" fill="{txt}" font-weight="700">Brain</text>
<text x="150" y="62" text-anchor="middle" font-family="sans-serif" font-size="9" fill="{sub}">{fmt(brain)} ng/mL</text>

<!-- Neck -->
<rect x="138" y="95" width="24" height="20" rx="6" fill="{outline}" opacity="0.4"/>

<!-- Shoulders + torso outline -->
<path d="M 90 115 Q 70 118 62 140 L 50 290 Q 52 298 60 298 L 65 298"
      fill="none" stroke="{outline}" stroke-width="1.5" opacity="0.35"/>
<path d="M 210 115 Q 230 118 238 140 L 250 290 Q 248 298 240 298 L 235 298"
      fill="none" stroke="{outline}" stroke-width="1.5" opacity="0.35"/>
<path d="M 90 115 L 210 115 L 220 310 L 80 310 Z"
      fill="{outline}" opacity="0.08" stroke="{outline}" stroke-width="1.5"/>

<!-- Left Lung -->
<ellipse cx="122" cy="178" rx="26" ry="38" fill="#0f4c75" opacity="0.7" stroke="{outline}" stroke-width="1"/>
<text x="122" y="175" text-anchor="middle" font-family="sans-serif" font-size="8" fill="{txt}">L. Lung</text>

<!-- Right Lung -->
<ellipse cx="178" cy="178" rx="26" ry="38" fill="#0f4c75" opacity="0.7" stroke="{outline}" stroke-width="1"/>
<text x="178" y="175" text-anchor="middle" font-family="sans-serif" font-size="8" fill="{txt}">R. Lung</text>

<!-- Heart (plasma proxy) -->
<ellipse cx="150" cy="162" rx="20" ry="22" fill="{pc}" opacity="0.9" stroke="{outline}" stroke-width="1.5"/>
<text x="150" y="158" text-anchor="middle" font-family="sans-serif" font-size="9" fill="{txt}" font-weight="600">Heart</text>
<text x="150" y="171" text-anchor="middle" font-family="sans-serif" font-size="8" fill="{sub}">{fmt(plasma)}</text>

<!-- Liver -->
<ellipse cx="138" cy="237" rx="37" ry="26" fill="{lc}" opacity="0.9" stroke="{outline}" stroke-width="1.5"/>
<text x="138" y="234" text-anchor="middle" font-family="sans-serif" font-size="11" fill="{txt}" font-weight="700">Liver</text>
<text x="138" y="248" text-anchor="middle" font-family="sans-serif" font-size="8" fill="{sub}">{fmt(liver)} ng/mL</text>

<!-- Left Kidney -->
<ellipse cx="112" cy="276" rx="18" ry="24" fill="{kc}" opacity="0.9" stroke="{outline}" stroke-width="1.5"/>
<text x="112" y="273" text-anchor="middle" font-family="sans-serif" font-size="7.5" fill="{txt}">L. Kidney</text>
<text x="112" y="285" text-anchor="middle" font-family="sans-serif" font-size="7" fill="{sub}">{fmt(kidney)}</text>

<!-- Right Kidney -->
<ellipse cx="188" cy="276" rx="18" ry="24" fill="{kc}" opacity="0.9" stroke="{outline}" stroke-width="1.5"/>
<text x="188" y="273" text-anchor="middle" font-family="sans-serif" font-size="7.5" fill="{txt}">R. Kidney</text>
<text x="188" y="285" text-anchor="middle" font-family="sans-serif" font-size="7" fill="{sub}">{fmt(kidney)}</text>

<!-- Left Arm -->
<rect x="55" y="118" width="32" height="148" rx="16" fill="{mc}" opacity="0.75" stroke="{outline}" stroke-width="1.2"/>
<text x="71" y="208" text-anchor="middle" font-family="sans-serif" font-size="8" fill="{txt}"
      transform="rotate(-90 71 208)">Muscle</text>

<!-- Right Arm -->
<rect x="213" y="118" width="32" height="148" rx="16" fill="{mc}" opacity="0.75" stroke="{outline}" stroke-width="1.2"/>
<text x="229" y="208" text-anchor="middle" font-family="sans-serif" font-size="8" fill="{txt}"
      transform="rotate(90 229 208)">Muscle</text>

<!-- Pelvis band -->
<rect x="85" y="308" width="130" height="24" rx="8" fill="{outline}" opacity="0.15" stroke="{outline}" stroke-width="1"/>

<!-- Left Leg -->
<rect x="88" y="330" width="50" height="200" rx="24" fill="{mc}" opacity="0.75" stroke="{outline}" stroke-width="1.2"/>

<!-- Right Leg -->
<rect x="162" y="330" width="50" height="200" rx="24" fill="{mc}" opacity="0.75" stroke="{outline}" stroke-width="1.2"/>

<!-- Plasma legend label -->
<text x="220" y="162" font-family="sans-serif" font-size="9" fill="{pc}" font-weight="700">Plasma</text>
<text x="220" y="174" font-family="sans-serif" font-size="8" fill="{sub}">{fmt(plasma)} ng/mL</text>

<!-- Muscle concentration label -->
<text x="150" y="435" text-anchor="middle" font-family="sans-serif" font-size="8" fill="{sub}">{fmt(muscle)} ng/mL</text>

<!-- Blood vessel lines (animated) -->
<line x1="150" y1="184" x2="150" y2="213" stroke="{pc}" stroke-width="1.5" stroke-dasharray="4 2" opacity="0.7"/>
<line x1="138" y1="213" x2="138" y2="213" stroke="{pc}" stroke-width="1" opacity="0.5"/>
<line x1="150" y1="184" x2="112" y2="252" stroke="{pc}" stroke-width="1" stroke-dasharray="3 2" opacity="0.45"/>
<line x1="150" y1="184" x2="188" y2="252" stroke="{pc}" stroke-width="1" stroke-dasharray="3 2" opacity="0.45"/>
</svg>"""


# ── App init ───────────────────────────────────────────────────────────────────

app = dash.Dash(
    __name__,
    suppress_callback_exceptions=True,
    assets_folder="assets",
    title="NeuroRepurpose · CipherQ",
    meta_tags=[
        {"name": "viewport", "content": "width=device-width, initial-scale=1"},
        {"charset": "UTF-8"},
    ],
)
server = app.server


# ── SVG icons ─────────────────────────────────────────────────────────────────

def _icon(name: str) -> html.Span:
    paths = {
        "dashboard": "M3 13h8V3H3v10zm0 8h8v-6H3v6zm10 0h8V11h-8v10zm0-18v6h8V3h-8z",
        "search":    "M21 21l-6-6m2-5a7 7 0 11-14 0 7 7 0 0114 0z",
        "flask":     "M9 3v11l-5 7h16l-5-7V3M9 3h6",
        "graph":     "M4 6h4M4 12h7M4 18h4M16 6l4 4-4 4M20 10H9",
        "database":  "M4 7c0-1.1 3.6-2 8-2s8 .9 8 2m-16 0v10c0 1.1 3.6 2 8 2s8-.9 8-2V7m-16 5c0 1.1 3.6 2 8 2s8-.9 8-2",
        "sun":       "M12 3v1m0 16v1m9-9h-1M4 12H3m15.36-6.36-.71.71M6.34 17.66l-.7.7M17.66 17.66l.7.7M6.34 6.34l-.7-.7M12 8a4 4 0 000 8 4 4 0 000-8z",
        "moon":      "M21 12.79A9 9 0 1 1 11.21 3 7 7 0 0 0 21 12.79z",
        "arrow-left":"M19 12H5m7-7-7 7 7 7",
        "download":  "M4 16v1a3 3 0 003 3h10a3 3 0 003-3v-1m-4-4l-4 4m0 0l-4-4m4 4V4",
        "external":  "M18 13v6a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2V8a2 2 0 0 1 2-2h6m4-3h6v6m-11 5L21 3",
    }
    d = paths.get(name, "")
    return html.Span(
        html.Svg(
            html.Path(d=d, strokeLinecap="round", strokeLinejoin="round"),
            xmlns="http://www.w3.org/2000/svg",
            viewBox="0 0 24 24",
            fill="none",
            stroke="currentColor",
            strokeWidth="2",
            style={"width": "16px", "height": "16px"},
        ),
        style={"display": "inline-flex", "alignItems": "center"},
    )


# ── Sidebar ────────────────────────────────────────────────────────────────────

def _sidebar(active: str, dark: bool) -> html.Div:
    nav = [
        ("dashboard", "Dashboard",      "/",         "dashboard"),
        ("discover",  "Discover Drugs",  "/discover", "search"),
        ("graph",     "Knowledge Graph", "/graph",    "graph"),
        ("database",  "Data Explorer",   "/database", "database"),
    ]
    items = []
    for key, label, href, icon in nav:
        cls = "nav-item" + (" active" if active == key else "")
        items.append(
            html.A(
                [_icon(icon), label],
                href=href, className=cls, id=f"nav-{key}",
            )
        )

    theme_icon = "sun" if dark else "moon"
    theme_label = "Light Mode" if dark else "Dark Mode"

    return html.Div(
        id="sidebar",
        className="dark" if dark else "light",
        children=[
            # Logo
            html.Div(
                [
                    html.Div([
                        html.Span("Neuro", style={"color": "var(--cyan)"}),
                        "Repurpose",
                    ], className="sidebar-logo-wordmark"),
                    html.Div("CipherQ Intelligence Platform", className="sidebar-logo-sub"),
                ],
                className="sidebar-logo",
            ),
            # Nav
            html.Nav(items, className="sidebar-nav"),
            # Footer
            html.Div(
                [
                    html.Button(
                        [_icon(theme_icon), theme_label],
                        id="theme-btn",
                        className="theme-toggle",
                    ),
                    html.Div(
                        [
                            html.Div(className=f"db-dot {'online' if DB_OK else 'offline'}"),
                            f"DB {'connected' if DB_OK else 'offline'}",
                        ],
                        className="db-status",
                    ),
                ],
                className="sidebar-footer",
            ),
        ],
    )


# ── Page: Dashboard ────────────────────────────────────────────────────────────

def _page_dashboard(dark: bool) -> html.Div:
    stats = get_stats() if DB_OK else {}

    stat_cards = []
    stat_defs = [
        ("compounds",      "Compounds",       "#0ea5e9"),
        ("targets",        "Targets",         "#10b981"),
        ("indications",    "Indications",     "#8b5cf6"),
        ("mesh_diseases",  "MeSH Diseases",   "#f59e0b"),
    ]
    for key, label, color in stat_defs:
        val = stats.get(key, 0)
        stat_cards.append(
            html.Div(
                [
                    html.Div(f"{val:,}", className="stat-value"),
                    html.Div(label, className="stat-label"),
                ],
                className="stat-card",
                style={"--accent": color},
            )
        )

    return html.Div([
        html.Div([
            html.Div([
                html.Span("Drug Repurposing ", className="hero-title"),
                html.Span("Intelligence", style={"color": "var(--cyan)"}),
            ], style={"fontSize": "2rem", "fontWeight": "900", "letterSpacing": "-0.04em"}),
            html.P(
                "Search any neurological or psychiatric condition to discover repurposing candidates "
                "with multi-signal evidence: indication data, target overlap, in-vitro activity, and knowledge graph paths.",
                className="hero-desc",
            ),
            html.Div(
                [
                    dcc.Input(id="hero-search", type="text", placeholder="e.g. Alzheimer Disease, Parkinson, MS…",
                              debounce=False, className="search-wrap",
                              style={"width": "440px", "padding": "0.75rem 1rem",
                                     "background": "var(--raised)", "border": "1px solid var(--border)",
                                     "borderRadius": "8px", "color": "var(--text)", "fontSize": "0.9rem",
                                     "outline": "none"}),
                    html.Button([_icon("search"), " Search"],
                                id="hero-search-btn", n_clicks=0, className="btn btn-primary"),
                ],
                style={"display": "flex", "gap": "0.75rem", "marginTop": "1.5rem", "alignItems": "center"},
            ),
        ], className="hero"),

        # Stats grid
        html.Div(stat_cards, className="grid-4", style={"marginBottom": "1.5rem"}),

        # Quick links
        html.Div([
            html.Div("Quick Start", style={"fontWeight": "700", "fontSize": "0.85rem",
                                            "color": "var(--muted)", "textTransform": "uppercase",
                                            "letterSpacing": "0.07em", "marginBottom": "0.75rem"}),
            html.Div(
                [
                    html.A("Alzheimer Disease", href="/discover?q=alzheimer+disease",
                           style={"padding": "0.4rem 0.8rem", "background": "var(--raised)",
                                  "border": "1px solid var(--border)", "borderRadius": "6px",
                                  "color": "var(--cyan)", "textDecoration": "none", "fontSize": "0.8rem"}),
                    html.A("Parkinson Disease", href="/discover?q=parkinson+disease",
                           style={"padding": "0.4rem 0.8rem", "background": "var(--raised)",
                                  "border": "1px solid var(--border)", "borderRadius": "6px",
                                  "color": "var(--cyan)", "textDecoration": "none", "fontSize": "0.8rem"}),
                    html.A("Multiple Sclerosis", href="/discover?q=multiple+sclerosis",
                           style={"padding": "0.4rem 0.8rem", "background": "var(--raised)",
                                  "border": "1px solid var(--border)", "borderRadius": "6px",
                                  "color": "var(--cyan)", "textDecoration": "none", "fontSize": "0.8rem"}),
                    html.A("ALS", href="/discover?q=als",
                           style={"padding": "0.4rem 0.8rem", "background": "var(--raised)",
                                  "border": "1px solid var(--border)", "borderRadius": "6px",
                                  "color": "var(--cyan)", "textDecoration": "none", "fontSize": "0.8rem"}),
                ],
                style={"display": "flex", "flexWrap": "wrap", "gap": "0.5rem"},
            ),
        ], className="card"),
    ])


# ── Page: Discover ─────────────────────────────────────────────────────────────

def _page_discover(results: list, disease_name: str, dark: bool) -> html.Div:

    def compound_card(c: dict, idx: int) -> html.Div:
        score = c.get("score", 0)
        sc = _score_color(score)
        phase = int(c.get("max_phase") or 0)
        breakdown = c.get("score_breakdown", {})
        mech = (c.get("mechanisms") or "")[:60]
        targ = (c.get("targets") or "")[:50]

        score_mini_bars = []
        for bk, bl, bcol in [("indication_score", "Indication", "#0ea5e9"),
                               ("target_score",    "Target",     "#8b5cf6"),
                               ("activity_score",  "Activity",   "#10b981"),
                               ("network_score",   "Network",    "#f59e0b")]:
            bv = breakdown.get(bk, 0) * 100
            score_mini_bars.append(html.Div([
                html.Div(bl, className="score-row-label"),
                html.Div(html.Div(
                    style={"width": f"{bv:.0f}%", "height": "3px",
                           "background": bcol, "borderRadius": "2px"}
                ), className="score-row-track"),
            ], className="score-row"))

        return html.Div(
            id={"type": "compound-card", "index": idx},
            n_clicks=0,
            className="compound-card",
            children=[
                html.Div([
                    html.Div([
                        html.Div(c.get("name", "Unknown"), className="compound-name"),
                        html.Div(c.get("chembl_id", ""), className="compound-chembl"),
                    ]),
                    html.Div([
                        html.Div(f"{score:.2f}", className="score-value", style={"color": sc}),
                        html.Div("Score", className="score-label"),
                    ], className="score-badge"),
                ], className="compound-card-header"),

                html.Div(score, id={"type": "score-bar-data", "index": idx},
                         style={"display": "none"}),

                html.Div([
                    html.Span(_phase_label(phase), className=f"phase-pill {_phase_class(phase)}"),
                    *(([html.Span(mech, className="chip violet")] if mech else [])),
                    *(([html.Span(targ, className="chip cyan")] if targ else [])),
                ], className="compound-meta"),

                html.Div(score_mini_bars, className="score-grid", style={"marginTop": "0.75rem"}),
            ],
        )

    content_area = html.Div(id="results-content", children=[
        html.Div([
            dcc.Input(id="discover-search", type="text",
                      value=disease_name or "",
                      placeholder="e.g. Alzheimer Disease, Parkinson, MS…",
                      debounce=False,
                      style={"flex": "1", "padding": "0.65rem 1rem",
                             "background": "var(--card)", "border": "1px solid var(--border)",
                             "borderRadius": "8px", "color": "var(--text)", "fontSize": "0.9rem",
                             "outline": "none"}),
            html.Button([_icon("search"), " Search"], id="discover-btn",
                        n_clicks=0, className="btn btn-primary"),
        ], style={"display": "flex", "gap": "0.75rem", "marginBottom": "1.5rem", "alignItems": "center"}),

        html.Div(id="search-status"),

        html.Div(
            [compound_card(c, i) for i, c in enumerate(results)] if results
            else [html.Div([
                html.P("Search for a disease to discover repurposing candidates.",
                       style={"color": "var(--muted)", "textAlign": "center", "padding": "3rem 0"}),
            ])],
            className="grid-auto",
            id="compound-grid",
        ),
    ])

    return html.Div([
        html.Div([
            html.Div(["Discover · ",
                      html.Span(disease_name or "Search a disease", style={"color": "var(--cyan)"})],
                     className="page-title"),
            html.Div(
                f"Found {len(results)} candidates" if results else "Enter a disease to begin",
                className="page-subtitle"),
        ], className="page-header"),
        content_area,
    ])


# ── Analysis helpers ───────────────────────────────────────────────────────────

def _tab_properties(c: dict, dark: bool) -> html.Div:
    smiles = c.get("smiles", "")
    rows = [
        ("ChEMBL ID",      c.get("chembl_id", "—"), False),
        ("Molecular Weight", f"{c.get('mw', '—')} g/mol" if c.get('mw') else "—", False),
        ("LogP (alogp)",   f"{c.get('alogp', '—')}", c.get('alogp', 999) < 5),
        ("PSA",            f"{c.get('psa', '—')} Å²" if c.get('psa') else "—", c.get('psa', 999) < 140),
        ("HBA",            str(c.get('hba', '—')), False),
        ("HBD",            str(c.get('hbd', '—')), False),
        ("Ro5 Violations", str(c.get('ro5_violations', '—')), int(c.get('ro5_violations') or 0) == 0),
        ("Max Phase",      _phase_label(c.get('max_phase')), False),
        ("SMILES",         (smiles[:60] + "…") if smiles and len(smiles) > 60 else (smiles or "—"), False),
    ]
    kv_rows = [
        html.Div([
            html.Span(k, className="kv-key"),
            html.Span(v, className=f"kv-val {'kv-pass' if ok else ''}"),
        ], className="kv-row")
        for k, v, ok in rows
    ]

    # RDKit descriptors
    extra = []
    if RDKIT_OK and smiles:
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                qed = round(Descriptors.qed(mol), 3)
                tpsa = round(Descriptors.TPSA(mol), 1)
                rot = Descriptors.NumRotatableBonds(mol)
                aromatic = Descriptors.NumAromaticRings(mol)
                extra = [
                    html.Div([html.Span("QED Score", className="kv-key"),
                               html.Span(str(qed), className=f"kv-val {'kv-pass' if qed > 0.5 else ''}"),
                               ], className="kv-row"),
                    html.Div([html.Span("TPSA (RDKit)", className="kv-key"),
                               html.Span(f"{tpsa} Å²", className="kv-val")], className="kv-row"),
                    html.Div([html.Span("Rotatable Bonds", className="kv-key"),
                               html.Span(str(rot), className="kv-val")], className="kv-row"),
                    html.Div([html.Span("Aromatic Rings", className="kv-key"),
                               html.Span(str(aromatic), className="kv-val")], className="kv-row"),
                ]
        except Exception:
            pass

    return html.Div([
        html.Div(kv_rows + extra, className="card"),
    ])


def _tab_bioactivity(cid: int, dark: bool) -> html.Div:
    targets = get_compound_targets(cid)
    activities = get_compound_activities(cid)

    tgt_section = html.Div()
    if targets:
        rows = [html.Tr([
            html.Td(t.get("gene_symbol") or "—"),
            html.Td(t.get("name", "")[:45]),
            html.Td(t.get("mechanism", "")[:40]),
            html.Td(t.get("action_type", "")[:20]),
            html.Td(str(t.get("confidence") or "—")),
        ]) for t in targets[:20]]
        tgt_section = html.Div([
            html.Div("Targets & Mechanisms", className="card-title mb-05"),
            html.Table([
                html.Thead(html.Tr([
                    html.Th("Gene"), html.Th("Target"), html.Th("Mechanism"),
                    html.Th("Action"), html.Th("Confidence"),
                ])),
                html.Tbody(rows),
            ], className="data-table"),
        ], className="card mb-1")

    act_fig = _empty_fig("No activity data.", dark)
    if activities:
        df = pd.DataFrame(activities)
        df = df.dropna(subset=["pchembl_value"])
        if not df.empty:
            df["label"] = df.get("gene_symbol", df.get("target_name", "")).fillna("").astype(str)
            df = df.sort_values("pchembl_value", ascending=False).head(20)
            act_fig = px.bar(
                df, x="pchembl_value", y="label", orientation="h",
                labels={"pchembl_value": "pChEMBL Value (-log₁₀ IC₅₀)", "label": "Target"},
                color="pchembl_value", color_continuous_scale=[[0, "#1e3a5f"], [0.5, "#0ea5e9"], [1, "#38bdf8"]],
            )
            act_fig.update_layout(**_plotly_theme(dark), height=max(250, len(df) * 28))
            act_fig.update_traces(marker_line_width=0)

    return html.Div([
        tgt_section,
        html.Div([
            html.Div("Binding Activity (pChEMBL)", className="card-title mb-05"),
            dcc.Graph(figure=act_fig, config={"displayModeBar": False}),
        ], className="card"),
    ])


def _tab_indications(cid: int, dark: bool) -> html.Div:
    inds = get_compound_indications(cid)
    if not inds:
        return html.Div(html.Div("No indication data.", style={"color": "var(--muted)"}), className="card")
    rows = [html.Tr([
        html.Td(i.get("disease", "—")),
        html.Td(i.get("mesh_id", "—"), style={"fontFamily": "monospace", "fontSize": "0.78rem"}),
        html.Td(_phase_label(i.get("max_phase"))),
        html.Td(i.get("source", "—")),
    ]) for i in inds[:30]]
    return html.Div([
        html.Table([
            html.Thead(html.Tr([html.Th("Disease"), html.Th("MeSH ID"), html.Th("Phase"), html.Th("Source")])),
            html.Tbody(rows),
        ], className="data-table"),
    ], className="card")


def _tab_clinical(drug_name: str, disease_name: str) -> html.Div:
    trials = _fetch_trials(drug_name, disease_name)
    papers = _fetch_papers(drug_name, disease_name)

    trial_rows = []
    for t in trials:
        status_color = {"COMPLETED": "var(--emerald)", "RECRUITING": "var(--cyan)",
                        "TERMINATED": "var(--rose)"}.get(t.get("status", ""), "var(--muted)")
        trial_rows.append(html.Tr([
            html.Td(html.A(t["nct_id"], href=t["url"], target="_blank",
                           style={"color": "var(--cyan)", "textDecoration": "none", "fontSize": "0.75rem",
                                  "fontFamily": "monospace"})),
            html.Td(t["title"], style={"fontSize": "0.78rem"}),
            html.Td(t.get("phase", "—"), style={"fontSize": "0.75rem"}),
            html.Td(t.get("status", "—"), style={"color": status_color, "fontSize": "0.75rem",
                                                   "fontWeight": "600"}),
        ]))

    paper_rows = []
    for p in papers:
        paper_rows.append(html.Tr([
            html.Td(html.A(p["pmid"], href=p["url"], target="_blank",
                           style={"color": "var(--cyan)", "textDecoration": "none",
                                  "fontSize": "0.75rem", "fontFamily": "monospace"})),
            html.Td(p["title"], style={"fontSize": "0.78rem"}),
            html.Td(p["journal"], style={"fontSize": "0.75rem", "color": "var(--muted)"}),
            html.Td(p["year"],    style={"fontSize": "0.75rem", "color": "var(--muted)"}),
        ]))

    trials_section = html.Div([
        html.Div([
            html.Div("Clinical Trials", className="card-title"),
            html.Span(f"{len(trials)} found", className="chip cyan"),
        ], className="card-header"),
        html.Table([
            html.Thead(html.Tr([html.Th("NCT ID"), html.Th("Title"), html.Th("Phase"), html.Th("Status")])),
            html.Tbody(trial_rows if trial_rows else [html.Tr([html.Td("No trials found.", colSpan=4, style={"color": "var(--muted)"})])]),
        ], className="data-table"),
    ], className="card mb-1")

    papers_section = html.Div([
        html.Div([
            html.Div("PubMed Literature", className="card-title"),
            html.Span(f"{len(papers)} found", className="chip green"),
        ], className="card-header"),
        html.Table([
            html.Thead(html.Tr([html.Th("PMID"), html.Th("Title"), html.Th("Journal"), html.Th("Year")])),
            html.Tbody(paper_rows if paper_rows else [html.Tr([html.Td("No papers found.", colSpan=4, style={"color": "var(--muted)"})])]),
        ], className="data-table"),
    ], className="card")

    return html.Div([trials_section, papers_section])


def _tab_pbpk(c: dict, dark: bool) -> html.Div:
    name = c.get("name", "Drug")
    mw   = float(c.get("mw") or 350)
    logp = float(c.get("alogp") or 2.5)

    return html.Div([
        html.Div([
            # Controls
            html.Div([
                html.Div([
                    html.Label("Dose (mg)"),
                    dcc.Input(id="pbpk-dose", type="number", value=100, min=1, max=2000, step=1,
                              style={"width": "100%"}),
                ]),
                html.Div([
                    html.Label("Route"),
                    dcc.Dropdown(id="pbpk-route",
                                 options=[{"label": "Oral", "value": "oral"},
                                          {"label": "IV", "value": "iv"}],
                                 value="oral", clearable=False,
                                 style={"background": "var(--card)", "color": "var(--text)"}),
                ]),
                html.Div([
                    html.Label("Duration (h)"),
                    dcc.Input(id="pbpk-hours", type="number", value=24, min=1, max=168, step=1,
                              style={"width": "100%"}),
                ]),
                html.Button("Run Simulation", id="pbpk-run-btn", n_clicks=0,
                             className="btn btn-primary", style={"marginTop": "1.4rem", "width": "100%"}),
            ], style={"display": "grid", "gridTemplateColumns": "1fr 1fr 1fr 1fr", "gap": "1rem",
                      "alignItems": "end", "marginBottom": "1.25rem"}),
        ], className="card mb-1"),

        html.Div(id="pbpk-result-area"),

        # Hidden stores
        dcc.Store(id="pbpk-compound-data", data={"name": name, "mw": mw, "logp": logp}),
    ])


def _tab_docking(c: dict, dark: bool) -> html.Div:
    name = c.get("name", "Drug")
    smiles = c.get("smiles", "")

    target_options = [
        {"label": "BACE1 (Alzheimer)", "value": "BACE1"},
        {"label": "AChE (Alzheimer/Parkinson)", "value": "ACHE"},
        {"label": "LRRK2 (Parkinson)", "value": "LRRK2"},
        {"label": "MAO-B (Parkinson)", "value": "MAOB"},
        {"label": "PPARG (Metabolic/Neuro)", "value": "PPARG"},
        {"label": "GSK3B (Alzheimer/Bipolar)", "value": "GSK3B"},
        {"label": "CDK5 (Neurodegeneration)", "value": "CDK5"},
        {"label": "HDAC1 (Epigenetics)", "value": "HDAC1"},
    ]

    api_status = html.Div(
        "NVIDIA DiffDock API ready" if DOCK_OK else
        "NVIDIA API key not set — set NVIDIA_API_KEY to enable DiffDock. AutoDock Vina available as fallback.",
        className=f"alert {'alert-success' if DOCK_OK else 'alert-warning'}",
        style={"marginBottom": "1rem"},
    )

    return html.Div([
        api_status,
        html.Div([
            html.Div([
                html.Label("Target Protein"),
                dcc.Dropdown(id="dock-target", options=target_options,
                             value="BACE1", clearable=False,
                             style={"background": "var(--card)", "color": "var(--text)"}),
            ]),
            html.Div([
                html.Label("Ligand SMILES"),
                dcc.Input(id="dock-smiles", type="text", value=smiles,
                          placeholder="SMILES string",
                          style={"width": "100%"}),
            ], style={"gridColumn": "span 2"}),
            html.Button("Run Docking", id="dock-run-btn", n_clicks=0,
                         className="btn btn-primary", style={"marginTop": "1.4rem"}),
        ], style={"display": "grid", "gridTemplateColumns": "200px 1fr 140px", "gap": "1rem",
                  "alignItems": "end", "marginBottom": "1.25rem"}, className="card"),

        html.Div(id="dock-result-area"),
        dcc.Store(id="dock-compound-data", data={"name": name, "smiles": smiles}),
    ])


def _tab_quantum(c: dict, dark: bool) -> html.Div:
    name   = c.get("name", "Drug")
    smiles = c.get("smiles", "")

    if not QUANTUM_OK:
        return html.Div(html.Div("Quantum optimization not available.", className="alert alert-warning"), className="card")
    if not smiles:
        return html.Div(html.Div("No SMILES available for this compound.", className="alert alert-warning"), className="card")

    try:
        result = run_quantum_optimization(name, smiles)
        if not result or not result.get("success"):
            return html.Div(html.Div("Optimization failed.", className="alert alert-error"), className="card")

        rows = [
            ("Optimization Status", "Optimized" if result.get("optimized") else "Failed",
             result.get("optimized", False)),
            ("Force Field",          result.get("force_field", "—"), True),
            ("Energy (initial)",     f"{result.get('initial_energy', '—'):.3f} kcal/mol" if isinstance(result.get('initial_energy'), (int, float)) else "—", False),
            ("Energy (optimized)",   f"{result.get('optimized_energy', '—'):.3f} kcal/mol" if isinstance(result.get('optimized_energy'), (int, float)) else "—", False),
            ("Energy Reduction",     f"{result.get('energy_delta', 0):.3f} kcal/mol", True),
            ("HOMO-LUMO Gap",        f"{result.get('homo_lumo_gap', '—'):.2f} eV" if isinstance(result.get('homo_lumo_gap'), (int, float)) else "—", True),
            ("Dipole Moment",        f"{result.get('dipole_moment', '—'):.2f} D" if isinstance(result.get('dipole_moment'), (int, float)) else "—", False),
            ("Polarizability",       f"{result.get('polarizability', '—'):.1f} Å³" if isinstance(result.get('polarizability'), (int, float)) else "—", False),
            ("QED Score",            f"{result.get('qed', '—'):.3f}" if isinstance(result.get('qed'), (int, float)) else "—",
             float(result.get('qed') or 0) > 0.5),
        ]

        kv = [
            html.Div([
                html.Span(k, className="kv-key"),
                html.Span(v, className=f"kv-val {'kv-pass' if ok else ''}"),
            ], className="kv-row")
            for k, v, ok in rows
        ]

        # Conformer visualization note
        note = html.Div(
            "3D conformer optimized using UFF force field. Electronic properties estimated via Gasteiger charges.",
            className="text-muted small", style={"marginTop": "0.5rem"},
        )

        return html.Div([html.Div(kv + [note], className="card")])

    except Exception as e:
        return html.Div(html.Div(f"Optimization error: {e}", className="alert alert-error"), className="card")


def _tab_3d(c: dict, dark: bool) -> html.Div:
    smiles = c.get("smiles", "")
    if not smiles:
        return html.Div(html.Div("No SMILES for this compound.", className="alert alert-warning"), className="card")
    html_src = _generate_3d_html(smiles, dark)
    return html.Div([
        html.Div([
            html.Div("3D Molecular Structure", className="card-title"),
            html.Span("Stick + VDW Surface · Cyan carbon coloring", className="text-muted small"),
        ], className="card-header"),
        html.Div(
            html.Iframe(srcDoc=html_src, style={"width": "100%", "height": "460px", "border": "none"}),
            className="viewer-3d-wrap",
        ),
    ], className="card")


# ── Page: Analysis ─────────────────────────────────────────────────────────────

def _page_analysis(compound: dict, mesh_ids: list, disease_name: str, dark: bool) -> html.Div:
    if not compound:
        return html.Div([
            html.Div("No compound selected.", style={"color": "var(--muted)", "padding": "2rem"}),
            html.A("Back to Discover", href="/discover", className="btn btn-secondary"),
        ])

    name  = compound.get("name", "Unknown")
    cid   = compound.get("id")
    score = compound.get("score", 0)
    sc    = _score_color(score)
    breakdown = compound.get("score_breakdown", {})

    # Score breakdown chart
    bd_fig = go.Figure()
    if breakdown:
        labels = ["Indication", "Target", "Activity", "Network"]
        values = [breakdown.get("indication_score", 0),
                  breakdown.get("target_score", 0),
                  breakdown.get("activity_score", 0),
                  breakdown.get("network_score", 0)]
        colors = ["#0ea5e9", "#8b5cf6", "#10b981", "#f59e0b"]
        bd_fig = go.Figure(go.Bar(x=labels, y=[v * 100 for v in values],
                                   marker_color=colors, marker_line_width=0))
        bd_fig.update_layout(**_plotly_theme(dark), height=180,
                             yaxis_title="Score %", showlegend=False)
        bd_fig.update_layout(margin=dict(l=40, r=10, t=10, b=30))

    header = html.Div([
        html.A([_icon("arrow-left"), " Back"], href="/discover",
               className="btn btn-ghost", style={"marginBottom": "1rem"}),
        html.Div([
            html.Div([
                html.Div(name, className="drug-header-name"),
                html.Div(compound.get("chembl_id", ""), className="drug-header-id"),
                html.Div([
                    html.Span(_phase_label(compound.get("max_phase")),
                              className=f"phase-pill {_phase_class(compound.get('max_phase'))}"),
                    html.Span(disease_name or "", className="chip cyan"),
                ], style={"display": "flex", "gap": "0.5rem", "marginTop": "0.5rem", "flexWrap": "wrap"}),
            ], className="drug-header-info"),
            html.Div([
                html.Div(f"{score:.2f}", className="drug-score-big", style={"color": sc}),
                html.Div("Repurposing Score", className="drug-score-label"),
                html.Div(f"Phase bonus: +{breakdown.get('phase_bonus', 0):.3f}",
                          style={"fontSize": "0.68rem", "color": "var(--muted)", "marginTop": "2px"}),
            ], className="drug-header-score"),
        ], className="drug-header"),

        html.Div([
            html.Div("Score Breakdown", className="card-title",
                     style={"marginBottom": "0.5rem", "fontSize": "0.72rem"}),
            dcc.Graph(figure=bd_fig, config={"displayModeBar": False}),
        ], className="card", style={"marginBottom": "1rem"}),
    ])

    # Tabs
    tab_ids = ["properties", "bioactivity", "indications", "clinical",
               "pbpk", "docking", "quantum", "3d-structure"]
    tab_labels = ["Properties", "Bioactivity", "Indications", "Clinical Evidence",
                  "PBPK Simulation", "Docking", "Quantum/Optimization", "3D Structure"]

    tab_bar = html.Div([
        html.Button(label, id=f"tab-btn-{tid}", className=f"tab-btn {'active' if tid == 'properties' else ''}",
                    n_clicks=0)
        for label, tid in zip(tab_labels, tab_ids)
    ], className="tab-bar")

    tab_content = html.Div(id="analysis-tab-content",
                            children=_tab_properties(compound, dark))

    return html.Div([header, tab_bar, tab_content,
                     dcc.Store(id="analysis-compound", data=compound),
                     dcc.Store(id="analysis-mesh-ids", data=mesh_ids or []),
                     dcc.Store(id="analysis-disease-name", data=disease_name or ""),
                     dcc.Store(id="active-tab", data="properties")])


# ── Page: Knowledge Graph ──────────────────────────────────────────────────────

def _page_graph(mesh_ids: list, disease_name: str, dark: bool) -> html.Div:
    elements, legend = build_graph(mesh_ids=mesh_ids or None, max_compounds=25, max_targets=35)

    legend_items = [
        html.Div([
            html.Div(style={"width": "10px", "height": "10px", "borderRadius": "50%",
                            "background": color, "flexShrink": "0"}),
            html.Span(kind, style={"fontSize": "0.72rem", "color": "var(--muted)"}),
        ], style={"display": "flex", "gap": "0.4rem", "alignItems": "center"})
        for kind, color in legend.items()
    ]

    cyto_graph = cyto.Cytoscape(
        id="cyto-graph",
        elements=elements,
        layout={"name": "cose", "animate": False, "randomize": True,
                "nodeRepulsion": 8000, "idealEdgeLength": 100, "gravity": 0.4},
        style={"width": "100%", "height": "640px"},
        stylesheet=CYTO_STYLESHEET,
        responsive=True,
    ) if elements else html.Div(
        "No graph data — search a disease first on the Discover page.",
        style={"color": "var(--muted)", "textAlign": "center", "padding": "4rem"}
    )

    node_info = html.Div(id="graph-node-info",
                         style={"minHeight": "60px", "color": "var(--muted)", "fontSize": "0.8rem"})

    return html.Div([
        html.Div([
            html.Div(["Knowledge Graph",
                      html.Span(" · BioCypher", style={"color": "var(--cyan)", "fontSize": "0.9rem"})],
                     className="page-title"),
            html.Div(f"{len([e for e in elements if 'source' not in e.get('data', {})])} nodes · "
                     f"{len([e for e in elements if 'source' in e.get('data', {})])} edges"
                     if elements else "Search a disease to populate the graph",
                     className="page-subtitle"),
        ], className="page-header"),

        html.Div([
            html.Div([
                html.Div([
                    html.Span("Legend: ", style={"fontSize": "0.72rem", "color": "var(--muted)"}),
                    *legend_items,
                ], style={"display": "flex", "flexWrap": "wrap", "gap": "0.8rem", "alignItems": "center"}),
                html.Div([
                    html.Button("Cose Layout", id="layout-cose", n_clicks=0, className="btn btn-ghost",
                                style={"fontSize": "0.75rem", "padding": "0.3rem 0.7rem"}),
                    html.Button("Circle Layout", id="layout-circle", n_clicks=0, className="btn btn-ghost",
                                style={"fontSize": "0.75rem", "padding": "0.3rem 0.7rem"}),
                    html.Button("Grid Layout", id="layout-grid", n_clicks=0, className="btn btn-ghost",
                                style={"fontSize": "0.75rem", "padding": "0.3rem 0.7rem"}),
                ], style={"display": "flex", "gap": "0.4rem"}),
            ], style={"display": "flex", "justifyContent": "space-between", "alignItems": "center",
                      "marginBottom": "0.75rem"}),
            html.Div(cyto_graph, className="cyto-wrap"),
            node_info,
        ], className="card"),
    ])


# ── Page: Database ─────────────────────────────────────────────────────────────

def _page_database(dark: bool) -> html.Div:
    stats = get_stats() if DB_OK else {}
    stat_rows = [html.Tr([html.Td(k.replace("_", " ").title()), html.Td(f"{v:,}")])
                 for k, v in stats.items()]

    return html.Div([
        html.Div([
            html.Div("Data Explorer", className="page-title"),
            html.Div("Browse compounds, targets, diseases and statistics.", className="page-subtitle"),
        ], className="page-header"),

        html.Div([
            html.Div([
                dcc.Input(id="db-search-input", type="text",
                          placeholder="Search by compound name, mechanism, or target…",
                          debounce=True,
                          style={"flex": "1", "padding": "0.55rem 0.9rem",
                                 "background": "var(--card)", "border": "1px solid var(--border)",
                                 "borderRadius": "8px", "color": "var(--text)", "fontSize": "0.85rem",
                                 "outline": "none"}),
                html.Button("Search", id="db-search-btn", n_clicks=0,
                             className="btn btn-primary", style={"padding": "0.55rem 1rem"}),
            ], style={"display": "flex", "gap": "0.6rem", "marginBottom": "1rem"}),
            html.Div(id="db-results"),
        ], className="card mb-1"),

        html.Div([
            html.Div([
                html.Div("Database Statistics", className="card-title"),
            ], className="card-header"),
            html.Table([
                html.Thead(html.Tr([html.Th("Table"), html.Th("Records")])),
                html.Tbody(stat_rows if stat_rows else [html.Tr([html.Td("DB offline", colSpan=2)])]),
            ], className="data-table"),
        ], className="card"),
    ])


# ── App layout ─────────────────────────────────────────────────────────────────

app.layout = html.Div(
    id="app-root",
    className="dark",
    children=[
        dcc.Location(id="url", refresh=False),
        dcc.Store(id="theme-store", data="dark"),
        dcc.Store(id="compound-store", data=None),
        dcc.Store(id="results-store", data=None),
        dcc.Store(id="disease-store", data={"name": "", "mesh_ids": []}),

        # Sidebar placeholder (updated by theme callback)
        html.Div(id="sidebar-wrap"),

        # Main content
        html.Div(id="main-content"),
    ],
)


# ── Callbacks ──────────────────────────────────────────────────────────────────

@app.callback(
    Output("sidebar-wrap", "children"),
    Output("main-content", "children"),
    Input("url", "pathname"),
    Input("url", "search"),
    Input("theme-store", "data"),
    State("compound-store", "data"),
    State("results-store", "data"),
    State("disease-store", "data"),
)
def route(pathname, search, theme, compound, results, disease_store):
    dark = (theme == "dark")
    disease_name = (disease_store or {}).get("name", "") if disease_store else ""
    mesh_ids     = (disease_store or {}).get("mesh_ids", []) if disease_store else []

    # Parse query string for pre-filled search
    query_from_url = ""
    if search and "q=" in search:
        query_from_url = search.split("q=")[-1].replace("+", " ").replace("%20", " ")

    p = (pathname or "/").rstrip("/") or "/"

    if p == "/" or p == "":
        page_key = "dashboard"
        content = _page_dashboard(dark)
    elif p == "/discover":
        page_key = "discover"
        # If URL has a query, auto-run search
        if query_from_url and not results:
            resolved, expanded, compounds = _resolve_and_fetch(query_from_url)
            results = compounds
            disease_name = resolved[0]["heading"] if resolved else query_from_url
        content = _page_discover(results or [], disease_name, dark)
    elif p == "/analysis":
        page_key = "discover"  # keep nav on discover
        content = _page_analysis(compound, mesh_ids, disease_name, dark)
    elif p == "/graph":
        page_key = "graph"
        content = _page_graph(mesh_ids, disease_name, dark)
    elif p == "/database":
        page_key = "database"
        content = _page_database(dark)
    else:
        page_key = "dashboard"
        content = _page_dashboard(dark)

    sidebar = _sidebar(page_key, dark)
    return sidebar, content


@app.callback(
    Output("theme-store", "data"),
    Output("app-root", "className"),
    Input("theme-btn", "n_clicks"),
    State("theme-store", "data"),
    prevent_initial_call=True,
)
def toggle_theme(n, theme):
    new = "light" if theme == "dark" else "dark"
    return new, new


@app.callback(
    Output("results-store", "data"),
    Output("disease-store", "data"),
    Output("url", "pathname"),
    Output("url", "search"),
    Input("hero-search-btn", "n_clicks"),
    Input("discover-btn", "n_clicks"),
    State("hero-search", "value"),
    State("discover-search", "value"),
    prevent_initial_call=True,
)
def run_search(n1, n2, hero_q, disc_q):
    trigger = ctx.triggered_id
    query = hero_q if trigger == "hero-search-btn" else disc_q
    if not query or not query.strip():
        return no_update, no_update, no_update, no_update
    query = query.strip()
    resolved, expanded, compounds = _resolve_and_fetch(query)
    disease_name = resolved[0]["heading"] if resolved else query
    disease_data = {"name": disease_name, "mesh_ids": expanded}
    return compounds, disease_data, "/discover", ""


@app.callback(
    Output("compound-store", "data"),
    Output("url", "pathname"),
    Input({"type": "compound-card", "index": ALL}, "n_clicks"),
    State("results-store", "data"),
    prevent_initial_call=True,
)
def select_compound(n_clicks_list, results):
    if not any(n_clicks_list) or not results:
        return no_update, no_update
    triggered = ctx.triggered_id
    if triggered is None:
        return no_update, no_update
    idx = triggered["index"]
    if idx < len(results):
        return results[idx], "/analysis"
    return no_update, no_update


@app.callback(
    Output("analysis-tab-content", "children"),
    Output("active-tab", "data"),
    Input("tab-btn-properties", "n_clicks"),
    Input("tab-btn-bioactivity", "n_clicks"),
    Input("tab-btn-indications", "n_clicks"),
    Input("tab-btn-clinical", "n_clicks"),
    Input("tab-btn-pbpk", "n_clicks"),
    Input("tab-btn-docking", "n_clicks"),
    Input("tab-btn-quantum", "n_clicks"),
    Input("tab-btn-3d-structure", "n_clicks"),
    State("analysis-compound", "data"),
    State("analysis-mesh-ids", "data"),
    State("analysis-disease-name", "data"),
    State("theme-store", "data"),
    prevent_initial_call=True,
)
def switch_tab(p1, p2, p3, p4, p5, p6, p7, p8,
               compound, mesh_ids, disease_name, theme):
    dark = (theme == "dark")
    tid = ctx.triggered_id
    if not tid or not compound:
        return no_update, no_update

    cid = compound.get("id")

    if tid == "tab-btn-properties":     return _tab_properties(compound, dark), "properties"
    if tid == "tab-btn-bioactivity":    return _tab_bioactivity(cid, dark), "bioactivity"
    if tid == "tab-btn-indications":    return _tab_indications(cid, dark), "indications"
    if tid == "tab-btn-clinical":       return _tab_clinical(compound.get("name",""), disease_name), "clinical"
    if tid == "tab-btn-pbpk":           return _tab_pbpk(compound, dark), "pbpk"
    if tid == "tab-btn-docking":        return _tab_docking(compound, dark), "docking"
    if tid == "tab-btn-quantum":        return _tab_quantum(compound, dark), "quantum"
    if tid == "tab-btn-3d-structure":   return _tab_3d(compound, dark), "3d-structure"
    return no_update, no_update


@app.callback(
    Output("pbpk-result-area", "children"),
    Input("pbpk-run-btn", "n_clicks"),
    State("pbpk-compound-data", "data"),
    State("pbpk-dose", "value"),
    State("pbpk-route", "value"),
    State("pbpk-hours", "value"),
    State("theme-store", "data"),
    prevent_initial_call=True,
)
def run_pbpk(n, comp_data, dose, route, hours, theme):
    if not n or not PBPK_OK:
        return html.Div("PBPK module not available.", className="alert alert-warning")
    dark = (theme == "dark")
    name = (comp_data or {}).get("name", "Drug")
    mw   = float((comp_data or {}).get("mw", 350))
    logp = float((comp_data or {}).get("logp", 2.5))
    dose  = float(dose or 100)
    hours = float(hours or 24)

    try:
        sim = PBPKSimulator()
        res = sim.simulate_drug_exposure(
            drug_name=name, molecular_weight=mw, logp=logp,
            dose_mg=dose, route=route or "oral", duration_hours=hours,
        )
        if not res.get("success"):
            return html.Div(f"Simulation failed.", className="alert alert-error")

        t  = res["time_hours"]
        pl = res["plasma_concentration_ng_ml"]
        li = res["liver_concentration_ng_ml"]
        br = res["brain_concentration_ng_ml"]
        tg = res["target_tissue_concentration_ng_ml"]
        pk = res.get("pk_metrics", {})

        # Concentration-time chart
        fig = go.Figure()
        colors_map = {"Plasma": "#0ea5e9", "Liver": "#10b981",
                      "Brain": "#8b5cf6", "Target Tissue": "#f59e0b"}
        for label, data in [("Plasma", pl), ("Liver", li), ("Brain", br), ("Target Tissue", tg)]:
            fig.add_trace(go.Scatter(x=t, y=data, name=label, mode="lines",
                                      line=dict(color=colors_map[label], width=2.5)))
        fig.update_layout(**_plotly_theme(dark), height=320,
                          xaxis_title="Time (hours)", yaxis_title="Concentration (ng/mL)",
                          legend=dict(orientation="h", yanchor="bottom", y=1.01, xanchor="left", x=0))

        # Snapshot for body diagram (at peak)
        peak_idx = int(np.argmax(pl))
        body_concs = {"plasma": pl[peak_idx], "liver": li[peak_idx],
                      "brain": br[peak_idx], "kidney": tg[peak_idx] * 0.7,
                      "muscle": tg[peak_idx] * 0.4}

        body_svg = _human_body_svg(body_concs, dark)

        # PK metrics
        pk_rows = [
            ("Cmax (plasma)",    f"{pk.get('cmax_plasma', '—'):.2f} ng/mL" if isinstance(pk.get('cmax_plasma'), (int,float)) else "—"),
            ("Tmax",             f"{pk.get('tmax', '—'):.1f} h" if isinstance(pk.get('tmax'), (int,float)) else "—"),
            ("AUC",              f"{pk.get('auc', '—'):.1f} ng·h/mL" if isinstance(pk.get('auc'), (int,float)) else "—"),
            ("Half-life (t½)",   f"{pk.get('half_life', '—'):.1f} h" if isinstance(pk.get('half_life'), (int,float)) else "—"),
            ("Brain/Plasma Ratio", f"{pk.get('brain_plasma_ratio', '—'):.3f}" if isinstance(pk.get('brain_plasma_ratio'), (int,float)) else "—"),
            ("Clearance",        f"{pk.get('clearance', '—'):.3f} L/h/kg" if isinstance(pk.get('clearance'), (int,float)) else "—"),
        ]

        return html.Div([
            html.Div([
                # Human body + PK metrics
                html.Div([
                    html.Div("Body Distribution (at Cmax)", className="card-title mb-05"),
                    html.Div(html.Div(
                        dangerouslySetInnerHTML={"__html": body_svg},
                        style={"display": "flex", "justifyContent": "center"},
                    ), className="pbpk-body-wrap"),
                ], className="card"),
                html.Div([
                    html.Div("PK Metrics", className="card-title mb-05"),
                    *[html.Div([
                        html.Span(k, className="kv-key"),
                        html.Span(v, className="kv-val"),
                    ], className="kv-row") for k, v in pk_rows],
                ], className="card"),
            ], style={"display": "grid", "gridTemplateColumns": "240px 1fr", "gap": "1rem",
                      "marginBottom": "1rem"}),

            html.Div([
                html.Div("Concentration–Time Profile", className="card-title mb-05"),
                dcc.Graph(figure=fig, config={"displayModeBar": True}),
            ], className="card"),
        ])

    except Exception as e:
        return html.Div(f"Simulation error: {e}", className="alert alert-error")


@app.callback(
    Output("dock-result-area", "children"),
    Input("dock-run-btn", "n_clicks"),
    State("dock-compound-data", "data"),
    State("dock-target", "value"),
    State("dock-smiles", "value"),
    State("theme-store", "data"),
    prevent_initial_call=True,
)
def run_docking(n, comp_data, target_name, smiles, theme):
    if not n:
        return no_update
    dark = (theme == "dark")
    name = (comp_data or {}).get("name", "Drug")

    if not _dock_svc:
        return html.Div("Docking service not available.", className="alert alert-error")

    try:
        result = _dock_svc.perform_docking(
            drug_name=name,
            target_name=target_name,
            ligand_smiles=smiles or (comp_data or {}).get("smiles"),
        )

        if not result.get("success"):
            method = result.get("docking_method", "Docking")
            return html.Div([
                html.Div(f"Docking failed: {result.get('error', 'Unknown error')}",
                          className="alert alert-error"),
                html.Div(f"Method attempted: {method}", className="text-muted small",
                          style={"marginTop": "0.4rem"}),
            ])

        poses = result.get("poses", [])
        confs = result.get("confidence_scores", [])
        affs  = result.get("binding_affinities", [])
        method = result.get("docking_method", "Docking")

        pose_cards = []
        for i, (conf, aff) in enumerate(zip(confs[:8], affs[:8] if affs else [None]*8)):
            aff_str = f"{aff:.2f} kcal/mol" if aff is not None else "—"
            pose_cards.append(html.Div([
                html.Div(f"Pose {i+1}", className="dock-pose-rank"),
                html.Div([
                    html.Span(f"Confidence: {conf:.3f}" if conf is not None else "—",
                               style={"color": "var(--cyan)", "fontSize": "0.82rem"}),
                    html.Span(" · ", style={"color": "var(--muted)"}),
                    html.Span(aff_str, className="dock-pose-score"),
                ]),
            ], className="dock-pose-card"))

        # 3D viewer for best pose SDF
        viewer_section = html.Div()
        if poses and RDKIT_OK and PY3DMOL_OK:
            try:
                best_pose_sdf = poses[0] if isinstance(poses[0], str) else ""
                if best_pose_sdf:
                    bg = "#020817" if dark else "#f8fafc"
                    view = py3Dmol.view(width=760, height=420, linked=False)
                    view.addModel(best_pose_sdf, "sdf")
                    view.setStyle({}, {"stick": {"colorscheme": "cyanCarbon", "radius": 0.14}})
                    view.setBackgroundColor(bg)
                    view.zoomTo()
                    viewer_section = html.Div([
                        html.Div("Best Docking Pose", className="card-title mb-05"),
                        html.Iframe(srcDoc=view._make_html(),
                                    style={"width": "100%", "height": "420px", "border": "none"}),
                    ], className="card mt-1")
            except Exception:
                pass

        return html.Div([
            html.Div([
                html.Span(f"Method: {method}",
                           className=f"chip {'green' if not result.get('fallback') else 'amber'}"),
                html.Span(f"{len(poses)} poses generated", className="chip cyan"),
                html.Span(f"Target: {target_name}", className="chip violet"),
            ], style={"display": "flex", "flexWrap": "wrap", "gap": "0.4rem", "marginBottom": "1rem"}),
            html.Div(pose_cards, style={"display": "flex", "flexDirection": "column", "gap": "0.4rem"},
                     className="card mb-1"),
            viewer_section,
        ])

    except Exception as e:
        return html.Div(f"Docking error: {e}", className="alert alert-error")


@app.callback(
    Output("db-results", "children"),
    Input("db-search-btn", "n_clicks"),
    Input("db-search-input", "value"),
    State("theme-store", "data"),
    prevent_initial_call=True,
)
def db_search_cb(n, q, theme):
    dark = (theme == "dark")
    if not q or not q.strip():
        return html.Div()
    compounds = db_search(q.strip(), limit=40)
    if not compounds:
        return html.Div("No results.", className="text-muted")
    df = pd.DataFrame(compounds)
    keep = [c for c in ["chembl_id", "name", "max_phase", "mw", "alogp", "psa", "mechanisms", "targets"]
            if c in df.columns]
    return html.Div([
        html.Div(f"{len(compounds)} results", className="text-muted small mb-05"),
        html.Div(
            html.Table(
                [html.Thead(html.Tr([html.Th(h) for h in keep]))] +
                [html.Tbody([html.Tr([html.Td(str(row.get(k, "—"))[:60]) for k in keep])
                              for row in compounds[:30]])],
                className="data-table",
            ),
            style={"overflowX": "auto"},
        ),
    ])


@app.callback(
    Output("cyto-graph", "layout"),
    Input("layout-cose", "n_clicks"),
    Input("layout-circle", "n_clicks"),
    Input("layout-grid", "n_clicks"),
    prevent_initial_call=True,
)
def change_layout(c1, c2, c3):
    tid = ctx.triggered_id
    names = {"layout-cose": "cose", "layout-circle": "circle", "layout-grid": "grid"}
    return {"name": names.get(tid, "cose"), "animate": True}


@app.callback(
    Output("graph-node-info", "children"),
    Input("cyto-graph", "tapNodeData"),
    prevent_initial_call=True,
)
def show_node_info(data):
    if not data:
        return no_update
    kind = data.get("kind", "")
    label = data.get("full_label", data.get("label", ""))
    parts = [html.Strong(label), f" ({kind})"]
    if data.get("chembl_id"):
        parts += [" · ChEMBL: ", html.Code(data["chembl_id"])]
    if data.get("max_phase"):
        parts += [" · ", _phase_label(data["max_phase"])]
    if data.get("gene_symbol"):
        parts += [" · Gene: ", html.Code(data["gene_symbol"])]
    if data.get("mesh_id"):
        parts += [" · MeSH: ", html.Code(data["mesh_id"])]
    return html.Div(parts, style={"padding": "0.5rem 0", "fontSize": "0.82rem"})


# ── Entry point ────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    app.run(debug=False, host="127.0.0.1", port=8050)
