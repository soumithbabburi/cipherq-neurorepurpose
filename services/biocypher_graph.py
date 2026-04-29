"""
BioCypher Knowledge Graph Adapter
Pulls compound-target-disease data from PostgreSQL and formats it
as nodes + edges compatible with dash-cytoscape and BioCypher conventions.
"""
import logging
from typing import Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

try:
    from services.neuro_db_service import _query
    _DB_OK = True
except Exception:
    _DB_OK = False
    def _query(*a, **k): return []

try:
    import biocypher
    _BC_OK = True
except ImportError:
    _BC_OK = False

# ── Node color palette ─────────────────────────────────────────────────────────
_COLORS = {
    "Compound":  {"bg": "#0ea5e9", "border": "#0284c7", "text": "#fff"},  # sky
    "Target":    {"bg": "#10b981", "border": "#059669", "text": "#fff"},  # emerald
    "Disease":   {"bg": "#f59e0b", "border": "#d97706", "text": "#fff"},  # amber
    "Mechanism": {"bg": "#8b5cf6", "border": "#7c3aed", "text": "#fff"},  # violet
    "Pathway":   {"bg": "#f43f5e", "border": "#e11d48", "text": "#fff"},  # rose
    "Gene":      {"bg": "#06b6d4", "border": "#0891b2", "text": "#fff"},  # cyan
}


# ── Cytoscape element builders ─────────────────────────────────────────────────

def _node(id_: str, label: str, kind: str, data: dict = None) -> dict:
    col = _COLORS.get(kind, {"bg": "#334155", "border": "#475569", "text": "#fff"})
    return {
        "data": {
            "id": id_,
            "label": label[:28] + ("…" if len(label) > 28 else ""),
            "full_label": label,
            "kind": kind,
            "color": col["bg"],
            "border": col["border"],
            **(data or {}),
        }
    }


def _edge(source: str, target: str, label: str, weight: float = 1.0) -> dict:
    return {
        "data": {
            "source": source,
            "target": target,
            "label": label,
            "weight": weight,
        }
    }


# ── Main graph builder ─────────────────────────────────────────────────────────

def build_graph(
    mesh_ids: List[str] = None,
    compound_id: int = None,
    max_compounds: int = 20,
    max_targets: int = 30,
) -> Tuple[List[dict], dict]:
    """
    Build a Cytoscape graph for the given disease or compound context.

    Returns (elements, legend) where:
      elements = list of {data:{...}} node/edge dicts
      legend   = {kind: color} mapping for the UI legend
    """
    elements: List[dict] = []
    seen_nodes: set = set()
    seen_edges: set = set()

    def add_node(n):
        nid = n["data"]["id"]
        if nid not in seen_nodes:
            seen_nodes.add(nid)
            elements.append(n)

    def add_edge(e):
        key = (e["data"]["source"], e["data"]["target"], e["data"]["label"])
        if key not in seen_edges:
            seen_edges.add(key)
            elements.append(e)

    # ── Disease nodes ──────────────────────────────────────────────────────────
    if mesh_ids:
        rows = _query(
            "SELECT mesh_id, heading FROM mesh_diseases "
            "WHERE mesh_id = ANY(%s::text[]) LIMIT 10",
            (mesh_ids,),
        )
        for r in rows:
            add_node(_node(
                f"dis_{r['mesh_id']}", r["heading"], "Disease",
                {"mesh_id": r["mesh_id"]},
            ))

    # ── Compound nodes ─────────────────────────────────────────────────────────
    if mesh_ids:
        comp_rows = _query(
            """
            SELECT DISTINCT c.id, c.chembl_id, c.name, c.max_phase,
                   i.mesh_id AS d_mesh
            FROM compounds c
            JOIN indications i ON i.compound_id = c.id
            WHERE i.mesh_id = ANY(%s::text[])
              AND c.max_phase >= 2
            ORDER BY c.max_phase DESC NULLS LAST
            LIMIT %s
            """,
            (mesh_ids, max_compounds),
        )
        for r in comp_rows:
            cid = f"cmp_{r['id']}"
            add_node(_node(cid, r["name"], "Compound", {
                "chembl_id": r["chembl_id"],
                "max_phase": r["max_phase"],
                "db_id": r["id"],
            }))
            dmid = f"dis_{r['d_mesh']}"
            if dmid in seen_nodes:
                add_edge(_edge(cid, dmid, "treats", 1.5))

    elif compound_id:
        r = _query("SELECT id,chembl_id,name,max_phase FROM compounds WHERE id=%s", (compound_id,))
        if r:
            r = r[0]
            add_node(_node(f"cmp_{r['id']}", r["name"], "Compound", {
                "chembl_id": r["chembl_id"], "max_phase": r["max_phase"], "db_id": r["id"],
            }))
            ind_rows = _query(
                "SELECT i.mesh_id, md.heading FROM indications i "
                "LEFT JOIN mesh_diseases md ON md.mesh_id = i.mesh_id "
                "WHERE i.compound_id = %s LIMIT 8",
                (compound_id,),
            )
            for ir in ind_rows:
                did = f"dis_{ir['mesh_id']}"
                label = ir.get("heading") or ir["mesh_id"]
                add_node(_node(did, label, "Disease", {"mesh_id": ir["mesh_id"]}))
                add_edge(_edge(f"cmp_{compound_id}", did, "treats", 1.5))

    # ── Target nodes + compound-target edges ──────────────────────────────────
    all_cmp_db_ids = [
        n["data"]["db_id"] for n in elements
        if n["data"].get("kind") == "Compound" and "db_id" in n["data"]
    ]
    if all_cmp_db_ids:
        tgt_rows = _query(
            """
            SELECT DISTINCT m.compound_id, t.id AS target_id,
                   t.name AS tname, t.gene_symbol, m.mechanism
            FROM mechanisms m
            JOIN targets t ON t.id = m.target_id
            WHERE m.compound_id = ANY(%s::int[])
            ORDER BY t.name
            LIMIT %s
            """,
            (all_cmp_db_ids, max_targets),
        )
        for r in tgt_rows:
            tid = f"tgt_{r['target_id']}"
            label = r.get("gene_symbol") or r["tname"]
            add_node(_node(tid, label, "Target", {
                "gene_symbol": r["gene_symbol"],
                "full_name": r["tname"],
            }))
            cid = f"cmp_{r['compound_id']}"
            if cid in seen_nodes:
                add_edge(_edge(cid, tid, r.get("mechanism") or "binds", 1.0))

    legend = {kind: col["bg"] for kind, col in _COLORS.items() if kind in
              {n["data"]["kind"] for n in elements if "kind" in n.get("data", {})}}
    return elements, legend


# ── Cytoscape stylesheet ───────────────────────────────────────────────────────

CYTO_STYLESHEET = [
    {
        "selector": "node",
        "style": {
            "background-color": "data(color)",
            "border-color": "data(border)",
            "border-width": 2,
            "label": "data(label)",
            "color": "#f8fafc",
            "font-size": "11px",
            "text-valign": "center",
            "text-halign": "center",
            "text-wrap": "wrap",
            "text-max-width": "80px",
            "width": 48,
            "height": 48,
        },
    },
    {
        "selector": "node[kind='Disease']",
        "style": {"width": 64, "height": 64, "font-size": "12px", "font-weight": "bold"},
    },
    {
        "selector": "node[kind='Compound']",
        "style": {"shape": "round-rectangle"},
    },
    {
        "selector": "node:selected",
        "style": {
            "border-width": 4,
            "border-color": "#f0f9ff",
            "box-shadow": "0 0 12px rgba(14,165,233,0.8)",
        },
    },
    {
        "selector": "edge",
        "style": {
            "curve-style": "bezier",
            "line-color": "#334155",
            "target-arrow-color": "#475569",
            "target-arrow-shape": "triangle",
            "arrow-scale": 1.2,
            "width": "mapData(weight, 0.5, 2, 1, 3)",
            "label": "data(label)",
            "font-size": "9px",
            "color": "#94a3b8",
            "text-background-color": "#0f172a",
            "text-background-opacity": 0.7,
            "text-background-padding": "2px",
        },
    },
    {
        "selector": "edge[label='treats']",
        "style": {"line-color": "#f59e0b", "target-arrow-color": "#f59e0b", "width": 2.5},
    },
]
