"""
Hetionet-backed, entity-agnostic knowledge graph for the /graph explorer.
════════════════════════════════════════════════════════════════════════════════
Hetionet (55,926 nodes / 2.26M edges, in the platform DB) is ONE unified graph that
already contains every entity kind — Compound, Gene, Disease, Pathway, Anatomy,
Side-effect, Pharmacologic-class… — so a drug, gene, pathway OR disease query all
resolve the SAME way, and click-to-expand is just a 1-hop neighbour query. Node ids
encode the kind ('Compound::…','Gene::ABCB1','Disease::DOID:…','Pathway::…').

Returns cytoscape-style elements {nodes with kind/color/size, edges with relation
label + weight} so the existing /graph front-end can render them directly.
"""
from __future__ import annotations

import logging
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)

# metaedge code → human relation label. Codes read source-kind / relation / target-kind
# (CbG = Compound-binds-Gene, DaG = Disease-associates-Gene, GpPW = Gene-participates-Pathway…)
_REL = {
    "CtD": "treats", "CpD": "palliates", "CbG": "binds", "CuG": "up-regulates",
    "CdG": "down-regulates", "DaG": "associated", "DuG": "up-regulates",
    "DdG": "down-regulates", "DlA": "localizes to", "GpPW": "participates in",
    "GpBP": "participates in", "GpMF": "has function", "GpCC": "in component",
    "GiG": "interacts", "Gr>G": "regulates", "GcG": "covaries", "AeG": "expresses",
    "AuG": "over-expresses", "AdG": "under-expresses", "CcSE": "causes",
    "CrC": "resembles", "DrD": "resembles", "PCiC": "in class",
}
# The relations worth showing (skip the huge noisy ones — GpBP 559k, AeG 526k — unless a
# user deliberately expands a gene; then they still come through 1-hop, capped per relation).
_CORE = {"CtD", "CpD", "CbG", "CuG", "CdG", "DaG", "DuG", "DdG", "DlA",
         "GpPW", "GiG", "PCiC", "CcSE", "GpMF"}
# relation → edge weight (importance), for front-end sizing / particle speed
_REL_W = {"CtD": 1.0, "CpD": 0.85, "DaG": 0.8, "CbG": 0.78, "CuG": 0.6, "CdG": 0.6,
          "DuG": 0.55, "DdG": 0.55, "GpPW": 0.55, "GiG": 0.45, "DlA": 0.5,
          "PCiC": 0.5, "CcSE": 0.35, "GpMF": 0.4}
_COLOR = {"Compound": "#3b82f6", "Gene": "#22d3ee", "Disease": "#f59e0b",
          "Pathway": "#a78bfa", "Anatomy": "#34d399", "Symptom": "#fb7185",
          "Side Effect": "#f472b6", "Pharmacologic Class": "#818cf8",
          "Biological Process": "#c084fc", "Molecular Function": "#67e8f9",
          "Cellular Component": "#5eead4"}


def _q(sql, params=()):
    from services.neuro_db_service import _query
    return _query(sql, params)


def _kind_of(node_id: str) -> str:
    return node_id.split("::", 1)[0] if "::" in node_id else "Unknown"


def resolve_entity(query: str) -> Optional[Dict]:
    """Find the Hetionet node best matching a free-text query — ANY kind. Exact
    (case-insensitive) name first, then prefix, then contains; preferring the kinds a
    user most likely means (Compound → Disease → Pathway → Gene)."""
    q = (query or "").strip()
    if not q:
        return None
    order = "CASE kind WHEN 'Compound' THEN 0 WHEN 'Disease' THEN 1 WHEN 'Pathway' THEN 2 WHEN 'Gene' THEN 3 ELSE 9 END"
    for cond, param in ((f"lower(name)=lower(%s)", q),
                        ("name ILIKE %s", q + "%"),
                        ("name ILIKE %s", "%" + q + "%")):
        rows = _q(f"SELECT id,name,kind FROM hetionet_nodes WHERE {cond} ORDER BY {order}, length(name) LIMIT 6", (param,))
        if rows:
            return rows[0]
    return None


def _neighbours(node_id: str, rels=None, cap_per_rel: int = 9) -> List[Dict]:
    """1-hop edges touching node_id over the given relations, capped per relation so a
    hub gene doesn't dump 500 pathway edges. UNION of two indexed lookups (out/in)."""
    rels = list(rels or _CORE)
    rows = _q(
        """
        SELECT source_id, target_id, metaedge FROM (
          SELECT source_id, target_id, metaedge FROM hetionet_edges
            WHERE source_id = %s AND metaedge = ANY(%s)
          UNION ALL
          SELECT source_id, target_id, metaedge FROM hetionet_edges
            WHERE target_id = %s AND metaedge = ANY(%s)
        ) e
        """, (node_id, rels, node_id, rels))
    # cap per metaedge (keep a balanced spread of relation types)
    by_rel: Dict[str, List[Dict]] = {}
    for r in rows:
        by_rel.setdefault(r["metaedge"], []).append(r)
    out: List[Dict] = []
    for me, es in by_rel.items():
        out.extend(es[:cap_per_rel])
    return out


def _elements(anchor: Dict, edges: List[Dict], existing_ids=None):
    """Assemble cytoscape-style elements from an anchor + its edges."""
    node_ids = {anchor["id"]}
    for e in edges:
        node_ids.add(e["source_id"]); node_ids.add(e["target_id"])
    names = {}
    if node_ids:
        rows = _q("SELECT id,name,kind FROM hetionet_nodes WHERE id = ANY(%s)", (list(node_ids),))
        names = {r["id"]: r for r in rows}
    els = []
    existing = existing_ids or set()
    for nid in node_ids:
        if nid in existing:
            continue
        meta = names.get(nid, {})
        kind = meta.get("kind") or _kind_of(nid)
        is_anchor = (nid == anchor["id"])
        els.append({"data": {
            "id": nid, "label": (meta.get("name") or nid.split("::")[-1])[:40],
            "full_label": meta.get("name") or nid, "kind": kind,
            "color": _COLOR.get(kind, "#94a3b8"),
            "size": 70 if is_anchor else 34, "anchor": is_anchor,
            "hetionet_id": nid}})
    for e in edges:
        me = e["metaedge"]
        els.append({"data": {
            "source": e["source_id"], "target": e["target_id"],
            "label": _REL.get(me, me), "metaedge": me,
            "weight": _REL_W.get(me, 0.4)}})
    return els


def build(query: str, limit: int = 70) -> Dict:
    """Entity-agnostic graph around whatever the query resolves to (drug/gene/pathway/
    disease). 1-hop neighbourhood — click a node to expand further."""
    anchor = resolve_entity(query)
    if not anchor:
        return {"elements": [], "legend": {}, "anchor": None}
    edges = _neighbours(anchor["id"])
    els = _elements(anchor, edges)
    kinds = {el["data"]["kind"] for el in els if "kind" in el.get("data", {})}
    legend = {k: _COLOR.get(k, "#94a3b8") for k in kinds}
    return {"elements": els, "legend": legend, "hetionet": True,
            "anchor": {"id": anchor["id"], "name": anchor["name"], "kind": anchor["kind"]}}


def expand(node_id: str, existing_ids=None) -> Dict:
    """1-hop neighbours of a node the user clicked — for graph traversal. Returns only
    the NEW elements not already on the canvas."""
    if not node_id:
        return {"elements": []}
    row = _q("SELECT id,name,kind FROM hetionet_nodes WHERE id=%s LIMIT 1", (node_id,))
    if not row:
        return {"elements": []}
    anchor = row[0]
    edges = _neighbours(node_id)
    els = _elements(anchor, edges, existing_ids=set(existing_ids or []))
    return {"elements": els}


def available() -> bool:
    try:
        return bool(_q("SELECT 1 FROM hetionet_nodes LIMIT 1"))
    except Exception:
        return False
