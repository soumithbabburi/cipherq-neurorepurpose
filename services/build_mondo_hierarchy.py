"""
Build the MONDO hierarchy edge table from the local MONDO ontology JSON.
═══════════════════════════════════════════════════════════════════════════════
MeSH has orphan rows (e.g. Churg-Strauss/EGPA with no tree numbers or parents), so
subtype expansion misses subtypes that are legitimately under a parent. MONDO is a
complete, curated disease ontology with proper is_a edges — so we extract its
child→parent graph once into `mondo_edges` inside disease_value.db (which already
holds MONDO id ↔ label ↔ MeSH xref), giving the platform a robust hierarchy backbone.

Run once (and whenever mondo-simple.json is refreshed):
    python -m services.build_mondo_hierarchy
"""
from __future__ import annotations

import json
import re
import sqlite3
from pathlib import Path

_REF = Path(__file__).parent.parent / "data" / "disease_reference"
_DB = _REF / "disease_value.db"
_JSON = _REF / "mondo-simple.json"


def _mid(uri: str):
    m = re.search(r"MONDO_(\d+)", uri or "")
    return f"MONDO:{m.group(1)}" if m else None


def build() -> int:
    if not _JSON.exists():
        raise FileNotFoundError(f"MONDO JSON not found: {_JSON}")
    data = json.load(open(_JSON, encoding="utf-8"))
    graph = data["graphs"][0]

    edges = []
    for e in graph.get("edges", []):
        if e.get("pred") != "is_a":
            continue
        child, parent = _mid(e.get("sub")), _mid(e.get("obj"))
        if child and parent and child != parent:
            edges.append((child, parent))

    # Labels for MONDO ids the diseases table may not carry (obsolete/rare) — so a
    # descendant always has a name even if it isn't in the value table.
    labels = []
    for n in graph.get("nodes", []):
        mid = _mid(n.get("id"))
        lbl = (n.get("lbl") or "").strip()
        if mid and lbl and not lbl.lower().startswith("obsolete"):
            labels.append((mid, lbl))

    con = sqlite3.connect(_DB)
    cur = con.cursor()
    cur.execute("DROP TABLE IF EXISTS mondo_edges")
    cur.execute("CREATE TABLE mondo_edges (child TEXT, parent TEXT)")
    cur.executemany("INSERT INTO mondo_edges VALUES (?, ?)", edges)
    cur.execute("CREATE INDEX ix_me_child ON mondo_edges(child)")
    cur.execute("CREATE INDEX ix_me_parent ON mondo_edges(parent)")

    cur.execute("DROP TABLE IF EXISTS mondo_labels")
    cur.execute("CREATE TABLE mondo_labels (mondo_id TEXT PRIMARY KEY, label TEXT)")
    cur.executemany("INSERT OR REPLACE INTO mondo_labels VALUES (?, ?)", labels)

    con.commit()
    con.close()
    return len(edges)


if __name__ == "__main__":
    n = build()
    print(f"mondo_edges built: {n} is_a edges")
