"""
Build a local training dataset from PrimeKG (Harvard Dataverse doi:10.7910/DVN/IXA7BM).
════════════════════════════════════════════════════════════════════════════════
PrimeKG is the open medical knowledge graph TxGNN is built on: 129,375 nodes,
4,050,249 relationships, 17,080 diseases. This turns the raw dump (nodes.tab +
edges.csv) into compact artifacts our predictor can train on and our platform can
resolve queries against:

  data/primekg/triples.tsv      h_index<TAB>relation<TAB>t_index   (int indices)
  data/primekg/pkg_nodes.json   {index: [type, ext_id, name, source]}
  data/primekg/pkg_index.json   resolution maps: drug (DrugBank id + name) -> index,
                                 disease name -> index, gene symbol -> index; plus the
                                 relation vocabulary with counts.

Node ids are already usable without an external map: gene nodes carry the HGNC symbol
as node_name, drug nodes carry the DrugBank id as node_id, disease nodes carry the
readable MONDO(_grouped) name. So we train on PrimeKG's own node indices and resolve our
queries in by DrugBank id / drug name / disease name / gene symbol.

Run:  python -m services.build_primekg
"""
from __future__ import annotations

import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

_DIR = Path(__file__).parent.parent / "data" / "primekg"
_NODES = _DIR / "nodes.tab"
_EDGES = _DIR / "edges.csv"


def _unq(s: str) -> str:
    s = (s or "").strip()
    if len(s) >= 2 and s[0] == '"' and s[-1] == '"':
        s = s[1:-1]
    return s


def load_nodes():
    """index -> [node_type, node_id, node_name, node_source]. nodes.tab is TAB-separated
    with quoted values: node_index, node_id, node_type, node_name, node_source."""
    nodes = {}
    with open(_NODES, encoding="utf-8") as f:
        header = f.readline()
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            idx = _unq(parts[0])
            try:
                idx = int(idx)
            except ValueError:
                continue
            nodes[idx] = [_unq(parts[2]), _unq(parts[1]), _unq(parts[3]), _unq(parts[4])]
    return nodes


def _sniff_edge_columns(path: Path):
    """Return (relation_col, x_col, y_col) indices from the edges.csv header."""
    with open(path, encoding="utf-8") as f:
        header = f.readline().strip()
    cols = [c.strip().strip('"').lower() for c in header.split(",")]
    def find(*names):
        for n in names:
            if n in cols:
                return cols.index(n)
        return None
    rel = find("relation", "display_relation", "edge_type")
    x = find("x_index", "x_idx", "source", "head")
    y = find("y_index", "y_idx", "target", "tail")
    return rel, x, y, cols


def build():
    if not _EDGES.exists():
        sys.exit(f"edges.csv not found at {_EDGES} — download it first (Phase 0).")
    print("Loading nodes.tab ...")
    nodes = load_nodes()
    print(f"  {len(nodes):,} nodes")

    rel_c, x_c, y_c, cols = _sniff_edge_columns(_EDGES)
    print(f"edges.csv columns: {cols}")
    if None in (rel_c, x_c, y_c):
        sys.exit(f"Could not locate relation/x/y columns in {cols}")

    rel_counts = Counter()
    # drug-disease relation breakdown (the ones that matter: indication/contraindication/off-label)
    dd_rel = Counter()
    type_of = {i: n[0] for i, n in nodes.items()}
    n_edges = 0
    triples_path = _DIR / "triples.tsv"
    with open(_EDGES, encoding="utf-8") as f, open(triples_path, "w", encoding="utf-8") as out:
        r = csv.reader(f)
        next(r, None)  # header
        for row in r:
            if len(row) <= max(rel_c, x_c, y_c):
                continue
            rel = row[rel_c].strip().strip('"')
            try:
                xi = int(row[x_c]); yi = int(row[y_c])
            except ValueError:
                continue
            rel_counts[rel] += 1
            out.write(f"{xi}\t{rel}\t{yi}\n")
            tx, ty = type_of.get(xi, ""), type_of.get(yi, "")
            if {tx, ty} == {"drug", "disease"}:
                dd_rel[rel] += 1
            n_edges += 1
    print(f"  wrote {n_edges:,} triples -> {triples_path.name}")

    # Resolution maps for query-time lookup into PrimeKG's node space.
    drug_by_dbid, drug_by_name, disease_by_name, gene_by_symbol = {}, {}, {}, {}
    for idx, (ntype, nid, name, src) in nodes.items():
        low = name.lower().strip()
        if ntype == "drug":
            if nid:
                drug_by_dbid[nid.upper()] = idx
            if low:
                drug_by_name.setdefault(low, idx)
        elif ntype == "disease":
            if low:
                disease_by_name.setdefault(low, idx)
        elif ntype == "gene/protein":
            if name:
                gene_by_symbol.setdefault(name.upper(), idx)

    json.dump({int(i): v for i, v in nodes.items()},
              open(_DIR / "pkg_nodes.json", "w", encoding="utf-8"))
    json.dump({
        "relations": dict(rel_counts.most_common()),
        "drug_disease_relations": dict(dd_rel.most_common()),
        "drug_by_dbid": drug_by_dbid,
        "drug_by_name": drug_by_name,
        "disease_by_name": disease_by_name,
        "gene_by_symbol": gene_by_symbol,
    }, open(_DIR / "pkg_index.json", "w", encoding="utf-8"))

    print("\n=== node types ===")
    nt = Counter(n[0] for n in nodes.values())
    for k, v in nt.most_common():
        print(f"  {v:>7,}  {k}")
    print("\n=== relation types (all) ===")
    for k, v in rel_counts.most_common():
        print(f"  {v:>9,}  {k}")
    print("\n=== DRUG-DISEASE relations (the repurposing / contraindication labels) ===")
    for k, v in dd_rel.most_common():
        print(f"  {v:>7,}  {k}")
    print(f"\nresolution maps: {len(drug_by_dbid):,} drugs(DBID) / {len(drug_by_name):,} drugs(name) / "
          f"{len(disease_by_name):,} diseases / {len(gene_by_symbol):,} genes")


if __name__ == "__main__":
    build()
