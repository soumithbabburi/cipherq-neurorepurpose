"""
DRKG subgraph builder  —  the modern-KG upgrade (replaces 2016 Hetionet's 213-disease
ceiling with DRKG's ~5.9k diseases / ~24k genes / DrugBank compound space)
═══════════════════════════════════════════════════════════════════════════════════
DRKG (Drug Repurposing Knowledge Graph, gnn4dr/DRKG) unifies 6 sources — DrugBank,
Hetionet, GNBR, STRING, IntAct, DGIdb — into 5.87M edges. We keep the repurposing-
relevant subgraph (drug-gene, gene-disease, drug-disease, gene-pathway, gene-gene PPI)
and drop the bulk that does not help drug->disease link prediction (compound-compound
similarity, anatomy expression, GO terms, side-effects). The result is a ~2.9M-triple
graph on the SAME order of magnitude as the Hetionet embedding we already train, so the
numpy DistMult trainer (services/kg_embedding.DistMultKGE) handles it on CPU.

Entity namespaces (so downstream can map our IDs in):
  Compound::DB#####      DrugBank   (map from chembl_id via UniChem, or name via Hetionet)
  Disease::MESH:D######  MeSH       (map from our disease mesh_id directly)
  Disease::DOID:######   Disease Ontology
  Gene::<entrez>         NCBI gene  (map from gene symbol via mygene)

Outputs (data/drkg/):
  drkg_graph.npz         triples(int32 Nx3), nodes(str), rels(str)
  drkg_names.json        {node_id: display_name} for Compound(DB)+Disease(MESH/DOID)
"""
from __future__ import annotations

import json
import logging
from pathlib import Path

import numpy as np

logger = logging.getLogger(__name__)
log = lambda m: (logger.info(m), print(m, flush=True))

_DRKG_DIR = Path(__file__).parent.parent / "data" / "drkg"
_TSV = _DRKG_DIR / "drkg.tsv"
_GRAPH = _DRKG_DIR / "drkg_graph.npz"
_NAMES = _DRKG_DIR / "drkg_names.json"

# Keep an edge only if its {head_type, tail_type} is one of these (order-free).
_KEEP_PAIRS = {
    frozenset({"Compound", "Gene"}),
    frozenset({"Gene", "Disease"}),
    frozenset({"Compound", "Disease"}),
    frozenset({"Gene", "Pathway"}),
    frozenset({"Gene"}),            # Gene-Gene PPI
}


def _etype(entity: str) -> str:
    # "Compound::DB00294" -> "Compound"; "Gene::2157" -> "Gene"
    return entity.split("::", 1)[0]


def build(cap_gene_gene: int | None = None) -> None:
    """Parse DRKG, keep the repurposing subgraph, write int-id triples + name map."""
    if not _TSV.exists():
        raise FileNotFoundError(f"DRKG edges not found at {_TSV} — run the download first")

    node2id: dict[str, int] = {}
    rel2id: dict[str, int] = {}
    other: list[tuple[int, int, int]] = []   # non-PPI repurposing triples (kept whole)
    gene_gene: list[tuple[int, int, int]] = []  # PPI (random-subsampled if capped)
    seen = 0

    def nid(x: str) -> int:
        i = node2id.get(x)
        if i is None:
            i = node2id[x] = len(node2id)
        return i

    def rid(x: str) -> int:
        i = rel2id.get(x)
        if i is None:
            i = rel2id[x] = len(rel2id)
        return i

    with _TSV.open("r", encoding="utf-8") as fh:
        for line in fh:
            seen += 1
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 3:
                continue
            h, r, t = parts
            pair = frozenset({_etype(h), _etype(t)})
            if pair not in _KEEP_PAIRS:
                continue
            trip = (nid(h), rid(r), nid(t))
            (gene_gene if pair == frozenset({"Gene"}) else other).append(trip)
            if seen % 1_000_000 == 0:
                log(f"  scanned {seen:,}  kept {len(other) + len(gene_gene):,}  "
                    f"(nodes {len(node2id):,})")

    # PPI is 82% of edges and only a secondary signal for drug->disease link prediction —
    # random-subsample it so DistMult trains in minutes without losing gene-gene structure.
    if cap_gene_gene is not None and len(gene_gene) > cap_gene_gene:
        rng = np.random.default_rng(42)
        pick = rng.choice(len(gene_gene), size=cap_gene_gene, replace=False)
        gene_gene = [gene_gene[i] for i in pick]
    gg_count = len(gene_gene)
    tri = np.asarray(other + gene_gene, dtype=np.int32)
    nodes = np.empty(len(node2id), dtype=object)
    for name, i in node2id.items():
        nodes[i] = name
    rels = np.empty(len(rel2id), dtype=object)
    for name, i in rel2id.items():
        rels[i] = name

    np.savez_compressed(_GRAPH, triples=tri,
                        nodes=nodes.astype("U64"), rels=rels.astype("U64"))
    log(f"kept {len(tri):,} triples / {len(node2id):,} nodes / {len(rel2id)} relations "
        f"(gene-gene {gg_count:,}) -> {_GRAPH.name}")

    _build_names(list(node2id.keys()))


def _build_names(node_keys: list[str]) -> None:
    """Display names for Compound(DrugBank) + Disease(MeSH/DOID) nodes, from our DB."""
    names: dict[str, str] = {}
    try:
        import psycopg2
        from config import db_params
        conn = psycopg2.connect(**db_params())
        cur = conn.cursor()
        # DrugBank + DOID names live in hetionet_nodes (already DRKG-namespaced).
        cur.execute("SELECT id, name FROM hetionet_nodes WHERE kind IN ('Compound','Disease')")
        for nid_, nm in cur.fetchall():
            if nid_ and nm:
                names[nid_] = nm
        # MeSH disease names from the mesh_diseases table (mesh_id -> label).
        cur.execute("SELECT column_name FROM information_schema.columns WHERE table_name='mesh_diseases'")
        cols = [r[0] for r in cur.fetchall()]
        idc = next((c for c in ("mesh_id", "descriptor_ui", "ui", "id") if c in cols), None)
        nmc = next((c for c in ("heading", "name", "label", "descriptor_name", "term") if c in cols), None)
        if idc and nmc:
            cur.execute(f"SELECT {idc}, {nmc} FROM mesh_diseases")
            for mid, nm in cur.fetchall():
                if mid and nm:
                    names[f"Disease::MESH:{mid}"] = nm
        conn.close()
    except Exception as e:  # pragma: no cover - name map is best-effort
        log(f"  name map: DB lookup failed ({e}); names limited to Hetionet subset")

    covered = sum(1 for k in node_keys if k in names)
    _NAMES.write_text(json.dumps(names), encoding="utf-8")
    log(f"names: {len(names):,} labels ({covered:,} DRKG nodes covered) -> {_NAMES.name}")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    build()
