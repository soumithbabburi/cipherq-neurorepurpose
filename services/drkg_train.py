"""
Train DRKG DistMult embeddings  —  the modern-KG link-prediction model.
═══════════════════════════════════════════════════════════════════════════════════
Rebuilds the repurposing subgraph with the PPI capped for tractable CPU training, then
fits the numpy DistMult (services/kg_embedding.DistMultKGE) so a (Compound, treats,
Disease) triple can be scored across DRKG's ~5.9k diseases (vs Hetionet's 213).

Output: data/drkg/drkg_embeddings.npz  {E, R, nodes, rels, dim}
"""
from __future__ import annotations

import time
from pathlib import Path

import numpy as np

from services.drkg_build import build
from services.kg_embedding import DistMultKGE

_DRKG = Path(__file__).parent.parent / "data" / "drkg"
_GRAPH = _DRKG / "drkg_graph.npz"
_EMB = _DRKG / "drkg_embeddings.npz"


def main(cap_gene_gene: int = 400_000, dim: int = 64, epochs: int = 20) -> None:
    print(f"[drkg] rebuilding subgraph (PPI cap {cap_gene_gene:,})…", flush=True)
    build(cap_gene_gene=cap_gene_gene)

    g = np.load(_GRAPH, allow_pickle=True)
    tri = g["triples"].astype(np.int64)
    nodes, rels = g["nodes"], g["rels"]
    nN, nR = len(nodes), len(rels)
    print(f"[drkg] training DistMult: {len(tri):,} triples / {nN:,} nodes / {nR} rels / "
          f"dim {dim} / {epochs} epochs", flush=True)

    t0 = time.time()
    model = DistMultKGE(nN, nR, dim=dim)
    model.fit(tri, epochs=epochs, log=lambda s: print(f"  {s}", flush=True))
    print(f"[drkg] trained in {(time.time() - t0) / 60:.1f} min", flush=True)

    np.savez_compressed(_EMB, E=model.E.astype(np.float32), R=model.R.astype(np.float32),
                        nodes=nodes.astype("U64"), rels=rels.astype("U64"),
                        dim=np.int32(dim))
    print(f"[drkg] saved -> {_EMB}", flush=True)


if __name__ == "__main__":
    main()
