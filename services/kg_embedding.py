"""
Knowledge-graph embedding (DistMult) — repurposing link prediction
═══════════════════════════════════════════════════════════════════════════
Learns vector embeddings for every node and relation in the Hetionet v1.0 graph
so that a (Compound, treats, Disease) triple can be *scored* — capturing the
indirect, multi-hop biology that the direct target-overlap score misses.

Why DistMult, in numpy: the host has no torch/pykeen (and a 2 GB DL install is
unjustified here). DistMult is a standard, competitive KGE whose score is a
simple trilinear product, so it trains on this graph (≈55k nodes, 2.25M triples)
on CPU with vectorised mini-batch Adagrad in a few minutes.

  score(h, r, t) = Σ_k  E[h]_k · R[r]_k · E[t]_k

Training: logistic loss with negative sampling (corrupt the tail). Sparse
per-batch Adagrad updates (only the rows touched in a batch), so each step costs
O(batch), not O(|nodes|).

This module only builds/trains/scores the model. The leakage-safe evaluation
(held-out treat edges) lives in validation/validate_kg_model.py.
"""
from __future__ import annotations

import time
import logging
import numpy as np
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

HETIONET_SOURCE = "hetionet_v1.0"     # the coherent single-namespace subgraph
TREATS = "CtD"                         # Compound-treats-Disease (the target relation)


def _sigmoid(x: np.ndarray) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-np.clip(x, -30, 30)))


def load_hetionet_triples(source: str = HETIONET_SOURCE
                          ) -> Tuple[np.ndarray, Dict[str, int], Dict[str, int], List[str]]:
    """Load the graph as integer triples (head, rel, tail) from the local DB.

    Returns (triples Nx3 int32, node2id, rel2id, id2node)."""
    import psycopg2
    from config import db_params

    node2id: Dict[str, int] = {}
    rel2id: Dict[str, int] = {}
    id2node: List[str] = []

    def nid(name: str) -> int:
        i = node2id.get(name)
        if i is None:
            i = len(id2node)
            node2id[name] = i
            id2node.append(name)
        return i

    def rid(name: str) -> int:
        i = rel2id.get(name)
        if i is None:
            i = len(rel2id)
            rel2id[name] = i
        return i

    conn = psycopg2.connect(**db_params())
    H: List[int] = []
    R: List[int] = []
    T: List[int] = []
    try:
        cur = conn.cursor(name="kg_stream")     # server-side cursor (streams)
        cur.itersize = 50000
        cur.execute("SELECT source_id, metaedge, target_id FROM hetionet_edges "
                    "WHERE source = %s", (source,))
        for s, m, t in cur:
            H.append(nid(s)); R.append(rid(m)); T.append(nid(t))
        cur.close()
    finally:
        conn.close()

    triples = np.empty((len(H), 3), dtype=np.int32)
    triples[:, 0] = H
    triples[:, 1] = R
    triples[:, 2] = T
    return triples, node2id, rel2id, id2node


class DistMultKGE:
    def __init__(self, n_nodes: int, n_rels: int, dim: int = 64, seed: int = 42):
        rng = np.random.default_rng(seed)
        scale = 6.0 / np.sqrt(dim)
        self.E = ((rng.random((n_nodes, dim), dtype=np.float64) * 2 - 1) * scale)
        self.R = ((rng.random((n_rels, dim), dtype=np.float64) * 2 - 1) * scale)
        self.dim = dim
        self.n_nodes = n_nodes
        self.n_rels = n_rels

    def score(self, h, r, t) -> np.ndarray:
        h = np.asarray(h); r = np.asarray(r); t = np.asarray(t)
        return np.sum(self.E[h] * self.R[r] * self.E[t], axis=-1)

    def fit(self, triples: np.ndarray, epochs: int = 25, batch: int = 2048,
            n_neg: int = 5, lr: float = 0.1, reg: float = 1e-6,
            seed: int = 0, log: Callable = print) -> "DistMultKGE":
        rng = np.random.default_rng(seed)
        H, Rr, T = triples[:, 0], triples[:, 1], triples[:, 2]
        N = len(triples)
        GE = np.full_like(self.E, 1e-8)     # Adagrad accumulators
        GR = np.full_like(self.R, 1e-8)

        for ep in range(epochs):
            t0 = time.time()
            perm = rng.permutation(N)
            tot, nb = 0.0, 0
            for start in range(0, N, batch):
                idx = perm[start:start + batch]
                h, r, t = H[idx], Rr[idx], T[idx]
                B = len(h)
                neg_t = rng.integers(0, self.n_nodes, size=B * n_neg)

                # combined positive + negative samples
                hh = np.concatenate([h, np.repeat(h, n_neg)])
                rr = np.concatenate([r, np.repeat(r, n_neg)])
                tt = np.concatenate([t, neg_t])
                yy = np.concatenate([np.ones(B), np.zeros(B * n_neg)])

                eh, er, et = self.E[hh], self.R[rr], self.E[tt]
                s = np.sum(eh * er * et, axis=1)
                p = _sigmoid(s)
                d = (p - yy)                       # dLoss/dscore  (M,)
                tot += float(-np.mean(yy * np.log(p + 1e-12) +
                                      (1 - yy) * np.log(1 - p + 1e-12)))
                nb += 1

                dcol = d[:, None]
                gE_h = dcol * (er * et)
                gE_t = dcol * (er * eh)
                gR_ = dcol * (eh * et)

                # sparse Adagrad on touched rows only → O(batch) per step
                allnodes = np.concatenate([hh, tt])
                uniqn, invn = np.unique(allnodes, return_inverse=True)
                M = len(hh)
                gradE = np.zeros((len(uniqn), self.dim))
                np.add.at(gradE, invn[:M], gE_h)
                np.add.at(gradE, invn[M:], gE_t)
                gradE = gradE / B + reg * self.E[uniqn]
                GE[uniqn] += gradE * gradE
                self.E[uniqn] -= lr * gradE / np.sqrt(GE[uniqn])

                uniqr, invr = np.unique(rr, return_inverse=True)
                gradR = np.zeros((len(uniqr), self.dim))
                np.add.at(gradR, invr, gR_)
                gradR = gradR / B + reg * self.R[uniqr]
                GR[uniqr] += gradR * gradR
                self.R[uniqr] -= lr * gradR / np.sqrt(GR[uniqr])

            log(f"  epoch {ep + 1:>2}/{epochs}  loss={tot / max(1, nb):.4f}  "
                f"({time.time() - t0:.1f}s)")
        return self

    def save(self, path, node2id: Dict[str, int], rel2id: Dict[str, int]):
        np.savez_compressed(path, E=self.E.astype(np.float32), R=self.R.astype(np.float32),
                            nodes=np.array(list(node2id.keys())),
                            node_ids=np.array(list(node2id.values()), dtype=np.int32),
                            rels=np.array(list(rel2id.keys())),
                            rel_ids=np.array(list(rel2id.values()), dtype=np.int32),
                            dim=self.dim)

    @classmethod
    def load(cls, path) -> Tuple["DistMultKGE", Dict[str, int], Dict[str, int]]:
        d = np.load(path, allow_pickle=True)
        m = cls(d["E"].shape[0], d["R"].shape[0], int(d["dim"]))
        m.E = d["E"].astype(np.float64)
        m.R = d["R"].astype(np.float64)
        node2id = {n: int(i) for n, i in zip(d["nodes"], d["node_ids"])}
        rel2id = {n: int(i) for n, i in zip(d["rels"], d["rel_ids"])}
        return m, node2id, rel2id
