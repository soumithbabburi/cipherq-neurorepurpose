"""
Metapath / Degree-Weighted Path Count (DWPC) features  —  Rephetio-style
═══════════════════════════════════════════════════════════════════════════
Builds the feature family that actually works on Hetionet (Himmelstein et al.,
eLife 2017): for a (Compound, Disease) pair, count the typed multi-hop paths
between them along each metapath, DOWN-WEIGHTING paths through high-degree hub
nodes so promiscuous genes/popular drugs don't dominate.

DWPC, in matrix form: for each metaedge build a degree-weighted sparse matrix
  W = diag(rowdeg^-w) · A · diag(coldeg^-w)        (w = 0.4)
then a metapath's DWPC = the chained product of its (possibly transposed) W's,
yielding a Compound×Disease matrix. Pure scipy.sparse — CPU-feasible.

LEAKAGE: we use ONLY biology metapaths (none traverse the treat edge CtD/CpD),
so predicting "treats" from them is leakage-free by construction — no edge
masking needed for the repoDB decision gate.

This module builds the features; validate_metapath.py does the leakage-controlled
evaluation against repoDB.
"""
from __future__ import annotations

import logging
import numpy as np
import scipy.sparse as sp
from typing import Dict, List, Tuple

logger = logging.getLogger(__name__)

HETIONET_SOURCE = "hetionet_v1.0"
DAMP = 0.4   # DWPC degree-damping exponent (Himmelstein 2017)

# Metaedges we need, with (source_kind, target_kind).
_METAEDGES = {
    "CbG": ("Compound", "Gene"), "CuG": ("Compound", "Gene"), "CdG": ("Compound", "Gene"),
    "DaG": ("Disease", "Gene"),  "DuG": ("Disease", "Gene"),  "DdG": ("Disease", "Gene"),
    "GiG": ("Gene", "Gene"),     "Gr>G": ("Gene", "Gene"),
    "GpPW": ("Gene", "Pathway"), "CrC": ("Compound", "Compound"),
}

# Biology metapaths C→D (no CtD/CpD ⇒ leakage-free). Each step = (metaedge, transpose?).
METAPATHS: Dict[str, List[Tuple[str, bool]]] = {
    "CbGaD":      [("CbG", False), ("DaG", True)],                              # target is a disease gene
    "CdGuD":      [("CdG", False), ("DuG", True)],                              # down-reg a disease-up gene
    "CuGdD":      [("CuG", False), ("DdG", True)],                              # up-reg a disease-down gene
    "CbGiGaD":    [("CbG", False), ("GiG", False), ("DaG", True)],             # target interacts w/ disease gene
    "CbGrGaD":    [("CbG", False), ("Gr>G", False), ("DaG", True)],            # target regulates disease gene
    "CbG_PW_GaD": [("CbG", False), ("GpPW", False), ("GpPW", True), ("DaG", True)],  # shared pathway
    "CrCbGaD":    [("CrC", False), ("CbG", False), ("DaG", True)],             # similar drug hits disease gene
    "CbGaD_via_DuG": [("CbG", False), ("DuG", True)],                          # target up-regulated in disease
    "CbGaD_via_DdG": [("CbG", False), ("DdG", True)],                          # target down-regulated in disease
}


class MetapathFeatures:
    def __init__(self):
        self.kind_index: Dict[str, Dict[str, int]] = {}   # kind -> {node_id: row}
        self.kind_ids: Dict[str, List[str]] = {}
        self.W: Dict[str, sp.csr_matrix] = {}             # metaedge -> degree-weighted matrix
        self.dwpc: Dict[str, np.ndarray] = {}             # metapath -> Compound×Disease dense
        self.cmp_index: Dict[str, int] = {}               # compound node id -> col in dwpc rows
        self.dis_index: Dict[str, int] = {}

    # ── build ────────────────────────────────────────────────────────────────
    def build(self, log=print):
        import psycopg2
        from config import db_params
        conn = psycopg2.connect(**db_params())
        cur = conn.cursor()

        # node indexes per kind
        cur.execute("SELECT id, kind FROM hetionet_nodes")
        for nid, kind in cur.fetchall():
            self.kind_index.setdefault(kind, {})
            idx = self.kind_index[kind]
            if nid not in idx:
                idx[nid] = len(idx)
        for kind, idx in self.kind_index.items():
            ids = [None] * len(idx)
            for nid, i in idx.items():
                ids[i] = nid
            self.kind_ids[kind] = ids
        log(f"  node kinds: " + ", ".join(f"{k}={len(v):,}" for k, v in self.kind_index.items()
                                           if k in ('Compound', 'Gene', 'Disease', 'Pathway')))

        # build degree-weighted matrix per metaedge
        for me, (sk, tk) in _METAEDGES.items():
            cur.execute("SELECT source_id, target_id FROM hetionet_edges "
                        "WHERE metaedge = %s AND source = %s", (me, HETIONET_SOURCE))
            rows, cols = [], []
            si, ti = self.kind_index.get(sk, {}), self.kind_index.get(tk, {})
            for s, t in cur.fetchall():
                if s in si and t in ti:
                    rows.append(si[s]); cols.append(ti[t])
            A = sp.csr_matrix((np.ones(len(rows)), (rows, cols)),
                              shape=(len(si), len(ti)), dtype=np.float64)
            self.W[me] = self._degree_weight(A)
            log(f"  {me:5} {A.shape} nnz={A.nnz:,}")
        conn.close()

        # compound / disease output indexes (rows/cols of every DWPC matrix)
        self.cmp_index = self.kind_index["Compound"]
        self.dis_index = self.kind_index["Disease"]

        # precompute each metapath's Compound×Disease DWPC (dense, small)
        for name, steps in METAPATHS.items():
            try:
                self.dwpc[name] = self._chain(steps)
            except Exception as e:
                log(f"  metapath {name} skipped: {e}")
        log(f"  computed {len(self.dwpc)} metapath DWPC matrices")
        return self

    @staticmethod
    def _degree_weight(A: sp.csr_matrix) -> sp.csr_matrix:
        rowdeg = np.asarray(A.sum(1)).ravel()
        coldeg = np.asarray(A.sum(0)).ravel()
        rw = np.zeros_like(rowdeg); nz = rowdeg > 0; rw[nz] = rowdeg[nz] ** (-DAMP)
        cw = np.zeros_like(coldeg); nz = coldeg > 0; cw[nz] = coldeg[nz] ** (-DAMP)
        return (sp.diags(rw) @ A @ sp.diags(cw)).tocsr()

    def _oriented(self, metaedge: str, transpose: bool) -> sp.csr_matrix:
        W = self.W[metaedge]
        return W.T.tocsr() if transpose else W

    def _chain(self, steps: List[Tuple[str, bool]]) -> np.ndarray:
        M = self._oriented(*steps[0])
        for me, tr in steps[1:]:
            M = M @ self._oriented(me, tr)
        # result is Compound × Disease (sparse) → dense
        return np.asarray(M.todense())

    # ── feature lookup ─────────────────────────────────────────────────────────
    def features(self, compound_node: str, disease_node: str) -> Dict[str, float]:
        ci = self.cmp_index.get(compound_node)
        di = self.dis_index.get(disease_node)
        if ci is None or di is None:
            return {}
        return {name: float(mat[ci, di]) for name, mat in self.dwpc.items()}

    def feature_names(self) -> List[str]:
        return list(self.dwpc.keys())


_singleton: "MetapathFeatures" = None


def get_features_engine(log=print) -> "MetapathFeatures":
    global _singleton
    if _singleton is None:
        _singleton = MetapathFeatures().build(log=log)
    return _singleton
