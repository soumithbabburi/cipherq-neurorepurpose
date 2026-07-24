"""
PrimeKG universe ranker — a real message-passing GNN (R-GCN, CPU).
════════════════════════════════════════════════════════════════════════════════════
The shallow treats classifier over DistMult-embedding pair features is at its ceiling
(validation/experiment_primekg_ranker_v2: negative sampling REGRESSES it). The missing
ingredient is graph STRUCTURE — a drug's indication should be inferable from its targets'
neighbourhoods, which only message passing captures. This trains a relational GCN that
does exactly that, on CPU, with no torch-geometric dependency.

Design (chosen for CPU tractability + correctness):
  * Input node features = the pretrained DistMult embeddings E (frozen) — a huge head start,
    so the GNN only has to LEARN THE STRUCTURE, not embeddings from scratch.
  * R-GCN with the 30 PrimeKG relations bucketed into ~8 semantic GROUPS (full 30-matrix
    R-GCN is memory-heavy; grouping stays relational and fits in RAM). Row-normalised
    (mean) aggregation per group + a self transform. 2 layers.
  * Direction-aware DistMult decoder on a learned `indication` relation vector (init from
    R[indication]); trained with train-drug indications as positives and their
    CONTRAINDICATIONS as hard negatives + random-disease easy negatives.
  * ZERO-SHOT protocol: the held-out (compound-disjoint, seed 42) test drugs' treatment
    edges are removed from BOTH the message-passing graph and supervision. Their
    drug->target edges REMAIN, so the model must generalise through shared biology.

Output:
  data/primekg/pkg_gnn.pt   {H (refined embeddings), r_ind, meta}
Run:  .venv312/Scripts/python.exe -m services.train_primekg_gnn
"""
from __future__ import annotations

import json
import random
import time
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn

_D = Path(__file__).parent.parent / "data" / "primekg"

# 30 PrimeKG relations -> semantic groups (keeps R-GCN tractable on CPU)
_GROUP = {
    "protein_protein": "ppi",
    "drug_protein": "drug_target",
    "disease_protein": "disease_gene",
    "contraindication": "treat", "indication": "treat", "off-label use": "treat",
    "drug_drug": "drug_drug",
    "drug_effect": "drug_effect",
    "phenotype_protein": "phenotype", "phenotype_phenotype": "phenotype",
    "disease_phenotype_negative": "phenotype", "disease_phenotype_positive": "phenotype",
    "disease_disease": "disease_gene",
    "bioprocess_bioprocess": "function", "molfunc_molfunc": "function",
    "cellcomp_cellcomp": "function", "molfunc_protein": "function",
    "cellcomp_protein": "function", "bioprocess_protein": "function",
    "pathway_pathway": "function", "pathway_protein": "function",
    "exposure_protein": "context", "exposure_disease": "context",
    "exposure_exposure": "context", "exposure_bioprocess": "context",
    "exposure_molfunc": "context", "exposure_cellcomp": "context",
    "anatomy_anatomy": "context", "anatomy_protein_present": "context",
    "anatomy_protein_absent": "context",
}
_TREAT_RELS = {"contraindication", "indication", "off-label use"}


def _reproduce_split(labels, seed=42, holdout=0.15):
    rng = random.Random(seed)
    drugs = [int(d) for d, r in labels.items() if r.get("indication")]
    rng.shuffle(drugs)
    return set(drugs[: int(len(drugs) * holdout)])


def _row_norm_adj(h, t, n):
    """Row-normalised (mean-aggregation) sparse adjacency from head->tail index arrays."""
    idx = torch.tensor(np.vstack([t, h]), dtype=torch.long)      # A[t,h]: message h->t
    deg = torch.zeros(n).scatter_add_(0, torch.tensor(t, dtype=torch.long),
                                      torch.ones(len(t)))
    val = 1.0 / deg[torch.tensor(t, dtype=torch.long)].clamp(min=1.0)
    return torch.sparse_coo_tensor(idx, val, (n, n)).coalesce()


class RGCNLayer(nn.Module):
    def __init__(self, dim, n_groups):
        super().__init__()
        self.Wg = nn.Parameter(torch.empty(n_groups, dim, dim))
        self.Ws = nn.Parameter(torch.empty(dim, dim))
        nn.init.xavier_uniform_(self.Wg)
        nn.init.xavier_uniform_(self.Ws)

    def forward(self, H, adjs):
        out = H @ self.Ws
        for gi, A in enumerate(adjs):
            out = out + torch.sparse.mm(A, H @ self.Wg[gi])
        return out


class RGCN(nn.Module):
    def __init__(self, E, n_groups, r_ind_init, layers=2):
        super().__init__()
        self.register_buffer("X", torch.tensor(E, dtype=torch.float32))  # frozen features
        dim = E.shape[1]
        self.layers = nn.ModuleList([RGCNLayer(dim, n_groups) for _ in range(layers)])
        self.act = nn.ReLU()
        self.r_ind = nn.Parameter(torch.tensor(r_ind_init, dtype=torch.float32))

    def encode(self, adjs):
        # RESIDUAL message passing: each layer adds a structural CORRECTION on top of the
        # pretrained DistMult geometry rather than replacing it. Without this, 2 layers of
        # mean aggregation over PrimeKG's huge-degree generic groups (ppi/function/context)
        # over-smooth the drug/disease-specific signal and ranking collapses below raw
        # DistMult (measured: R@10 0.039). The skip keeps the informative input geometry.
        H = self.X
        for lyr in self.layers:
            H = H + self.act(lyr(H, adjs))
        return H

    def score(self, H, d_idx, z_idx):
        return (H[d_idx] * self.r_ind * H[z_idx]).sum(-1)


def main(seed=42, epochs=60, lr=0.01, neg_per_pos=8, layers=2):
    torch.manual_seed(seed)
    rng = random.Random(seed)
    t0 = time.time()

    emb = np.load(_D / "pkg_embeddings.npz"); E = emb["E"]; R = emb["R"]
    nodes = json.load(open(_D / "pkg_nodes.json"))
    labels = json.load(open(_D / "pkg_labels.json"))
    rels = json.load(open(_D / "pkg_rels.json"))
    n = E.shape[0]
    disease_ids = np.array([int(i) for i, v in nodes.items() if v[0] == "disease"], dtype=np.int64)
    test_drugs = _reproduce_split(labels)
    print(f"[gnn] {n} nodes | {len(disease_ids)} diseases | {len(test_drugs)} held-out drugs "
          f"| torch {torch.__version__} {torch.get_num_threads()} threads", flush=True)

    # ── load triples, drop test-drug treatment edges, bucket into groups ──────────
    import pandas as pd
    df = pd.read_csv(_D / "triples.tsv", sep="\t", header=None, names=["h", "r", "t"],
                     dtype={"h": np.int64, "t": np.int64, "r": str})
    print(f"[gnn] {len(df):,} triples loaded ({time.time()-t0:.0f}s)", flush=True)
    is_treat = df["r"].isin(_TREAT_RELS).to_numpy()
    h_in_test = df["h"].isin(test_drugs).to_numpy()
    t_in_test = df["t"].isin(test_drugs).to_numpy()
    drop = is_treat & (h_in_test | t_in_test)          # zero-shot: hide test-drug treatments
    df = df[~drop]
    print(f"[gnn] dropped {int(drop.sum()):,} test-drug treatment edges (zero-shot)", flush=True)
    df["g"] = df["r"].map(_GROUP)
    groups = sorted(df["g"].dropna().unique())
    adjs = []
    for g in groups:
        sub = df[df["g"] == g]
        adjs.append(_row_norm_adj(sub["h"].to_numpy(), sub["t"].to_numpy(), n))
    print(f"[gnn] {len(groups)} relation groups: {groups} ({time.time()-t0:.0f}s)", flush=True)

    # ── supervision: train-drug indications (pos) + contraindications (hard neg) ──
    train_pos = []      # (drug, disease)
    train_hard = []     # (drug, disease) contraindications
    dset = set(int(z) for z in disease_ids)
    for d, r in labels.items():
        di = int(d)
        if di in test_drugs:
            continue
        for z in r.get("indication", []):
            if z in dset:
                train_pos.append((di, z))
        for z in r.get("contraindication", []):
            if z in dset:
                train_hard.append((di, z))
    train_pos = np.array(train_pos, dtype=np.int64)
    train_hard = np.array(train_hard, dtype=np.int64)
    print(f"[gnn] supervision: {len(train_pos)} pos / {len(train_hard)} hard-neg contra", flush=True)

    model = RGCN(E, len(groups), R[rels["indication"]], layers=layers)
    opt = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=1e-5)
    bce = nn.BCEWithLogitsLoss()

    pos_d = torch.tensor(train_pos[:, 0]); pos_z = torch.tensor(train_pos[:, 1])
    hard_d = torch.tensor(train_hard[:, 0]); hard_z = torch.tensor(train_hard[:, 1])
    dis_t = torch.tensor(disease_ids)

    model.train()
    for ep in range(1, epochs + 1):
        opt.zero_grad()
        H = model.encode(adjs)
        # random easy negatives: same drugs as positives, random diseases
        neg_d = pos_d.repeat(neg_per_pos)
        neg_z = dis_t[torch.randint(0, len(dis_t), (len(pos_d) * neg_per_pos,))]
        d_all = torch.cat([pos_d, hard_d, neg_d])
        z_all = torch.cat([pos_z, hard_z, neg_z])
        y_all = torch.cat([torch.ones(len(pos_d)),
                           torch.zeros(len(hard_d)),
                           torch.zeros(len(neg_d))])
        logits = model.score(H, d_all, z_all)
        loss = bce(logits, y_all)
        loss.backward()
        opt.step()
        if ep == 1 or ep % 10 == 0:
            print(f"[gnn] epoch {ep:3d}  loss {loss.item():.4f}  ({time.time()-t0:.0f}s)", flush=True)

    model.eval()
    with torch.no_grad():
        H = model.encode(adjs).detach()
    torch.save({"H": H, "r_ind": model.r_ind.detach(), "groups": groups,
                "meta": {"epochs": epochs, "layers": layers, "neg_per_pos": neg_per_pos,
                         "seed": seed}}, _D / "pkg_gnn.pt")
    print(f"[gnn] saved {_D / 'pkg_gnn.pt'} ({time.time()-t0:.0f}s)", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
