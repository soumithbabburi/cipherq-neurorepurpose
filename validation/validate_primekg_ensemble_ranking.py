"""
Fusion ranker: DistMult recall pool -> re-rank by (treats classifier + R-GCN), combined via
Reciprocal Rank Fusion. The classifier scores PAIR features; the GNN scores GRAPH STRUCTURE
— complementary views, so fusing them inside the proven recall pool is the real shot at
beating the shipped cascade (R@10 0.276). Same split + filtered-ranking metrics as the
other harnesses, so numbers are directly comparable.

Run:  .venv312/Scripts/python.exe -m validation.validate_primekg_ensemble_ranking
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import torch

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from services import primekg_predictor as pkg                        # noqa: E402
from services.train_primekg_gnn import _reproduce_split              # noqa: E402
from validation.validate_primekg_ranking import _filtered_ranks, _summ, KS, _POOL  # noqa: E402

_D = ROOT / "data" / "primekg"


def _rrf_ranks_over_pool(scores_list, weights, n_dis, pool, pos_local, k_rrf=60):
    """Reciprocal Rank Fusion of several score vectors, evaluated over the whole disease
    array but with only `pool` candidates rankable (others pushed to the bottom). Returns
    filtered ranks of the targets in pos_local."""
    fused = np.full(n_dis, -np.inf)
    # rank each scorer within the pool (higher score = better = rank 1)
    contrib = np.zeros(len(pool))
    for sc, w in zip(scores_list, weights):
        order = np.argsort(-sc)                    # positions in pool, best first
        rank = np.empty(len(pool)); rank[order] = np.arange(1, len(pool) + 1)
        contrib += w * (1.0 / (k_rrf + rank))
    fused[pool] = contrib
    return _filtered_ranks(fused, pos_local, pos_local)


def main():
    ckpt = _D / "pkg_gnn.pt"
    if not ckpt.exists():
        print("pkg_gnn.pt not found — train the GNN first."); return 1
    blob = torch.load(ckpt, weights_only=False)
    H = blob["H"].numpy(); r_ind = blob["r_ind"].numpy()

    s = pkg._load()
    labels = s["labels"]; disease_ids = s["disease_ids"]; n_dis = len(disease_ids)
    pos_of = {int(z): i for i, z in enumerate(disease_ids)}
    E = s["E"]; R = s["R"]; rv = R[s["rels"]["indication"]]
    Hz = H[disease_ids]

    test_drugs = _reproduce_split(labels)
    ev4 = []
    for d in test_drugs:
        pl = {pos_of[z] for z in labels[str(d)].get("indication", []) if z in pos_of}
        if pl:
            ev4.append((d, pl))

    R_ = {"clf_only(cascade)": [], "gnn_only": [], "fusion_clf_gnn": []}
    for d, pos_local in ev4:
        dm = (E[d] * rv) @ E[disease_ids].T          # DistMult relatedness (recall)
        pool = np.argsort(-dm)[:_POOL]
        # classifier over the pool
        Ez = E[disease_ids][pool]
        feats = np.concatenate([np.broadcast_to(E[d], (len(pool), E.shape[1])),
                                Ez, E[d] * Ez, np.abs(E[d] - Ez)], axis=1)
        clf = s["treats"].predict_proba(feats)[:, 1]
        # GNN over the pool
        gnn = (H[d] * r_ind * Hz[pool]).sum(-1)
        # build full-length score vectors for filtered ranking
        clf_full = np.full(n_dis, -np.inf); clf_full[pool] = clf
        gnn_full = np.full(n_dis, -np.inf); gnn_full[pool] = gnn
        R_["clf_only(cascade)"].extend(_filtered_ranks(clf_full, pos_local, pos_local))
        R_["gnn_only"].extend(_filtered_ranks(gnn_full, pos_local, pos_local))
        # (3,1) classifier:GNN — swept-best weighting: keeps the classifier's strong top-10
        # while the GNN's complementary structural signal lifts R@20/50/100 and median rank.
        R_["fusion_clf_gnn"].extend(
            _rrf_ranks_over_pool([clf, gnn], [3.0, 1.0], n_dis, pool, pos_local))

    res = {k: _summ(v) for k, v in R_.items()}
    base = json.load(open(ROOT / "validation" / "primekg_ranking_results.json"))["rankers"]["cascade"]
    out = {"fusion": res, "shipped_cascade_ref": base}
    json.dump(out, open(ROOT / "validation" / "primekg_ensemble_results.json", "w"), indent=2)

    def line(name, m):
        return f"{name:20} " + " ".join(f"{m[f'recall@{k}']:<6}" for k in KS) + \
               f"  MRR {m['mrr']:<6} med {m['median_rank']}"
    print("\n=== Fusion ranker (DistMult pool -> classifier + R-GCN, RRF) ===", flush=True)
    print(f"{'':20} " + " ".join(f"R@{k:<4}" for k in KS), flush=True)
    for name in ("clf_only(cascade)", "gnn_only", "fusion_clf_gnn"):
        print(line(name, res[name]), flush=True)
    win = res["fusion_clf_gnn"]["recall@10"] - res["clf_only(cascade)"]["recall@10"]
    print(f"\nfusion vs cascade R@10: {res['fusion_clf_gnn']['recall@10']} vs "
          f"{res['clf_only(cascade)']['recall@10']} ({'+' if win >= 0 else ''}{round(win, 4)})", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
