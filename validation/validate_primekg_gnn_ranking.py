"""
Evaluate the R-GCN universe ranker on the held-out ranking protocol and compare to the
shipped cascade baseline. Uses the SAME compound-disjoint split and filtered-ranking
metrics as validate_primekg_ranking, so the numbers are directly comparable.

Run:  .venv312/Scripts/python.exe -m validation.validate_primekg_gnn_ranking
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

from services import primekg_predictor as pkg                       # noqa: E402
from services.train_primekg_gnn import _reproduce_split             # noqa: E402
from validation.validate_primekg_ranking import _filtered_ranks, _summ, KS  # noqa: E402

_D = ROOT / "data" / "primekg"


def main():
    ckpt = _D / "pkg_gnn.pt"
    if not ckpt.exists():
        print("pkg_gnn.pt not found — train it first (services.train_primekg_gnn)."); return 1
    blob = torch.load(ckpt, weights_only=False)
    H = blob["H"]; r_ind = blob["r_ind"]
    s = pkg._load()
    labels = s["labels"]; disease_ids = s["disease_ids"]
    pos_of = {int(z): i for i, z in enumerate(disease_ids)}
    dis_t = torch.tensor(disease_ids)

    test_drugs = _reproduce_split(labels)
    ev4 = []
    for d in test_drugs:
        pl = {pos_of[z] for z in labels[str(d)].get("indication", []) if z in pos_of}
        if pl:
            ev4.append((d, pl))

    Hz = H[dis_t]                                   # (n_dis, dim)
    ranks = []
    with torch.no_grad():
        for d, pos_local in ev4:
            scores = (H[d] * r_ind * Hz).sum(-1).numpy()
            ranks.extend(_filtered_ranks(scores, pos_local, pos_local))
    gnn = _summ(ranks)

    base = json.load(open(ROOT / "validation" / "primekg_ranking_results.json"))["rankers"]
    out = {"gnn": gnn, "baseline_cascade": base["cascade"],
           "baseline_distmult": base["distmult"], "meta": blob.get("meta", {})}
    json.dump(out, open(ROOT / "validation" / "primekg_gnn_ranking_results.json", "w"), indent=2)

    def line(name, m):
        return f"{name:18} " + " ".join(f"{m[f'recall@{k}']:<6}" for k in KS) + \
               f"  MRR {m['mrr']:<6} med {m['median_rank']}"
    print("\n=== R-GCN universe ranker vs baselines (held-out, all diseases) ===", flush=True)
    print(f"{'':18} " + " ".join(f"R@{k:<4}" for k in KS), flush=True)
    print(line("distmult", base["distmult"]), flush=True)
    print(line("cascade (shipped)", base["cascade"]), flush=True)
    print(line("R-GCN (new)", gnn), flush=True)
    lift = gnn["recall@10"] - base["cascade"]["recall@10"]
    print(f"\nR@10 vs shipped cascade: {gnn['recall@10']} vs {base['cascade']['recall@10']} "
          f"({'+' if lift >= 0 else ''}{round(lift, 4)})", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
