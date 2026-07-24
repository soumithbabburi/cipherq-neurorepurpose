"""
Train a DistMult link-prediction model on PrimeKG, and validate it (Phase 2).
════════════════════════════════════════════════════════════════════════════════
Reuses the proven numpy DistMult (services.kg_embedding.DistMultKGE) that already backs
the DRKG predictor. Trains on the full PrimeKG triple set (17,080 diseases) with a
held-out slice of `indication` edges, then measures how well the model recovers those
held-out true indications ABOVE two kinds of negative — random drug/disease pairs and the
drug's real `contraindication` edges. That second AUC is the honest test: can the model
tell a true indication from a genuine contraindication for the same drug?

PrimeKG node_index is used directly as the embedding row (no remap), so the predictor can
resolve a query drug/disease -> node_index (via pkg_index.json) -> embedding.

Output:
  data/primekg/pkg_embeddings.npz   {E, R, dim}
  data/primekg/pkg_rels.json        {relation_name: rel_id}
  data/primekg/pkg_eval.json        recovery-AUC metrics

Run:  python -m services.train_primekg
"""
from __future__ import annotations

import json
import random
from pathlib import Path

import numpy as np

from services.kg_embedding import DistMultKGE

_D = Path(__file__).parent.parent / "data" / "primekg"
_TRIP = _D / "triples.tsv"


def _auc(pos, neg) -> float:
    """Ranking AUC = P(random positive scores above random negative)."""
    pos = np.asarray(pos, dtype=np.float64); neg = np.asarray(neg, dtype=np.float64)
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    # rank-sum (Mann-Whitney) formulation, robust for large n
    allv = np.concatenate([pos, neg])
    order = allv.argsort()
    ranks = np.empty_like(order, dtype=np.float64)
    ranks[order] = np.arange(1, len(allv) + 1)
    r_pos = ranks[:len(pos)].sum()
    return float((r_pos - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg)))


def main(dim: int = 64, epochs: int = 20, holdout: float = 0.10, seed: int = 42):
    rng = random.Random(seed)
    print(f"[primekg] reading {_TRIP.name} ...", flush=True)
    rel2id: dict = {}
    H, R, T = [], [], []
    ind_test = []          # held-out (drug_idx, disease_idx) indication pairs
    contra = {}            # drug_idx -> set(disease_idx) contraindications (for eval negatives)
    n_ind = 0
    with open(_TRIP, encoding="utf-8") as f:
        for line in f:
            h, r, t = line.rstrip("\n").split("\t")
            h = int(h); t = int(t)
            rid = rel2id.setdefault(r, len(rel2id))
            if r == "indication":
                n_ind += 1
                if rng.random() < holdout:
                    ind_test.append((h, t))
                    continue                    # exclude held-out from training
            elif r == "contraindication":
                contra.setdefault(h, set()).add(t)
                contra.setdefault(t, set()).add(h)
            H.append(h); R.append(rid); T.append(t)

    tri = np.column_stack([np.array(H, dtype=np.int64),
                           np.array(R, dtype=np.int64),
                           np.array(T, dtype=np.int64)])
    nN = int(tri[:, [0, 2]].max()) + 1
    nR = len(rel2id)
    print(f"[primekg] {len(tri):,} training triples / {nN:,} nodes / {nR} rels / "
          f"indication edges {n_ind:,} (held out {len(ind_test):,}) / dim {dim} / {epochs} ep",
          flush=True)

    model = DistMultKGE(nN, nR, dim=dim)
    model.fit(tri, epochs=epochs, log=lambda s: print(f"  {s}", flush=True))

    np.savez_compressed(_D / "pkg_embeddings.npz",
                        E=model.E.astype(np.float32), R=model.R.astype(np.float32),
                        dim=np.int32(dim))
    json.dump(rel2id, open(_D / "pkg_rels.json", "w"))
    print(f"[primekg] saved embeddings + rels", flush=True)

    # ── Validation: recover held-out indications above negatives ──────────────
    ind_r = rel2id["indication"]
    E, Rmat = model.E, model.R
    def sc(h, t):
        return float(np.sum(E[h] * Rmat[ind_r] * E[t]))

    # disease node pool (for random negatives) from pkg_nodes.json
    nodes = json.load(open(_D / "pkg_nodes.json"))
    disease_ids = [int(i) for i, v in nodes.items() if v[0] == "disease"]

    pos_scores, rand_neg, contra_neg = [], [], []
    for (h, t) in ind_test:
        pos_scores.append(sc(h, t))
        # random disease negative for this drug
        for _ in range(1):
            rd = disease_ids[rng.randrange(len(disease_ids))]
            rand_neg.append(sc(h, rd))
        # a real contraindication for this drug, if any (the hard negative)
        cset = contra.get(h)
        if cset:
            cd = rng.choice(list(cset))
            contra_neg.append(sc(h, cd))

    auc_rand = _auc(pos_scores, rand_neg)
    auc_contra = _auc(pos_scores, contra_neg)
    metrics = {
        "held_out_indications": len(ind_test),
        "auc_indication_vs_random": round(auc_rand, 4),
        "auc_indication_vs_contraindication": round(auc_contra, 4),
        "n_contra_negatives": len(contra_neg),
        "dim": dim, "epochs": epochs,
    }
    json.dump(metrics, open(_D / "pkg_eval.json", "w"), indent=2)
    print("\n=== PrimeKG DistMult recovery ===", flush=True)
    print(f"  held-out indications:               {len(ind_test):,}", flush=True)
    print(f"  AUC  indication vs RANDOM disease:   {auc_rand:.4f}", flush=True)
    print(f"  AUC  indication vs CONTRAINDICATION: {auc_contra:.4f}   (the honest test)", flush=True)


if __name__ == "__main__":
    main()
