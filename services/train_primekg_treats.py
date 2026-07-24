"""
Direction-aware TREATS classifier over PrimeKG (Phase A — the precise generator).
════════════════════════════════════════════════════════════════════════════════
The unsupervised DistMult on PrimeKG captured drug-disease RELATEDNESS but not DIRECTION
(indication vs contraindication AUC 0.57) — it surfaced a drug's toxicities as "candidates".
This fixes that the way TxGNN does conceptually: SUPERVISED metric learning that separates
indications from contraindications, using PrimeKG's real contraindication edges as HARD
negatives (and random diseases as easy negatives), on top of the DistMult node embeddings
as features. No GPU/DGL needed.

Honest evaluation: COMPOUND-DISJOINT split (whole drugs held out), so the AUC measures
generalisation to UNSEEN drugs (zero-shot-ish), and the headline metric is
indication-vs-CONTRAINDICATION (the hard task the plain model failed), not vs random.

Output:
  data/primekg/pkg_treats.pkl    {model, feature_spec, meta}
  data/primekg/pkg_treats_eval.json
Run:  python -m services.train_primekg_treats
"""
from __future__ import annotations

import json
import random
from pathlib import Path

import numpy as np

_D = Path(__file__).parent.parent / "data" / "primekg"


def _feat(E, d_idx, z_idx):
    """Pair features from node embeddings: [Ed, Ez, Ed*Ez, |Ed-Ez|]."""
    ed = E[d_idx]; ez = E[z_idx]
    return np.concatenate([ed, ez, ed * ez, np.abs(ed - ez)], axis=-1)


def main(seed: int = 42, easy_neg_per_pos: int = 2, holdout_drugs: float = 0.15):
    rng = random.Random(seed)
    np.random.seed(seed)
    emb = np.load(_D / "pkg_embeddings.npz"); E = emb["E"]
    labels = json.load(open(_D / "pkg_labels.json"))
    nodes = json.load(open(_D / "pkg_nodes.json"))
    disease_ids = [int(i) for i, v in nodes.items() if v[0] == "disease"]

    drugs = [int(d) for d, rels in labels.items() if rels.get("indication")]
    rng.shuffle(drugs)
    n_hold = int(len(drugs) * holdout_drugs)
    test_drugs = set(drugs[:n_hold]); train_drugs = set(drugs[n_hold:])
    print(f"[treats] {len(drugs)} drugs w/ indications -> train {len(train_drugs)} / test {len(test_drugs)}",
          flush=True)

    def build(drug_set):
        X, y, meta = [], [], []
        for d in drug_set:
            rels = labels[str(d)]
            pos = rels.get("indication", [])
            hard = rels.get("contraindication", [])
            for z in pos:
                X.append(_feat(E, d, z)); y.append(1); meta.append((d, z, "ind"))
            for z in hard:
                X.append(_feat(E, d, z)); y.append(0); meta.append((d, z, "contra"))
            # easy random negatives
            for _ in range(len(pos) * easy_neg_per_pos):
                z = disease_ids[rng.randrange(len(disease_ids))]
                X.append(_feat(E, d, z)); y.append(0); meta.append((d, z, "rand"))
        return np.array(X, dtype=np.float32), np.array(y), meta

    Xtr, ytr, _ = build(train_drugs)
    Xte, yte, mte = build(test_drugs)
    print(f"[treats] train {Xtr.shape}  test {Xte.shape}  (pos rate tr={ytr.mean():.3f})", flush=True)

    from sklearn.ensemble import HistGradientBoostingClassifier
    from sklearn.metrics import roc_auc_score
    clf = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.08,
                                         max_depth=6, l2_regularization=1.0,
                                         random_state=seed)
    clf.fit(Xtr, ytr)

    p_te = clf.predict_proba(Xte)[:, 1]
    # overall (ind vs all-neg)
    auc_all = roc_auc_score(yte, p_te)
    # honest hard task: ind vs CONTRA only (drop random negatives)
    idx_ic = [i for i, m in enumerate(mte) if m[2] in ("ind", "contra")]
    y_ic = np.array([1 if mte[i][2] == "ind" else 0 for i in idx_ic])
    auc_ic = roc_auc_score(y_ic, p_te[idx_ic]) if len(set(y_ic)) == 2 else float("nan")
    # ind vs random only
    idx_ir = [i for i, m in enumerate(mte) if m[2] in ("ind", "rand")]
    y_ir = np.array([1 if mte[i][2] == "ind" else 0 for i in idx_ir])
    auc_ir = roc_auc_score(y_ir, p_te[idx_ir]) if len(set(y_ir)) == 2 else float("nan")

    import joblib
    joblib.dump({"model": clf, "dim": int(E.shape[1]),
                 "feature": "[Ed,Ez,Ed*Ez,|Ed-Ez|]"}, _D / "pkg_treats.pkl")
    metrics = {
        "compound_disjoint": True,
        "test_drugs": len(test_drugs),
        "auc_ind_vs_all_neg": round(float(auc_all), 4),
        "auc_ind_vs_contraindication": round(float(auc_ic), 4),
        "auc_ind_vs_random": round(float(auc_ir), 4),
        "baseline_distmult_ind_vs_contra": 0.5703,
        "n_test_pairs": int(len(yte)),
    }
    json.dump(metrics, open(_D / "pkg_treats_eval.json", "w"), indent=2)
    print("\n=== TREATS classifier (compound-disjoint / zero-shot to unseen drugs) ===", flush=True)
    print(f"  AUC  indication vs CONTRAINDICATION : {auc_ic:.4f}   (was 0.5703 for plain DistMult)", flush=True)
    print(f"  AUC  indication vs random           : {auc_ir:.4f}", flush=True)
    print(f"  AUC  indication vs all negatives     : {auc_all:.4f}", flush=True)


if __name__ == "__main__":
    main()
