"""
EXPERIMENT: can broader + hard negative sampling improve the PrimeKG universe ranker
without regressing the direction-aware AUC? (No GPU/GNN required.)
════════════════════════════════════════════════════════════════════════════════════
Baseline (validate_primekg_ranking): shipped cascade R@10 0.276, median rank 36. The
treats classifier is trained with only 2 random negatives per positive + contraindications
as hard negatives — it barely sees the 17k-disease space it must rank, so its far tail is
noisy. Hypothesis: exposing it to (a) many more random negatives and (b) HARD negatives
mined from high DistMult relatedness that are NOT labeled indications (the exact
false-positive shape) will sharpen the tail and lift ranking recall, while keeping the
ind-vs-contra discrimination that the contraindication gate depends on.

Trains to a SEPARATE artifact (pkg_treats_v2.pkl) — never overwrites the shipped model.
Reports ind-vs-contra AUC (must stay ~0.98) AND held-out ranking vs the shipped cascade.
Run:  .venv312/Scripts/python.exe -m validation.experiment_primekg_ranker_v2

RESULT (2026-07-24): HYPOTHESIS REFUTED. Broad + hard negatives REGRESSED everything —
cascade R@10 0.276 -> 0.021, ind-vs-contra AUC 0.984 -> 0.911. A broad-random-only variant
(no relatedness-mined hard negatives) failed the same way (R@10 0.015, AUC 0.928), so it is
the class imbalance itself, not the mining, that breaks a tree-based classifier here: the
DistMult relatedness features that drive recall correlate with true indications, so training
against "related-but-unlabeled" diseases as negatives suppresses the very signal ranking
needs. The shipped 2-negative config is near this shallow model's ceiling. CONCLUSION: a
better universe ranker needs GRAPH STRUCTURE (message passing / a real GNN), not more
negative engineering on DistMult-embedding features — which is exactly what the roadmap
predicted. Kept as a documented negative result; the shipped model is unchanged.
"""
from __future__ import annotations

import json
import random
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from services import primekg_predictor as pkg                      # noqa: E402
from validation.validate_primekg_ranking import (                  # noqa: E402
    _reproduce_split, _distmult_all, _filtered_ranks, _summ, KS, _POOL)

_D = ROOT / "data" / "primekg"


def _feat_rows(E, d_idx, z_idxs):
    ed = E[d_idx]; Ez = E[np.asarray(z_idxs)]
    n = len(z_idxs)
    return np.concatenate(
        [np.broadcast_to(ed, (n, ed.shape[0])), Ez, ed * Ez, np.abs(ed - Ez)], axis=1)


def _build(E, s, drug_set, labels, disease_ids, rng,
           easy_neg=15, hard_neg=10, hard_pool=200):
    """Rows for a drug set. Positives = indications; negatives = contraindications (hard,
    curated) + DistMult-mined hard negatives (look related, not labeled) + random."""
    dis_arr = disease_ids
    pos_of = {int(z): i for i, z in enumerate(dis_arr)}
    X, y, meta = [], [], []
    for d in drug_set:
        rels = labels[str(d)]
        pos = [z for z in rels.get("indication", []) if z in pos_of]
        if not pos:
            continue
        contra = [z for z in rels.get("contraindication", []) if z in pos_of]
        offl = set(rels.get("off-label use", []))
        exclude = set(pos) | set(contra) | offl
        # positives
        for z in pos:
            X.append(_feat_rows(E, d, [pos_of[z]])[0]); y.append(1); meta.append("ind")
        # curated hard negatives (contraindications)
        for z in contra:
            X.append(_feat_rows(E, d, [pos_of[z]])[0]); y.append(0); meta.append("contra")
        # DistMult-mined hard negatives: most-related diseases that are NOT labeled indications
        dm = _distmult_all(s, d, dis_arr)
        top = np.argsort(-dm)[:hard_pool]
        mined = [j for j in top if int(dis_arr[j]) not in exclude][: hard_neg * max(1, len(pos))]
        for j in mined:
            X.append(_feat_rows(E, d, [int(j)])[0]); y.append(0); meta.append("hard")
        # easy random negatives
        for _ in range(len(pos) * easy_neg):
            j = rng.randrange(len(dis_arr))
            if int(dis_arr[j]) in exclude:
                continue
            X.append(_feat_rows(E, d, [int(j)])[0]); y.append(0); meta.append("rand")
    return np.asarray(X, dtype=np.float32), np.asarray(y), meta


def _rank_eval(clf, s, disease_ids, ev4):
    """Held-out ranking with an arbitrary classifier (mirrors the harness, incl. cascade)."""
    E = s["E"]; n_dis = len(disease_ids)
    R = {"cascade_v2": [], "treats_full_v2": []}
    for d, pos_local in ev4:
        Ez = E[disease_ids]
        feats = np.concatenate(
            [np.broadcast_to(E[d], (n_dis, E.shape[1])), Ez, E[d] * Ez, np.abs(E[d] - Ez)], axis=1)
        tf = clf.predict_proba(feats)[:, 1]
        dm = _distmult_all(s, d, disease_ids)
        pool = np.argsort(-dm)[:_POOL]
        casc = np.full(n_dis, -np.inf); casc[pool] = tf[pool]
        R["cascade_v2"].extend(_filtered_ranks(casc, pos_local, pos_local))
        R["treats_full_v2"].extend(_filtered_ranks(tf, pos_local, pos_local))
    return {k: _summ(v) for k, v in R.items()}


def main(seed: int = 42):
    from sklearn.ensemble import HistGradientBoostingClassifier
    from sklearn.metrics import roc_auc_score
    s = pkg._load()
    if s is None or s.get("treats") is None:
        print("PrimeKG artifacts not present."); return 1
    E = s["E"]; labels = s["labels"]; disease_ids = s["disease_ids"]
    pos_of = {int(z): i for i, z in enumerate(disease_ids)}
    rng = random.Random(seed)

    test_drugs = _reproduce_split(labels)
    train_drugs = set(int(d) for d, r in labels.items() if r.get("indication")) - test_drugs
    ev4 = []
    for d in test_drugs:
        pl = {pos_of[z] for z in labels[str(d)].get("indication", []) if z in pos_of}
        if pl:
            ev4.append((d, pl))

    t0 = time.time()
    print(f"[v2] building training rows (broad + hard negatives)...", flush=True)
    Xtr, ytr, mtr = _build(E, s, train_drugs, labels, disease_ids, rng)
    print(f"[v2] train {Xtr.shape} pos_rate={ytr.mean():.3f} "
          f"(contra={mtr.count('contra')} hard={mtr.count('hard')} rand={mtr.count('rand')}) "
          f"{time.time()-t0:.0f}s", flush=True)

    clf = HistGradientBoostingClassifier(max_iter=400, learning_rate=0.08, max_depth=6,
                                         l2_regularization=1.0, random_state=seed)
    clf.fit(Xtr, ytr)
    print(f"[v2] trained {time.time()-t0:.0f}s", flush=True)

    # direction AUC on held-out drugs (ind vs contra) — must not regress
    Xte, yte, mte = _build(E, s, test_drugs, labels, disease_ids, rng, easy_neg=2, hard_neg=0)
    pte = clf.predict_proba(Xte)[:, 1]
    ic = [i for i, m in enumerate(mte) if m in ("ind", "contra")]
    yic = np.array([1 if mte[i] == "ind" else 0 for i in ic])
    auc_ic = float(roc_auc_score(yic, pte[ic])) if len(set(yic)) == 2 else float("nan")

    print(f"[v2] ranking eval over {len(ev4)} held-out drugs...", flush=True)
    rank = _rank_eval(clf, s, disease_ids, ev4)

    base = json.load(open(ROOT / "validation" / "primekg_ranking_results.json"))["rankers"]
    out = {
        "auc_ind_vs_contra_v2": round(auc_ic, 4),
        "auc_ind_vs_contra_shipped": 0.9839,
        "baseline_cascade": base["cascade"],
        "v2_cascade": rank["cascade_v2"],
        "v2_treats_full": rank["treats_full_v2"],
        "train_rows": int(Xtr.shape[0]),
        "elapsed_s": round(time.time() - t0, 1),
    }
    json.dump(out, open(ROOT / "validation" / "primekg_ranker_v2_results.json", "w"), indent=2)

    def line(name, m):
        return f"{name:16} " + " ".join(f"{m[f'recall@{k}']:<5}" for k in KS) + \
               f"  MRR {m['mrr']:<6} med {m['median_rank']}"
    print("\n=== v2 negative-sampling experiment (held-out ranking) ===", flush=True)
    print(f"ind-vs-contra AUC:  shipped 0.9839   ->   v2 {auc_ic:.4f}", flush=True)
    print(f"{'':16} " + " ".join(f"R@{k:<4}" for k in KS), flush=True)
    print(line("baseline_cascade", base["cascade"]), flush=True)
    print(line("v2_cascade", rank["cascade_v2"]), flush=True)
    print(line("v2_treats_full", rank["treats_full_v2"]), flush=True)

    win = rank["cascade_v2"]["recall@10"] > base["cascade"]["recall@10"] and auc_ic >= 0.97
    print(f"\nVERDICT: {'v2 beats baseline (candidate for promotion)' if win else 'no clear win — keep shipped model'}",
          flush=True)
    if win:
        import joblib
        joblib.dump({"model": clf, "dim": int(E.shape[1]),
                     "feature": "[Ed,Ez,Ed*Ez,|Ed-Ez|]"}, _D / "pkg_treats_v2.pkl")
        print(f"wrote {_D / 'pkg_treats_v2.pkl'} (NOT yet promoted to pkg_treats.pkl)", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
