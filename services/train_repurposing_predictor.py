"""
Train the repurposing predictor  →  data/repurposing_predictor.pkl
═══════════════════════════════════════════════════════════════════════════════
A REAL, validated probability-of-treatment model — the rigorous answer to "how do
we decide a number." Trained on the Hetionet gold standard (755 curated
Compound-treats-Disease edges) using the leakage-free DWPC metapath features
(Rephetio; Himmelstein 2017). None of the features traverse the treat edge, so
predicting "treats" from them is leakage-free by construction.

Honest evaluation — two numbers, both reported:
  • stratified 5-fold AUC   — standard cross-validation over pairs.
  • compound-disjoint AUC   — hold out ENTIRE drugs; test only on unseen drugs.
                              This is the number that matters for repurposing: can
                              it rank indications for a drug it has never trained on?

If the held-out AUC clears the bar, the fitted model + a data-derived probability
cutoff (max-Youden on out-of-fold predictions) are saved for services to consume.

Run:  .venv312\\Scripts\\python.exe -m services.train_repurposing_predictor
"""
from __future__ import annotations

import json
import logging
from pathlib import Path

import numpy as np

logging.basicConfig(level=logging.INFO, format="%(message)s")
log = logging.getLogger("train_predictor").info

_ROOT = Path(__file__).parent.parent
_OUT = _ROOT / "data" / "repurposing_predictor.pkl"
_METRICS = _ROOT / "data" / "repurposing_predictor_metrics.json"

NEG_RATIO = 10          # negatives per positive
SEED = 0


def _load_labels(eng):
    import psycopg2
    from config import db_params
    conn = psycopg2.connect(**db_params()); cur = conn.cursor()
    cur.execute("SELECT source_id, target_id FROM hetionet_edges "
                "WHERE metaedge='CtD' AND source='hetionet_v1.0'")
    ci, di = eng.cmp_index, eng.dis_index
    pos = [(ci[c], di[d]) for c, d in cur.fetchall() if c in ci and d in di]
    conn.close()
    return pos


def _build_xy(eng, pos_idx):
    feat_names = eng.feature_names()
    mats = [eng.dwpc[n] for n in feat_names]
    n_cmp, n_dis = mats[0].shape
    pos_set = set(pos_idx)
    rng = np.random.default_rng(SEED)
    neg_idx = []
    target_neg = len(pos_idx) * NEG_RATIO
    while len(neg_idx) < target_neg:
        c = int(rng.integers(n_cmp)); d = int(rng.integers(n_dis))
        if (c, d) not in pos_set:
            neg_idx.append((c, d))
    pairs = pos_idx + neg_idx
    X = np.array([[m[c, d] for m in mats] for c, d in pairs], dtype=float)
    y = np.array([1] * len(pos_idx) + [0] * len(neg_idx), dtype=int)
    groups = np.array([c for c, d in pairs])          # compound index → group
    X = np.log1p(X)                                    # DWPC is heavily skewed
    return X, y, groups, feat_names


def _cv_auc(X, y, groups):
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import StratifiedKFold, GroupKFold
    from sklearn.metrics import roc_auc_score
    from sklearn.preprocessing import StandardScaler

    def _fit_eval(splitter, split_args):
        aucs, oof_p, oof_y = [], [], []
        for tr, te in splitter.split(*split_args):
            sc = StandardScaler().fit(X[tr])
            clf = LogisticRegression(max_iter=2000, class_weight="balanced",
                                     C=1.0).fit(sc.transform(X[tr]), y[tr])
            p = clf.predict_proba(sc.transform(X[te]))[:, 1]
            aucs.append(roc_auc_score(y[te], p))
            oof_p.extend(p); oof_y.extend(y[te])
        return float(np.mean(aucs)), float(np.std(aucs)), np.array(oof_p), np.array(oof_y)

    skf = StratifiedKFold(5, shuffle=True, random_state=SEED)
    strat_auc, strat_sd, oof_p, oof_y = _fit_eval(skf, (X, y))
    gkf = GroupKFold(5)
    grp_auc, grp_sd, _, _ = _fit_eval(gkf, (X, y, groups))
    return {"stratified_auc": round(strat_auc, 4), "stratified_sd": round(strat_sd, 4),
            "compound_disjoint_auc": round(grp_auc, 4), "compound_disjoint_sd": round(grp_sd, 4)}, oof_p, oof_y


def _youden_cutoff(p, y):
    """Probability threshold maximising (TPR − FPR) on out-of-fold predictions."""
    order = np.argsort(-p)
    P, N = int(y.sum()), int((1 - y).sum())
    tp = fp = 0; best_j = -1; best_t = 0.5
    for i in order:
        if y[i] == 1: tp += 1
        else: fp += 1
        j = tp / max(1, P) - fp / max(1, N)
        if j > best_j:
            best_j = j; best_t = float(p[i])
    return round(best_t, 4), round(best_j, 4)


def main():
    from services.metapath_features import get_features_engine
    log("Building DWPC metapath features (Hetionet)...")
    eng = get_features_engine(log=log)
    pos_idx = _load_labels(eng)
    log(f"Gold-standard positives (CtD edges in-graph): {len(pos_idx)}")
    if len(pos_idx) < 50:
        log("Too few positives to train — aborting."); return

    X, y, groups, feat_names = _build_xy(eng, pos_idx)
    log(f"Dataset: {len(y)} pairs ({int(y.sum())} pos, {int((1-y).sum())} neg), "
        f"{X.shape[1]} DWPC features")

    metrics, oof_p, oof_y = _cv_auc(X, y, groups)
    cutoff, jstat = _youden_cutoff(oof_p, oof_y)
    metrics["youden_prob_cutoff"] = cutoff
    metrics["youden_j"] = jstat
    metrics["n_pos"] = int(y.sum()); metrics["n_neg"] = int((1 - y).sum())
    metrics["features"] = feat_names

    log("\n" + "=" * 60)
    log(f"stratified 5-fold AUC      : {metrics['stratified_auc']} ± {metrics['stratified_sd']}")
    log(f"compound-disjoint AUC      : {metrics['compound_disjoint_auc']} ± {metrics['compound_disjoint_sd']}  (generalises to NEW drugs)")
    log(f"Youden probability cutoff  : {cutoff}  (J={jstat})")

    # Fit the final model on ALL data and persist (only if it actually works).
    ok = metrics["compound_disjoint_auc"] >= 0.65
    if ok:
        import joblib
        from sklearn.linear_model import LogisticRegression
        from sklearn.preprocessing import StandardScaler
        scaler = StandardScaler().fit(X)
        clf = LogisticRegression(max_iter=2000, class_weight="balanced", C=1.0)
        clf.fit(scaler.transform(X), y)
        joblib.dump({"model": clf, "scaler": scaler, "features": feat_names,
                     "metrics": metrics, "neg_ratio": NEG_RATIO}, _OUT)
        log(f"\nSaved predictor → {_OUT}")
    else:
        log(f"\ncompound-disjoint AUC {metrics['compound_disjoint_auc']} < 0.65 — NOT saving "
            "(would not be a trustworthy predictor).")

    _METRICS.write_text(json.dumps(metrics, indent=2), encoding="utf-8")
    log(f"Wrote metrics → {_METRICS}")


if __name__ == "__main__":
    main()
