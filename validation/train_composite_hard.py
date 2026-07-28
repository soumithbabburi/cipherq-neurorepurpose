"""
Phase 1 of the learned composite (validation/LEARNED_COMPOSITE_DESIGN.md):
train + evaluate on the HARD task, DISEASE-DISJOINT, and decide.

Loads the Phase-0 feature matrix (data/composite_dataset.json: 8 non-leaky evidence
features + approved/failed label + disease group + the hand-composite score). Trains a
regularized logistic and a GBM, evaluates them with out-of-fold predictions under a
DISEASE-disjoint GroupKFold (the honest generalization test), calibrates isotonically
out-of-fold, and bootstraps a 95% CI on the AUROC gap vs the hand composite.

SHIP GATE: adopt the learned composite only if its disease-disjoint OOF AUROC beats the
hand composite by >= +0.01 AND the bootstrap CI on the gap excludes zero. Otherwise keep
the hand weights (an honest, likely-possible outcome given the small hard-label set).

DOES NOT touch data/composite_model.pkl (the production hook) — that wiring is Phase 2,
which must also reconcile the feature set (8 non-leaky here vs the hook's 11) and keep
prior-art additive. Writes validation/composite_hard_results.json; saves the fitted model
to data/composite_model_hard.pkl (a NON-production path) only if the gate passes.

Run:  .venv312/Scripts/python.exe -m validation.train_composite_hard
"""
import sys
import json
import datetime
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))

DATASET = ROOT / "data" / "composite_dataset.json"
OUT_METRICS = HERE / "composite_hard_results.json"
OUT_MODEL = ROOT / "data" / "composite_model_hard.pkl"   # NON-production path
SHIP_MARGIN = 0.01
SEED = 0


def _auroc(y, p):
    from sklearn.metrics import roc_auc_score
    y = np.asarray(y); p = np.asarray(p)
    if len(np.unique(y)) < 2:
        return float("nan")
    return float(roc_auc_score(y, p))


def run():
    from sklearn.linear_model import LogisticRegression
    from sklearn.ensemble import GradientBoostingClassifier
    from sklearn.model_selection import GroupKFold
    from sklearn.preprocessing import StandardScaler
    from sklearn.isotonic import IsotonicRegression

    d = json.loads(DATASET.read_text())
    rows = d["rows"]
    feats = d["features"]
    X = np.array([r["features"] for r in rows], dtype=float)
    y = np.array([int(r["label"]) for r in rows])
    groups = np.array([r["disease"] for r in rows])
    hand = np.array([float(r.get("hand_composite", 0.0)) for r in rows])
    n, npos = len(y), int(y.sum())
    ndis = len(set(groups))
    print(f"[p1] {n} rows ({npos} approved / {n-npos} failed) across {ndis} diseases", flush=True)

    n_splits = min(5, ndis)
    gkf = GroupKFold(n_splits=n_splits)

    def oof(model_fn):
        pred = np.full(n, np.nan)
        for tr, te in gkf.split(X, y, groups):
            sc = StandardScaler().fit(X[tr])
            m = model_fn().fit(sc.transform(X[tr]), y[tr])
            pred[te] = m.predict_proba(sc.transform(X[te]))[:, 1]
        return pred

    lr_oof = oof(lambda: LogisticRegression(max_iter=2000, class_weight="balanced", C=1.0))
    gb_oof = oof(lambda: GradientBoostingClassifier(random_state=SEED, max_depth=3, n_estimators=150))

    au_hand = _auroc(y, hand)
    au_lr = _auroc(y, lr_oof)
    au_gb = _auroc(y, gb_oof)
    best_name, best_oof, best_au = max(
        [("logistic", lr_oof, au_lr), ("gbm", gb_oof, au_gb)], key=lambda t: t[2])

    # bootstrap CI on (best_learned - hand) AUROC gap, paired over pairs
    rng = np.random.default_rng(SEED)
    idx = np.arange(n)
    diffs = []
    for _ in range(2000):
        bi = rng.choice(idx, size=n, replace=True)
        if len(np.unique(y[bi])) < 2:
            continue
        a_l = _auroc(y[bi], best_oof[bi])
        a_h = _auroc(y[bi], hand[bi])
        if not (np.isnan(a_l) or np.isnan(a_h)):
            diffs.append(a_l - a_h)
    lo, hi = (round(float(np.percentile(diffs, 2.5)), 4),
              round(float(np.percentile(diffs, 97.5)), 4)) if diffs else (None, None)
    gap = round(best_au - au_hand, 4)

    # out-of-fold isotonic calibration quality (ECE) of the best model
    iso = IsotonicRegression(out_of_bounds="clip").fit(best_oof, y)
    cal = iso.predict(best_oof)
    bins = np.linspace(0, 1, 11)
    ece = 0.0
    for b in range(10):
        m = (cal >= bins[b]) & (cal < bins[b + 1])
        if m.any():
            ece += (m.mean()) * abs(cal[m].mean() - y[m].mean())

    # interpretable logistic weights (full-data fit, for reporting)
    sc = StandardScaler().fit(X)
    lr_full = LogisticRegression(max_iter=2000, class_weight="balanced").fit(sc.transform(X), y)
    weights = dict(zip(feats, [round(float(w), 3) for w in lr_full.coef_[0]]))

    ship = bool(gap >= SHIP_MARGIN and lo is not None and lo > 0)
    res = {
        "run_at": datetime.datetime.now().isoformat(timespec="seconds"),
        "n": n, "n_pos": npos, "n_diseases": ndis, "cv": f"disease-disjoint GroupKFold({n_splits})",
        "auroc_hand_composite": round(au_hand, 4),
        "auroc_learned_logistic": round(au_lr, 4),
        "auroc_learned_gbm": round(au_gb, 4),
        "best_model": best_name, "best_auroc": round(best_au, 4),
        "gap_vs_hand": gap, "gap_ci95": [lo, hi],
        "calibration_ece_isotonic_oof": round(float(ece), 4),
        "logistic_weights": weights,
        "ship": ship, "ship_margin": SHIP_MARGIN,
        "verdict": ("SHIP: learned composite beats hand on disease-disjoint approved-vs-failed."
                    if ship else
                    "DO NOT SHIP: learned composite does not clear +0.01 with a CI excluding zero "
                    "on the hard task; keep the validated hand composite."),
    }
    OUT_METRICS.write_text(json.dumps(res, indent=2), encoding="utf-8")
    print(json.dumps(res, indent=2), flush=True)

    if ship:
        import joblib
        joblib.dump({"model": lr_full if best_name == "logistic" else
                     GradientBoostingClassifier(random_state=SEED, max_depth=3, n_estimators=150)
                     .fit(sc.transform(X), y),
                     "scaler": sc, "calibrator": iso, "features": feats, "metrics": res}, OUT_MODEL)
        print(f"[p1] gate PASSED -> saved {OUT_MODEL} (NON-production; Phase 2 wires it)", flush=True)
    else:
        print("[p1] gate FAILED -> nothing saved; hand composite stands.", flush=True)
    return res


if __name__ == "__main__":
    run()
