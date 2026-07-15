"""
Learn the composite weights  →  data/composite_model.pkl   (gap #4, the capstone)
═══════════════════════════════════════════════════════════════════════════════
The evidence composite was a HAND-WEIGHTED sum (0.25 target, 0.20 pathway…) — never
validated, and we proved it ranks RepoDB Failed above Approved (AUC 0.41). This
replaces the hand weights with a supervised model fitted on the Hetionet CtD gold
standard, using the composite's own sub-signals as features, with a COMPOUND-DISJOINT
split (train on some drugs, test on unseen ones — the honest number).

DWPC plausibility is deliberately EXCLUDED as a feature: it is itself trained on CtD
(AUC 0.98), so including it would just relearn DWPC. We fit the ORTHOGONAL evidence
layer (target/pathway/PPI/clinical/network/directional/proliferation/signature/
direction) and keep DWPC as its own axis. Reports learned AUC vs the hand-weighted
composite on the same held-out set.

Run:  .venv312\\Scripts\\python.exe -m services.train_composite [n_pos]
"""
from __future__ import annotations

import json
import random
import sys
from pathlib import Path

import numpy as np

_ROOT = Path(__file__).parent.parent
_OUT = _ROOT / "data" / "composite_model.pkl"
_METRICS = _ROOT / "data" / "composite_model_metrics.json"

FEATURES = ["target", "pathway", "ppi", "clinical", "indication", "regulatory",
            "network", "directional", "proliferation", "signature", "direction"]
SEED = 0


def _feats(r: dict) -> list:
    s = r.get("scores", {}) or {}
    sb = r.get("score_breakdown", {}) or {}
    return [
        float(s.get("target", 0) or 0), float(s.get("pathway", 0) or 0),
        float(s.get("ppi", 0) or 0), float(s.get("clinical", 0) or 0),
        float(s.get("indication", 0) or 0), float(s.get("regulatory", 0) or 0),
        float(sb.get("network_score", 0) or 0),
        float((r.get("directional_evidence", {}) or {}).get("signal", 0) or 0),
        float((r.get("proliferation", {}) or {}).get("score", 0) or 0),
        float((r.get("signature_reversal", {}) or {}).get("connectivity", 0) or 0),
        float((r.get("direction", {}) or {}).get("factor", 1.0) or 1.0),
    ]


def _labelled_pairs(n_pos: int):
    """Hetionet CtD positives + random negatives, as (drug_name, disease_name, label,
    compound_id)."""
    import psycopg2
    from config import db_params
    conn = psycopg2.connect(**db_params()); cur = conn.cursor()
    cur.execute("SELECT c.name AS drug, d.name AS dis, e.source_id "
                "FROM hetionet_edges e "
                "JOIN hetionet_nodes c ON c.id=e.source_id "
                "JOIN hetionet_nodes d ON d.id=e.target_id "
                "WHERE e.metaedge='CtD' AND e.source='hetionet_v1.0'")
    ctd = [(r[0], r[1], r[2]) for r in cur.fetchall() if r[0] and r[1]]
    cur.execute("SELECT id, name FROM hetionet_nodes WHERE kind='Compound' AND name IS NOT NULL")
    comps = cur.fetchall()
    cur.execute("SELECT name FROM hetionet_nodes WHERE kind='Disease' AND name IS NOT NULL")
    diss = [r[0] for r in cur.fetchall()]
    conn.close()
    rng = random.Random(SEED); rng.shuffle(ctd)
    pos = [(d, z, 1, cid) for d, z, cid in ctd[:n_pos]]
    pos_set = {(d.lower(), z.lower()) for d, z, _ in ctd}
    neg = []
    while len(neg) < n_pos * 3:
        cid, cname = rng.choice(comps); z = rng.choice(diss)
        if (cname.lower(), z.lower()) not in pos_set:
            neg.append((cname, z, 0, cid))
    return pos + neg


def main(n_pos: int = 250):
    from services.reverse_repurposing import canonical_pair_score, resolve_drug
    pairs = _labelled_pairs(n_pos)
    print(f"scoring {len(pairs)} pairs ({n_pos} CtD positives + {n_pos*3} random)…", flush=True)

    X, y, groups, hand = [], [], [], []
    resolved_cid = {}
    for k, (drug, dis, label, cid) in enumerate(pairs, 1):
        try:
            # resolve drug → chembl_id once (so direction/proliferation features fire
            # exactly as they do in production)
            if drug not in resolved_cid:
                resolved_cid[drug] = (resolve_drug(drug) or {}).get("chembl_id", "")
            chembl = resolved_cid[drug]
            r = canonical_pair_score(chembl, dis, drug_name=drug)
            if r.get("error"):
                continue
            X.append(_feats(r)); y.append(label); groups.append(cid)
            hand.append(float(r.get("composite_score", 0.0)))
        except Exception:
            pass
        if k % 50 == 0:
            print(f"  {k}/{len(pairs)} scored, {len(y)} usable", flush=True)

    X = np.array(X); y = np.array(y); groups = np.array(groups); hand = np.array(hand)
    print(f"dataset: {len(y)} pairs ({int(y.sum())} pos)", flush=True)

    from sklearn.linear_model import LogisticRegression
    from sklearn.ensemble import GradientBoostingClassifier
    from sklearn.model_selection import GroupKFold
    from sklearn.metrics import roc_auc_score
    from sklearn.preprocessing import StandardScaler

    def _cv(model_fn):
        gkf = GroupKFold(5); aucs = []
        for tr, te in gkf.split(X, y, groups):
            sc = StandardScaler().fit(X[tr])
            m = model_fn().fit(sc.transform(X[tr]), y[tr])
            p = m.predict_proba(sc.transform(X[te]))[:, 1]
            aucs.append(roc_auc_score(y[te], p))
        return round(float(np.mean(aucs)), 4), round(float(np.std(aucs)), 4)

    hand_auc = round(float(roc_auc_score(y, hand)), 4)               # baseline
    lr_auc, lr_sd = _cv(lambda: LogisticRegression(max_iter=2000, class_weight="balanced"))
    gb_auc, gb_sd = _cv(lambda: GradientBoostingClassifier(random_state=0, max_depth=3, n_estimators=150))

    # Fit final logistic (interpretable weights) on all data
    scaler = StandardScaler().fit(X)
    clf = LogisticRegression(max_iter=2000, class_weight="balanced").fit(scaler.transform(X), y)
    weights = dict(zip(FEATURES, [round(float(w), 3) for w in clf.coef_[0]]))

    metrics = {"n": int(len(y)), "n_pos": int(y.sum()),
               "auc_hand_weighted_composite": hand_auc,
               "auc_learned_logistic_compound_disjoint": lr_auc, "lr_sd": lr_sd,
               "auc_learned_gbm_compound_disjoint": gb_auc, "gb_sd": gb_sd,
               "learned_weights": weights, "features": FEATURES}
    _METRICS.write_text(json.dumps(metrics, indent=2))
    print("\n" + "=" * 60)
    print(f"hand-weighted composite AUC   : {hand_auc}   (the un-validated baseline)")
    print(f"learned logistic  (cmpd-disj) : {lr_auc} ± {lr_sd}")
    print(f"learned GBM       (cmpd-disj) : {gb_auc} ± {gb_sd}")
    print("learned logistic weights:", weights)

    if lr_auc >= hand_auc:            # only ship if it actually beats the hand weights
        import joblib
        joblib.dump({"model": clf, "scaler": scaler, "features": FEATURES,
                     "metrics": metrics}, _OUT)
        print(f"\nSaved learned composite → {_OUT}")
    else:
        print("\nLearned model did not beat the hand weights — NOT saving.")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 250)
