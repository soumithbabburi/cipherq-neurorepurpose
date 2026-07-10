"""
External validation of the trained predictor against RepoDB's REAL failures.
═══════════════════════════════════════════════════════════════════════════════
The predictor was trained on Hetionet CtD with RANDOM negatives (easy). The honest
test is whether it still ranks RepoDB Approved pairs above RepoDB Failed pairs —
drugs that were actually TRIED in that indication and failed (hard negatives).

Maps RepoDB drugbank_id → Hetionet Compound::<DBid> (exact) and ind_name → Hetionet
Disease by name. Only Hetionet-covered pairs can be scored (213 diseases), so
coverage is partial by design — reported honestly.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

_ROOT = Path(__file__).parent.parent
_REPODB = _ROOT / "data" / "external" / "repodb" / "full.csv"
_MODEL = _ROOT / "data" / "repurposing_predictor.pkl"
_FAIL = {"Terminated", "Withdrawn", "Suspended"}


def _auc(scores, labels):
    pairs = sorted(zip(scores, labels))
    ranks = {}; i = 0
    while i < len(pairs):
        j = i
        while j < len(pairs) and pairs[j][0] == pairs[i][0]:
            j += 1
        r = (i + j - 1) / 2.0 + 1
        for k in range(i, j):
            ranks[k] = r
        i = j
    P = sum(labels); N = len(labels) - P
    if P == 0 or N == 0:
        return None
    sp = sum(ranks[k] for k, (s, l) in enumerate(pairs) if l == 1)
    return round((sp - P * (P + 1) / 2.0) / (P * N), 4)


def main():
    import joblib
    from services.metapath_features import get_features_engine
    bundle = joblib.load(_MODEL)
    clf, scaler, feats = bundle["model"], bundle["scaler"], bundle["features"]

    eng = get_features_engine(log=lambda *a, **k: None)
    mats = [eng.dwpc[f] for f in feats]

    # disease name → node index
    import psycopg2
    from config import db_params
    conn = psycopg2.connect(**db_params()); cur = conn.cursor()
    cur.execute("SELECT id, name FROM hetionet_nodes WHERE kind='Disease'")
    name2node = {}
    for nid, nm in cur.fetchall():
        if nm:
            name2node[nm.strip().lower()] = nid
    conn.close()

    def resolve_disease(ind_name):
        n = (ind_name or "").strip().lower()
        if n in name2node:
            return name2node[n]
        # loose substring either direction
        for nm, nid in name2node.items():
            if n and (n in nm or nm in n) and abs(len(n) - len(nm)) < 12:
                return nid
        return None

    rows = list(csv.DictReader(open(_REPODB, encoding="utf-8", errors="ignore")))
    ci, di = eng.cmp_index, eng.dis_index
    scored, labels, n_seen, n_resolved = [], [], 0, 0
    seen = set()
    for r in rows:
        st = r.get("status"); db = r.get("drugbank_id"); ind = r.get("ind_name")
        if st != "Approved" and st not in _FAIL:
            continue
        key = (db, ind, st)
        if not db or key in seen:
            continue
        seen.add(key); n_seen += 1
        cnode = f"Compound::{db}"
        dnode = resolve_disease(ind)
        if cnode not in ci or not dnode or dnode not in di:
            continue
        c, d = ci[cnode], di[dnode]
        x = np.log1p(np.array([[m[c, d] for m in mats]], dtype=float))
        p = float(clf.predict_proba(scaler.transform(x))[0, 1])
        scored.append(p); labels.append(1 if st == "Approved" else 0)
        n_resolved += 1

    pos = [s for s, l in zip(scored, labels) if l == 1]
    neg = [s for s, l in zip(scored, labels) if l == 0]
    auc = _auc(scored, labels)
    out = {
        "n_seen": n_seen, "n_resolved": n_resolved,
        "n_approved": len(pos), "n_failed": len(neg),
        "external_auc_approved_vs_failed": auc,
        "prob_approved_mean": round(float(np.mean(pos)), 4) if pos else None,
        "prob_approved_median": round(float(np.median(pos)), 4) if pos else None,
        "prob_failed_mean": round(float(np.mean(neg)), 4) if neg else None,
        "prob_failed_median": round(float(np.median(neg)), 4) if neg else None,
    }
    (_ROOT / "data" / "repodb_external_metrics.json").write_text(json.dumps(out, indent=2))
    print("=" * 60)
    print(f"RepoDB coverage: {n_resolved}/{n_seen} pairs mapped into Hetionet")
    print(f"  approved={len(pos)}  failed={len(neg)}")
    print(f"EXTERNAL AUC (Approved vs Failed, hard negatives): {auc}")
    print(f"  P(treats) approved: median {out['prob_approved_median']}, mean {out['prob_approved_mean']}")
    print(f"  P(treats) failed  : median {out['prob_failed_median']}, mean {out['prob_failed_mean']}")


if __name__ == "__main__":
    main()
