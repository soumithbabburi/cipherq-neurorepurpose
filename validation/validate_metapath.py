"""
Metapath/DWPC validation on repoDB — the DECISION GATE
═══════════════════════════════════════════════════════════════════════════
Same-pairs head-to-head on the EXTERNAL repoDB benchmark: do Rephetio-style
metapath (DWPC) features beat the current mechanistic score (~0.73)?

For every repoDB (drug, disease) pair that maps to a Hetionet node AND has drug
targets, we compute BOTH:
  • mechanistic score (the live engine's biology score, reused verbatim)
  • metapath DWPC features (biology-only ⇒ leakage-free)
then compare on identical pairs:
  - each metapath alone (UNSUPERVISED — the fair comparison vs the mechanistic
    score, which is also unsupervised)
  - a grouped-CV logistic regression over all metapaths (the Rephetio model;
    supervised, so reported separately)
  - rank-average ensemble of mechanistic + metapath-LR
  - label-shuffle negative control

Decision: integrate only if metapath clearly beats the mechanistic score here.

Usage:  python validation/validate_metapath.py [max_diseases]
"""
import sys
import json
import datetime
import numpy as np
from pathlib import Path
from collections import defaultdict

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))

from validation.validate_predictions import (
    load_repodb, _resolve_drug_targets, mech_score, auroc, bedroc, DEFAULT_REPODB)
from validation.validate_concordance import get_engine


def run(max_diseases: int = 80):
    log_lines = []

    def log(m=""):
        print(m, flush=True); log_lines.append(str(m))

    log("=" * 72)
    log("METAPATH/DWPC VALIDATION on repoDB  -  decision gate vs mechanistic 0.73")
    log(f"Run: {datetime.datetime.now().isoformat(timespec='seconds')}")
    log("=" * 72)

    from services.disease_ontology import resolve_disease, get_ppi_network
    from services.signature_engine import drug_signature
    from services import kg_score
    from services.metapath_features import get_features_engine

    repo = load_repodb(DEFAULT_REPODB)
    all_drugs = sorted({d for dis in repo.values() for d in dis})
    drug_genes, drug_actions, drug_chembl = _resolve_drug_targets(get_engine(), all_drugs, "fallback")
    mapped = set(drug_genes)
    drug_sig = {n: drug_signature([{"gene": g, "action": a} for g, a in drug_actions[n].items()])
                for n in mapped}
    log(f"\nrepoDB drugs mapped to ChEMBL targets: {len(mapped):,}")

    if not kg_score.is_available():
        log("ERROR: kg_score crosswalk unavailable."); sys.exit(1)
    log("Building metapath DWPC engine…")
    mp = get_features_engine(log=log)

    # diseases mappable to Hetionet, OT-resolvable, with both classes
    cand = []
    for dis, drugs in repo.items():
        dnode = kg_score.resolve_disease(dis)
        if not dnode:
            continue
        labeled = {d: y for d, y in drugs.items() if d in mapped}
        npos = sum(1 for y in labeled.values() if y == 1)
        nneg = len(labeled) - npos
        if npos >= 1 and nneg >= 1 and len(labeled) >= 5:
            cand.append((dis, dnode, len(labeled)))
    cand.sort(key=lambda x: -x[2])
    cand = cand[:max_diseases]
    log(f"Hetionet-mappable evaluable diseases: {len(cand)}")

    rows = []   # (mech, featvec-dict, label, disease)
    for dis, dnode, _ in cand:
        dinfo = resolve_disease(dis) or {}
        dgenes = [t["gene_symbol"] for t in dinfo.get("targets", [])[:40]]
        if len(dgenes) < 5:
            continue
        dpath = dinfo.get("pathways", [])
        ppi = get_ppi_network(dgenes[:15])
        for d, y in repo[dis].items():
            if d not in mapped:
                continue
            cnode = kg_score.resolve_compound(name=d)
            if not cnode:
                continue
            feats = mp.features(cnode, dnode)
            if not feats:
                continue
            g = sorted(drug_genes[d])
            mech = mech_score(g, drug_sig[d], dgenes, dpath, ppi, "aware")
            rows.append((mech, feats, y, dis))

    n_pos = sum(1 for r in rows if r[2] == 1)
    log(f"\nMetapath-mappable pairs: {len(rows)} (approved {n_pos} / failed {len(rows) - n_pos})")
    if len(rows) < 40 or n_pos == 0 or n_pos == len(rows):
        log("Too few mappable pairs for a reliable gate."); sys.exit(1)

    fnames = mp.feature_names()
    X = np.array([[r[1].get(f, 0.0) for f in fnames] for r in rows])
    Xl = np.log1p(X)
    y = np.array([r[2] for r in rows])
    mech = np.array([r[0] for r in rows])
    groups = [r[3] for r in rows]

    au_mech = auroc(list(zip(mech.tolist(), y.tolist())))
    per_mp = {}
    for i, f in enumerate(fnames):
        per_mp[f] = round(auroc(list(zip(X[:, i].tolist(), y.tolist()))) or 0.5, 4)

    # grouped-CV logistic regression over all metapaths (Rephetio model)
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import GroupKFold
    from sklearn.preprocessing import StandardScaler
    from scipy.stats import rankdata

    proba = np.zeros(len(y))
    gkf = GroupKFold(n_splits=min(5, len(set(groups))))
    for tr, te in gkf.split(Xl, y, groups):
        sc = StandardScaler().fit(Xl[tr])
        lr = LogisticRegression(max_iter=2000, C=1.0, class_weight="balanced")
        lr.fit(sc.transform(Xl[tr]), y[tr])
        proba[te] = lr.predict_proba(sc.transform(Xl[te]))[:, 1]
    au_lr = auroc(list(zip(proba.tolist(), y.tolist())))

    ens = rankdata(mech) + rankdata(proba)
    au_ens = auroc(list(zip(ens.tolist(), y.tolist())))

    # negative control: shuffle labels, redo CV
    rng = np.random.default_rng(7)
    ysh = y.copy(); rng.shuffle(ysh)
    proba_sh = np.zeros(len(y))
    for tr, te in gkf.split(Xl, ysh, groups):
        sc = StandardScaler().fit(Xl[tr])
        lr = LogisticRegression(max_iter=2000, C=1.0, class_weight="balanced").fit(sc.transform(Xl[tr]), ysh[tr])
        proba_sh[te] = lr.predict_proba(sc.transform(Xl[te]))[:, 1]
    au_ctrl = auroc(list(zip(proba_sh.tolist(), ysh.tolist())))

    # full-data coefficients (interpretability)
    sc = StandardScaler().fit(Xl)
    lr_full = LogisticRegression(max_iter=2000, C=1.0, class_weight="balanced").fit(sc.transform(Xl), y)
    coefs = sorted(zip(fnames, lr_full.coef_[0].tolist()), key=lambda x: -abs(x[1]))

    log("\nPer-metapath UNSUPERVISED AUROC (fair vs mechanistic):")
    for f, a in sorted(per_mp.items(), key=lambda x: -x[1]):
        log(f"  {f:16} {a}")
    log(f"\nMechanistic (same pairs):        {au_mech:.4f}")
    log(f"Metapath CV-LR (Rephetio model): {au_lr:.4f}")
    log(f"Ensemble mech + metapath-LR:     {au_ens:.4f}")
    log(f"Negative control (shuffled CV):  {au_ctrl:.4f}  (expect ~0.5)")
    log("Top metapath coefficients:")
    for f, c in coefs[:6]:
        log(f"  {f:16} {c:+.3f}")

    best_mp = max(per_mp.values())
    lr_lift = round((au_lr or 0) - (au_mech or 0), 4)
    ens_lift = round((au_ens or 0) - (au_mech or 0), 4)
    verdict = ("INTEGRATE — metapath model clearly beats the mechanistic score."
               if (au_lr or 0) - (au_mech or 0) > 0.02 or ens_lift > 0.02 else
               "DO NOT integrate — metapath does not clearly beat the mechanistic score on the external benchmark.")
    sev = "INFO" if "INTEGRATE —" in verdict else "WARN"

    result = {
        "run_at": datetime.datetime.now().isoformat(timespec="seconds"),
        "benchmark": "repoDB external, Hetionet-mappable subset (same-pairs head-to-head)",
        "mappable_pairs": len(rows),
        "approved": n_pos, "failed": len(rows) - n_pos,
        "per_metapath_unsupervised_auroc": per_mp,
        "best_single_metapath_auroc": round(best_mp, 4),
        "mechanistic_auroc": round(au_mech, 4),
        "metapath_cv_lr_auroc": round(au_lr, 4),
        "ensemble_auroc": round(au_ens, 4),
        "cv_lr_lift_over_mech": lr_lift,
        "ensemble_lift_over_mech": ens_lift,
        "negative_control_auroc": round(au_ctrl, 4),
        "top_coefficients": [{"metapath": f, "coef": round(c, 3)} for f, c in coefs[:8]],
        "verdict": verdict,
        "findings": [{
            "id": "MP-01", "severity": sev,
            "title": "Metapath/DWPC vs mechanistic on external repoDB (decision gate)",
            "detail": f"On {len(rows)} Hetionet-mappable repoDB pairs: best single metapath "
                      f"{best_mp:.3f} (unsupervised), mechanistic {au_mech:.3f}, metapath CV-LR {au_lr:.3f} "
                      f"(lift {lr_lift}), mech+metapath ensemble {au_ens:.3f} (lift {ens_lift}); "
                      f"negative control {au_ctrl:.3f}. {verdict}",
            "capa": "Integrate the metapath signal only on a positive, reproducible external lift; "
                    "otherwise keep the mechanistic score (simpler) and revisit features."},
        ],
    }
    (HERE / "metapath_results.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    (HERE / "metapath_run.log").write_text("\n".join(log_lines), encoding="utf-8")
    log(f"\n  {verdict}")
    log(f"  Wrote: {HERE / 'metapath_results.json'}")
    return result


if __name__ == "__main__":
    md = int(sys.argv[1]) if len(sys.argv) > 1 else 80
    try:
        run(md)
    except Exception:
        import traceback; traceback.print_exc(); sys.exit(1)
