"""
EXPERIMENT — does the signature-reversal bonus add ranking signal on repoDB?
═══════════════════════════════════════════════════════════════════════════════
The scoring audit flagged that signature_reversal (CMap/CREEDS connectivity, a
+0.15 additive bonus when reversing / x0.85 when mimicking) has NEVER been
measured against the external gold standard. This measures it the repo's way, on
its COVERED subset (like experiment_plausibility_rerank did for DWPC):

  On repoDB (approved vs FAILED), does the CMap connectivity separate approvals
  from failures, and does blending it into the mechanistic score improve AUROC?

Ship criterion (repo culture): keep/strengthen the bonus only if a blend beats the
mechanistic score on the covered subset by >= +0.01 AUROC; else keep it as a
display/mechanism flag, not a ranking driver. Read-only. Writes
validation/reversal_ablation_results.json.
"""
import sys, json, datetime
from pathlib import Path
from scipy.stats import rankdata

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))
from validation.validate_predictions import (
    DEFAULT_REPODB, load_repodb, _resolve_drug_targets, mech_score, auroc)


def _blend_auroc(mech, sig, labels, w):
    rm, rs = rankdata(mech), rankdata(sig)
    blended = (1.0 - w) * rm + w * rs
    return auroc(list(zip(blended.tolist(), labels)))


def run(max_diseases=60):
    from validation.validate_concordance import get_engine
    from services.disease_ontology import resolve_disease, get_ppi_network
    from services.signature_engine import drug_signature
    from services.signature_reversal import reversal_score

    repo = load_repodb(str(DEFAULT_REPODB))
    all_drugs = sorted({d for dis in repo.values() for d in dis})
    dg, da, _dc = _resolve_drug_targets(get_engine(), all_drugs, "fallback")
    mapped = set(dg)
    drug_sig = {n: drug_signature([{"gene": g, "action": a} for g, a in da[n].items()]) for n in mapped}

    cand = []
    for dis, drugs in repo.items():
        lab = {d: y for d, y in drugs.items() if d in mapped}
        np_, nn = sum(1 for y in lab.values() if y == 1), sum(1 for y in lab.values() if y == 0)
        if np_ >= 1 and nn >= 1 and (np_ + nn) >= 5:
            cand.append((dis, np_ + nn))
    cand.sort(key=lambda x: -x[1]); cand = cand[:max_diseases]

    cov = []   # (mech, connectivity, label) on covered pairs
    n_all = 0
    for dis, _n in cand:
        info = resolve_disease(dis) or {}
        dgenes = [t["gene_symbol"] for t in info.get("targets", [])[:40]]
        if len(dgenes) < 5:
            continue
        dpath = info.get("pathways", [])
        ppi = get_ppi_network(dgenes[:15]) if dgenes else {}
        for d, y in [(d, y) for d, y in repo[dis].items() if d in mapped]:
            n_all += 1
            rv = reversal_score(d, dis)
            if not rv.get("covered"):
                continue
            g = sorted(dg[d])
            base = mech_score(g, drug_sig[d], dgenes, dpath, ppi, "aware")
            cov.append((base, float(rv.get("connectivity") or 0.0), y))

    def fmt(x):
        return None if x is None else round(x, 4)

    res = {"run_at": datetime.datetime.now().isoformat(timespec="seconds"),
           "evaluable_pairs": n_all, "covered_pairs": len(cov)}
    if len(cov) >= 20:
        mech = [m for m, _, _ in cov]; sig = [s for _, s, _ in cov]; ys = [y for _, _, y in cov]
        au_m = auroc(list(zip(mech, ys)))
        au_s = auroc(list(zip(sig, ys)))
        sweep = {f"w_{w:.2f}": fmt(_blend_auroc(mech, sig, ys, w)) for w in (0.0, 0.25, 0.5)}
        best_w = max((0.0, 0.25, 0.5), key=lambda w: (_blend_auroc(mech, sig, ys, w) or 0))
        best = _blend_auroc(mech, sig, ys, best_w)
        lift = None if (best is None or au_m is None) else round(best - au_m, 4)
        res.update({
            "approved": sum(ys), "failed": len(ys) - sum(ys),
            "mechanism_auroc": fmt(au_m), "connectivity_only_auroc": fmt(au_s),
            "blend_sweep": sweep, "best_weight": best_w, "best_blend_auroc": fmt(best), "best_lift": lift,
            "verdict": ("WIRE/STRENGTHEN: connectivity adds ranking signal on the covered subset."
                        if (lift is not None and lift >= 0.01 and best_w > 0) else
                        "KEEP AS DISPLAY/FLAG: connectivity does not beat mechanism by the >=0.01 bar "
                        "on the covered subset; not a ranking driver."),
        })
    else:
        res["verdict"] = f"Too few covered pairs ({len(cov)}) to judge."

    (HERE / "reversal_ablation_results.json").write_text(json.dumps(res, indent=2), encoding="utf-8")
    print(json.dumps(res, indent=2))
    return res


if __name__ == "__main__":
    run(int(sys.argv[1]) if len(sys.argv) > 1 else 60)
