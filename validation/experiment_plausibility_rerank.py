"""
EXPERIMENT — does the validated DWPC plausibility model improve the RANKING?
═══════════════════════════════════════════════════════════════════════════
The DWPC/metapath plausibility model (services/repurposing_predictor.py) reports
AUC 0.98 — but that is treats-vs-RANDOM (random non-edges as negatives). It is
computed on every reverse-screen pair and shown as a `plausibility` field, but it
does NOT enter the ranking sort key. This experiment asks the honest question,
on the SAME external gold standard the mechanistic score is validated against:

  On repoDB (Approved vs tried-and-FAILED — hard negatives), does blending the
  DWPC P(treats) into the mechanistic score improve AUROC over mechanism alone?

This mirrors the existing KG-embedding ensemble experiment in validate_predictions
(finding KGE-04, which said "do NOT wire in" because it degraded the external
benchmark). Same method, different signal: the DWPC plausibility model.

Ship criterion (matches the repo's culture — see train_composite.py, KGE-04):
  wire the blend into the ranking ONLY if the ensemble AUROC beats mechanism
  alone on this external benchmark by a non-trivial margin on the covered subset.
Otherwise report the negative result and leave the ranking unchanged.

Read-only. Reuses validate_predictions' machinery so the mechanism score and the
pair selection are IDENTICAL to the shipped validator. Writes
validation/plausibility_rerank_results.json.

Usage:  python validation/experiment_plausibility_rerank.py [max_diseases]
"""
import sys
import json
import datetime
import random
from pathlib import Path

import numpy as np
from scipy.stats import rankdata

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))

from validation.validate_predictions import (
    DEFAULT_REPODB, load_repodb, _resolve_drug_targets, mech_score, auroc,
    average_precision, enrichment_factor,
)


def _blend_auroc(mech, plaus, labels, w):
    """Rank-average blend with weight w on plausibility (0..1). w=0 -> mech only."""
    rm = rankdata(mech)
    rp = rankdata(plaus)
    blended = (1.0 - w) * rm + w * rp
    return auroc(list(zip(blended.tolist(), labels)))


def run(max_diseases=60):
    log_lines = []

    def log(msg=""):
        print(msg)
        log_lines.append(str(msg))

    log("=" * 72)
    log("EXPERIMENT — DWPC plausibility as a RANK signal on repoDB (Approved vs failed)")
    log(f"Run: {datetime.datetime.now().isoformat(timespec='seconds')}")
    log("=" * 72)

    from validation.validate_concordance import get_engine
    from services.disease_ontology import resolve_disease, get_ppi_network
    from services.signature_engine import drug_signature
    from services.repurposing_predictor import plausibility

    repo = load_repodb(str(DEFAULT_REPODB))
    all_drugs = sorted({d for dis in repo.values() for d in dis})
    eng = get_engine()
    drug_genes, drug_actions, drug_chembl = _resolve_drug_targets(eng, all_drugs, "fallback")
    mapped = set(drug_genes)
    log(f"repoDB drugs mapped to ChEMBL w/ targets: {len(mapped):,}/{len(all_drugs):,}")

    drug_sig = {name: drug_signature([{"gene": g, "action": a}
                                      for g, a in drug_actions[name].items()])
                for name in mapped}

    cand_dis = []
    for dis, drugs in repo.items():
        labeled = {d: y for d, y in drugs.items() if d in mapped}
        npos = sum(1 for y in labeled.values() if y == 1)
        nneg = sum(1 for y in labeled.values() if y == 0)
        if npos >= 1 and nneg >= 1 and (npos + nneg) >= 5:
            cand_dis.append((dis, npos + nneg))
    cand_dis.sort(key=lambda x: -x[1])
    cand_dis = cand_dis[:max_diseases]

    # Collect, per pair: mechanism score, plausibility P (or None), label
    all_pairs = []            # (mech, label) — all evaluable pairs (full benchmark)
    cov = []                  # (mech, plaus, label) — plausibility-COVERED subset only
    resolved_ok = 0
    for dis, _tot in cand_dis:
        dinfo = resolve_disease(dis) or {}
        dgenes = [t["gene_symbol"] for t in dinfo.get("targets", [])[:40]]
        dpath = dinfo.get("pathways", [])
        if len(dgenes) < 5:
            continue
        resolved_ok += 1
        ppi = get_ppi_network(dgenes[:15]) if dgenes else {}
        for d, y in [(d, y) for d, y in repo[dis].items() if d in mapped]:
            g = sorted(drug_genes[d])
            sa = mech_score(g, drug_sig[d], dgenes, dpath, ppi, "aware")
            all_pairs.append((sa, y))
            pr = plausibility(d, dis)
            if pr is not None and pr.get("probability") is not None:
                cov.append((sa, float(pr["probability"]), y))

    n_pos = sum(1 for _, y in all_pairs if y == 1)
    log(f"Diseases resolved: {resolved_ok} | evaluable pairs: {len(all_pairs):,} "
        f"(approved {n_pos:,} / failed {len(all_pairs) - n_pos:,})")
    log(f"Plausibility-COVERED pairs (both signals present): {len(cov):,} "
        f"({100.0 * len(cov) / max(1, len(all_pairs)):.0f}% of evaluable)")

    def fmt(x):
        return None if x is None else round(x, 4)

    results = {
        "run_at": datetime.datetime.now().isoformat(timespec="seconds"),
        "gold_standard": "repoDB Approved vs failed (Terminated/Withdrawn/Suspended) — hard negatives",
        "note": "DWPC plausibility (repurposing_predictor.py, AUC 0.98 vs RANDOM) tested as a "
                "rank signal on the harder external approved-vs-failed benchmark, on the covered subset.",
        "full_benchmark": {
            "pairs": len(all_pairs),
            "mechanistic_auroc": fmt(auroc(all_pairs)),
        },
    }

    if len(cov) >= 20:
        mech = [m for m, _, _ in cov]
        plaus = [p for _, p, _ in cov]
        ys = [y for _, _, y in cov]
        cov_pos = sum(ys)
        au_m = auroc(list(zip(mech, ys)))
        au_p = auroc(list(zip(plaus, ys)))
        # rank-average ensemble (identical method to KGE-04) + a weight sweep
        sweep = {f"w_plaus_{w:.2f}": fmt(_blend_auroc(mech, plaus, ys, w))
                 for w in (0.0, 0.25, 0.5, 0.75, 1.0)}
        au_e = _blend_auroc(mech, plaus, ys, 0.5)     # the canonical 50/50 rank-average
        best_w = max((0.0, 0.25, 0.5, 0.75, 1.0),
                     key=lambda w: (_blend_auroc(mech, plaus, ys, w) or 0))
        best_e = _blend_auroc(mech, plaus, ys, best_w)
        lift = None if (au_e is None or au_m is None) else round(au_e - au_m, 4)
        best_lift = None if (best_e is None or au_m is None) else round(best_e - au_m, 4)
        results["covered_subset"] = {
            "pairs": len(cov),
            "approved": cov_pos,
            "failed": len(cov) - cov_pos,
            "mechanistic_auroc": fmt(au_m),
            "plausibility_only_auroc": fmt(au_p),
            "ensemble_50_50_auroc": fmt(au_e),
            "ensemble_lift_over_mechanistic": lift,
            "weight_sweep": sweep,
            "best_weight_on_plausibility": best_w,
            "best_ensemble_auroc": fmt(best_e),
            "best_lift": best_lift,
        }
        SHIP_MARGIN = 0.01
        ship = (best_lift is not None and best_lift >= SHIP_MARGIN and best_w > 0.0)
        results["verdict"] = {
            "ship": ship,
            "ship_margin_required": SHIP_MARGIN,
            "recommendation": (
                f"WIRE IN: blend plausibility into the rank at weight ~{best_w:.2f} "
                f"(best lift {best_lift} over mechanism, covered subset)."
                if ship else
                "DO NOT WIRE IN: plausibility does not beat mechanism alone on the external "
                "approved-vs-failed benchmark by the required margin. Leave the ranking "
                "unchanged; keep plausibility as a display/recall field (like KGE-04)."),
        }
        log("\nCOVERED SUBSET (both signals present):")
        log(f"  mechanism      AUROC {fmt(au_m)}")
        log(f"  plausibility   AUROC {fmt(au_p)}  (its 0.98-vs-random does NOT define this task)")
        log(f"  ensemble 50/50 AUROC {fmt(au_e)}  (lift {lift})")
        log(f"  weight sweep: {sweep}")
        log(f"  best weight on plausibility: {best_w} -> AUROC {fmt(best_e)} (lift {best_lift})")
        log(f"\n  VERDICT: {'SHIP' if ship else 'DO NOT SHIP'} — {results['verdict']['recommendation']}")
    else:
        results["covered_subset"] = {"pairs": len(cov),
                                     "note": "too few covered pairs (<20) to judge"}
        results["verdict"] = {"ship": False,
                              "recommendation": "Insufficient coverage to judge; do not change ranking."}
        log(f"\n  Too few covered pairs ({len(cov)}) to judge — no ranking change.")

    (HERE / "plausibility_rerank_results.json").write_text(
        json.dumps(results, indent=2), encoding="utf-8")
    (HERE / "plausibility_rerank_run.log").write_text("\n".join(log_lines), encoding="utf-8")
    log(f"\n  Wrote: {HERE / 'plausibility_rerank_results.json'}")
    return results


if __name__ == "__main__":
    md = int(sys.argv[1]) if len(sys.argv) > 1 else 60
    try:
        run(md)
    except Exception:
        import traceback
        traceback.print_exc()
        sys.exit(1)
