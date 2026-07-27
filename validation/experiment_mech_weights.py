"""
EXPERIMENT — should the mechanism composite be re-weighted toward target overlap?
═══════════════════════════════════════════════════════════════════════════════
predictions_results.json shows target-overlap ONLY (0.7464) beating the full
target+pathway+ppi triad (0.7347) on repoDB, with pathway dilutive. But a +0.012
AUROC gap on ~1.5k pairs is close to noise. This experiment settles it properly:

  * computes the three mechanism sub-scores (target, pathway-blind, ppi) ONCE per
    repoDB pair, then sweeps weightings as cheap linear combos;
  * reports overall AUROC AND per-disease mean AUROC (a re-weight that lifts the
    pooled number but hurts per-disease ranking is NOT a win);
  * bootstraps a 95% CI on the AUROC GAP vs the current weighting, so we only ship
    a change whose improvement is real, not sampling noise.

Ship criterion: adopt a weighting only if its AUROC gap vs current has a bootstrap
95% CI lower bound >= 0 (non-inferior) AND a positive point estimate AND per-disease
mean AUROC does not regress. Read-only. Writes validation/mech_weights_results.json.
"""
import sys
import json
import datetime
from pathlib import Path
from collections import defaultdict

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))

from validation.validate_predictions import (
    DEFAULT_REPODB, load_repodb, _resolve_drug_targets, auroc,
)

# weightings over (target, pathway, ppi)
_VARIANTS = {
    "current_triad":   (0.25, 0.20, 0.20),
    "target_only":     (1.00, 0.00, 0.00),
    "no_pathway":      (0.25, 0.00, 0.20),
    "target_heavy":    (0.45, 0.10, 0.20),
    "target_heavy_lowpw": (0.50, 0.05, 0.15),
    "mostly_target":   (0.60, 0.10, 0.15),
}


def _combo(tpq, w):
    ws = sum(w) or 1.0
    return (w[0] * tpq[0] + w[1] * tpq[1] + w[2] * tpq[2]) / ws


def run(max_diseases=60, n_boot=2000):
    from validation.validate_concordance import get_engine
    from services.disease_ontology import resolve_disease, get_ppi_network
    from services.repurposing_engine import (
        _score_target_overlap, _score_pathway_overlap, _score_ppi_network)

    repo = load_repodb(str(DEFAULT_REPODB))
    all_drugs = sorted({d for dis in repo.values() for d in dis})
    eng = get_engine()
    drug_genes, _da, _dc = _resolve_drug_targets(eng, all_drugs, "fallback")
    mapped = set(drug_genes)

    cand = []
    for dis, drugs in repo.items():
        lab = {d: y for d, y in drugs.items() if d in mapped}
        npos = sum(1 for y in lab.values() if y == 1)
        nneg = sum(1 for y in lab.values() if y == 0)
        if npos >= 1 and nneg >= 1 and (npos + nneg) >= 5:
            cand.append((dis, npos + nneg))
    cand.sort(key=lambda x: -x[1])
    cand = cand[:max_diseases]

    # per pair: (t, p_blind, q) sub-scores + label + disease index
    tpq, labels, dis_idx = [], [], []
    per_dis = defaultdict(list)
    resolved = 0
    for di, (dis, _n) in enumerate(cand):
        info = resolve_disease(dis) or {}
        dgenes = [t["gene_symbol"] for t in info.get("targets", [])[:40]]
        dpath = info.get("pathways", [])
        if len(dgenes) < 5:
            continue
        resolved += 1
        ppi = get_ppi_network(dgenes[:15]) if dgenes else {}
        for d, y in [(d, y) for d, y in repo[dis].items() if d in mapped]:
            g = sorted(drug_genes[d])
            t = _score_target_overlap(g, dgenes)
            p = _score_pathway_overlap(g, dpath, drug_signature=None)   # direction-blind
            q = _score_ppi_network(g, ppi)
            tpq.append((t, p, q)); labels.append(y); dis_idx.append(di)
            per_dis[di].append(((t, p, q), y))

    n = len(labels)
    tpq_arr = tpq
    ys = labels

    def variant_scores(w):
        return [_combo(x, w) for x in tpq_arr]

    def overall_auroc(w):
        return auroc(list(zip(variant_scores(w), ys)))

    def per_disease_mean(w):
        aus = []
        for di, pairs in per_dis.items():
            sc = [(_combo(x, w), y) for x, y in pairs]
            a = auroc(sc)
            if a is not None:
                aus.append(a)
        return round(sum(aus) / len(aus), 4) if aus else None

    base_w = _VARIANTS["current_triad"]
    base_pooled = variant_scores(base_w)

    rng = np.random.default_rng(42)
    idx = np.arange(n)
    boot_idx = [rng.choice(idx, size=n, replace=True) for _ in range(n_boot)]

    def boot_gap_ci(w):
        """95% CI on AUROC(variant) - AUROC(current) via paired bootstrap over pairs."""
        vs = variant_scores(w)
        diffs = []
        for bi in boot_idx:
            av = auroc([(vs[i], ys[i]) for i in bi])
            ab = auroc([(base_pooled[i], ys[i]) for i in bi])
            if av is not None and ab is not None:
                diffs.append(av - ab)
        if not diffs:
            return None, None
        lo, hi = np.percentile(diffs, [2.5, 97.5])
        return round(float(lo), 4), round(float(hi), 4)

    results = {
        "run_at": datetime.datetime.now().isoformat(timespec="seconds"),
        "pairs": n, "diseases": resolved, "n_boot": n_boot,
        "baseline": "current_triad (0.25/0.20/0.20)",
        "variants": {},
    }
    base_overall = overall_auroc(base_w)
    for name, w in _VARIANTS.items():
        ov = overall_auroc(w)
        pdm = per_disease_mean(w)
        lo, hi = boot_gap_ci(w) if name != "current_triad" else (0.0, 0.0)
        gap = None if (ov is None or base_overall is None) else round(ov - base_overall, 4)
        results["variants"][name] = {
            "weights_target_pathway_ppi": list(w),
            "overall_auroc": None if ov is None else round(ov, 4),
            "per_disease_mean_auroc": pdm,
            "gap_vs_current": gap,
            "gap_ci95": [lo, hi],
        }

    base_pdm = results["variants"]["current_triad"]["per_disease_mean_auroc"]
    # ship decision: positive point gap, CI lower bound >= 0, per-disease not worse
    winners = []
    for name, v in results["variants"].items():
        if name == "current_triad":
            continue
        lo = v["gap_ci95"][0]
        if (v["gap_vs_current"] or 0) > 0 and lo is not None and lo >= 0 \
           and v["per_disease_mean_auroc"] is not None and v["per_disease_mean_auroc"] >= base_pdm - 0.005:
            winners.append((name, v["gap_vs_current"], lo))
    winners.sort(key=lambda x: -x[1])
    results["winners_meeting_ship_criterion"] = [
        {"variant": w[0], "gap": w[1], "ci_low": w[2]} for w in winners]
    results["recommendation"] = (
        f"Re-weight toward '{winners[0][0]}' (gap +{winners[0][1]}, CI low {winners[0][2]})."
        if winners else
        "No weighting robustly beats the current triad (CI includes 0 or per-disease regresses). "
        "Keep current weights; the target-only edge is within sampling noise.")

    (HERE / "mech_weights_results.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
    print(json.dumps({k: results["variants"][k] for k in results["variants"]}, indent=2))
    print("\nRECOMMENDATION:", results["recommendation"])
    return results


if __name__ == "__main__":
    md = int(sys.argv[1]) if len(sys.argv) > 1 else 60
    try:
        run(md)
    except Exception:
        import traceback; traceback.print_exc(); sys.exit(1)
