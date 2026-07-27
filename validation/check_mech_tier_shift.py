"""
Tier re-check for the 2026-07-27 mechanism re-weight (pathway dropped).
Recomputes the engine's mechanism component under the OLD (0.25/0.20/0.20) and
NEW (0.36/0/0.29) weights across the repoDB pairs, using the engine's exact
renorm (mass 0.65, denom floor 0.45), and reports how the score distribution and
the tier-boundary crossings (0.40 Promising, 0.60 Strong) shift. A re-weight that
improves ranking but blows up the tier populations would NOT ship.
"""
import sys, json, datetime
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))
from validation.validate_predictions import DEFAULT_REPODB, load_repodb, _resolve_drug_targets, auroc

OLD = {"target": 0.25, "pathway": 0.20, "ppi": 0.20}
NEW = {"target": 0.36, "pathway": 0.00, "ppi": 0.29}
FLOOR, MASS = 0.45, 0.65


def mech_component(W, t, p, q):
    dims = {"target": t, "pathway": p, "ppi": q}
    active = [k for k, v in dims.items() if v > 1e-9]
    if not active:
        return 0.0
    w_active = sum(W[k] for k in active)
    mech_raw = sum(dims[k] * W[k] for k in active)
    return (mech_raw / max(FLOOR, w_active)) * MASS


def run(max_diseases=60):
    from validation.validate_concordance import get_engine
    from services.disease_ontology import resolve_disease, get_ppi_network
    from services.repurposing_engine import _score_target_overlap, _score_pathway_overlap, _score_ppi_network

    repo = load_repodb(str(DEFAULT_REPODB))
    all_drugs = sorted({d for dis in repo.values() for d in dis})
    dg, _da, _dc = _resolve_drug_targets(get_engine(), all_drugs, "fallback")
    mapped = set(dg)
    cand = []
    for dis, drugs in repo.items():
        lab = {d: y for d, y in drugs.items() if d in mapped}
        np_, nn = sum(1 for y in lab.values() if y == 1), sum(1 for y in lab.values() if y == 0)
        if np_ >= 1 and nn >= 1 and (np_ + nn) >= 5:
            cand.append((dis, np_ + nn))
    cand.sort(key=lambda x: -x[1]); cand = cand[:max_diseases]

    old_s, new_s, ys = [], [], []
    for dis, _n in cand:
        info = resolve_disease(dis) or {}
        dgenes = [t["gene_symbol"] for t in info.get("targets", [])[:40]]
        if len(dgenes) < 5:
            continue
        dpath = info.get("pathways", [])
        ppi = get_ppi_network(dgenes[:15]) if dgenes else {}
        for d, y in [(d, y) for d, y in repo[dis].items() if d in mapped]:
            g = sorted(dg[d])
            t = _score_target_overlap(g, dgenes)
            p = _score_pathway_overlap(g, dpath, drug_signature=None)
            q = _score_ppi_network(g, ppi)
            old_s.append(mech_component(OLD, t, p, q))
            new_s.append(mech_component(NEW, t, p, q))
            ys.append(y)

    old_s, new_s = np.array(old_s), np.array(new_s)
    n = len(ys)

    def frac(a, thr):
        return round(float((a >= thr).mean()), 4)

    out = {
        "run_at": datetime.datetime.now().isoformat(timespec="seconds"),
        "pairs": n,
        "mech_component_auroc": {"old": round(auroc(list(zip(old_s.tolist(), ys))), 4),
                                 "new": round(auroc(list(zip(new_s.tolist(), ys))), 4)},
        "mean": {"old": round(float(old_s.mean()), 4), "new": round(float(new_s.mean()), 4)},
        "frac_ge_0.40_promising": {"old": frac(old_s, 0.40), "new": frac(new_s, 0.40)},
        "frac_ge_0.60_strong": {"old": frac(old_s, 0.60), "new": frac(new_s, 0.60)},
        "max": {"old": round(float(old_s.max()), 4), "new": round(float(new_s.max()), 4)},
    }
    (HERE / "mech_tier_shift.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))
    return out


if __name__ == "__main__":
    run(int(sys.argv[1]) if len(sys.argv) > 1 else 60)
