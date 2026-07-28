"""
Phase 0 of the learned composite (validation/LEARNED_COMPOSITE_DESIGN.md):
build + cache the LEAKAGE-FREE feature matrix for repoDB approved-vs-failed.

Per (drug, disease) pair it runs the real engine (canonical_pair_score) and records the
8 NON-LEAKY evidence features + the approved/failed label + the disease group (for the
disease-disjoint CV). Clinical/indication/regulatory are excluded (leaky); DWPC and
penalties stay outside the learned model.

The full-engine call hits per-pair services (proliferation -> mygene/PubMed, network ->
local DB); the disk cache warms across pairs. Progress is written INCREMENTALLY to
data/composite_dataset.json so a stall never loses work (resumes by skipping done pairs).

Run:  .venv312/Scripts/python.exe -m validation.build_composite_dataset [max_diseases]
"""
import sys
import json
import datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))

from validation.validate_predictions import DEFAULT_REPODB, load_repodb
from validation.validate_concordance import get_engine
from validation.validate_predictions import _resolve_drug_targets

OUT = ROOT / "data" / "composite_dataset.json"
FEATURES = ["target", "pathway", "ppi", "network", "directional",
            "proliferation", "signature", "direction"]   # NON-LEAKY only


def _feats(r: dict):
    s = r.get("scores", {}) or {}
    sb = r.get("score_breakdown", {}) or {}
    return [
        float(s.get("target", 0) or 0),
        float(s.get("pathway", 0) or 0),
        float(s.get("ppi", 0) or 0),
        float(sb.get("network_score", 0) or 0),
        float((r.get("directional_evidence", {}) or {}).get("signal", 0) or 0),
        float((r.get("proliferation", {}) or {}).get("score", 0) or 0),
        float((r.get("signature_reversal", {}) or {}).get("connectivity", 0) or 0),
        float((r.get("direction", {}) or {}).get("factor", 1.0) or 1.0),
    ]


def run(max_diseases: int = 80):
    from services.reverse_repurposing import canonical_pair_score, resolve_drug

    repo = load_repodb(str(DEFAULT_REPODB))
    all_drugs = sorted({d for dis in repo.values() for d in dis})
    eng = get_engine()
    drug_genes, _da, drug_chembl = _resolve_drug_targets(eng, all_drugs, "fallback")
    mapped = set(drug_genes)

    # diseases with BOTH classes (needed for a valid per-disease AUROC fold), most pairs first
    cand = []
    for dis, drugs in repo.items():
        lab = {d: y for d, y in drugs.items() if d in mapped}
        npos = sum(1 for y in lab.values() if y == 1)
        nneg = sum(1 for y in lab.values() if y == 0)
        if npos >= 1 and nneg >= 1 and (npos + nneg) >= 5:
            cand.append((dis, npos + nneg))
    cand.sort(key=lambda x: -x[1])
    cand = cand[:max_diseases]

    # resume: load any already-built rows, skip their pairs
    done = {}
    rows = []
    if OUT.exists():
        try:
            prev = json.loads(OUT.read_text())
            rows = prev.get("rows", [])
            done = {(r["drug"], r["disease"]) for r in rows}
        except Exception:
            rows, done = [], set()

    total_pairs = sum(1 for dis, _ in cand for d in repo[dis] if d in mapped)
    print(f"[ds] {len(cand)} diseases | ~{total_pairs} pairs | {len(rows)} already done", flush=True)

    def _save():
        OUT.write_text(json.dumps({
            "built_at": None, "features": FEATURES,
            "gold_standard": "repoDB approved(1) vs failed(0); non-leaky evidence features",
            "n_rows": len(rows), "rows": rows}, indent=1), encoding="utf-8")

    n = 0
    for dis, _tot in cand:
        for d, y in [(d, y) for d, y in repo[dis].items() if d in mapped]:
            if (d, dis) in done:
                continue
            n += 1
            try:
                chembl = drug_chembl.get(d, "")
                r = canonical_pair_score(chembl, dis, drug_name=d)
                if r.get("error"):
                    continue
                rows.append({"drug": d, "disease": dis, "label": int(y),
                             "features": _feats(r),
                             "hand_composite": float(r.get("composite_score", 0.0) or 0.0)})
            except Exception as e:
                if n <= 3:
                    print(f"  [warn] {d}/{dis}: {e}", flush=True)
            if n % 25 == 0:
                _save()
                print(f"  built {len(rows)} rows ({n} attempted)…", flush=True)
    _save()
    npos = sum(1 for r in rows if r["label"] == 1)
    print(f"[ds] DONE: {len(rows)} rows ({npos} approved / {len(rows)-npos} failed) -> {OUT}", flush=True)
    print(f"[ds] wrote {OUT}", flush=True)


if __name__ == "__main__":
    run(int(sys.argv[1]) if len(sys.argv) > 1 else 80)
