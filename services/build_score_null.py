"""
Build the null score distribution for calibration (P1).
Scores RANDOM drug-disease pairs (mostly unrelated = background) and writes their
composite scores to data/score_null.json incrementally, so a partial run is still
usable. Run offline:  .venv312\\Scripts\\python.exe -m services.build_score_null [N]
"""
import csv, json, random, sys, time
from pathlib import Path

_ROOT = Path(__file__).parent.parent
NULL_FILE = _ROOT / "data" / "score_null.json"


def main(n: int = 120, seed: int = 3):
    from services.reverse_repurposing import canonical_pair_score
    rows = list(csv.DictReader(
        open(_ROOT / "data/external/repodb/full.csv", encoding="utf-8", errors="ignore")))
    drugs = sorted({r["drug_name"] for r in rows if r["drug_name"] and r["drug_name"] != "NA"})
    diseases = sorted({r["ind_name"] for r in rows if r["ind_name"] and r["ind_name"] != "NA"})
    rng = random.Random(seed)
    scores = []
    if NULL_FILE.exists():
        try:
            scores = json.loads(NULL_FILE.read_text()).get("scores", [])
        except Exception:
            scores = []
    t0 = time.time()
    for i in range(n):
        drug = rng.choice(drugs)
        dis = rng.choice(diseases)
        try:
            res = canonical_pair_score("", dis, drug_name=drug)
            s = res.get("composite_score")
            if res.get("error") or s is None:
                continue
            scores.append(round(float(s), 4))
            if len(scores) % 5 == 0:
                _write(scores)
                print(f"  {len(scores)} pairs  (last {drug}->{dis} = {s:.3f}, "
                      f"{time.time()-t0:.0f}s)", flush=True)
        except Exception as e:
            print(f"  ERR {drug}->{dis}: {e}", flush=True)
    _write(scores)
    print(f"DONE: {len(scores)} null scores written to {NULL_FILE}", flush=True)


def _write(scores):
    ss = sorted(scores)
    NULL_FILE.write_text(json.dumps({
        "scores": scores,
        "n": len(ss),
        "sorted": ss,
        "mean": round(sum(ss) / len(ss), 4) if ss else 0,
    }))


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 120)
