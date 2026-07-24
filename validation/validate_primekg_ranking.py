"""
PrimeKG UNIVERSE-RANKING validation (the honest baseline before a better generator).
════════════════════════════════════════════════════════════════════════════════════
Every existing PrimeKG metric is a PAIRWISE AUC (indication vs contra / vs random on
labeled pairs). That measures discrimination on a curated slice — it does NOT measure the
thing a universe generator must do: "given a drug, rank ALL 17,080 diseases and put its
true indications near the top." The tail is exactly where top_diseases_for_drug was flagged
as noisy (metformin -> common cold), and it was never quantified. This harness quantifies it.

Protocol (honest, reproducible):
  * COMPOUND-DISJOINT split reproduced bit-for-bit from train_primekg_treats.py
    (seed 42, 15% of indication-bearing drugs held out) so we only ever rank drugs the
    treats classifier NEVER trained on (true zero-shot).
  * For each held-out drug, score every disease and rank its true indications.
  * FILTERED ranking (KG-completion convention): a drug's OTHER known indications are
    removed from the candidate pool when ranking a target, so co-indications don't
    penalise each other.
  * Metrics per ranker: recall@{10,20,50,100}, MRR, median rank, mean rank.

Rankers compared:
  A. distmult    — DistMult relatedness only (the recall stage alone; the pre-direction baseline)
  B. cascade     — the SHIPPED top_diseases_for_drug: DistMult top-400 pool, treats re-rank
  C. treats_full — the treats classifier scored over ALL 17,080 diseases (no pool cap)

Output:
  validation/primekg_ranking_results.json
Run:  .venv312/Scripts/python.exe -m validation.validate_primekg_ranking
"""
from __future__ import annotations

import json
import random
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from services import primekg_predictor as pkg  # noqa: E402

KS = (10, 20, 50, 100)
_POOL = 400   # must match primekg_predictor.top_diseases_for_drug


def _reproduce_split(labels, seed=42, holdout_drugs=0.15):
    """Exact compound-disjoint split from train_primekg_treats.main()."""
    rng = random.Random(seed)
    drugs = [int(d) for d, rels in labels.items() if rels.get("indication")]
    rng.shuffle(drugs)
    n_hold = int(len(drugs) * holdout_drugs)
    return set(drugs[:n_hold])          # test drugs (never seen by the classifier)


def _treats_all(s, d_idx, disease_ids):
    """P(treats) for (drug, every disease) — vectorised classifier scores."""
    E = s["E"]; ed = E[d_idx]; Ez = E[disease_ids]
    n = len(disease_ids)
    feats = np.concatenate(
        [np.broadcast_to(ed, (n, ed.shape[0])), Ez, ed * Ez, np.abs(ed - Ez)], axis=1)
    return s["treats"].predict_proba(feats)[:, 1]


def _distmult_all(s, d_idx, disease_ids):
    rv = s["R"][s["rels"]["indication"]]
    return (s["E"][d_idx] * rv) @ s["E"][disease_ids].T


def _filtered_ranks(scores, pos_local, all_pos_local):
    """Ranks (1-based) of each target in pos_local, after removing the drug's OTHER known
    positives from the candidate set. `scores`/indices are in disease-array-local space."""
    ranks = []
    other = all_pos_local - set(pos_local)
    keep = np.ones(len(scores), dtype=bool)
    for j in other:
        keep[j] = False
    kept_scores = scores[keep]
    # map each target's score to its rank among kept candidates
    for j in pos_local:
        # target must itself be kept (it is: it's not in `other`)
        r = int(np.sum(kept_scores > scores[j]) + 1)
        ranks.append(r)
    return ranks


def _summ(ranks):
    ranks = np.array(ranks, dtype=np.float64)
    out = {f"recall@{k}": round(float(np.mean(ranks <= k)), 4) for k in KS}
    out["mrr"] = round(float(np.mean(1.0 / ranks)), 4)
    out["median_rank"] = int(np.median(ranks))
    out["mean_rank"] = round(float(np.mean(ranks)), 1)
    out["n_pairs"] = int(len(ranks))
    return out


def main():
    s = pkg._load()
    if s is None or s.get("treats") is None:
        print("PrimeKG artifacts / treats classifier not present — cannot evaluate.")
        return 1
    labels = s["labels"]
    disease_ids = s["disease_ids"]
    n_dis = len(disease_ids)
    # disease node-index -> local position in the disease_ids array
    pos_of = {int(z): i for i, z in enumerate(disease_ids)}

    test_drugs = _reproduce_split(labels)
    # keep test drugs that actually have >=1 indication resolvable into the disease array
    ev4 = []
    for d in test_drugs:
        pos = [pos_of[z] for z in labels[str(d)].get("indication", []) if z in pos_of]
        if pos:
            ev4.append((d, set(pos)))
    print(f"[rank-eval] {n_dis} diseases | {len(ev4)} held-out test drugs with indications", flush=True)

    R = {name: [] for name in ("distmult", "cascade", "treats_full")}
    t0 = time.time()
    for n, (d, pos_local) in enumerate(ev4, 1):
        dm = _distmult_all(s, d, disease_ids)          # (n_dis,)
        tf = _treats_all(s, d, disease_ids)            # (n_dis,)
        # cascade: re-rank only the DistMult top-_POOL pool with the classifier; everything
        # outside the pool is unreachable -> assign a rank worse than any in-pool item.
        pool = np.argsort(-dm)[:_POOL]
        casc = np.full(n_dis, -np.inf)
        casc[pool] = tf[pool]
        R["distmult"].extend(_filtered_ranks(dm, pos_local, pos_local))
        R["cascade"].extend(_filtered_ranks(casc, pos_local, pos_local))
        R["treats_full"].extend(_filtered_ranks(tf, pos_local, pos_local))
        if n % 25 == 0:
            print(f"  {n}/{len(ev4)} drugs  ({time.time()-t0:.0f}s)", flush=True)

    results = {name: _summ(r) for name, r in R.items()}
    meta = {
        "protocol": "compound-disjoint (seed 42, 15% holdout), filtered ranking over all diseases",
        "n_diseases": int(n_dis),
        "n_test_drugs": len(ev4),
        "elapsed_s": round(time.time() - t0, 1),
    }
    out = {"meta": meta, "rankers": results}
    outp = ROOT / "validation" / "primekg_ranking_results.json"
    json.dump(out, open(outp, "w"), indent=2)

    print("\n=== PrimeKG universe-ranking (held-out drugs, all diseases) ===", flush=True)
    hdr = f"{'ranker':12} " + " ".join(f"R@{k:<4}" for k in KS) + "  MRR    medRank"
    print(hdr, flush=True)
    for name in ("distmult", "cascade", "treats_full"):
        m = results[name]
        row = f"{name:12} " + " ".join(f"{m[f'recall@{k}']:<5}" for k in KS)
        row += f"  {m['mrr']:<6} {m['median_rank']}"
        print(row, flush=True)
    print(f"\nwrote {outp}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
