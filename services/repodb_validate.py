"""
RepoDB validation  —  does the repurposing scorer separate SUCCESS from FAILURE?
═══════════════════════════════════════════════════════════════════════════════
Scores a stratified sample of RepoDB pairs — Approved (label 1) vs Failed
(Terminated/Withdrawn/Suspended, label 0) — and reports ROC-AUC + the per-class
score distributions. This is the honest test the platform was missing: if the
scorer cannot rank real repurposing successes above real failures, the tier
cutoffs are meaningless and scoring must be fixed before a threshold is set.

Two scores per pair:
  • composite     — the full canonical score. For an ON-LABEL Approved pair the
                    indication/clinical dims fire, so this is expected to be
                    optimistic (mild outcome leakage) — reported for reference.
  • mechanistic   — target/pathway/PPI/network only, the part that GENERALISES to
                    a novel candidate (no trial/indication signal). This is the
                    fair test and the number the cutoff work should rely on.

Cache writes are disabled so a parallel run cannot corrupt reverse_cache.json.

Run:  .venv312\\Scripts\\python.exe -m services.repodb_validate [n_per_class]
"""
from __future__ import annotations

import csv
import json
import random
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

_ROOT = Path(__file__).parent.parent
_REPODB = _ROOT / "data" / "external" / "repodb" / "full.csv"
_OUT = _ROOT / "data" / "repodb_validation.json"

_FAIL = {"Terminated", "Withdrawn", "Suspended"}


def _sample(n_per_class: int, seed: int = 7):
    rows = list(csv.DictReader(open(_REPODB, encoding="utf-8", errors="ignore")))
    pos, neg = [], []
    seen = set()
    for r in rows:
        drug, ind, st = r.get("drug_name"), r.get("ind_name"), r.get("status")
        if not drug or not ind or drug == "NA" or ind == "NA":
            continue
        key = (drug.lower(), ind.lower())
        if key in seen:
            continue
        seen.add(key)
        if st == "Approved":
            pos.append((drug, ind))
        elif st in _FAIL:
            neg.append((drug, ind))
    rng = random.Random(seed)
    rng.shuffle(pos); rng.shuffle(neg)
    # oversample ~2x target to survive resolution failures
    return pos[:n_per_class * 2], neg[:n_per_class * 2]


def _mech_score(sr: dict) -> float:
    """Mechanistic-only score (no indication/clinical/regulatory leakage).
    Renormalised target/pathway/ppi + the KG network bonus."""
    s = sr.get("scores", {}) or {}
    net = (sr.get("score_breakdown", {}) or {}).get("network_score", 0.0) or 0.0
    mech = 0.45 * s.get("target", 0) + 0.30 * s.get("pathway", 0) + 0.25 * s.get("ppi", 0)
    return round(min(1.0, mech + 0.15 * float(net)), 4)


def main(n_per_class: int = 120):
    # Protect the live cache: disable writes for this run.
    from services import reverse_repurposing as rr
    rr._cache_put = lambda *a, **k: None
    from services.reverse_repurposing import canonical_pair_score

    pos_pairs, neg_pairs = _sample(n_per_class)
    jobs = [(d, i, 1) for d, i in pos_pairs] + [(d, i, 0) for d, i in neg_pairs]
    print(f"scoring up to {len(jobs)} pairs "
          f"({len(pos_pairs)} approved-pool, {len(neg_pairs)} failed-pool), "
          f"target {n_per_class}/class...", flush=True)

    results = []
    t0 = time.time()

    def _score(job):
        drug, ind, label = job
        try:
            r = canonical_pair_score("", ind, drug_name=drug)
            if r.get("error"):
                return None
            # require BOTH sides resolved for a meaningful mechanistic comparison
            if not r.get("drug_genes"):
                return None
            return {"drug": drug, "ind": ind, "label": label,
                    "composite": float(r.get("composite_score", 0.0)),
                    "mech": _mech_score(r),
                    "scores": r.get("scores", {})}
        except Exception:
            return None

    done = 0
    with ThreadPoolExecutor(max_workers=8) as ex:
        futs = {ex.submit(_score, j): j for j in jobs}
        for fut in as_completed(futs):
            done += 1
            row = fut.result()
            if row:
                results.append(row)
            if done % 20 == 0:
                np_ = sum(1 for r in results if r["label"] == 1)
                nn_ = sum(1 for r in results if r["label"] == 0)
                print(f"  {done}/{len(jobs)} scored, {len(results)} resolved "
                      f"({np_} approved, {nn_} failed), {time.time()-t0:.0f}s", flush=True)

    _report(results)


def _auc(scores_labels):
    """ROC-AUC via the Mann-Whitney U statistic (no sklearn dependency)."""
    pos = [s for s, l in scores_labels if l == 1]
    neg = [s for s, l in scores_labels if l == 0]
    if not pos or not neg:
        return None
    # rank-based AUC
    allv = sorted([(s, l) for s, l in scores_labels])
    ranks = {}
    i = 0
    while i < len(allv):
        j = i
        while j < len(allv) and allv[j][0] == allv[i][0]:
            j += 1
        avg_rank = (i + j - 1) / 2.0 + 1
        for k in range(i, j):
            ranks[k] = avg_rank
        i = j
    sum_pos = sum(ranks[k] for k, (s, l) in enumerate(allv) if l == 1)
    n_p, n_n = len(pos), len(neg)
    auc = (sum_pos - n_p * (n_p + 1) / 2.0) / (n_p * n_n)
    return round(auc, 4)


def _stats(vals):
    if not vals:
        return {}
    v = sorted(vals)
    n = len(v)
    return {"n": n, "mean": round(sum(v) / n, 4),
            "median": round(v[n // 2], 4),
            "p25": round(v[n // 4], 4), "p75": round(v[3 * n // 4], 4)}


def _report(results):
    pos = [r for r in results if r["label"] == 1]
    neg = [r for r in results if r["label"] == 0]
    out = {
        "n_resolved": len(results), "n_approved": len(pos), "n_failed": len(neg),
        "auc_composite": _auc([(r["composite"], r["label"]) for r in results]),
        "auc_mechanistic": _auc([(r["mech"], r["label"]) for r in results]),
        "composite_approved": _stats([r["composite"] for r in pos]),
        "composite_failed": _stats([r["composite"] for r in neg]),
        "mech_approved": _stats([r["mech"] for r in pos]),
        "mech_failed": _stats([r["mech"] for r in neg]),
        "rows": results,
    }
    _OUT.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print("\n" + "=" * 64)
    print(f"RESOLVED {len(results)} pairs  ({len(pos)} approved, {len(neg)} failed)")
    print(f"AUC  composite (leaky, ref) : {out['auc_composite']}")
    print(f"AUC  mechanistic (fair test): {out['auc_mechanistic']}")
    print(f"mech  approved: {out['mech_approved']}")
    print(f"mech  failed  : {out['mech_failed']}")
    print(f"comp  approved: {out['composite_approved']}")
    print(f"comp  failed  : {out['composite_failed']}")
    print(f"\nwrote {_OUT}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 120)
