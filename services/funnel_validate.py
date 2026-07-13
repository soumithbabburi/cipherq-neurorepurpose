"""
Does the potency + window funnel recover the SUCCESS signal the graph model lacks?
═══════════════════════════════════════════════════════════════════════════════
The DWPC plausibility model scores AUC 0.42 on RepoDB Approved-vs-Failed (it ranks
plausibility, not success). This experiment asks the payoff question: does adding
the physical funnel — potency (IC50 vs the disease target) and therapeutic window
(free Cmax / IC50) — separate real successes from real failures better than the
graph alone?

For each RepoDB pair: disease → Open Targets driver genes → the drug's best measured
pChEMBL against them (potency), then PBPK free-Cmax margin (window). Reports AUC for
each signal and for the stacked funnel, on the covered subset. Honest about coverage.
"""
from __future__ import annotations

import csv
import json
import random
from functools import lru_cache
from pathlib import Path

import numpy as np

_ROOT = Path(__file__).parent.parent
_REPODB = _ROOT / "data" / "external" / "repodb" / "full.csv"
_FAIL = {"Terminated", "Withdrawn", "Suspended"}


def _auc(scores, labels):
    xs = [(s, l) for s, l in zip(scores, labels) if s is not None]
    if not xs:
        return None, 0
    xs.sort()
    ranks = {}; i = 0
    while i < len(xs):
        j = i
        while j < len(xs) and xs[j][0] == xs[i][0]:
            j += 1
        r = (i + j - 1) / 2.0 + 1
        for k in range(i, j):
            ranks[k] = r
        i = j
    P = sum(l for _, l in xs); N = len(xs) - P
    if P == 0 or N == 0:
        return None, len(xs)
    sp = sum(ranks[k] for k, (s, l) in enumerate(xs) if l == 1)
    return round((sp - P * (P + 1) / 2.0) / (P * N), 4), len(xs)


@lru_cache(maxsize=2048)
def _disease_genes(name: str):
    try:
        from services.disease_ontology import resolve_disease
        d = resolve_disease(name) or {}
        return tuple(t["gene_symbol"] for t in d.get("targets", [])[:40] if t.get("gene_symbol"))
    except Exception:
        return tuple()


def _sample(n_per_class, seed=11):
    rows = list(csv.DictReader(open(_REPODB, encoding="utf-8", errors="ignore")))
    pos, neg, seen = [], [], set()
    for r in rows:
        d, ind, st, cid = r.get("drug_name"), r.get("ind_name"), r.get("status"), r.get("drugbank_id")
        if not d or not ind or d == "NA":
            continue
        k = (d.lower(), ind.lower())
        if k in seen:
            continue
        seen.add(k)
        if st == "Approved":
            pos.append((d, ind))
        elif st in _FAIL:
            neg.append((d, ind))
    rng = random.Random(seed); rng.shuffle(pos); rng.shuffle(neg)
    return pos[:n_per_class], neg[:n_per_class]


def main(n_per_class=180):
    from services.lead_viability import potency, therapeutic_window
    from services.repurposing_predictor import plausibility

    pos, neg = _sample(n_per_class)
    jobs = [(d, i, 1) for d, i in pos] + [(d, i, 0) for d, i in neg]
    rows = []
    for k, (drug, ind, label) in enumerate(jobs, 1):
        genes = list(_disease_genes(ind))
        pot = potency(drug, genes) if genes else {"covered": False}
        plaus = plausibility(drug, ind)
        rec = {"drug": drug, "ind": ind, "label": label,
               "plaus": (plaus or {}).get("probability"),
               "pchembl": pot.get("pchembl") if pot.get("covered") else None,
               "psource": pot.get("source") if pot.get("covered") else None,
               "margin": None}
        if pot.get("covered") and pot.get("ic50_nm"):
            win = therapeutic_window(drug, pot["ic50_nm"])
            rec["margin"] = win.get("margin") if win.get("covered") else None
        rows.append(rec)
        if k % 25 == 0:
            nc = sum(1 for r in rows if r["pchembl"] is not None)
            print(f"  {k}/{len(jobs)} — potency-covered {nc}", flush=True)

    lab = [r["label"] for r in rows]
    # individual signals
    auc_plaus, n_plaus = _auc([r["plaus"] for r in rows], lab)
    auc_pot, n_pot = _auc([r["pchembl"] for r in rows], lab)
    auc_win, n_win = _auc([r["margin"] for r in rows], lab)

    # stacked funnel score on the potency-covered subset: normalise + combine
    def _funnel(r):
        if r["pchembl"] is None:
            return None
        pl = r["plaus"] if r["plaus"] is not None else 0.15
        pot_s = min(1.0, max(0.0, (r["pchembl"] - 4.0) / 5.0))     # pChEMBL 4→0, 9→1
        win_s = 1.0 if r["margin"] is None else min(1.0, (r["margin"]) / 10.0)
        return pl * pot_s * win_s
    auc_funnel, n_funnel = _auc([_funnel(r) for r in rows], lab)

    # potency+window only (no graph), covered subset
    def _pw(r):
        if r["pchembl"] is None:
            return None
        pot_s = min(1.0, max(0.0, (r["pchembl"] - 4.0) / 5.0))
        win_s = 1.0 if r["margin"] is None else min(1.0, (r["margin"]) / 10.0)
        return pot_s * win_s
    auc_pw, n_pw = _auc([_pw(r) for r in rows], lab)

    # AUC broken out BY potency source — the key insight: disease-specific potency
    # carries the success signal; generic fallback potency does not (it only supports
    # the sloppy-binder elimination). Reported honestly so we know what to trust.
    by_source = {}
    for src in ("disease-target", "mechanism-target", "best-measured"):
        sub = [(r["pchembl"], r["label"]) for r in rows if r.get("psource") == src]
        a, n = _auc([s for s, _ in sub], [l for _, l in sub])
        by_source[src] = {"auc": a, "n": n}
    # disease-target OR mechanism-target (pharmacologically specific), excluding generic
    spec = [(r["pchembl"], r["label"]) for r in rows
            if r.get("psource") in ("disease-target", "mechanism-target")]
    auc_spec, n_spec = _auc([s for s, _ in spec], [l for _, l in spec])

    out = {
        "n_pairs": len(rows),
        "n_approved": sum(lab), "n_failed": len(lab) - sum(lab),
        "coverage_potency": n_pot, "coverage_window": n_win,
        "auc_plausibility_only": auc_plaus, "n_plausibility": n_plaus,
        "auc_potency_only": auc_pot,
        "auc_window_only": auc_win,
        "auc_potency_x_window": auc_pw,
        "auc_full_funnel": auc_funnel,
        "auc_potency_by_source": by_source,
        "auc_potency_specific": {"auc": auc_spec, "n": n_spec},
    }
    (_ROOT / "data" / "funnel_validation.json").write_text(json.dumps(out, indent=2))
    print("\n" + "=" * 64)
    print(f"pairs={len(rows)}  approved={sum(lab)}  failed={len(lab)-sum(lab)}")
    print(f"potency coverage: {n_pot}/{len(rows)}   window coverage: {n_win}/{len(rows)}")
    print("-- AUC vs Approved/Failed (higher = better success separation) --")
    print(f"  plausibility (graph) only : {auc_plaus}   (n={n_plaus})   [baseline ~0.42]")
    print(f"  potency (IC50) only       : {auc_pot}")
    print(f"  window (Cmax/IC50) only   : {auc_win}")
    print(f"  potency × window          : {auc_pw}   (n={n_pw})")
    print(f"  FULL funnel (plaus×pot×win): {auc_funnel}   (n={n_funnel})")
    print("-- potency AUC BY SOURCE (the key insight) --")
    for src, v in by_source.items():
        print(f"  {src:<16}: AUC {v['auc']}   (n={v['n']})")
    print(f"  disease+mechanism (specific): AUC {auc_spec}   (n={n_spec})")


if __name__ == "__main__":
    import sys
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 180)
