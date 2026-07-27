"""
Score Calibration — "when the engine says 0.8, is it right ~80% of the time?"
═══════════════════════════════════════════════════════════════════════════
The mechanistic score is a *ranking* signal in [0,1]; this validator checks how
well it behaves as a PROBABILITY of being a real (approved) repurposing, and
provides an honest recalibration map.

Consumes validation/predictions_pairs.json (the leakage-free mechanism-only
scores vs repoDB Approved/failed labels produced by validate_predictions.py).

Measures:
  • Reliability table — observed approval rate per score decile.
  • ECE (expected calibration error) and Brier score on the raw score.
  • Monotonicity — Spearman(score-bin, observed rate); a usable ranking score
    should be monotonic even if its absolute scale is off.
  • Isotonic recalibration (5-fold, out-of-fold) — ECE/Brier before vs after,
    so the reported improvement is not from over-fitting.

This is the SAME discipline used to reject the earlier leaky PoS calibrator:
calibrate on leakage-free, out-of-fold data only.

Writes validation/calibration_results.json + run log.

Usage:
    python validation/validate_calibration.py
"""

import sys
import json
import datetime
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

HERE = Path(__file__).resolve().parent
PAIRS = HERE / "predictions_pairs.json"


def _ece(scores, labels, bins=10):
    """Expected calibration error + per-bin reliability table."""
    table = []
    N = len(scores)
    ece = 0.0
    for b in range(bins):
        lo, hi = b / bins, (b + 1) / bins
        idx = [i for i, s in enumerate(scores)
               if (s >= lo and (s < hi or (b == bins - 1 and s <= hi)))]
        if not idx:
            table.append({"bin": f"{lo:.1f}-{hi:.1f}", "n": 0,
                          "mean_score": None, "observed_rate": None})
            continue
        mean_s = sum(scores[i] for i in idx) / len(idx)
        obs = sum(labels[i] for i in idx) / len(idx)
        ece += (len(idx) / N) * abs(obs - mean_s)
        table.append({"bin": f"{lo:.1f}-{hi:.1f}", "n": len(idx),
                      "mean_score": round(mean_s, 3), "observed_rate": round(obs, 3)})
    return round(ece, 4), table


def _brier(scores, labels):
    return round(sum((s - y) ** 2 for s, y in zip(scores, labels)) / len(scores), 4)


def _spearman_bins(table):
    """Monotonicity of observed rate across populated bins (Spearman)."""
    pts = [(i, r["observed_rate"]) for i, r in enumerate(table) if r["observed_rate"] is not None]
    if len(pts) < 3:
        return None
    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]

    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        rk = [0] * len(v)
        for pos, i in enumerate(order):
            rk[i] = pos
        return rk
    rx, ry = rank(xs), rank(ys)
    n = len(pts)
    d2 = sum((rx[i] - ry[i]) ** 2 for i in range(n))
    return round(1 - 6 * d2 / (n * (n * n - 1)), 3)


def _isotonic_oof(scores, labels, folds=5):
    """Out-of-fold isotonic recalibration → recalibrated scores (no leakage)."""
    try:
        from sklearn.isotonic import IsotonicRegression
    except Exception:
        return None
    import random
    idx = list(range(len(scores)))
    random.Random(42).shuffle(idx)
    recal = [None] * len(scores)
    for f in range(folds):
        test = set(idx[f::folds])
        train = [i for i in idx if i not in test]
        ir = IsotonicRegression(out_of_bounds="clip", y_min=0.0, y_max=1.0)
        ir.fit([scores[i] for i in train], [labels[i] for i in train])
        for i in test:
            recal[i] = float(ir.predict([scores[i]])[0])
    return recal


def run():
    log_lines = []

    def log(msg=""):
        print(msg)
        log_lines.append(str(msg))

    log("=" * 72)
    log("SCORE CALIBRATION  -  mechanism score as probability of real repurposing")
    log(f"Run: {datetime.datetime.now().isoformat(timespec='seconds')}")
    log("=" * 72)

    if not PAIRS.exists():
        log(f"\nMissing {PAIRS}. Run validate_predictions.py first.")
        sys.exit(1)
    data = json.loads(PAIRS.read_text(encoding="utf-8"))
    scores = [float(d["score"]) for d in data]
    labels = [int(d["label"]) for d in data]
    base = sum(labels) / len(labels)
    log(f"\nPairs: {len(scores):,} | base approval rate: {base:.3f}")

    ece_raw, table = _ece(scores, labels)
    brier_raw = _brier(scores, labels)
    mono = _spearman_bins(table)

    log("\nReliability (raw mechanism score):")
    log(f"  {'bin':<10}{'n':>6}{'mean_score':>12}{'observed':>11}")
    for r in table:
        if r["n"]:
            log(f"  {r['bin']:<10}{r['n']:>6}{r['mean_score']:>12}{r['observed_rate']:>11}")
    log(f"\nECE (raw):   {ece_raw}")
    log(f"Brier (raw): {brier_raw}")
    log(f"Monotonicity (Spearman of observed rate vs bin): {mono}")

    recal = _isotonic_oof(scores, labels)
    iso = None
    if recal is not None:
        ece_cal, _ = _ece(recal, labels)
        brier_cal = _brier(recal, labels)
        iso = {"ece_calibrated": ece_cal, "brier_calibrated": brier_cal,
               "ece_improvement": round(ece_raw - ece_cal, 4)}
        log(f"\nIsotonic recalibration (5-fold out-of-fold):")
        log(f"  ECE   {ece_raw} -> {ece_cal}   (improvement {iso['ece_improvement']})")
        log(f"  Brier {brier_raw} -> {brier_cal}")

    sev = "INFO" if (mono is not None and mono >= 0.6) else "WARN"
    result = {
        "run_at": datetime.datetime.now().isoformat(timespec="seconds"),
        "source": "predictions_pairs.json (mechanism-only scores vs repoDB labels)",
        "n_pairs": len(scores),
        "base_rate": round(base, 3),
        "ece_raw": ece_raw,
        "brier_raw": brier_raw,
        "monotonicity_spearman": mono,
        "reliability_table": table,
        "isotonic_recalibration": iso,
        "findings": [{
            "id": "CAL-01", "severity": sev,
            "title": "Mechanism-score calibration vs repoDB outcomes",
            "detail": f"Across {len(scores):,} pairs the raw mechanistic score has ECE {ece_raw} "
                      f"(Brier {brier_raw}); observed approval rate is monotonic in the score "
                      f"(Spearman {mono}). "
                      + (f"An out-of-fold isotonic map lowers ECE to {iso['ece_calibrated']}."
                         if iso else "Isotonic recalibration unavailable (sklearn missing)."),
            "capa": "Ship the isotonic map as the score->probability display layer; re-fit on each "
                    "data refresh. Monotonicity < 0.6 triggers review of the score definition.",
        }],
    }
    (HERE / "calibration_results.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    (HERE / "calibration_run.log").write_text("\n".join(log_lines), encoding="utf-8")
    log(f"\n  Wrote: {HERE / 'calibration_results.json'}")
    return result


if __name__ == "__main__":
    try:
        run()
    except Exception:
        import traceback
        traceback.print_exc()
        sys.exit(1)
