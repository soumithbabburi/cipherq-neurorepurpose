# SOP-003 — Predictive (Results) Validation

**Standard Operating Procedure**
**Document ID:** SOP-003
**Version:** 1.0
**Effective date:** 2026-06-29
**Owner:** Solix Technologies — Data Team
**Applies to:** The CipherQ / RepurposeIQ repurposing engine — i.e. whether the
ranked drug→indication results it produces are *correct*, not just whether the
underlying data is clean.

> Companion to **SOP-001/002** (data validation). Those certify the *inputs*; this
> certifies the *output*. A reviewer who accepts the data still asks "but do your
> predictions recover real repurposing?" — SOP-003 answers that with a number.

---

## 1. Purpose

Data validation cannot prove the engine works; only a benchmark against known
outcomes can. This SOP measures whether the engine's mechanistic score separates
**real (approved) repurposings** from **plausible-but-failed** ones, on an external
gold standard, **with retrospective leakage controlled**.

It answers the question that actually gates a repurposing deal:
*"Would your platform have ranked known repurposing successes above failures, using
only biology it could have known beforehand?"*

## 2. Scope

- **In scope:** the engine's mechanistic ranking (target overlap + direction-aware
  pathway + PPI proximity) and the calibration of its score to a probability.
- **Out of scope:** clinical-trial and literature features — deliberately **removed**
  for this test because, in a backward-looking benchmark, they encode the label
  (leakage). They remain in the live product, where they are legitimate.

## 3. Gold standard & negatives

| Item | Source |
|---|---|
| Positives | repoDB status = **Approved** drug-indication pairs |
| Hard negatives | repoDB status = **Terminated / Withdrawn / Suspended** (tried, failed) |
| Reference | Brown & Patel, *A standard database for drug repositioning*, **Sci Data** 2017; repoDB v2.x |
| Drug mapping | repoDB drug name → ChEMBL (local `chembl_33`); drugs without a curated mechanism are excluded and the coverage % reported |
| Disease biology | Open Targets gene/pathway sets; STRING PPI (cached) |

repoDB file staged at `data/external/repodb/full.csv` (downloaded by the operator,
mirroring the GtoPdb handling in SOP for concordance).

## 4. Standards applied

| Framework | Role |
|---|---|
| **GAMP 5** (ISPE) | `validate_predictions.py` / `validate_calibration.py` are the executable OQ/PQ test protocol |
| **DAMA-DMBOK** | metrics framework (discrimination, early recognition, calibration) |
| Leakage control | clinical/literature features ablated; out-of-fold calibration only |

## 5. Procedure (read-only)

```bash
python validation/validate_predictions.py [repodb.csv] [max_diseases]
python validation/validate_calibration.py
```

| Step | Check | Acceptance criterion |
|---|---|---|
| 1 | Map repoDB drugs → ChEMBL + targets/actions | coverage reported; <40% triggers review |
| 2 | Resolve disease biology (genes/pathways/PPI), cached | ≥5 disease genes or disease skipped |
| 3 | Mechanism-only score per pair (no clinical/literature) | leakage features confirmed absent |
| 4 | **Discrimination**: AUROC, AUPRC, BEDROC(α=20) | AUROC ≥ 0.65 INFO · 0.55–0.65 WARN · <0.55 FAIL |
| 5 | **Baselines**: random, drug-popularity, target-only | engine must beat random + popularity |
| 6 | **Ablation**: direction-aware vs blind vs excluded pathway | report Δ; if Δ ≤ 0 persistently, revert term |
| 7 | **Negative control**: shuffle labels | AUROC ≈ 0.5 |
| 8 | **Calibration**: ECE, Brier, monotonicity; OOF isotonic | report raw vs recalibrated; monotonicity ≥ 0.6 |

## 6. Handling findings

Logged in **CAPA-LOG-predictive.md** (IDs `PRED-`, `CAL-`). Severity: **HIGH**
(customer-facing claim wrong) · **MEDIUM** · **LOW** · **INFO** (by-design /
control). A negative finding (e.g. a feature that does not help) is recorded and
acted on, not hidden.

## 7. Frequency

- On every data refresh (ChEMBL / Hetionet / Open Targets) and repoDB update.
- Before any customer-facing *predictive-performance* claim.
- Quarterly standing review; track AUROC over time for drift.

## 8. Records

| Record | Location |
|---|---|
| Executable protocols | `validation/validate_predictions.py`, `validation/validate_calibration.py` |
| Results | `validation/predictions_results.json`, `validation/calibration_results.json` |
| Scored-pairs sidecar | `validation/predictions_pairs.json` |
| Run transcripts | `validation/predictions_run.log`, `validation/calibration_run.log` |
| CAPA log | `validation/CAPA-LOG-predictive.md` |
| Unified scorecard | `DATA_QUALITY.md` §3b |

## 9. Honest boundary (state to customers)

A retrospective benchmark proves **historical recovery on a fixed snapshot**. It is
**not** a guarantee of prospective performance, and it measures the *mechanistic*
component only. Phrase as "retrospectively validated on repoDB with leakage
controlled" — never "clinically validated" or "FDA-validated."

## 10. Revision history

| Version | Date | Author | Change |
|---|---|---|---|
| 1.0 | 2026-06-29 | Solix Data Team | Initial issue; first predictive validation vs repoDB (AUROC 0.71 mechanism-only) + score calibration |
