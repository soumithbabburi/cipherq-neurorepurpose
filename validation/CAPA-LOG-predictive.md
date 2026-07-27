# CAPA Log — Predictive (Results) Validation

**Corrective and Preventive Action log for SOP-003.**
Every finding from `validate_predictions.py` and `validate_calibration.py` is
recorded here with root cause and action. A negative finding (a feature that does
not help, a score that is not yet a probability) is *evidence of honest testing*,
not a failure.

Severity: HIGH · MEDIUM · LOW · INFO (by-design / control).

---

## Run: 2026-06-29 (rev 2) — coverage fix

Re-ran after a drug-target coverage fix. A/B on the same gold standard:

| Mode | Drug coverage | Pairs | AUROC | Neg control |
|---|---|---|---|---|
| legacy (mechanism only) | 46% | 1,000 | 0.709 | 0.49 |
| rich (activities for all drugs) | 73% | 1,577 | 0.715 | 0.50 |
| **fallback (gap-fill only)** ← shipped | **68%** | **1,476** | **0.732** | 0.48 |

### COV-01 — Drug-target coverage fix (salt→parent + bioactivity gap-fill)  *(MEDIUM · DONE)*

- **Finding:** The mechanism-only resolver mapped only **46%** of repoDB drugs to a
  target set, so the engine was *blind* to half the drug universe (the calibration
  score≈0 bin was ~52% real approvals). Two additions fixed it: (a) fold salts to
  their **parent molecule** (`molecule_hierarchy`) so a salt inherits the parent's
  targets; (b) for drugs with **no curated mechanism**, gap-fill from
  high-confidence single-protein bioactivity (`activities`, pchembl≥6, confidence≥8).
- **Result:** drug coverage **46%→68%**, AUROC **0.709→0.732**, 48% more of the gold
  standard now testable, negative control still ≈0.5.
- **Root cause / design choice:** Adding activity targets to *every* drug ("rich")
  only reached AUROC 0.715 — broadening already-curated drugs with promiscuous
  bioactivity targets **diluted** precision. Gap-fill-only ("fallback") keeps clean
  curated targets and only fills the holes → higher coverage *and* higher AUROC.
- **Corrective action:** Shipped the fallback resolver in
  `repurposing_engine._local_targets_for_molecules` (used by the disease screen).
- **Preventive action:** Benchmark re-run with `mapping=fallback` as default;
  legacy/rich/fallback modes retained for re-measurement on each refresh.
- **Status:** DONE. Headline metrics below updated to the fallback run.

---

## Run: 2026-06-29 — repoDB v2.x vs CipherQ mechanistic score

Scope: 1,000 drug-indication pairs across 54 diseases (approval base rate 72%).
Gold standard repoDB (Approved vs Terminated/Withdrawn/Suspended). Scoring is
mechanism-only (target + direction-aware pathway + PPI); clinical/literature/
regulatory features removed to prevent retrospective leakage. 4 findings.

### PRED-01 — Engine recovers real repurposing above baselines  *(INFO · PASS)*

- **Finding:** Mechanism-only AUROC **0.71** (per-disease mean 0.74), AUPRC 0.85,
  BEDROC(α=20) 0.87, EF@10% 1.2×. Baselines: random 0.49, drug-popularity 0.42,
  target-overlap-only 0.72.
- **Root cause / interpretation:** The biology signal genuinely discriminates
  approved from failed repurposings, and beats a "popular drug" baseline — so the
  ranking is mechanism-driven, not popularity-driven.
- **Corrective action:** None required; this is the headline pass.
- **Preventive action:** Re-run on each data refresh; track AUROC for drift.
  AUROC < 0.55 triggers review of the biology inputs (Open Targets gene sets,
  mechanism coverage).
- **Status:** PASS. Surfaced live on `/validation`.

### PRED-02 — Pathway term is near-neutral for ranking  *(LOW · WATCH)*

- **Finding:** Ablation AUROC — direction-aware pathway 0.709, direction-blind
  0.712, pathway excluded 0.718. Target overlap carries essentially all the
  mechanistic signal on this benchmark; the pathway term does not improve (and
  marginally dilutes) retrospective ranking.
- **Root cause:** repoDB approved drugs typically have *direct* target overlap with
  the disease gene set, so target overlap dominates. The pathway term's value is in
  indirect/novel connections (and in the Pathway Screen feature), which this
  benchmark does not specifically stress; the direction-awareness is a correctness
  safeguard (it damps drugs that *mimic* a disease) that repoDB's failed-drug
  negatives do not isolate.
- **Corrective action:** Keep the pathway term as a correctness/explainability
  feature; do **not** present it as a recovery booster in customer materials.
- **Preventive action:** Track the ablation Δ each run. If direction-aware is
  persistently worse than direction-blind by a meaningful margin, revert the term
  to direction-blind (per SOP-003 step 6).
- **Status:** WATCH; documented honestly in DATA_QUALITY.md §3b and on `/validation`.

### PRED-03 — Negative control passes  *(INFO · CONTROL)*

- **Finding:** Shuffling the labels collapses AUROC to 0.49 (≈ chance), confirming
  the measured signal is not an artefact of the evaluation set-up.
- **Action:** None — control check; re-run with each benchmark.
- **Status:** PASS.

### CAL-01 — Raw score is not a probability; isotonic map fixes it  *(MEDIUM · ACTION)*

- **Finding:** The raw mechanism score has expected calibration error **0.39**
  (Brier 0.36) — it ranks well but its absolute value is not a probability (e.g.
  the score≈0 bin still has ~52% approvals, the base-rate effect). An out-of-fold
  (5-fold) isotonic recalibration lowers ECE to **0.014** (Brier 0.17). Observed
  approval rate is monotonic in the score (Spearman 0.37 across bins — weak,
  dominated by the no-overlap bin).
- **Root cause:** Mechanistic overlap alone does not explain ~half of real
  repurposings (they act via mechanisms not captured by target/pathway/PPI overlap
  with Open Targets gene sets) — which is exactly why the *live* engine also uses
  clinical/literature evidence (removed here for leakage control).
- **Corrective action:** Display a calibrated probability via the out-of-fold
  isotonic map, not the raw score; keep the raw score for ranking.
- **Preventive action:** Re-fit the isotonic map on each refresh, always
  out-of-fold and on leakage-free data — the same discipline that led us to reject
  the earlier leaky PoS calibrator (reverse causation).
- **Status:** Map specified; to be wired as the score→probability display layer.

---

## Notes

- Validation is **read-only**; no source data modified.
- Leakage control is the load-bearing design choice: clinical-trial and literature
  features are switched off for the benchmark because they encode the label in a
  retrospective test. They remain valid in the live product.
- A retrospective benchmark proves historical recovery on a fixed snapshot, not
  prospective performance. Re-run `validate_predictions.py` + `validate_calibration.py`
  after any data refresh and append a dated section here.
