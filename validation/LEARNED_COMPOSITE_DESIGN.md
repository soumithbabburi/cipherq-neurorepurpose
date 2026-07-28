# Scope: a learned composite on the HARD labels (the capstone fix)

**Goal.** Replace the hand-weighted evidence composite with a supervised model trained on the
task that actually matters — repoDB **approved-vs-FAILED** — evaluated **disease-disjoint**, so it
(a) captures the target-overlap signal W1 couldn't ship (blocked on the mechanism-mass ceiling),
(b) does honest feature selection across the bonus stack W3 couldn't ablate piecemeal, and
(c) outputs a calibrated 0-1 probability that dissolves the tier-band problem.

Status: DESIGN ONLY. Nothing is trained or wired by this document.

## 1. Why the previous learned composite (services/train_composite.py) lost — and why that was the wrong test

- It trained AND evaluated on **Hetionet CtD positives vs RANDOM negatives** (train_composite.py:51-76)
  — an EASY task the hand composite was effectively designed around (AUC_hand 0.969). The learned
  model got 0.927 and was correctly not saved.
- But CtD-vs-random is not the decision task. The task that matters is **approved vs tried-and-FAILED
  repurposings** (hard negatives), where the hand composite is KNOWN to fail (AUROC ~0.41-0.42,
  data/repodb_external_metrics.json). The learned model was never trained or scored on it.
- Its learned weights (target +1.56, pathway -0.61 — data/composite_model_metrics.json) were actually
  CORRECT in direction (they independently reproduce our pathway-is-dilutive finding). The exercise
  measured the wrong benchmark, not a bad model.

## 2. The design

### 2.1 Labels
repoDB (data/external/repodb/full.csv): positives = status "Approved"; negatives = "Terminated"/
"Withdrawn"/"Suspended" (hard negatives). Same loader as validate_predictions.load_repodb.

### 2.2 Features — NON-LEAKY evidence only
A retrospective approved-vs-failed benchmark leaks through any feature that encodes phase/approval.
So the learned model uses ONLY the leakage-free evidence sub-signals:

    target, pathway, ppi, network, directional, proliferation, signature, direction    (8 features)

EXCLUDED and kept OUTSIDE the learned model:
- clinical, indication, regulatory — leaky (encode the label); remain ADDITIVE prior-art in the
  composite exactly as today (legitimate for forward scoring, invalid for retrospective eval).
- DWPC plausibility — circular (itself trained on treats); stays its own axis.
- All PENALTIES (safety, CCH, coverage, phantom, registry, trial-failure, harmful-direction,
  PrimeKG contraindication) — gates on evidence, remain outside (as composite_model.py already does).

This feature set is precisely what answers W3: the model assigns each bonus a learned weight, so
"does signature/proliferation/network help?" is decided by the fit, not by hand.

### 2.3 Cross-validation — DISEASE-disjoint (the honest generalization test)
GroupKFold grouped by DISEASE (not compound). This is the setting that exposed the KG zero-shot
collapse; it prevents the "wins on home data, loses externally" trap (KGE-04). Report per-fold AUROC
+ a bootstrap CI on the gap vs the hand composite. Optionally a nested compound-AND-disease-disjoint
check as a stress test.

### 2.4 Model + calibration
Logistic (interpretable, regularized, class_weight balanced) AND GBM; choose by disease-disjoint CV
AUROC. Then calibrate (isotonic/Platt, out-of-fold) so the output is a true P(0-1). The calibrated P
REPLACES the mechanism+bonus evidence term; prior-art (clinical/indication/regulatory) and penalties
stay additive/multiplicative on top, unchanged.

### 2.5 Tiers — the W1 problem DISSOLVES
Because the learned output is a genuine 0-1 probability (not a mechanism-mass-capped sum), the tier
bands are re-anchored on the P distribution with no 0.65 ceiling artifact. This is the clean home for
W1's target-overlap win that simple re-weighting could not calibrate.

## 3. Ship criterion (repo discipline)
SAVE + wire the learned composite ONLY if its disease-disjoint CV AUROC on repoDB approved-vs-failed
beats BOTH the current hand composite AND a target-heavy re-weighted hand variant by >= +0.01 with a
bootstrap 95% CI excluding zero. Otherwise keep hand weights and document the negative result (as the
previous attempt honestly did). A "no significant improvement" outcome is valid and likely-possible
given the small hard-label set.

## 4. Honest constraints and risks
- **Feature-generation cost is the long pole.** Building the training matrix runs the FULL engine
  (canonical_pair_score) per pair to get the bonus features, which hits per-pair external APIs
  (proliferation -> mygene/PubMed, safety not needed here, network -> local DB). ~1,476 pairs (or the
  ~574 fully-covered subset). This is a one-time OFFLINE run (hours, cache to data/); it is the main
  effort. Cache the feature matrix so training/iteration is fast afterwards.
- **Small hard-label set** (~422 approved / 152 failed covered) => limited statistical power; the CI
  may straddle zero => no ship. Honest and acceptable.
- **Restructuring risk.** Separating the learned-evidence term from the additive prior-art in the
  composite assembly (repurposing_engine.py:951-1057) touches validated behavior. Guard with the
  existing harnesses + tests/test_primekg_gate.py before wiring.
- **Overfitting.** Mitigated by disease-disjoint CV, small regularized feature set, bootstrap CI, and
  the hard ship bar.

## 5. Phased plan
- **Phase 0 (long pole):** build & cache the leakage-free feature matrix — per repoDB pair: the 8
  non-leaky features + approved/failed label + disease group. New validation/build_composite_dataset.py.
- **Phase 1:** train logistic + GBM, disease-disjoint CV, calibrate, bootstrap the gap vs hand.
  DECISION GATE: proceed only if it clears the ship criterion.
- **Phase 2 (only if Phase 1 wins):** restructure the composite (learned evidence + additive
  prior-art + penalties), re-anchor tiers on the calibrated P, regression-test (harnesses + gate
  tests), and wire behind the existing composite_model.learned_composite hook.

## 6. What this resolves
W1 (target re-weight, calibratable via the P distribution), W3 (bonus feature selection by the fit),
and the core "hand-weighted, un-learned composite" weakness — all with one measured, honest artifact,
IF it beats the baseline on the task that matters. If it doesn't, we keep the validated hand composite
and have proven, honestly, that the hand weights are not beatable on the available hard data.
