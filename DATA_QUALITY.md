# Data Quality Scorecard

**Unified data-validation summary for the shared data platform**
**Last updated:** 2026-06-16
**Standards:** ALCOA+ (FDA / WHO / EU GMP Annex 11), GAMP 5
**Owner:** Solix Technologies — Data Team

This document is the single, auditable summary of how our underlying data is
validated. It does **not** claim the data is perfect — it states, per dataset and
per quality dimension, what we measured and what we do about what we found.

> **Why one document for two platforms.** The POZ Platform and CipherQ read the
> **same** `chembl_33` database on `localhost:5433`. We therefore validate the
> shared data **once** and certify it for both. This is also why the two
> platforms cannot show conflicting numbers for the same drug: there is one
> certified source, not two independent cleanings.

---

## 1. Dataset inventory

| Dataset | Source (authority) | Version / date | Size | Used by | Validation |
|---|---|---|---|---|---|
| ChEMBL bioactivity (Ki/IC50/EC50/Kd) | ChEMBL 33 (EMBL-EBI) | ChEMBL_33 · 2023-05-31 | 4,104,166 rows | POZ + CipherQ | **SOP-002** ✓ (2026-06-16) |
| ChEMBL Cmax universe | ChEMBL 33 + FDA labels | ChEMBL_33 · 2023-05-31 | 20,212 drugs | POZ + CipherQ | **SOP-001** ✓ (2026-06-02) |
| Hetionet edges | Hetionet v1.0 (het.io) | v1.0 · loaded 2026-06-16 | 2,250,197 real + 11,466 ChEMBL-derived · 55,926 nodes | CipherQ | Provenance-tagged ✓; full validation pending |
| MeSH diseases | NLM MeSH | loaded subset | 1,938 diseases | CipherQ | Pending (module 3) |
| Indications | ChEMBL / curated | — | 733 | CipherQ | Pending |
| Mechanisms (action_type) | ChEMBL | ChEMBL_33 | 10,733 | CipherQ | **MECH-01** ✓ (2026-06-29) — 99.3% direction-concordant vs IUPHAR |
| Orange Book | FDA Orange Book | products 48,215 | 48,215 products | POZ + CipherQ | Pending |
| PDB structures (cache) | RCSB PDB | live-fetched | ~11 cached | CipherQ | Pending (per-fetch validated in pipeline) |

---

## 2. Scorecard — six dimensions × dataset

Legend: ✓ measured & within criterion · ⚠ measured, finding logged (see CAPA) ·
○ not yet validated (module pending).

| Dataset | Completeness | Accuracy | Consistency | Validity | Uniqueness | Timeliness |
|---|---|---|---|---|---|---|
| **ChEMBL bioactivity** | ✓ pchembl 64.9% / units 91.8% / conf 100% | ✓ 91% concordant vs IUPHAR/GtoPdb; noise ≤0.42 log (within 0.5 floor) | ✓ endpoints kept separate | ⚠ 58.3% high-conf; suspect rows flagged (VAL-B01) | ✓ dedup→median (raw kept) | ✓ ChEMBL_33 (2023-05-31) |
| **ChEMBL Cmax** | ✓ full accounting reconciles | ✓ spot-audit 100% (10 drugs) | ✓ salt→parent collapse | ✓ all values source-tagged | ✓ parent dedup | ✓ ChEMBL_33 |
| **Hetionet edges** | ✓ 2.25M edges, 24 metaedge types (KG-01) | ✓ integration of validated sources (KG-04) | ✓ 100% referential integrity, 0 self-loops (KG-02) | ✓ `source` tag 100% (hetionet_v1.0 vs chembl_derived) | ✓ 0 duplicate edge groups (KG-03) | ✓ v1.0 (loaded 2026-06-16) |
| **MeSH diseases** | ✓ 1,938 loaded | ✓ golden-set 100% correctness | ✓ round-trip heading 100% / synonym 99.7% | ✓ entry-term mapping | ✓ descriptor unique | ✓ NLM MeSH |
| **Mechanism direction** | ✓ 10,733 mechanisms | ✓ 99.3% direction-concordant vs IUPHAR (MECH-01) | ✓ signed −1/+1 mapping | ✓ action_type curated | ✓ per (drug,target) | ✓ ChEMBL_33 |
| **Orange Book** | ○ | ○ | ○ | ○ | ○ | ○ FDA edition TBD |

---

## 3. Headline result — ChEMBL affinity noise, measured

The most common pharma objection is *"even ChEMBL has noise and discrepancies."*
We measured it. Across high-confidence (assay confidence ≥ 8), validity-clean,
replicate-measured compound–target pairs:

| Endpoint | Pairs | Median spread | Mean spread | % pairs disagree > 1 log | vs 0.5-log floor |
|---|---|---|---|---|---|
| Ki | 37,796 | 0.02 log | 0.19 log | 3.3% | within |
| IC50 | 111,186 | 0.16 log | 0.30 log | 5.7% | within |
| EC50 | 17,065 | 0.15 log | 0.28 log | 4.3% | within |
| Kd | 8,380 | 0.19 log | 0.42 log | 16.4% | within (noisiest) |

Interpretation: our measured reproducibility is consistent with the published
~0.5-log experimental floor for public affinity data (Kramer et al. 2012). We do
not hide the disagreeing pairs — we **flag them, certify the median, and carry the
spread as an uncertainty band** so downstream docking/PBPK/scoring inherit the
uncertainty rather than a false-precise point value.

---

## 3b. Predictive validation — does the engine recover real repurposing? (2026-06-29)

The sections above validate the **data**. This validates the **results**. Gold
standard: **repoDB** (Brown & Patel, *Sci Data* 2017) — approved drug-indication
pairs as positives, pairs that reached trials and **failed** (Terminated /
Withdrawn / Suspended) as hard negatives.

**Leakage control (the honest part).** The live engine weights clinical-trial
signal (30%) and literature co-mention — in a retrospective test those *are* the
label. So the benchmark scores **biology only**: target overlap + direction-aware
pathway + PPI proximity, with clinical / literature / regulatory features removed.

| Metric | Result |
|---|---|
| AUROC (approved vs failed), mechanism-only | **0.73** (per-disease mean **0.73**) |
| BEDROC (α=20, early recognition) | 0.85 |
| Baselines | random 0.49 · drug-popularity 0.43 · target-only 0.75 |
| Negative control (label shuffle) | 0.48 (≈ chance) ✓ |
| Drug-target coverage | **68%** of repoDB drugs scored (was 46% mechanism-only — see COV-01) |
| Pairs / diseases evaluated | 1,476 pairs · 55 diseases (base rate 67%) |
| Score calibration | raw ECE 0.38 → **isotonic OOF ECE 0.013**; monotonic in score |

**Honest findings (logged as CAPA COV-/PRED-/CAL-):**
- **Coverage fix (COV-01, 2026-06-29):** adding salt→parent mapping and a
  high-confidence bioactivity *gap-fill* (only for drugs with no curated mechanism)
  lifted drug coverage 46%→68% and AUROC 0.71→0.73. Broadening *every* drug with
  activity targets instead diluted precision (0.715), so the engine gap-fills only.
- The mechanistic score genuinely separates real from failed repurposings, beats
  random/popularity baselines, and the label-shuffle control collapses to chance.
- **Metapath/DWPC features (MP-01, 2026-06-29) — tested, NOT integrated:** the
  Rephetio method (degree-weighted path counts over Hetionet; ~0.97 in the
  literature). In-graph it works (single metapaths 0.6–0.7 reconstructing
  Hetionet's own treats), but on the **external repoDB benchmark** it failed:
  best single metapath ≈ 0.53 (chance), full CV logistic-regression model 0.59,
  vs **mechanistic 0.75** on the same pairs; ensembling hurt (0.73). Same pattern
  as the KG embedding — in-graph ≠ external. Root cause: Hetionet (sparse, curated
  2017; many drugs have zero binding edges) is a weaker biology source than the
  live Open Targets + ChEMBL data the mechanistic score already uses. **Both
  graph-based levers are ruled out by the external benchmark.**
- **Genetics-weighted disease genes (GEN-01, 2026-06-29) — tested, neutral:**
  weighting disease genes by Open Targets genetic_association (Nelson 2015) scored
  AUROC 0.730 vs flat 0.735 (−0.005) — neutral for retrospective recovery. Useful
  side-result: the tiny gap confirms the headline 0.73 is *not* meaningfully
  inflated by the overall-score known-drug channel (a leakage reassurance). Genetics
  enrichment is now fetched per disease and available as an opt-in, leakage-safer
  target weighting; kept OFF as the engine default (no measured benefit).
- **KG-embedding (KGE-01/03/04, 2026-06-29) — tested, NOT integrated:** a DistMult
  embedding of the Hetionet graph predicts *Hetionet's own* held-out treat edges at
  AUROC 0.67–0.69, and an in-graph gene-overlap+KG ensemble looked complementary
  (+0.042). **But on the external repoDB benchmark it failed:** KG alone scored
  AUROC 0.54 (≈ chance) at discriminating real approved-vs-failed repurposings, and
  ensembling it *degraded* the mechanistic score (0.757 → 0.703, **−0.053**). The
  in-graph gain did not generalise. **Decision: the KG embedding is NOT wired into
  the engine** — it would be a regression. Two benchmarks (in-graph + external) and
  the external one is decisive. Revisit only with metapath/DWPC features. This is
  the validation framework preventing a regression, not a setback.
- **Target overlap carries most of the signal**; the pathway term (direction-aware
  vs blind vs excluded ≈ 0.71/0.71/0.72) is near-neutral *for ranking* — it is kept
  as a correctness/explainability feature, not a recovery booster.
- The **raw score is not a probability** (ECE 0.39); an out-of-fold isotonic map
  makes it one (ECE 0.014). Calibrated on leakage-free data only — the same
  discipline that led us to reject the earlier leaky PoS calibrator.
- A retrospective benchmark proves **historical recovery on a fixed snapshot**, not
  prospective performance. That boundary is stated, not hidden.

Protocols: `validation/validate_predictions.py`, `validation/validate_calibration.py`
(SOP-003). Live on `/validation`.

---

## 4. Open findings (see CAPA logs)

| ID | Sev | Title | Disposition |
|---|---|---|---|
| VAL-B01 | MEDIUM | ChEMBL-flagged suspect affinity rows (~3,300) | Excluded from certified tier; raw retained |
| NOISE-B01 | INFO | High-disagreement pairs (e.g. METHOTREXATE/DHFR ~5 log) | Median ± spread; flagged with N and range |
| OBS-B01 | LOW | Kd is the noisiest endpoint (16.4% > 1 log) | Prefer Ki/IC50 where available |
| ACC-01 (SOP-001) | MEDIUM | Out-of-range Cmax (unit errors) | Under review per SOP-001 CAPA |

Full detail: `validation/CAPA-LOG-bioactivity.md`, `validation/CAPA-LOG.md` (POZ).

---

## 5. Known gaps / roadmap (honest)

1. **Hetionet edges — RESOLVED 2026-06-16.** The real Hetionet v1.0 edge file had
   failed to download (`edges.sif.gz` was 133 bytes; `edges.tsv` 0 bytes) and the
   importer had silently fallen back to ~11k ChEMBL-derived edges (only CbG/CtD).
   Re-downloaded the genuine 12 MB file, loaded **2,250,197 edges across 24
   metaedge types**, and added a `source` column (`hetionet_v1.0` vs
   `chembl_derived`). The repurposing scorer now reads **only** `hetionet_v1.0`
   and follows the indirect **Compound→Gene←Disease** path (not just direct CtD),
   so network evidence is real, orthogonal to ChEMBL, and explainable (returns the
   connecting genes). Remaining: map drug names/synonyms to HetioNet via DrugBank
   IDs so more compounds resolve (e.g. "Aspirin" → "Acetylsalicylic acid").
2. **MeSH entity-resolution audit — DONE 2026-06-17.** Round-trip identity over a
   300-descriptor sample of 1,938 MeSH diseases: heading 100%, **synonym 99.7%**
   (the silent-killer check), alias 100% (after fixing the alias map — was 60%),
   ambiguity 0%. Protocol: `validation/validate_resolution.py`; live on `/validation`.
3. **Cross-database concordance — DONE 2026-06-17.** ChEMBL affinities vs the
   independently-curated IUPHAR/Guide to Pharmacology (GtoPdb v2026.2):
   **91.0% of 5,754 shared compound-target pairs agree within 1 log** (median
   absolute difference 0.03 log). By endpoint: Ki 92.5%, IC50 90.8%, EC50 87.3%,
   Kd 91.0%. Protocol: `validation/validate_concordance.py`; result surfaced live
   on the `/validation` page. This is the external "accuracy" number.
4. **Orange Book + PDB** — formal completeness/timeliness checks.
5. **Apply the suspect-row and median±spread filters in the certified extract**
   so consumers read certified values by default.

---

## 6. How to reproduce (read-only)

```bash
# DATA validation — SELECT only, never modifies data
python validation/validate_bioactivity.py    # completeness, validity, noise (section 3)
python validation/validate_concordance.py     # accuracy vs IUPHAR/GtoPdb (91%)
python validation/validate_resolution.py       # disease-term resolution (consistency + golden set)
python validation/validate_mechanisms.py       # mechanism DIRECTION vs IUPHAR (99.3%)
python validation/validate_kg_accuracy.py      # knowledge-graph integrity & provenance

# RESULTS validation (section 3b) — needs data/external/repodb/full.csv
python validation/validate_predictions.py      # retrospective recovery vs repoDB (AUROC)
python validation/validate_calibration.py      # score → probability calibration

# Cmax (POZ SOP-001)
python validation/validate_data.py   # in the POZ repo
```

Outputs: `validation/validation_results_bioactivity.json` + run log. Re-run after
any ChEMBL refresh and append a dated section to the CAPA log.

---

*Validation is additive and read-only. No source data is modified or deleted;
exclusions are flags, and raw rows are always retained. A pre-validation backup of
all data was taken 2026-06-16 (see `poz_cipherq_backups/`).*
