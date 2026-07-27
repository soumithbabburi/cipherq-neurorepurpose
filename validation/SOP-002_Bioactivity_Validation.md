# SOP-002 — Shared Bioactivity Data Validation (chembl_33)

**Standard Operating Procedure**
**Document ID:** SOP-002
**Version:** 1.0
**Effective date:** 2026-06-16
**Owner:** Solix Technologies — Data Team
**Applies to:** All ChEMBL bioactivity/affinity data (Ki, IC50, EC50, Kd) in the
shared `chembl_33` database (localhost:5433) consumed by **both** platforms:
the POZ Platform and CipherQ (docking, PBPK, repurposing scoring).

> Companion to **SOP-001** (which covers the Cmax universe). SOP-001 + SOP-002
> together validate the full shared dataset. Both platforms read the same
> `chembl_33` database, so this validation is performed **once** and certifies
> the data for both — eliminating cross-platform discrepancies by construction.

---

## 1. Purpose

Affinity data is the part of ChEMBL a pharma reviewer means when they say "it has
noise and discrepancies." This SOP does **not** claim the data is clean. It
**measures** the noise using ChEMBL's own quality fields and bounds it against
the published experimental reproducibility floor, producing documented evidence
an auditor can re-run.

It answers the three recurring partner questions, specifically for affinities:

1. *How noisy is your affinity data, exactly?* → measured intra-pair spread per endpoint.
2. *How do you handle the discrepancies?* → flag, aggregate to median ± spread, keep raw.
3. *How do I know nothing was thrown away?* → additive only; raw rows retained, exclusions documented.

## 2. Scope

- **In scope:** Every `activities` row with `standard_type` ∈ {Ki, IC50, EC50, Kd}
  in `chembl_33` (~4.1M rows). Each endpoint is validated in its **own lane** —
  IC50/Ki/EC50/Kd are never mixed.
- **Out of scope:** Cmax data (SOP-001); non-affinity endpoints (Potency, AC50,
  ADMET) — to be added in later revisions; ML-derived values (validated separately).

## 3. Standards Applied

| Framework | Role |
|---|---|
| **ALCOA+** (FDA / WHO / EU GMP Annex 11) | Data-integrity principle each check maps to |
| **GAMP 5** (ISPE) | `validate_bioactivity.py` is the executable OQ/PQ test protocol |

**ALCOA+** = **A**ttributable, **L**egible, **C**ontemporaneous, **O**riginal,
**A**ccurate, **+ C**omplete, **C**onsistent, **E**nduring, **A**vailable.

## 4. Reference Standards (authoritative)

| Purpose | Authority |
|---|---|
| Bioactivity, structures, targets | ChEMBL 33 (EMBL-EBI), released 2023-05-31 |
| Target-assignment confidence | ChEMBL `assays.confidence_score` (≥ 8 = single direct protein) |
| Suspect-value flags | ChEMBL `data_validity_comment`, `potential_duplicate`, `standard_relation` |
| Experimental noise floor | ~0.5 log units (pKi/pIC50); Kramer et al. 2012, *J Med Chem* |
| Cross-database accuracy | BindingDB, IUPHAR/Guide to Pharmacology |

## 5. Procedure

Run the validation protocol (read-only — SELECT only, never modifies the DB):

```bash
# credentials read from .env (DB_USER / DB_PASSWORD) or environment
python validation/validate_bioactivity.py
```

Produces `validation_results_bioactivity.json` and `validation_run_bioactivity.log`.
Checks, each mapped to an ALCOA+ principle and an acceptance criterion:

| Step | ALCOA+ principle | Check | Acceptance criterion |
|---|---|---|---|
| 0 | Available / Enduring / **Timeliness** | DB reachable; ChEMBL version + release date stamped | `activities` > 0; version recorded |
| 1 | — | Affinity universe counted **per endpoint** (no mixing) | Counts reproducible; endpoints separate |
| 2 | **Attributable** | Every measurement → assay → target (tid) | 0 with no assay link; untargeted rows documented |
| 3 | **Accurate / Valid** | ChEMBL native flags: confidence_score, data_validity_comment, censored relation, potential_duplicate, pchembl range | High-confidence share reported; flagged rows excluded from certified tier, not deleted |
| 4 | **Accurate** | **Reproducibility** — intra-pair log spread per endpoint; % pairs disagreeing > 1 log | Mean spread ≤ 0.5 log floor; disagreeing pairs flagged |
| 5 | **Consistent / Unique** | Duplicate measurements per pair → aggregate to median (raw retained) | Dedup ratio reported; aggregation rule documented |
| 6 | **Complete** | Field coverage: pchembl, units, confidence | Coverage reported per field |

## 6. Handling Findings

Every discrepancy is logged in **CAPA-LOG-bioactivity.md** with finding ID, root
cause, corrective action, and preventive action. **Finding an error is not a
failure — hiding one is.** A confessed, measured error rate is the evidence of a
validated system.

Severity: **HIGH** (clinical-tier customer-facing) · **MEDIUM** (accuracy,
limited exposure) · **LOW** (minor) · **INFO** (by-design, documented).

## 7. Frequency

- On every ChEMBL refresh / re-load.
- Before any customer-facing affinity claim or pitch-deck number.
- Quarterly standing data-integrity review.

## 8. Records

| Record | Location |
|---|---|
| Executable protocol | `validation/validate_bioactivity.py` |
| Machine-readable results + scorecard | `validation/validation_results_bioactivity.json` |
| Run transcript | `validation/validation_run_bioactivity.log` |
| Deviation / CAPA log | `validation/CAPA-LOG-bioactivity.md` |
| Unified scorecard (all datasets) | `DATA_QUALITY.md` (repo root) |

## 9. Revision History

| Version | Date | Author | Change |
|---|---|---|---|
| 1.0 | 2026-06-16 | Solix Data Team | Initial issue; first full validation of 4.1M Ki/IC50/EC50/Kd rows in chembl_33 |
