# CAPA Log — Bioactivity Data Validation (chembl_33)

**Corrective and Preventive Action log for SOP-002.**
Every finding from `validate_bioactivity.py` is recorded here with root cause,
corrective action, and preventive action. Confessed and documented findings are
the evidence that the data is validated — not a sign that it is bad.

Severity: HIGH · MEDIUM · LOW · INFO (by-design, documented).

---

## Run: 2026-06-16 — ChEMBL_33 (released 2023-05-31)

Scope: 4,104,166 affinity rows (Ki / IC50 / EC50 / Kd). 2 findings.

### VAL-B01 — ChEMBL-flagged suspect measurements present  *(MEDIUM · OPEN)*

- **Finding:** ~3,300 affinity rows carry a non-null `data_validity_comment`
  (ChEMBL's own "this value looks wrong" flag): Ki 1,802, IC50 1,365, EC50 80,
  Kd 66.
- **Root cause:** Inherited from the upstream source. ChEMBL curators flag values
  that fail their automated/manual sanity checks but retain them for completeness.
- **Corrective action:** Certified-tier affinity queries filter
  `data_validity_comment IS NULL`. Flagged rows are **excluded from the certified
  value, not deleted** — they remain in `activities` and stay queryable.
- **Preventive action:** The certified data layer applies this filter by default;
  any consumer wanting flagged rows must opt in explicitly.
- **Status:** Filter specified; to be enforced in the certified extract (Phase 3 build).

### NOISE-B01 — High-disagreement compound-target pairs identified  *(INFO · BY DESIGN)*

- **Finding:** Replicate measurements for some compound-target pairs disagree by
  more than 1 log. Share of pairs with > 1 log spread: Ki 3.27%, IC50 5.69%,
  EC50 4.33%, **Kd 16.38%**. Worst observed examples include
  METHOTREXATE / DHFR (pKi range 5.6–10.7, ~5 log) and SIROLIMUS / FKBP1A
  (5.3–9.7). Measured mean intra-pair spread per endpoint: Ki 0.19, IC50 0.30,
  EC50 0.28, Kd 0.42 log — all within the published ~0.5-log experimental floor.
- **Root cause:** Genuine inter-laboratory / inter-assay experimental variability
  in public affinity data — expected, not a defect of our pipeline.
- **Corrective action:** For multi-measured pairs the certified value is the
  **median**, with the spread carried as an **uncertainty band** (not a false-
  precise point). Disagreeing pairs (> 1 log) are surfaced with N and range.
- **Preventive action:** Affinities propagate to docking / PBPK / scoring as
  median ± spread so downstream outputs inherit the uncertainty instead of hiding it.
- **Status:** Aggregation + uncertainty-propagation rule documented; raw rows retained.

### OBS-B01 — Kd is the noisiest endpoint  *(LOW · WATCH)*

- **Observation:** Kd shows the highest disagreement (16.4% of pairs > 1 log;
  p95 spread 1.60 log) and the lowest pchembl coverage trend, vs Ki/IC50/EC50.
- **Action:** When multiple endpoints exist for a target, prefer Ki/IC50 over Kd
  for certified affinity unless Kd is the only available measure; note Kd's wider
  band where used. No correction to source data.
- **Status:** Documented as a preference rule for the certified layer.

---

## Notes

- Validation is **read-only**; the source `chembl_33` database was not modified.
- All exclusions are *flags*, not deletions. Raw `activities` rows are intact and
  can be re-included by changing the certified-layer filter.
- Re-run `python validation/validate_bioactivity.py` after any ChEMBL refresh to
  regenerate these numbers; append a new dated section here.
