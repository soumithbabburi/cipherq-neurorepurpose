# Requirements Traceability Matrix (RTM)

**Purpose (GAMP 5 P0 item 2):** map each user requirement (URS) to the software
that implements it (FS) and to the objective evidence that verifies it (a
committed test or a validation harness + its dated result artifact). This is the
single index a pharma quality reviewer uses to confirm that every claim the
platform makes is backed by a re-runnable check.

**Verification types:**
`auto` = automated test in `tests/` (runs in the CI gate, see .github/workflows/ci.yml).
`harness` = a `validation/validate_*.py` script writing a dated result JSON (needs the
local chembl_33 DB + Open Targets cache; re-run on each data refresh).
`experiment` = a one-off measured study documented as a result JSON.
`doc` = a controlled document.

Last updated: 2026-07-27.

| URS | Requirement | FS (implementation) | Verification | Evidence artifact | Type |
|---|---|---|---|---|---|
| URS-01 | Bioactivity data is field-complete and within the ~0.5-log published noise floor | ChEMBL 33 data layer | validate_bioactivity.py | validation/validation_results_bioactivity.json | harness |
| URS-02 | Potency values agree with an independent expert DB (IUPHAR) within 1 log | Cross-DB concordance check | validate_concordance.py (CONC-01) | validation/concordance_results.json | harness |
| URS-03 | Disease terms resolve to the correct MeSH concept (consistency + curated golden set) | services/disease_resolver.py | validate_resolution.py (RES-00, RES-01) | validation/resolution_results.json | harness |
| URS-04 | Mechanism direction (inhibit vs activate) matches IUPHAR | services/mechanism_direction.py | validate_mechanisms.py (MECH-01) | validation/mechanisms_results.json | harness |
| URS-05 | The knowledge graph is provenance-tagged with referential integrity, no duplicates | services/hetionet_graph.py, provenance | validate_kg_accuracy.py (KG-01..04) | validation/kg_results.json | harness |
| URS-06 | The engine recovers real (Approved) repurposings vs failed, leakage-free | services/repurposing_engine.py (mechanism dims) | validate_predictions.py (PRED-01..03) | validation/predictions_results.json | harness |
| URS-07 | A raw composite score is made interpretable (tiered; calibration reported honestly) | services/score_calibration.py | validate_calibration.py (CAL-01) | validation/calibration_results.json | harness |
| URS-08 | Genetics-weighting is only adopted if it measurably helps recovery | services/repurposing_engine.py | validate_predictions.py (GEN-01: tested, neutral, kept opt-in) | validation/predictions_results.json | harness |
| URS-09 | KG-embedding / metapath signals are integrated ONLY if they beat mechanism on the external benchmark | services/kg_embedding.py, metapath_features.py | validate_predictions.py (KGE-04), validate_metapath.py (MP-01) — both NOT integrated | validation/predictions_results.json, metapath_results.json | harness |
| URS-10 | The validated plausibility model is not wired into the rank unless it improves it | services/repurposing_predictor.py | experiment_plausibility_rerank.py — rejected (0.34 AUROC on approved-vs-failed) | validation/plausibility_rerank_results.json | experiment |
| URS-11 | A ground-truth contraindication gates (kills) the composite score | services/primekg_predictor.py, repurposing_scorer.py | test_contraindication_kills_composite_in_canonical_scorer; test_why_not_surfaces_contraindication | tests/test_primekg_gate.py | auto |
| URS-12 | A validated indication is NOT penalized by the contraindication gate | services/repurposing_scorer.py | test_indication_not_penalized_in_canonical_scorer; test_why_not_silent_on_indication | tests/test_primekg_gate.py | auto |
| URS-13 | A safety-terminated trial produces a NEGATIVE outcome signal (never neutral/positive); verdict ties are order-independent | services/endpoint_parser.py | test_endpoint_terminated_safety_is_negative; test_endpoint_verdict_tie_is_order_independent | tests/test_cipherq.py | auto |
| URS-14 | Merely-studied (not approved) pairs keep every safety/clinical/coverage guardrail | services/repurposing_engine.py (_own_therapy) | Code review; regression via URS-11/12 gate tests | commit 9147c18 | manual + code review |
| URS-15 | No personal credentials are hardcoded in source; DB params come from env | config.py | test_no_hardcoded_personal_username_in_source; test_db_params_shape | tests/test_cipherq.py | auto |
| URS-16 | Docking selects the correct structure and lands the ligand in the named pocket | services/real_pdb_fetcher.py, dock_engine.py | test_docking_lands_in_named_pocket; test_cif_only_entry_is_converted | tests/test_cipherq.py | auto (skips without network) |
| URS-17 | When no suitable structure exists, the platform degrades honestly (no fabricated pose) | services/docking_service.py | test_no_structure_degrades_honestly | tests/test_cipherq.py | auto |
| URS-18 | The disk cache round-trips values and honors TTL expiry | services/diskcache.py | test_diskcache_roundtrip; test_diskcache_ttl_expiry | tests/test_cipherq.py | auto |
| URS-19 | Quantum descriptors come from a real QM engine (GFN2-xTB), not a stub | services/qc_engine.py | test_qc_engine_aspirin | tests/test_cipherq.py | auto (skips without xTB env) |
| URS-20 | Every output is traceable to the exact code + data versions that produced it | services/provenance.py (build_stamp, attached to lineage) | build_stamp() resolves the running commit SHA + source snapshots | services/provenance.py | manual (verified 2026-07-27) |
| URS-21 | Every candidate shows its reviewable basis alongside any score (CDS criterion 4) | templates/*, services/evidence_dossier.py, provenance | Design invariant; positioning rationale | validation/REGULATORY_POSITIONING.md | doc |
| URS-22 | The platform makes no claim of clinical-success prediction | services/repurposing_predictor.py (scope note), UI copy | Stated limit: approved-vs-failed AUROC 0.42 | validation/predictions_results.json, repodb_external_metrics.json | harness + doc |
| URS-23 | Operator actions are recorded in a tamper-evident, time-stamped audit trail that survives rotation (Part 11) | services/audit_trail.py; flask_app after_request hook | test_audit_trail_chain_and_tamper; test_audit_trail_rotation_keeps_chain | tests/test_cipherq.py | auto |
| URS-24 | System access is limited to authenticated users with role-based authority (viewer/analyst/admin) (Part 11) | services/auth.py; flask_app before_request guard + /login + /admin/users | test_auth_password_and_roles; test_auth_role_policy_and_user_listing | tests/test_cipherq.py | auto |

## Coverage summary

* Automated (CI-gated) tests: 21 across `tests/test_cipherq.py` and
  `tests/test_primekg_gate.py`.
* Validation harnesses: 10, each writing a dated result JSON.
* Documented negative results (tested and deliberately NOT integrated): KGE-04
  (KG embedding), MP-01 (metapath/DWPC), URS-10 (plausibility re-rank). Recording
  what was rejected, and why, is part of the evidence.

## Honesty notes

* URS-14 (guardrail suppression) is verified by code review; a dedicated automated
  test needs the full DB-backed scoring path, so it is a recommended follow-up (the
  outcome-signal logic in URS-13 is now covered by CI-gated unit tests).
* Harness rows require the local database and cache; they are re-run on each data
  refresh, and each run appends a dated entry to the relevant CAPA log. They do not
  run in the CI gate (no DB in CI).
