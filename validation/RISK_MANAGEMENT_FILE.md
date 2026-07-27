# Risk Management File (ISO 14971 scaffold)

**Scope:** RepurposeIQ, positioned as a non-device research / decision-support
tool (see validation/REGULATORY_POSITIONING.md). ISO 14971 is written for medical
devices; a non-device research tool is not obligated to it, but pharma quality
teams expect a hazard analysis, and applying the method is good practice. This is
a working scaffold, not a completed, signed risk file; a completed file requires
a named risk owner and periodic review under a QMS.

Last updated: 2026-07-27.

## 1. Intended use and reasonably foreseeable misuse

Intended use: surface and rank drug-repurposing hypotheses for human experts,
with the reviewable basis always shown. Foreseeable misuse: treating a ranking or
an "actionable" flag as a clinical or investment decision WITHOUT reviewing the
basis, or as a prediction of clinical success (which the platform explicitly does
not provide, approved-vs-failed AUROC 0.42).

## 2. Harm model (indirect, for a research tool)

The platform does not touch a patient, so it causes no direct physical harm. The
relevant harms are indirect: a wrong or overstated output could (a) misdirect R&D
spend or portfolio decisions, (b) mislead a 505(b)(2) feasibility assessment, or
(c) erode a partner's trust if a number cannot be reproduced or traced. Severity
is therefore scored in terms of decision impact, not patient injury.

## 3. Risk estimation scale

Severity: Low (cosmetic / easily caught), Moderate (could bias one decision),
High (could systematically misprioritize or mislead a partner).
Probability before controls: Low / Medium / High.

## 4. Hazard analysis (grounded in observed hazards)

Each row below is a REAL hazard identified in the 2026-07-27 review, its control,
and the objective evidence of the control. This is the value of a hazard log:
findings are recorded with cause and mitigation, not hidden.

| ID | Hazard / hazardous situation | Sev | Prob (pre) | Risk control | Evidence | Residual |
|---|---|---|---|---|---|---|
| H-01 | Reflected XSS: a crafted URL executes script in a reviewer's browser (integrity/access) | High | Medium | request args rendered via `\|tojson`, not `\|safe`; hostile-payload render test | commit ac1dd12; VALIDATION_NARRATIVE | Low |
| H-02 | A safety-terminated trial scored as neutral or positive, understating a real safety signal | High | Medium | terminated_safety branch added to the outcome signal; CI regression test | commit 5b7732b; test_endpoint_terminated_safety_is_negative | Low |
| H-03 | A merely-studied (not approved) pair has its safety/coverage guardrails suppressed and over-ranks | High | Medium | `_own_therapy` made approval-only; studied pairs keep every gate | commit 9147c18 | Low-Med (add DB-backed regression test) |
| H-04 | "Actionable" flag fires for weak candidates because the calibrated percentile saturates | High | High | actionability gated on the raw well-separated bands, not the saturated percentile | commit 2f59caa; validated vs the real null | Low |
| H-05 | CNS-penetration penalty silently disabled by a positional-arg bug (peripheral drug looks brain-penetrant) | Moderate | Medium | `cns_mpo(smiles=...)` keyword call; verified it returns a real score | commit 5a976c4 | Low |
| H-06 | A validated-looking model degrades the ranking if wired in naively (plausibility is anti-predictive on approved-vs-failed) | High | Medium | measured before shipping; NOT wired into rank (0.34 AUROC) | validation/plausibility_rerank_results.json | Low |
| H-07 | Output cannot be reproduced or traced to the data/code that produced it | Moderate | Medium | reproducibility stamp (code SHA + data versions) on every lineage record | services/provenance.build_stamp | Low |
| H-08 | Unauthorized access to screens or records; unattributable actions | High | Medium (multi-user) | authentication + role authority + tamper-evident audit trail (gated) | services/auth.py, audit_trail.py; tests | Low (when AUTH_ENABLED) |
| H-09 | Silent tampering with the audit record | High | Low | hash-chained entries + rotation chain hand-off; verify() detects any edit | test_audit_trail_chain_and_tamper, _rotation_keeps_chain | Low |
| H-10 | A user relies primarily on a score without reviewing its basis (CDS criterion 4 / misuse) | High | Medium | basis (mechanism, evidence, provenance) always co-presented; explicit no-clinical-success claim | REGULATORY_POSITIONING URS-21/22 | Medium (enforce as a permanent UI invariant) |
| H-11 | Data staleness: a decision made on an out-of-date source snapshot | Moderate | Medium | provenance freshness + integrity scoring; snapshot dates grounded in real mtimes | services/provenance.py | Low-Med |

## 5. Residual risk and the honesty position

After controls, the highest residual risks are behavioral, not code: H-10 (a user
over-relying on a score) and H-03 (needs a DB-backed regression test). These are
managed by the permanent design invariant that the reviewable basis is always
shown, by the explicit no-clinical-success claim, and by the stated boundary that
retrospective validation is not a prospective guarantee. No residual risk is
mitigated by a claim the platform cannot support.

## 6. Open items

* Assign a named risk owner and a review cadence (QMS activity).
* Add a DB-backed regression test for H-03 (guardrail suppression).
* Make H-10's "basis always shown" an enforced, tested UI invariant, not only a
  design convention.
