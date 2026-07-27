# Regulatory Positioning and Compliance Roadmap

**Scope:** RepurposeIQ (drug-repurposing hypothesis and triage platform).
**Lane chosen:** Research / decision-support (non-device).
**Status of this document:** a positioning and gap-analysis working paper. It is
NOT a compliance certificate, an FDA determination, or legal/regulatory advice.
It states the intended use, assesses the non-device rationale against the
governing FDA guidance, and lists the concrete work still required. A formal
determination requires qualified regulatory counsel and, for any regulated
deployment, an external audit.

Last updated: 2026-07-27.

---

## 1. Intended Use Statement

RepurposeIQ is a research and decision-support tool for drug-repurposing
scientists, translational researchers, and business-development teams. It
surfaces and ranks candidate drug and indication pairings from public
biomedical data (ChEMBL, a curated knowledge graph, Open Targets, MeSH) and
explains the mechanistic and evidentiary basis for each candidate.

The platform:

* produces hypotheses and a triage ranking for human experts to review;
* always presents the reasoning behind a candidate (target and pathway
  overlap, network proximity, clinical and literature signals, provenance);
* keeps a human expert in the loop for every downstream decision.

The platform does NOT:

* diagnose, treat, cure, mitigate, or prevent disease in an individual patient;
* generate patient-specific output at the point of care;
* replace the judgment of a clinician or a regulatory reviewer;
* claim probability of clinical success (see the validated limit: approved
  versus failed AUROC 0.42 in validation/predictions_results.json).

## 2. Why the platform is a non-device in this lane

There are two independent arguments, in order of strength.

### 2.1 Primary argument: it is a research informatics tool, not patient care

A repurposing platform used for drug research and development, portfolio
selection, and 505(b)(2) feasibility does not make patient-specific clinical
decisions. Its output is a research hypothesis and a commercial or scientific
prioritization, not a recommendation about the care of an identified patient.
Software that supports drug discovery and development, rather than the
diagnosis or treatment of a patient, falls outside the medical-device
definition in section 201(h) of the FD&C Act. This is the cleanest basis for
the non-device position and should be the stated primary rationale.

### 2.2 Fallback argument: the Clinical Decision Support (CDS) exemption

If any output is ever used to support a health-care professional (HCP) in a
prevention, diagnosis, or treatment decision, the platform must instead satisfy
the CDS criteria in section 520(o)(1)(E) of the FD&C Act (added by the 21st
Century Cures Act) and interpreted in FDA's final guidance *Clinical Decision
Support Software* (September 2022). Software is NON-device CDS only if it meets
ALL FOUR criteria:

| # | Criterion (paraphrased) | Platform self-assessment | Gap |
|---|---|---|---|
| 1 | Does NOT acquire, process, or analyze a medical image or a signal from an IVD or a signal-acquisition system | MET. The platform ingests curated database records, not images or device signals. | None |
| 2 | Displays, analyzes, or prints medical information (patient info or peer-reviewed studies, guidelines, etc.) | MET. It analyzes public pharmacology and knowledge-graph data and cites sources. | None |
| 3 | Provides recommendations to an HCP about prevention, diagnosis, or treatment (rather than a specific directive) | MET in form: it ranks candidate hypotheses; it does not issue a directive to treat a patient. | Keep wording non-directive in UI copy |
| 4 | Enables the HCP to INDEPENDENTLY REVIEW the basis, so the HCP does not rely PRIMARILY on the software | Largely MET by design: every candidate shows its mechanism (genes, pathways, network), evidence tier, provenance, and honest limitations. | See gap list: the basis must be co-presented with every score, sources must be identifiable and current, and no output may be framed as time-critical or autonomous. |

The strongest and weakest criterion is #4. The platform's explainability and
provenance layers are the direct evidence for it. The design rule that protects
criterion 4 is: **never present a score, tier, or "actionable" flag without the
reviewable basis alongside it.** This is a permanent product constraint in this
lane, not a one-time task.

## 3. Frameworks that still apply (non-device does not mean nothing applies)

Non-device positioning removes the SaMD / 510(k) / IEC 62304 obligations. It
does NOT remove data-integrity and computer-system-validation expectations when
outputs feed a regulated (GxP) process at a partner, or when records are relied
upon:

* **21 CFR Part 11 (electronic records and signatures)** applies if records are
  submitted to FDA or used to satisfy a predicate GxP rule. Requires audit
  trails, access control, e-signatures, and validation.
* **ALCOA+ (WHO TRS 996 Annex 5)** data-integrity principles. Substantially
  addressed at the data layer already (see VALIDATION_NARRATIVE.md).
* **GAMP 5 (2nd ed.) computer-system validation.** The validation suite plus the
  SOPs and CAPA logs are the executable OQ/PQ evidence; a requirements
  traceability matrix and a CI re-validation gate are the remaining pieces.
* **ISO 14971 risk management** is good practice even for a non-device research
  tool and is expected by pharma quality teams.

## 4. Gap analysis

Legend: DONE / PARTIAL / MISSING.

| Item | Framework | Status | Notes |
|---|---|---|---|
| Data-layer validation (bioactivity, concordance, resolution, KG integrity) | ALCOA+, GAMP 5 | DONE | validation/validate_*.py + results |
| Predictive validation vs external gold standard (repoDB) | GAMP 5 | DONE | AUROC 0.73, negative controls collapse to 0.5 |
| SOPs and CAPA logs | GAMP 5, MHRA | DONE | validation/SOP-*, CAPA-LOG-* |
| Provenance and freshness scoring | ALCOA+ Attributable | PARTIAL | services/provenance.py + build_stamp() attached to every lineage() record; extend the stamp to all API result payloads as needed |
| Intended Use Statement and non-device rationale | FD&C 201(h), CDS guidance | DONE (this document) | Needs counsel review before external use |
| Reproducibility stamping (data version + code commit into each result) | ALCOA+ Original, Part 11 | DONE | services/provenance.build_stamp() resolves the running code SHA + data snapshot versions; attached to lineage records |
| Audit trail (who ran what, when, on which data and code version) | Part 11 | PARTIAL | Implemented: services/audit_trail.py (hash-chained, tamper-evident, with size-based rotation that keeps the chain continuous across read-only segments) wired into every request when AUTH_ENABLED. Remaining: append-only external (WORM) storage — deployment-side |
| Authentication and access control | Part 11 | PARTIAL | Implemented: services/auth.py (hashed passwords, 3 roles, login/logout, gated by AUTH_ENABLED default off) + a declarative authority policy enforced in before_request (viewer reads, analyst writes, admin manages users) + an admin user API. Remaining: SSO/IdP federation — deployment-side |
| Electronic signatures | Part 11 | PARTIAL | Implemented: services/esign.py (re-auth second component, meaning + signer + timestamp, recorded into the tamper-evident audit chain) + /api/sign, /api/signatures. Remaining: bind signing to specific record types in the UI (dossier/screen sign-off) |
| Risk management file | ISO 14971 | PARTIAL | Scaffold done: validation/RISK_MANAGEMENT_FILE.md (11 hazards grounded in real findings, controls, residual risk). Remaining: named risk owner + review cadence (QMS) |
| Requirements traceability matrix (URS to spec to test) | GAMP 5 | DONE | validation/TRACEABILITY_MATRIX.md — 24 requirements mapped to implementation + evidence |
| CI re-validation gate (block on regression) | GAMP 5 re-provable | DONE | .github/workflows/ci.yml runs the 22-test suite on push/PR |
| Criterion-4 UI guardrail (basis always shown with score) | CDS guidance | PARTIAL | Explainability exists; enforce as an invariant |

## 5. Prioritized roadmap

**P0 (do first; lane-independent, high credibility, low risk):**

1. Reproducibility stamping. Write the code commit SHA and the exact source
   versions (ChEMBL 33, the knowledge-graph release, Open Targets snapshot date)
   into every result and provenance record, so any output can be regenerated and
   traced. Directly serves ALCOA+ Original and Part 11.
2. Requirements traceability matrix. One table mapping each user requirement to
   its functional spec and the validation test that proves it. Mostly assembling
   what already exists.
3. CI re-validation gate. Run the pytest suite (and, where data is available, the
   validation harnesses) on every change and block on regression.

**P1 (needed before any multi-user or partner deployment) — IMPLEMENTED, gated:**

4. Audit trail. DONE (services/audit_trail.py): append-only, hash-chained,
   tamper-evident event log stamped with actor, action, timestamp, and code SHA;
   wired into every request when AUTH_ENABLED; size-based rotation keeps the chain
   continuous across read-only segments. Remaining hardening: append-only
   external (WORM) storage — deployment-side, not application code.
5. Authentication and role-based access control. DONE (services/auth.py):
   hashed passwords, three roles, login/logout, secret-key hardening, gated by
   AUTH_ENABLED (default off). Role authority is enforced by a declarative policy
   in before_request (viewer reads, analyst state-changing verbs, admin /admin/*)
   plus an admin user-management API. Remaining: SSO/IdP federation — deployment-side.

**P2 (needed only if records are formally approved in the app, or a partner's
QMS requires it) — IMPLEMENTED:**

6. Electronic signatures (Part 11 subpart C). DONE (services/esign.py): a signing
   mechanism that re-authenticates (second component), records the signer, meaning,
   and timestamp into the tamper-evident audit chain, and exposes /api/sign +
   /api/signatures. Remaining: bind signing to specific record types in the UI.
7. Risk management file (ISO 14971). DONE as a scaffold
   (validation/RISK_MANAGEMENT_FILE.md): 11 hazards grounded in the real 2026-07-27
   findings, each with a control and residual risk. Remaining: a named risk owner
   and a review cadence, which are QMS (P3) activities.

**P3 (organizational, not code):**

8. Formal CDS / non-device determination reviewed by regulatory counsel.
9. External CSV audit if a pharma partner's GxP process consumes the output.

## 6. What this document does not claim

It does not claim the platform is FDA-approved, FDA-cleared, FDA-compliant,
GAMP-certified, or Part-11-validated. It does not constitute a regulatory
determination. The accurate external phrasing remains the one already used in
VALIDATION_NARRATIVE.md: the platform is *built to align with* the standards used
for GxP data integrity and computerized-system validation, and is positioned as
a non-device research and decision-support tool. Any stronger claim requires a
quality-management system, qualified regulatory counsel, and, where applicable,
an external audit.

## 7. References

* FD&C Act section 201(h) (device definition) and section 520(o)(1)(E) (CDS
  carve-out, added by the 21st Century Cures Act, 2016).
* FDA. *Clinical Decision Support Software — Guidance for Industry and FDA
  Staff.* September 2022.
* FDA. *Policy for Device Software Functions and Mobile Medical Applications.*
* FDA. *Data Integrity and Compliance With Drug CGMP.* 2018.
* WHO. *Guidance on Good Data and Record Management Practices.* TRS 996 Annex 5,
  2016 (ALCOA+).
* ISPE. *GAMP 5: A Risk-Based Approach to Compliant GxP Computerized Systems*,
  2nd ed., 2022.
* ISO 14971:2019. *Medical devices — Application of risk management to medical
  devices.*
