# External Validation Package (pharma-partner facing)

**Purpose.** The evidence a prospective pharma partner needs to decide whether RepurposeIQ is
trustworthy enough to (a) start a collaboration / due-diligence, and (b) — later, and only with
partner-run validation — inform prioritization. It states plainly what is validated, what is not,
when the platform must NOT be relied upon, and the two partner-run studies that would move it from
"credible decision-support" to "trusted for a development decision."

Honest posture (unchanged): RepurposeIQ is a **non-device, clinically-informed decision-support
layer**, not a standalone source of truth and not a clinical-outcome predictor. Every number below
is produced by a read-only, re-runnable script.

Last updated: 2026-07-28. Status: skeleton — the retrospective case studies (Section 3) are filled
after the current learned-composite run; the partner-run studies (Section 4) require partner data.

## 1. What IS validated (and how to re-run it)

- **Data layer (ALCOA+ / GAMP 5 framing).** Provenance + freshness on every source; 91% concordance
  with an independent authority (IUPHAR) within 1 log; measurement noise inside the ~0.5-log
  published floor; 100% disease-resolution on a curated golden set; 100% KG provenance coverage.
  Re-run: `validation/validate_*.py` (each writes a dated result + appends a CAPA entry).
- **Retrospective model recovery.** On repoDB **approved-vs-FAILED** (hard negatives — beating
  tried-and-failed repurposings, not random pairs), the leakage-free mechanistic score separates
  approved from failed at **AUROC ~0.73**; random baseline ~0.49; label-shuffle control collapses to
  ~0.49. Re-run: `validation/validate_predictions.py`.
- **Documented negative results** (what we tested and did NOT ship): KG-embedding ensemble (KGE-04),
  DWPC plausibility re-rank, signature-reversal ranking, and — critically — the PrimeKG zero-shot
  ranker's leakage-free collapse (AUROC 0.22) and its fix (disease-similarity transfer, certified
  recall@10 0.39). These are in `validation/*.json`; disclosing them is part of the trust case.

## 2. What is NOT validated (state this first, before any claim)

- **No prospective validation.** Every result is retrospective recovery on a fixed data snapshot.
- **Cannot predict clinical SUCCESS.** approved-vs-failed success labels are ~unlearnable (AUROC ~0.42);
  the platform ranks mechanistic PLAUSIBILITY, not probability of benefit. A repurposed drug still
  inherits ~8-14% Phase-I->approval odds.
- **Not validated on a partner's internal assets** — recovery on public repoDB does not guarantee
  performance on a specific proprietary portfolio.

## 3. When the platform must NOT be relied upon (hard "do-not-use" list)

- As the SOLE basis for a target-selection, go/no-go, or investment decision — it is a triage layer,
  not governance.
- For a **truly novel disease with no similar treated disease** (zero-shot) — the KG generator is
  recall-only there and near-chance on direction; treat hits as unverified hypotheses.
- Where **KG/data coverage is thin** (the platform emits a coverage-honesty flag) — a data gap is
  neither a confident positive nor negative.
- For **direction-sensitive cases** without confirmed biology — loss-of-function/congenital diseases,
  or any pair where inhibit-vs-activate hasn't been checked against the disease's GoF/LoF.
- For **biologics against intracellular targets** or **docking of antibodies** (category errors —
  suppressed, not scored).
- Any reading of a **cold single deltaG** without structure provenance (correct UniProt identity,
  domain, holo/apo, named pocket), or an apo/homology-only dock, as an affinity claim.
- As a **clinical / point-of-care** tool — it addresses researchers/BD, not patient care.

## 4. Partner-run studies that would earn development-decision trust (protocols)

These CANNOT be run without partner data; the protocols below make them turnkey.

**4.1 Blinded retrospective benchmark on internal assets.**
- Partner supplies N (>=30) of their own drug x indication pairs with KNOWN outcomes (approved / failed
  / deprioritized), labels withheld from us.
- We score all pairs blind; partner unblinds and computes AUROC/recall@k of our ranking vs the known
  outcomes, and vs the partner's own historical prioritization at the time.
- Success criterion (pre-registered): our ranking's AUROC on approved-vs-failed exceeds the partner's
  historical triage by a partner-set margin, disease- and compound-disjoint.

**4.2 Prospective silent-run pilot.**
- For one quarter, the platform scores the partner's live incoming repurposing questions in parallel
  ("silent") with their normal process; no decisions are changed.
- Compare our top-k hits against what the partner's team independently prioritized; track concordance,
  the novel hits we surfaced that they missed, and the deprioritizations we flagged.
- Pre-register what "useful" means (e.g. >=1 actioned novel hit per quarter, or measurable triage-time
  reduction) BEFORE the run.

**4.3 Reproducible public case studies (Section, filled next).**
- 10-20 repoDB cases: known approved repurposings the platform should RECOVER (rank high) + known
  trial FAILURES it should DEPRIORITIZE (rank low). Each with the actual score, rank, and the
  mechanism/evidence shown — demonstrating recovery of successes and rejection of failures on public
  ground truth, as a transparency artifact (NOT a substitute for 4.1/4.2).

## 5. The honest one-paragraph summary for a partner

"RepurposeIQ is a rigorously governed, evidence-traceable repurposing decision-support platform:
provenance and re-runnable validation on the data layer, an honest retrospective recovery benchmark
on hard negatives (AUROC ~0.73), and openly disclosed limits and failure modes. It is strong enough
to start a collaboration and support due-diligence, and it works as a decision-support layer over
your existing triage. It is NOT yet validated prospectively or on your internal assets, and it does
not predict clinical success — so it should augment, not replace, your target-selection governance.
The path to development-decision trust is a blinded internal benchmark and a prospective silent-run
pilot, both specified above."
