# Data Validation — Explained

**For managers, partners, and pharma customers of RepurposeIQ and CompoundIQ (POZ).**
Last updated: 2026-06-17

This document explains, in plain language, **what data our platforms run on, how we
checked that the data is trustworthy, and why a pharma partner can rely on it.**
You do not need a data-science background to read it. Every number here is something
we actually measured with a script that anyone can re-run; every reference is real.

---

## TL;DR (the 30-second version)

We build on the best public pharmaceutical databases (ChEMBL, a curated biomedical
knowledge graph, IUPHAR, MeSH). Public data has a reputation for being messy, so instead of *claiming* it's
clean, we **measured how good it is** and published the results inside the product.
The three numbers that matter:

- **91%** — how often our drug-potency data agrees with a *second, independent*
  expert database (IUPHAR), within normal experimental error.
- **Within the ~0.5-log "noise floor"** — our data is no noisier than the known,
  published limit of how reproducible lab measurements are in the first place.
- **100%** — correctness of our disease-name matching, checked against a
  hand-curated list of known-correct answers.

And we did all of this using the **same validation standards (ALCOA+, GAMP 5) that
regulated pharma companies use on their own systems** — which is the language a
pharma quality team already speaks.

---

## 1. The problem we are solving

Pharma scientists are trained to distrust two things: **black boxes** (models that
can't explain themselves) and **public data** (which they assume is scraped and
noisy). If we just said "we use ChEMBL, it's gold-standard," a serious reviewer
would shrug — because even good public databases contain measurement noise,
duplicates, and mapping errors.

So we reframed the problem. We treat public data not as "free stuff we downloaded,"
but as a **raw material that we put through engineering quality controls** — the
same kind a pharma company applies to its own data. Then we *show our work*.

---

## 2. Plain-language glossary (the terms you'll be asked about)

**The data sources**

- **ChEMBL** — a large, manually-curated database from the European Bioinformatics
  Institute (EMBL-EBI) that records how strongly drugs bind to their biological
  targets. It is the pharmaceutical industry's standard reference for this data.
- **Curated biomedical knowledge graph** — a giant map connecting drugs, genes,
  and diseases with millions of relationships (e.g. "this drug affects this gene,"
  "this gene is involved in this disease"). It lets us reason about *why* a drug
  might work for a new disease, through biology. (We build on an established,
  peer-reviewed open knowledge graph; the integration and reasoning layer is ours.)
- **IUPHAR / Guide to Pharmacology (GtoPdb)** — a *separate*, independently
  expert-curated pharmacology database. We use it as an outside referee to
  cross-check ChEMBL.
- **MeSH (Medical Subject Headings)** — the U.S. National Library of Medicine's
  official controlled vocabulary of medical/disease terms. It's how we make sure
  "Parkinson's," "PD," and "paralysis agitans" all map to the *one* correct disease.

**The measurement terms**

- **Ki, IC50, EC50, Kd** — different ways of measuring how potent a drug is at its
  target (smaller number = more potent). We keep them in separate "lanes" because
  they are not directly interchangeable.
- **Log unit** — potency is measured on a logarithmic scale, so **1 log unit = a
  10× difference**. When we say two measurements "agree within 1 log," we mean they
  are within an order of magnitude — the normal tolerance for this kind of lab data.
- **The ~0.5-log noise floor** — published research (Kramer 2012; see references)
  showed that even careful labs, measuring the *same* drug, disagree by about 0.5
  log on average. That's the natural noise floor of the field. Our job isn't to beat
  physics — it's to show our data sits *within* that known floor, not above it.

**The validation standards (this is what makes it "industry-level")**

- **ALCOA+** — the data-integrity standard used by the FDA, the UK's MHRA, and the
  WHO. It's an acronym for nine properties every trustworthy data record must have.
  In plain terms:
  - **A — Attributable:** you can trace every value back to its source.
  - **L — Legible:** it's readable and permanent.
  - **C — Contemporaneous:** recorded at the time, not back-filled.
  - **O — Original:** the raw record is kept, never overwritten.
  - **A — Accurate:** it's correct and error-checked.
  - **+ Complete:** nothing is silently dropped.
  - **Consistent:** it doesn't contradict itself.
  - **Enduring:** it survives over time.
  - **Available:** it can be retrieved when needed.
- **GAMP 5** — "Good Automated Manufacturing Practice," 5th edition, from ISPE. It's
  the pharma industry's playbook for **validating computer systems**. Its core idea
  is that validation must be a *repeatable test you can run again to re-prove the
  system* — not a one-time sign-off. Our validation scripts ARE that repeatable test.
  (In GAMP language these are the "OQ/PQ" tests — Operational and Performance
  Qualification — which just means "prove it works, and prove it keeps working.")
- **CAPA — Corrective And Preventive Action.** A formal log where every problem found
  is written down with its cause and its fix. In a validated system, *finding and
  recording errors is expected* — a clean-looking dashboard with no findings is
  actually a red flag.
- **The six data-quality dimensions** (from the DAMA data-management body of
  knowledge) — the standard checklist for "is this data good?": **completeness,
  accuracy, consistency, validity, uniqueness, timeliness.** We measured all six.

---

## 3. How we actually did the validation (the method)

We followed four principles a regulated pharma data team would recognize:

1. **Validate once, at the data layer.** Both products — RepurposeIQ and CompoundIQ
   — read the *same* ChEMBL database underneath. So we validate that shared source
   **one time** and certify it for both. This is also why the two products can't
   show contradictory raw facts about the same drug: there is one validated source,
   not two separate cleanings.

2. **Make it read-only and re-runnable.** Our validation is a set of scripts that
   only *read* the data (they can never change it) and can be re-run at any time —
   for example, every time we refresh the data. "Validated" therefore stays
   *provable on demand*, which is exactly what GAMP 5 asks for.

3. **Never delete — only flag (the ALCOA+ "Original" rule).** When we find a
   low-quality record, we don't erase it. We *tag* it and exclude it from the
   "certified" set, but the raw record stays. This means every decision is
   reversible and auditable — an auditor can always see what we excluded and why.

4. **Disclose everything, including weaknesses (CAPA).** Every issue we find is
   logged openly with a root cause and a corrective action. We publish these inside
   the product. We'd rather a customer see a measured, honest weakness than a
   suspiciously perfect dashboard.

---

## 4. The five validations we ran — what, why, and the result

**1) Source provenance & versioning** — *Can every dataset be traced and dated?*
We stamped every dataset with its exact source and version, so nothing is anonymous
or undated.
→ ChEMBL 33 (released mid-2023): 2.4M molecules, 20.3M potency measurements;
the curated biomedical knowledge graph: 2.25M relationships across 24 types.

**2) Bioactivity (potency) validation** — *Is the drug-potency data complete and how
noisy is it?* We checked field-completeness, applied ChEMBL's own quality flags
(confidence score, "suspect value" flags, censored values), kept each measurement
type (Ki/IC50/EC50/Kd) separate, and **measured the noise** by comparing repeated
measurements of the same drug-target pair.
→ Across 4.1M measurements, our noise is a median of **0.02–0.19 log** and a mean of
**≤0.42 log** — i.e. **inside the published ~0.5-log floor.** We are not artificially
smoothing the data; it is genuinely as reproducible as the field allows.

**3) Cross-database concordance** — *Does an independent expert database agree with
ChEMBL?* We took every drug-target potency value that appears in *both* ChEMBL and
the independently-curated IUPHAR database and compared them.
→ **91.0% of 5,754 shared values agree within 1 log** (median difference just 0.03
log). This is the strongest evidence of all, because it's ChEMBL checked against an
*outside authority*, not against itself.

**4) Disease entity-resolution audit** — *When a user types a disease, do we map it
to the right concept?* This is the highest-risk silent error: if "RA" mapped to the
wrong disease, every downstream result would be wrong while *looking* fine. We tested
this two different ways, and we keep them honestly separate:
- **Consistency** (round-trip identity: name → ID → name): proves the matching is
  *stable and deterministic*. Result: heading **100%**, synonym **99.7%**.
- **Correctness** (against a hand-curated "golden set" of 34 tricky synonyms with
  *known-correct* answers, e.g. "Lou Gehrig's disease" → Amyotrophic Lateral
  Sclerosis): proves the matching is *actually clinically right*, not just
  consistent. Result: **100% (34/34).**
We report both because consistency alone could be "consistently wrong" — correctness
is what proves it's right.

**5) Knowledge-graph integrity** — *Is our biological reasoning based on real data?*
We confirmed the platform uses the genuine, full curated knowledge graph (2.25M relationships),
tagged by source, so the "drug → gene → disease" reasoning is real and independent of
ChEMBL, and every conclusion can show the specific genes that connect a drug to a
disease.

---

## 5. The results at a glance

| Validation | The question it answers | Result |
|---|---|---|
| Provenance & versioning | Is everything traceable and dated? | ChEMBL 33 (mid-2023, 2.4M molecules, 20.3M activities); curated knowledge graph (2.25M edges) |
| Bioactivity / noise | How noisy is the potency data? | Median 0.02–0.19 log, mean ≤0.42 log — **within the ~0.5-log published floor** |
| Cross-database concordance | Does an independent database agree? | **91.0%** of 5,754 shared pairs agree within 1 log (median diff 0.03 log) |
| Entity resolution | Do disease names map correctly? | Consistency 99.7–100%; **correctness 100%** on a curated golden set |
| Knowledge-graph integrity | Is the biology reasoning real? | 2.25M edges, **100% provenance-tagged, 100% referential integrity**, 0 duplicates |
| Mechanism direction | Is inhibit-vs-activate right? | **99.3%** direction agreement with IUPHAR (1,410 shared drug-target pairs) |
| **Predictive recovery** | **Does the engine find real repurposing?** | **AUROC 0.73** (mechanism-only, leakage-free) vs repoDB; random 0.49; 68% drug coverage |
| **Score calibration** | **Is a score a probability?** | Raw ECE 0.38 → **0.013** after honest out-of-fold isotonic recalibration |
| **KG-embedding** | **Does graph structure add signal?** | In-graph ensemble looked promising (+0.042), but on the external repoDB benchmark KG alone ≈ chance (0.54) and ensembling *hurt* (−0.053) → **tested and NOT integrated** (would regress) |
| **Metapath/DWPC** | **Does the Rephetio method help?** | In-graph 0.6–0.7, but external repoDB CV-model 0.59 vs mechanistic 0.75 on same pairs → **tested and NOT integrated**; Hetionet is a weaker source than live OT+ChEMBL |

---

## 6. Why a pharma partner can trust this — by platform

Both platforms stand on the **same validated ChEMBL data layer.** What differs is the
*reasoning* built on top — so the trust argument is shared at the bottom and specific
at the top:

- **RepurposeIQ (drug repurposing).** A repurposing suggestion rests on (a) potency
  data that is validated and 91%-concordant with an independent database; (b) biology
  from the *real* curated knowledge graph, where every suggestion shows the actual genes
  linking the drug to the disease (explainable, not a black box); and (c) a
  provenance tag on every candidate. The trust comes from being **auditable and
  externally corroborated**, not from a confidence score alone.

- **CompoundIQ / POZ (505(b)(2) reformulation screening).** The drug-exposure (Cmax)
  data is validated under its own protocol (SOP-001) with a full accounting of every
  drug, and the underlying potency data shares the same ChEMBL validation. So a
  screening candidate's pharmacokinetic and potency inputs are traceable,
  version-stamped, and quality-flagged.

The two platforms **share the data layer but use distinct algorithms**: RepurposeIQ
reasons over the *knowledge graph*; CompoundIQ reasons over *pharmacokinetic* (Cmax)
scaling. Same validated facts underneath; different, purpose-built logic on top.

---

## 7. Validating the predictions (not just the data)

Everything above validates the **data**. We also now validate the **results** —
whether the engine actually recovers real repurposing — against an external gold
standard, **repoDB** (a curated database of drug-indication pairs that were either
**approved** or **failed in trials**).

The honest, hard version of the test:

- **We test biology only.** The live engine also uses clinical-trial and literature
  signals; in a *backward-looking* test those would be the answer (data leakage). So
  for the benchmark we **switch them off** and score on mechanism alone — target
  overlap, pathway, and protein-network proximity.
- **We use hard negatives.** The negatives aren't random drug-disease pairs; they're
  repurposings that were actually *tried and failed*. Beating those is real signal.

Result: the mechanism-only score separates approved from failed repurposings with
**AUROC 0.73** (0.73 averaged within disease), versus 0.49 for random and 0.43 for a
"popular drug" baseline. A label-shuffle control collapses to 0.48 (chance), proving
the signal is real. The raw score is a good *ranking* but is not itself a probability
(calibration error 0.38); an honest out-of-fold recalibration turns it into one
(error 0.013).

We also *measure our own improvements*: a coverage fix (mapping drug salts to their
parent molecule and filling gaps from high-confidence bioactivity) raised the share
of drugs the engine can score from 46% to 68% and AUROC from 0.71 to 0.73 — and we
kept the version that helped (gap-fill only) over the one that didn't (broaden every
drug), because the benchmark said so.

We also disclose what *doesn't* help: target overlap carries most of the signal, and
the pathway term is near-neutral for ranking (we keep it for correctness and
explainability, not because it boosts recovery).

**The remaining honest limit.** A retrospective benchmark proves *historical
recovery on a fixed data snapshot* — it is not a guarantee of prospective
performance. We track that as the data refreshes. Being explicit about this boundary
is what makes the number believable.

---

## 8. How to explain it in 60 seconds

**To a pharma scientist / data-quality team:**
> "We treat public data as a raw material under engineering controls. We validated
> the shared ChEMBL layer under ALCOA+ and GAMP 5 — read-only, re-runnable, additive
> (we flag, never delete), with an open CAPA log. The headline: 91% concordance with
> IUPHAR as an independent referee, measured noise inside the published 0.5-log floor,
> and 100% correctness on a curated disease-mapping golden set. And it's all published
> live in the product, not just in a slide. We also validate the *results*: on the
> repoDB gold standard, our leakage-free mechanistic score recovers approved-vs-failed
> repurposings at AUROC 0.73 — with the clinical/literature features switched off so
> it can't cheat. Prospective performance is the stated boundary we still track."

**To your manager (business framing):**
> "Customers don't trust 'we use public data.' So we built a pharma-grade validation
> layer — using the exact standards (ALCOA+, GAMP 5) their own QA teams use — and we
> show the proof inside the product. The standout number is 91% agreement with an
> independent expert database. This turns 'trust us' into 'here's the measured
> evidence,' which is what unblocks the technical due-diligence stage of a deal."

---

## 9. References (all real and citable)

**Standards we built to**
- FDA. *Data Integrity and Compliance With Drug CGMP: Questions and Answers — Guidance for Industry.* 2018.
- MHRA. *'GXP' Data Integrity Guidance and Definitions.* 2018.
- WHO. *Guidance on Good Data and Record Management Practices.* WHO Technical Report Series No. 996, Annex 5, 2016. (Defines **ALCOA+**.)
- EU GMP **Annex 11**: Computerised Systems (EudraLex Vol. 4).
- ISPE. *GAMP 5: A Risk-Based Approach to Compliant GxP Computerized Systems*, 2nd ed., 2022.
- DAMA International. *DAMA-DMBOK: Data Management Body of Knowledge*, 2nd ed., 2017. (The six data-quality dimensions.)

**Data sources**
- Zdrazil B, et al. "The ChEMBL Database in 2023." *Nucleic Acids Research* 2024, 52(D1):D1180–D1192. (ChEMBL, EMBL-EBI. Note: the paper sits in the 2024 NAR database issue; the ChEMBL 33 data release itself dates to mid-2023 — the citation year is the journal volume, not a newer data version.)
- Harding SD, et al. "The IUPHAR/BPS Guide to PHARMACOLOGY in 2024." *Nucleic Acids Research* 2024, 52(D1):D1438–D1449.
- U.S. National Library of Medicine. *Medical Subject Headings (MeSH).*

**The experimental-noise benchmark (why ~0.5 log is the floor)**
- Kramer C, Kalliokoski T, Gedeck P, Vulpetti A. "The Experimental Uncertainty of Heterogeneous Public Ki Data." *Journal of Medicinal Chemistry* 2012, 55(11):5165–5173.
- Kalliokoski T, Kramer C, Vulpetti A, Gedeck P. "Comparability of Mixed IC50 Data — A Statistical Analysis." *PLoS ONE* 2013, 8(4):e61007.

---

## 10. How to reproduce (read-only)

```bash
# Data validation
python validation/validate_bioactivity.py     # completeness, validity, noise
python validation/validate_concordance.py      # accuracy vs IUPHAR/GtoPdb (91%)
python validation/validate_resolution.py       # disease-term resolution (consistency + golden set)
python validation/validate_mechanisms.py       # mechanism direction vs IUPHAR (99.3%)
python validation/validate_kg_accuracy.py      # knowledge-graph integrity & provenance
python validation/validate_data.py             # Cmax universe (CompoundIQ / POZ, SOP-001)

# Predictive (results) validation — needs data/external/repodb/full.csv
python validation/validate_predictions.py      # retrospective recovery vs repoDB (AUROC)
python validation/validate_calibration.py      # score → probability calibration
```

Each script only reads the database, writes a machine-readable results file plus a
run log, and can be re-run on every data refresh (append a dated entry to the CAPA
log each time). The live versions of these results are shown on the platform's
**Data Validation** page.

---

## Appendix A — The six standards, explained

Five of these are pharma regulatory / data-integrity standards; one (DAMA-DMBOK)
is a general data-management standard. They divide the work: some define *what*
trustworthy data is, one requires you to validate the *system*, one tells you
*how*, and one provides the *measurement framework*.

**1. FDA — Data Integrity and Compliance With Drug CGMP (2018).**
A U.S. Food and Drug Administration guidance. "CGMP" = Current Good Manufacturing
Practice. This Q&A spells out the FDA's data-integrity expectations and is what FDA
inspectors reference when judging whether a company's data can be trusted. *How we
use it:* its core demands — every value traceable to source (Attributable) and the
raw record preserved (Original) — are why our validation is read-only and additive:
we flag bad records, never delete, and always keep the raw data.

**2. MHRA — 'GXP' Data Integrity Guidance and Definitions (2018).**
The UK equivalent, from the Medicines and Healthcare products Regulatory Agency.
"GXP" is the umbrella for all Good Practice standards (GMP, GLP, GCP, …). It is
widely regarded as the clearest, most practical statement of data integrity and is
cited internationally. *How we use it:* its emphasis on the full data lifecycle and
a "measure-and-disclose" posture is what our CAPA discipline maps to — we record
every finding, weaknesses included, rather than hiding them.

**3. WHO TRS 996, Annex 5 (2016) — defines ALCOA+.**
World Health Organization, Technical Report Series No. 996, Annex 5: "Guidance on
good data and record management practices." The global, region-neutral document
that formally lays out the ALCOA+ principles. *How we use it:* this is the source
of the ALCOA+ framing for the whole validation — every check maps to a principle
(source-tagging → Attributable, raw retention → Original, concordance → Accurate,
field coverage → Complete).

**4. EU GMP Annex 11 — Computerised Systems.**
Part of EudraLex Volume 4 (the EU's GMP rulebook). Annex 11 specifically governs
computer systems in regulated pharma — validation, audit trails, data security,
electronic records. *How we use it:* our pipeline and platform are computerised
systems, so Annex 11 is why we version-stamp the data (ChEMBL 33), keep an audit
trail (CAPA log + per-run logs), and treat validation as a controlled, repeatable
activity.

**5. ISPE — GAMP 5, 2nd ed. (2022).**
From the International Society for Pharmaceutical Engineering. "GAMP" = Good
Automated Manufacturing Practice. Where Annex 11 says "you must validate computer
systems," GAMP 5 is the practical how-to: a risk-based methodology with the
qualification lifecycle (IQ/OQ/PQ — Installation / Operational / Performance
Qualification). *How we use it:* our validation scripts ARE the executable OQ/PQ
test — re-runnable proof that data quality holds, re-proven on every refresh.

**6. DAMA-DMBOK, 2nd ed. (2017) — data-quality dimensions.**
From DAMA International (Data Management Association); "DMBOK" = Data Management Body
of Knowledge. The standard, vendor-neutral reference for data management as a
discipline; it defines the formal data-quality dimensions. *How we use it:* our
scorecard is organized around its six dimensions — completeness, accuracy,
consistency, validity, uniqueness, timeliness — and each validation module maps to
one or more.

**How they fit together**
- WHO + MHRA + FDA → define *what* trustworthy data is (ALCOA+).
- EU Annex 11 → requires you to *validate the computer system*.
- GAMP 5 → tells you *how* (risk-based, re-runnable OQ/PQ).
- DAMA-DMBOK → gives the *measurement framework* (the six dimensions).

**Honesty note.** We built our methodology to *align with* these standards; we are
not formally certified or audited against them (that would require a full
quality-management system and an external auditor). The accurate phrasing for a
customer is "built to the same standards used for GxP data integrity and
computerized-system validation" — never "FDA-compliant" or "GAMP-certified."
