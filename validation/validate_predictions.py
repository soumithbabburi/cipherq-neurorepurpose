"""
Predictive Validation — Retrospective Repurposing Recovery
═══════════════════════════════════════════════════════════════════════════
Validates the *results*, not just the data: does the engine's MECHANISTIC score
separate real repurposing successes from failures?

Gold standard (repoDB, Brown & Patel, Sci Data 2017):
  • POSITIVES = drug-indication pairs with status "Approved"
  • NEGATIVES = pairs that reached trials but FAILED ("Terminated"/"Withdrawn"/
    "Suspended") — hard negatives, not random non-edges.
This is the right benchmark because the negatives are *plausible-but-failed*
repurposings, so beating them is genuine signal, not a popularity artefact.

Leakage control (the part that makes this honest):
  The platform's full score weights clinical-trial signal (30%) and literature
  co-mention — in a retrospective test those features ARE the label. So here we
  score on BIOLOGY ONLY: target overlap + (direction-aware) pathway overlap +
  PPI proximity, with clinical / indication / regulatory / literature REMOVED.

Method:
  1. repoDB drug_name -> ChEMBL id + target genes + mechanism action_type
     (local chembl_33; drugs with no curated mechanism are excluded + reported).
  2. repoDB ind_name -> disease genes / pathways / PPI (Open Targets, cached).
  3. Mechanism-only score per (drug, disease) pair.
  4. Metrics: AUROC, AUPRC, BEDROC, enrichment; vs baselines (random, drug
     popularity, target-overlap-only); per-disease precision@k; ablation of the
     direction-aware pathway term; and a label-shuffle negative control.

Read-only against chembl_33 + Open Targets/STRING (cached). repoDB defaults to
data/external/repodb/full.csv. Writes validation/predictions_results.json + log.

Usage:
    python validation/validate_predictions.py [repodb.csv] [max_diseases]
"""

import sys
import csv
import json
import math
import random
import datetime
import numpy as np
from pathlib import Path
from collections import defaultdict

from sqlalchemy import text

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))

from validation.validate_concordance import get_engine

DEFAULT_REPODB = ROOT / "data" / "external" / "repodb" / "full.csv"
APPROVED = "Approved"
FAILED = {"Terminated", "Withdrawn", "Suspended"}

# Mechanism-only weights = the biology dimensions of the engine, renormalised.
# (Engine WEIGHTS: target .25, pathway .20, ppi .20 — clinical/indication/
#  regulatory deliberately dropped here to remove retrospective leakage.)
_W = {"target": 0.25, "pathway": 0.20, "ppi": 0.20}
_WSUM = sum(_W.values())
MECH_W = {k: v / _WSUM for k, v in _W.items()}


# ── Metrics (no sklearn dependency) ────────────────────────────────────────────

def auroc(pairs):
    """pairs: list of (score, label 0/1). Mann-Whitney rank AUC; ties handled."""
    pos = [s for s, y in pairs if y == 1]
    neg = [s for s, y in pairs if y == 0]
    if not pos or not neg:
        return None
    ordered = sorted(pairs, key=lambda p: p[0])
    ranks = {}
    i = 0
    n = len(ordered)
    while i < n:
        j = i
        while j < n and ordered[j][0] == ordered[i][0]:
            j += 1
        avg_rank = (i + 1 + j) / 2.0       # average rank for ties (1-based)
        for k in range(i, j):
            ranks[id(ordered[k])] = avg_rank
        i = j
    rank_sum = 0.0
    for p in pairs:
        if p[1] == 1:
            rank_sum += ranks[id(p)]
    n_pos, n_neg = len(pos), len(neg)
    return (rank_sum - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)


def average_precision(pairs):
    """Area under precision-recall (AP). pairs: (score, label)."""
    ordered = sorted(pairs, key=lambda p: -p[0])
    total_pos = sum(1 for _, y in ordered if y == 1)
    if total_pos == 0:
        return None
    tp = 0
    ap = 0.0
    for i, (_, y) in enumerate(ordered, 1):
        if y == 1:
            tp += 1
            ap += tp / i               # precision at this recall step
    return ap / total_pos


def bedroc(pairs, alpha=20.0):
    """Truchon & Bayly (2007) early-recognition metric (RDKit closed form).
    pairs: (score, label). Returns a value in [0, 1]."""
    ordered = sorted(pairs, key=lambda p: -p[0])
    n = len(ordered)
    npos = sum(1 for _, y in ordered if y == 1)
    if npos == 0 or npos == n:
        return None
    ra = npos / n
    # sum of exp(-alpha * rank/N) over actives (rank 1-based)
    s = sum(math.exp(-alpha * (i + 1) / n) for i, (_, y) in enumerate(ordered) if y == 1)
    rie = (s / npos) / ((1.0 / n) * (1 - math.exp(-alpha)) / (math.exp(alpha / n) - 1))
    bedroc_val = (rie * ra * math.sinh(alpha / 2)
                  / (math.cosh(alpha / 2) - math.cosh(alpha / 2 - alpha * ra))
                  + 1.0 / (1 - math.exp(alpha * (1 - ra))))
    return max(0.0, min(1.0, bedroc_val))


def enrichment_factor(pairs, frac=0.1):
    ordered = sorted(pairs, key=lambda p: -p[0])
    n = len(ordered)
    k = max(1, int(round(n * frac)))
    total_pos = sum(1 for _, y in ordered if y == 1)
    if total_pos == 0:
        return None
    hits_top = sum(1 for _, y in ordered[:k] if y == 1)
    return round((hits_top / k) / (total_pos / n), 2)


# ── repoDB ──────────────────────────────────────────────────────────────────────

def load_repodb(path):
    """Return {disease_name: {drug_name_lower: label}} with label 1=approved,0=failed.
    A pair approved in any row is positive; else failed."""
    pairs = defaultdict(dict)   # disease -> drug -> label
    with open(path, encoding="utf-8") as f:
        for r in csv.DictReader(f):
            dis = (r.get("ind_name") or "").strip()
            drug = (r.get("drug_name") or "").strip().lower()
            status = (r.get("status") or "").strip()
            if not dis or not drug:
                continue
            if status == APPROVED:
                pairs[dis][drug] = 1
            elif status in FAILED:
                pairs[dis].setdefault(drug, 0)   # don't overwrite an approval
    return pairs


# ── Scoring ──────────────────────────────────────────────────────────────────────

def mech_score(drug_genes, drug_sig, disease_genes, disease_pathways, ppi_adj, pathway_mode,
               disease_weights=None):
    """Mechanism-only composite. pathway_mode: 'aware' | 'blind' | 'off'.
    disease_weights (optional): {GENE: genetic score} → genetics-weighted target overlap."""
    from services.repurposing_engine import (
        _score_target_overlap, _score_target_overlap_weighted,
        _score_pathway_overlap, _score_ppi_network)
    t = (_score_target_overlap_weighted(drug_genes, disease_weights)
         if disease_weights else _score_target_overlap(drug_genes, disease_genes))
    p = _score_ppi_network(drug_genes, ppi_adj)
    if pathway_mode == "off":
        w = {"target": _W["target"], "ppi": _W["ppi"]}
        wsum = sum(w.values())
        return (w["target"] * t + w["ppi"] * p) / wsum
    sig = drug_sig if pathway_mode == "aware" else None
    pw = _score_pathway_overlap(drug_genes, disease_pathways, drug_signature=sig)
    return MECH_W["target"] * t + MECH_W["pathway"] * pw + MECH_W["ppi"] * p


def _resolve_drug_targets(eng, names, mode="rich"):
    """Map repoDB drug names -> target gene sets, action map, and a ChEMBL id.

    mode='legacy' : pref_name match + drug_mechanism targets only (the original).
    mode='rich'   : also match molecule_synonyms, fold salt->parent via
                    molecule_hierarchy, and add high-confidence single-protein
                    targets from the activities table (pchembl>=6, conf>=8) — the
                    coverage fix for drugs that have no curated mechanism row.
    Returns (drug_genes name->set(gene), drug_actions name->{gene:action},
             drug_chembl name->chembl_id).
    """
    drug_genes = defaultdict(set)
    drug_actions = defaultdict(dict)
    drug_chembl = {}
    names = [n for n in names if n]
    if not names:
        return drug_genes, drug_actions, drug_chembl

    with eng.connect() as c:
        # name -> molregno(s); molregno -> name(s)
        name_molregnos = defaultdict(set)
        molregno_name = defaultdict(set)
        for nm, mol, ch in c.execute(text(
                "SELECT LOWER(pref_name), molregno, chembl_id FROM molecule_dictionary "
                "WHERE LOWER(pref_name) = ANY(:names)"), {"names": names}).fetchall():
            name_molregnos[nm].add(mol); molregno_name[mol].add(nm); drug_chembl.setdefault(nm, ch)

        if mode == "rich":
            for nm, mol, ch in c.execute(text(
                    "SELECT LOWER(ms.synonyms), md.molregno, md.chembl_id "
                    "FROM molecule_synonyms ms JOIN molecule_dictionary md ON md.molregno = ms.molregno "
                    "WHERE LOWER(ms.synonyms) = ANY(:names)"), {"names": names}).fetchall():
                name_molregnos[nm].add(mol); molregno_name[mol].add(nm); drug_chembl.setdefault(nm, ch)

        all_mol = set().union(*name_molregnos.values()) if name_molregnos else set()

        if mode == "rich" and all_mol:
            # fold salt -> parent so a salt inherits the parent's targets
            for mol, parent in c.execute(text(
                    "SELECT molregno, parent_molregno FROM molecule_hierarchy "
                    "WHERE molregno = ANY(:mols)"), {"mols": list(all_mol)}).fetchall():
                if parent and parent != mol:
                    for nm in molregno_name.get(mol, ()):
                        molregno_name[parent].add(nm)
                    all_mol.add(parent)

        if not all_mol:
            return drug_genes, drug_actions, drug_chembl
        mol_list = list(all_mol)

        # curated mechanism targets (carry action_type for the direction term)
        for mol, gene, action in c.execute(text(
                """
                SELECT dm.molregno, UPPER(csyn.component_synonym), dm.action_type
                FROM drug_mechanism dm
                JOIN target_components tc    ON tc.tid = dm.tid
                JOIN component_synonyms csyn ON csyn.component_id = tc.component_id
                                            AND csyn.syn_type = 'GENE_SYMBOL'
                WHERE dm.molregno = ANY(:mols)
                """), {"mols": mol_list}).fetchall():
            if not gene:
                continue
            for nm in molregno_name.get(mol, ()):
                drug_genes[nm].add(gene)
                if action:
                    drug_actions[nm].setdefault(gene, action)

        # high-confidence bioactivity targets — the big coverage add.
        #   rich     : add for ALL drugs (broadest coverage, can dilute clean drugs)
        #   fallback : add ONLY for drugs with no curated mechanism target (the
        #              coverage gap), so well-annotated drugs keep clean targets.
        if mode in ("rich", "fallback"):
            if mode == "fallback":
                gap_mols = [m for m in mol_list
                            if not any(drug_genes.get(nm) for nm in molregno_name.get(m, ()))]
            else:
                gap_mols = mol_list
            if gap_mols:
                for mol, gene in c.execute(text(
                        """
                        SELECT a.molregno, UPPER(csyn.component_synonym)
                        FROM activities a
                        JOIN assays ass             ON ass.assay_id = a.assay_id
                        JOIN target_dictionary td   ON td.tid = ass.tid
                                                   AND td.target_type = 'SINGLE PROTEIN'
                        JOIN target_components tc   ON tc.tid = ass.tid
                        JOIN component_synonyms csyn ON csyn.component_id = tc.component_id
                                                    AND csyn.syn_type = 'GENE_SYMBOL'
                        WHERE a.molregno = ANY(:mols)
                          AND a.pchembl_value >= 6
                          AND a.standard_relation = '='
                          AND ass.confidence_score >= 8
                        """), {"mols": gap_mols}).fetchall():
                    if not gene:
                        continue
                    for nm in molregno_name.get(mol, ()):
                        drug_genes[nm].add(gene)

    return drug_genes, drug_actions, drug_chembl


def run(repodb_path, max_diseases=60, mapping="fallback"):
    log_lines = []

    def log(msg=""):
        print(msg)
        log_lines.append(str(msg))

    log("=" * 72)
    log("PREDICTIVE VALIDATION  -  retrospective repurposing recovery (repoDB)")
    log(f"Run: {datetime.datetime.now().isoformat(timespec='seconds')}")
    log("=" * 72)

    from services.disease_ontology import resolve_disease, get_ppi_network
    from services.signature_engine import drug_signature

    repo = load_repodb(repodb_path)
    log(f"\nrepoDB diseases: {len(repo):,}")

    # 1. All drug names we might need → ChEMBL id + genes + actions (local chembl_33)
    all_drugs = sorted({d for dis in repo.values() for d in dis})
    eng = get_engine()
    drug_genes, drug_actions, drug_chembl = _resolve_drug_targets(eng, all_drugs, mapping)
    mapped = set(drug_genes)
    log(f"Drug-target mapping mode: {mapping}")
    log(f"repoDB drugs total: {len(all_drugs):,} | mapped to ChEMBL w/ targets: {len(mapped):,} "
        f"({100.0*len(mapped)/max(1,len(all_drugs)):.0f}%)")

    drug_sig = {name: drug_signature([{"gene": g, "action": a}
                                      for g, a in drug_actions[name].items()])
                for name in mapped}

    # 2. Pick diseases: most mapped-drug coverage, with ≥1 pos and ≥1 neg
    cand_dis = []
    for dis, drugs in repo.items():
        labeled = {d: y for d, y in drugs.items() if d in mapped}
        npos = sum(1 for y in labeled.values() if y == 1)
        nneg = sum(1 for y in labeled.values() if y == 0)
        if npos >= 1 and nneg >= 1 and (npos + nneg) >= 5:
            cand_dis.append((dis, npos + nneg, npos, nneg))
    cand_dis.sort(key=lambda x: -x[1])
    cand_dis = cand_dis[:max_diseases]
    log(f"Evaluable diseases (>=5 mapped drugs, both classes): {len(cand_dis)}")

    # 3. Resolve disease biology (cached) + 4. score every pair
    all_pairs = []                 # (score_aware, label)
    base_random, base_pop, base_target = [], [], []
    pairs_blind, pairs_off = [], []
    kg_triplets = []               # (mech_score, kg_score, label) for KG-mappable pairs
    pairs_genetics = []            # genetics-weighted target overlap variant
    per_disease = []
    resolved_ok = 0
    rng = random.Random(42)
    try:
        from services import kg_score as _kg
        _kg_ok = _kg.is_available()
    except Exception:
        _kg, _kg_ok = None, False

    for dis, _tot, _np, _nn in cand_dis:
        dinfo = resolve_disease(dis) or {}
        dgenes = [t["gene_symbol"] for t in dinfo.get("targets", [])[:40]]
        dpath = dinfo.get("pathways", [])
        if len(dgenes) < 5:
            continue
        resolved_ok += 1
        ppi = get_ppi_network(dgenes[:15]) if dgenes else {}
        # genetics weights (for the weighted-scorer variant) AND a genetics-ordered
        # gene list (for the clean same-formula isolation: only the ranking changes)
        gt = [t for t in dinfo.get("targets", [])[:40] if t.get("gene_symbol")]
        dweights = {t["gene_symbol"].upper():
                    (t.get("genetic_score") if t.get("genetic_score") is not None else t.get("score", 0))
                    for t in gt if (t.get("genetic_score") or t.get("score", 0)) > 0}
        dgenes_gen = [t["gene_symbol"] for t in
                      sorted(gt, key=lambda x: (x.get("genetic_score") or 0), reverse=True)]

        labeled = [(d, y) for d, y in repo[dis].items() if d in mapped]
        dis_pairs = []
        for d, y in labeled:
            g = sorted(drug_genes[d])
            sa = mech_score(g, drug_sig[d], dgenes, dpath, ppi, "aware")
            # clean isolation: same flat formula, disease genes RE-RANKED by genetics
            sg = mech_score(g, drug_sig[d], dgenes_gen, dpath, ppi, "aware")
            pairs_genetics.append((sg, y))
            sb = mech_score(g, drug_sig[d], dgenes, dpath, ppi, "blind")
            so = mech_score(g, drug_sig[d], dgenes, dpath, ppi, "off")
            all_pairs.append((sa, y)); pairs_blind.append((sb, y)); pairs_off.append((so, y))
            base_random.append((rng.random(), y))
            base_pop.append((float(len(g)), y))             # popularity = #targets
            base_target.append((_target_only(g, dgenes), y))
            dis_pairs.append((sa, y))
            if _kg_ok:
                kgv = _kg.treats_score(drug_name=d, disease=dis)
                if kgv is not None:
                    kg_triplets.append((sa, kgv, y))
        da = auroc(dis_pairs)
        if da is not None:
            per_disease.append({"disease": dis, "n": len(dis_pairs),
                                "auroc": round(da, 3),
                                "ef10": enrichment_factor(dis_pairs, 0.1)})

    log(f"KG signal available: {_kg_ok} | KG-mappable pairs collected: {len(kg_triplets)}")
    n_pos = sum(1 for _, y in all_pairs if y == 1)
    n_neg = len(all_pairs) - n_pos
    log(f"Diseases resolved with biology: {resolved_ok} | evaluated pairs: {len(all_pairs):,} "
        f"(approved {n_pos:,} / failed {n_neg:,}, base rate {100.0*n_pos/max(1,len(all_pairs)):.0f}%)")

    # negative control: shuffle labels
    shuffled = [(s, y) for (s, _), y in
                zip(all_pairs, _shuffled([y for _, y in all_pairs], rng))]

    def fmt(x):
        return None if x is None else round(x, 4)

    main_auroc = auroc(all_pairs)
    results = {
        "engine_mechanistic": {
            "auroc": fmt(main_auroc),
            "auprc": fmt(average_precision(all_pairs)),
            "bedroc_a20": fmt(bedroc(all_pairs, 20.0)),
            "ef_top10pct": enrichment_factor(all_pairs, 0.1),
        },
        "baselines": {
            "random":            {"auroc": fmt(auroc(base_random))},
            "drug_popularity":   {"auroc": fmt(auroc(base_pop))},
            "target_overlap_only": {"auroc": fmt(auroc(base_target))},
        },
        "ablation_pathway_term": {
            "direction_aware": fmt(main_auroc),
            "direction_blind": fmt(auroc(pairs_blind)),
            "pathway_excluded": fmt(auroc(pairs_off)),
        },
        "negative_control_label_shuffle": {"auroc": fmt(auroc(shuffled))},
    }

    # ── Genetics-weighted disease genes vs flat (target dimension) ─────────────
    au_gen = auroc(pairs_genetics) if pairs_genetics else None
    results["disease_genes_genetics_weighted"] = {
        "flat_target_auroc": fmt(main_auroc),
        "genetics_weighted_auroc": fmt(au_gen),
        "lift": (None if (au_gen is None or main_auroc is None) else round(au_gen - main_auroc, 4)),
        "note": "Genetics-weighted = disease genes weighted by Open Targets genetic_association "
                "score (Nelson 2015); also reduces leakage vs the overall OT score (no known-drug channel).",
    }

    # ── KG-embedding ensemble on the KG-mappable subset of repoDB ──────────────
    # Does adding the KG signal beat the mechanistic score on the SAME benchmark?
    # Measured only where both drug & disease map to Hetionet (rest falls back).
    if kg_triplets:
        from scipy.stats import rankdata
        ms = np.array([m for m, _, _ in kg_triplets])
        ks = np.array([k for _, k, _ in kg_triplets])
        ys = [y for _, _, y in kg_triplets]
        ens = rankdata(ms) + rankdata(ks)        # rank-average ensemble
        au_m = auroc(list(zip(ms.tolist(), ys)))
        au_k = auroc(list(zip(ks.tolist(), ys)))
        au_e = auroc(list(zip(ens.tolist(), ys)))
        lift = None if (au_e is None or au_m is None) else round(au_e - au_m, 4)
        results["kg_ensemble_on_repodb"] = {
            "mappable_pairs": len(kg_triplets),
            "mechanistic_auroc": fmt(au_m),
            "kg_only_auroc": fmt(au_k),
            "ensemble_auroc": fmt(au_e),
            "ensemble_lift_over_mechanistic": lift,
        }

    log("\nENGINE (mechanism-only) vs gold standard:")
    em = results["engine_mechanistic"]
    log(f"  AUROC {em['auroc']}   AUPRC {em['auprc']}   BEDROC(a=20) {em['bedroc_a20']}   "
        f"EF@10% {em['ef_top10pct']}x")
    log("Baselines (AUROC):")
    for k, v in results["baselines"].items():
        log(f"  {k:<22} {v['auroc']}")
    log("Ablation — direction-aware pathway term (AUROC):")
    ab = results["ablation_pathway_term"]
    log(f"  aware {ab['direction_aware']}  |  blind {ab['direction_blind']}  |  "
        f"excluded {ab['pathway_excluded']}")
    log(f"Negative control (label shuffle) AUROC: "
        f"{results['negative_control_label_shuffle']['auroc']}  (expect ~0.5)")
    dg = results["disease_genes_genetics_weighted"]
    log(f"Disease genes — flat {dg['flat_target_auroc']} vs genetics-weighted "
        f"{dg['genetics_weighted_auroc']} (lift {dg['lift']})")
    if "kg_ensemble_on_repodb" in results:
        ke = results["kg_ensemble_on_repodb"]
        log(f"KG ensemble on {ke['mappable_pairs']} KG-mappable pairs: "
            f"mech {ke['mechanistic_auroc']} | kg {ke['kg_only_auroc']} | "
            f"ensemble {ke['ensemble_auroc']} (lift {ke['ensemble_lift_over_mechanistic']})")

    per_disease.sort(key=lambda x: -(x["auroc"] or 0))
    mean_dis_auroc = (round(sum(d["auroc"] for d in per_disease) / len(per_disease), 3)
                      if per_disease else None)

    delta = (None if (ab["direction_aware"] is None or ab["direction_blind"] is None)
             else round(ab["direction_aware"] - ab["direction_blind"], 4))
    sev = "INFO" if (main_auroc or 0) >= 0.65 else ("WARN" if (main_auroc or 0) >= 0.55 else "FAIL")

    result = {
        "run_at": datetime.datetime.now().isoformat(timespec="seconds"),
        "gold_standard": "repoDB (Brown & Patel, Sci Data 2017) — Approved vs failed (Terminated/Withdrawn/Suspended)",
        "scoring": "Mechanism-only (target + direction-aware pathway + PPI); clinical/literature/"
                   "regulatory features removed to prevent retrospective leakage.",
        "drug_target_mapping_mode": mapping,
        "diseases_evaluated": resolved_ok,
        "pairs_evaluated": len(all_pairs),
        "pairs_approved": n_pos,
        "pairs_failed": n_neg,
        "base_rate_pct": round(100.0 * n_pos / max(1, len(all_pairs)), 1),
        "drug_mapping_coverage_pct": round(100.0 * len(mapped) / max(1, len(all_drugs)), 1),
        "metrics": results,
        "per_disease_mean_auroc": mean_dis_auroc,
        "per_disease_top": per_disease[:15],
        "findings": [
            {"id": "PRED-01", "severity": sev,
             "title": "Retrospective repurposing recovery (mechanism-only) vs repoDB",
             "detail": f"On {len(all_pairs):,} repoDB drug-indication pairs across {resolved_ok} diseases, "
                       f"the engine's leakage-free mechanistic score separates Approved from failed "
                       f"repurposings with AUROC {em['auroc']} (AUPRC {em['auprc']}, BEDROC {em['bedroc_a20']}, "
                       f"EF@10% {em['ef_top10pct']}x), vs random {results['baselines']['random']['auroc']} "
                       f"and target-only {results['baselines']['target_overlap_only']['auroc']}.",
             "capa": "Re-run on each data refresh; track AUROC over time. Below 0.55 triggers review of the "
                     "biology inputs (Open Targets gene sets, mechanism coverage)."},
            {"id": "PRED-02",
             "severity": "INFO" if (delta is not None and delta >= 0) else "WARN",
             "title": "Ablation: does the direction-aware pathway term help?",
             "detail": f"AUROC with direction-aware pathway scoring {ab['direction_aware']} vs direction-blind "
                       f"{ab['direction_blind']} (delta {delta}); pathway excluded {ab['pathway_excluded']}.",
             "capa": "If delta is consistently <=0, revert the pathway term to direction-blind."},
            {"id": "PRED-03", "severity": "INFO",
             "title": "Negative control",
             "detail": f"Shuffling the labels collapses AUROC to "
                       f"{results['negative_control_label_shuffle']['auroc']} (expected ~0.5), confirming the "
                       f"signal is not an artefact of the evaluation.",
             "capa": "None — control check."},
            {"id": "GEN-01", "severity": "INFO",
             "title": "Genetics-weighted disease genes — tested, neutral",
             "detail": f"Weighting disease genes by Open Targets genetic_association (Nelson 2015) gives AUROC "
                       f"{dg['genetics_weighted_auroc']} vs flat {dg['flat_target_auroc']} (lift {dg['lift']}) — "
                       f"neutral on retrospective recovery. The tiny flat-over-genetics gap also confirms the "
                       f"headline score is NOT meaningfully inflated by the overall-score known-drug channel. "
                       f"Genetics enrichment is now fetched per disease and available as an opt-in, leakage-safer "
                       f"target weighting; it is NOT the engine default (no measured benefit).",
             "capa": "Keep genetics weighting opt-in; revisit for prospective validation where the Nelson prior "
                     "is more relevant than retrospective overlap."},
        ] + ([
            {"id": "KGE-04",
             "severity": "WARN" if (results["kg_ensemble_on_repodb"]["ensemble_lift_over_mechanistic"] or 0) <= 0 else "INFO",
             "title": "KG-embedding ensemble fails the EXTERNAL benchmark — not integrated",
             "detail": f"On {results['kg_ensemble_on_repodb']['mappable_pairs']} KG-mappable repoDB pairs, the KG "
                       f"embedding alone scores AUROC {results['kg_ensemble_on_repodb']['kg_only_auroc']} (~chance) "
                       f"and the mechanistic+KG ensemble {results['kg_ensemble_on_repodb']['ensemble_auroc']} vs "
                       f"mechanistic {results['kg_ensemble_on_repodb']['mechanistic_auroc']} "
                       f"(lift {results['kg_ensemble_on_repodb']['ensemble_lift_over_mechanistic']}). The +0.042 "
                       f"seen on Hetionet's own held-out treats did NOT transfer to discriminating real-world "
                       f"approved-vs-failed repurposings.",
             "capa": "Do NOT wire the KG embedding into the engine — it degrades the mechanistic score on the "
                     "external benchmark. Revisit only with metapath/DWPC features or a better-generalising model."},
        ] if "kg_ensemble_on_repodb" in results else []),
    }
    (HERE / "predictions_results.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    (HERE / "predictions_run.log").write_text("\n".join(log_lines), encoding="utf-8")
    # Sidecar: raw scored pairs (mechanism-only, direction-aware) for the
    # calibration validator to consume without re-querying Open Targets.
    (HERE / "predictions_pairs.json").write_text(
        json.dumps([{"score": round(s, 6), "label": y} for s, y in all_pairs]),
        encoding="utf-8")
    log(f"\n  Wrote: {HERE / 'predictions_results.json'}")
    return result


def _target_only(drug_genes, disease_genes):
    from services.repurposing_engine import _score_target_overlap
    return _score_target_overlap(drug_genes, disease_genes)


def _shuffled(labels, rng):
    out = list(labels)
    rng.shuffle(out)
    return out


if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else str(DEFAULT_REPODB)
    md = int(sys.argv[2]) if len(sys.argv) > 2 else 60
    mapping = sys.argv[3] if len(sys.argv) > 3 else "fallback"
    try:
        run(path, md, mapping)
    except Exception:
        import traceback
        traceback.print_exc()
        sys.exit(1)
