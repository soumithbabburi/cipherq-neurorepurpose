"""
Repurposing-score audit — is the canonical score rigged?

Scores three classes of drug->disease pairs through the SAME entry point the UI uses
(canonical_pair_score) and reports:
  1. SEPARATION  — do true indications score >> plausible-wrong >> random pairs?
                   (a rigged/always-high score fails this; a sound one separates.)
  2. INFLATION   — how far the applicable-weight renorm + stacked bonuses lift a pair
                   above its plain weighted-mean of the 6 sub-scores.
  3. DOUBLE-COUNT — do the independent bonuses (network / proliferation / directional)
                   fire off the SAME genes, crediting one biology multiple times?
"""
import json
import statistics as st
from services.reverse_repurposing import resolve_drug, indication_phase_for
from services.repurposing_engine import (
    score_compound_for_disease, _actions_for_molecules, WEIGHTS)
from services.disease_ontology import resolve_disease, get_ppi_network

TRUE_POS = [
    ("imatinib",     "Chronic Myeloid Leukemia"),
    ("adalimumab",   "Rheumatoid Arthritis"),
    ("metformin",    "Type 2 Diabetes Mellitus"),
    ("atorvastatin", "Hypercholesterolemia"),
    ("sildenafil",   "Erectile Dysfunction"),
    ("donepezil",    "Alzheimer Disease"),
    ("methotrexate", "Psoriasis"),
    ("gefitinib",    "Non-Small Cell Lung Carcinoma"),
]
HARD_NEG = [   # real drug, plausible-sounding but NOT its indication (different organ/axis)
    ("imatinib",     "Alzheimer Disease"),
    ("adalimumab",   "Type 2 Diabetes Mellitus"),
    ("metformin",    "Rheumatoid Arthritis"),
    ("atorvastatin", "Schizophrenia"),
    ("sildenafil",   "Alzheimer Disease"),
    ("donepezil",    "Essential Hypertension"),
    ("methotrexate", "Type 2 Diabetes Mellitus"),
    ("gefitinib",    "Major Depressive Disorder"),
]
RANDOM = [
    ("loratadine",  "Parkinson Disease"),
    ("omeprazole",  "Multiple Sclerosis"),
    ("ibuprofen",   "Glaucoma"),
    ("warfarin",    "Asthma"),
    ("amoxicillin", "Osteoporosis"),
    ("cetirizine",  "Epilepsy"),
]
# The blind spot the old harness never tested: the drug WAS TRIED in this disease but
# is NOT approved (studied / failed). These must NOT read "Strong / established therapy"
# on prior clinical art alone — a zero-mechanism failed pair (imatinib->PAH) especially.
STUDIED_FAILED = [
    ("imatinib",    "Pulmonary Hypertension"),          # IMPRES: failed efficacy + safety
    ("imatinib",    "Idiopathic Pulmonary Fibrosis"),   # trialed, not approved
    ("mepolizumab", "Eosinophilic Esophagitis"),        # Phase 2/3, missed symptom endpoint
    ("metformin",   "Breast Cancer"),                   # heavily studied, not approved
    ("sildenafil",  "Idiopathic Pulmonary Fibrosis"),   # STEP-IPF, no benefit
]


def _genes(x):
    """Best-effort gene list from a bonus sub-dict."""
    if not isinstance(x, dict):
        return set()
    for k in ("genes", "network_genes", "path_genes", "shared_genes"):
        v = x.get(k)
        if isinstance(v, list) and v:
            return {str(g).upper() for g in v}
    return set()


def _score_pair(drug, dis):
    """Replicate canonical_pair_score's input assembly, then call the engine directly
    so we get the FULL breakdown (mech_plausibility, network, all penalty multipliers)."""
    info = resolve_drug(drug) or {}
    drug_genes = info.get("targets", []) or []
    chembl_id = info.get("chembl_id", "") or ""
    max_phase = info.get("max_phase", 0) or 0
    indications = "; ".join(k["name"] for k in info.get("known_indications", []))
    name = info.get("name", drug) or drug
    dinfo = resolve_disease(dis) or {}
    disease_genes = [t["gene_symbol"] for t in dinfo.get("targets", [])[:40] if t.get("gene_symbol")]
    disease_pathways = dinfo.get("pathways", [])
    ppi = get_ppi_network(disease_genes[:15]) if disease_genes else {}
    # per-indication development phase (drug's phase in THIS disease) — mirrors what the
    # reverse engine passes, so prior clinical art is phase-scaled (approved vs studied).
    ind_phase = indication_phase_for(dis, info.get("known_indications", []))
    compound = {"chembl_id": chembl_id, "name": name, "max_phase": max_phase,
                "indications": indications, "targets": ";".join(drug_genes)}
    drug_actions = None
    if chembl_id:
        try:
            drug_actions = _actions_for_molecules([chembl_id]).get(chembl_id)
        except Exception:
            pass
    disease_weights = {t["gene_symbol"].upper():
                       (t.get("quality_score") or t.get("genetic_score") or t.get("score", 0.0))
                       for t in dinfo.get("targets", []) if t.get("gene_symbol")}
    return score_compound_for_disease(
        compound, dis, disease_genes, disease_pathways, ppi, drug_genes,
        drug_actions=drug_actions, disease_gene_weights=disease_weights or None,
        indication_phase=ind_phase), drug_genes


def audit(pair):
    drug, dis = pair
    try:
        sr, drug_genes = _score_pair(drug, dis)
    except Exception as e:
        return {"drug": drug, "disease": dis, "error": str(e), "composite": 0.0,
                "plain_wmean": 0.0, "inflation": 0.0, "penalty_factor": 1.0,
                "tier": None, "subscores": {}, "bonuses_fired": {}, "n_bonuses": 0,
                "shared_bonus_genes": []}
    s = sr.get("scores", {}) or {}
    final = float(sr.get("composite_score") or 0.0)
    mech = float(sr.get("mechanistic_plausibility") or 0.0)   # post renorm+bonuses, PRE penalty
    plain = sum(float(s.get(k, 0.0)) * WEIGHTS[k] for k in WEIGHTS)
    bd = sr.get("score_breakdown", {}) or {}
    prolif = sr.get("proliferation", {}) or {}
    rev = sr.get("signature_reversal", {}) or {}
    direc = sr.get("directional_evidence", {}) or {}
    bonuses = {
        "network":     float(bd.get("network_score") or 0.0),
        "prolif":      float(prolif.get("score") or 0.0) if prolif.get("match") else 0.0,
        "reversal":    float(rev.get("score") or 0.0) if rev.get("direction") == "reversing" else 0.0,
        "directional": float(direc.get("signal") or 0.0),
    }
    fired = {k: v for k, v in bonuses.items() if v > 0}
    # penalty multipliers actually applied to the composite
    pens = {
        "safety":      float((sr.get("safety") or {}).get("multiplier", 1.0)),
        "coverage":    float((sr.get("coverage") or {}).get("multiplier", 1.0)),
        "clinical":    float((sr.get("clinical_constraints") or {}).get("multiplier", 1.0)),
        "delivery":    float((sr.get("delivery") or {}).get("multiplier", 1.0)),
        "trial_fail":  float((sr.get("trial_failure") or {}).get("multiplier", 1.0)),
        "appropriate": float((sr.get("appropriateness") or {}).get("factor", 1.0)),
    }
    active_pens = {k: round(v, 3) for k, v in pens.items() if v < 0.999}
    # gene sets each bonus exposes → double-count detection
    genesets = {"network": _genes(bd), "directional": _genes(direc), "prolif": _genes(prolif)}
    genesets = {k: v for k, v in genesets.items() if v}
    shared, keys = set(), list(genesets)
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            shared |= (genesets[keys[i]] & genesets[keys[j]])
    return {
        "drug": drug, "disease": dis,
        "composite": round(final, 4),
        "mech_plaus": round(mech, 4),
        "plain_wmean": round(plain, 4),
        "inflation": round(mech - plain, 4),            # renorm+bonuses lift over naive mean
        "penalty_factor": round(final / mech, 3) if mech > 1e-6 else 1.0,
        "tier": (sr.get("calibration") or {}).get("tier"),
        "act_tier": (sr.get("actionability") or {}).get("tier"),
        "subscores": {k: round(float(s.get(k, 0.0)), 3) for k in WEIGHTS},
        "bonuses_fired": {k: round(v, 3) for k, v in fired.items()},
        "penalties": active_pens,
        "n_bonuses": len(fired),
        "shared_bonus_genes": sorted(shared),
    }


def summarize(name, rows):
    comps = [r["composite"] for r in rows]
    infl = [r["inflation"] for r in rows]
    print(f"\n=== {name}  (n={len(rows)}) ===")
    print(f"  composite: mean={st.mean(comps):.3f} median={st.median(comps):.3f} "
          f"min={min(comps):.3f} max={max(comps):.3f}")
    print(f"  inflation (mech_plaus - plain wmean): mean={st.mean(infl):+.3f} max={max(infl):+.3f}")
    for r in rows:
        dc = f"  DBL-COUNT={r['shared_bonus_genes']}" if r["shared_bonus_genes"] else ""
        print(f"  {r['composite']:.3f} [{str(r['tier'])[:9]:9}|{str(r.get('act_tier'))[:12]:12}] "
              f"{r['drug']}->{r['disease'][:22]:22} "
              f"plain={r['plain_wmean']:.2f} mech={r['mech_plaus']:.2f} pen={r['penalty_factor']:.2f} "
              f"b={r['bonuses_fired']} p={r.get('penalties',{})}{dc}")
    return comps


def auc(pos, neg):
    """Fraction of (pos,neg) pairs correctly ordered — a simple ranking AUC."""
    if not pos or not neg:
        return None
    wins = sum((p > n) + 0.5 * (p == n) for p in pos for n in neg)
    return wins / (len(pos) * len(neg))


if __name__ == "__main__":
    print("Scoring pairs through canonical_pair_score (live)...")
    tp = [audit(p) for p in TRUE_POS]
    hn = [audit(p) for p in HARD_NEG]
    rn = [audit(p) for p in RANDOM]
    sf = [audit(p) for p in STUDIED_FAILED]
    c_tp = summarize("TRUE POSITIVES", tp)
    c_hn = summarize("HARD NEGATIVES (plausible-wrong)", hn)
    c_sf = summarize("STUDIED / FAILED (tried, not approved)", sf)
    c_rn = summarize("RANDOM PAIRS", rn)
    print("\n=== SEPARATION (ranking AUC) ===")
    print(f"  TP vs RANDOM         : {auc(c_tp, c_rn)}")
    print(f"  TP vs HARD-NEG       : {auc(c_tp, c_hn)}")
    print(f"  TP vs STUDIED/FAILED : {auc(c_tp, c_sf)}")
    print(f"  HARD-NEG vs RANDOM   : {auc(c_hn, c_rn)}")
    # the F2 trust check: a studied/failed pair must NOT read 'Strong' on prior art alone
    strong_sf = [r for r in sf if str(r.get("tier")) == "Strong"]
    print(f"\n=== F2 CHECK: studied/failed pairs reading 'Strong': {len(strong_sf)}/{len(sf)} ===")
    for r in strong_sf:
        print(f"    !! {r['drug']}->{r['disease']}  {r['composite']}  subs={r['subscores']}")
    all_rows = tp + hn + rn
    dc_rows = [r for r in all_rows if r["shared_bonus_genes"]]
    multi = [r for r in all_rows if r["n_bonuses"] >= 2]
    print("\n=== DOUBLE-COUNTING ===")
    print(f"  pairs with >=2 bonuses firing: {len(multi)}/{len(all_rows)}")
    print(f"  pairs where bonuses share genes: {len(dc_rows)}/{len(all_rows)}")
    for r in dc_rows:
        print(f"    {r['drug']}->{r['disease'][:24]}: {list(r['bonuses_fired'])} share {r['shared_bonus_genes']}")
    json.dump(all_rows, open("scratch_audit_rows.json", "w"), indent=2)
    print("\nwrote scratch_audit_rows.json")
