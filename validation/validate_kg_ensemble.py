"""
KG-Embedding complementarity test — does it ADD to gene overlap?
═══════════════════════════════════════════════════════════════════════════
Standalone DistMult underperforms the mechanistic score. The real question for
keeping the KG work is: does the embedding carry signal ORTHOGONAL to direct
gene overlap, so that an ensemble beats overlap alone?

Self-contained on Hetionet (no cross-ontology mapping): on the SAME held-out 20%
of treat edges used by validate_kg_model.py, compare three scorers on identical
positives + sampled negatives:
  • gene-overlap  : |genes(compound) ∩ genes(disease)| from Hetionet biology
                    edges (CbG/CuG/CdG ∩ DaG/DuG/DdG) — the Hetionet-native
                    analog of the mechanistic target-overlap score.
  • kg-embedding  : DistMult treats-score from data/kg_embeddings.npz.
  • ensemble      : rank-average of the two.

Leakage-safe: the embedding was trained without the held-out treats; gene overlap
uses only biology edges (never the treat label).

Writes validation/kg_ensemble_results.json + log.
"""
import sys
import json
import datetime
import numpy as np
from pathlib import Path
from collections import defaultdict

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))

from services.kg_embedding import load_hetionet_triples, DistMultKGE, TREATS
from validation.validate_predictions import auroc

COMPOUND_GENE = {"CbG", "CuG", "CdG"}
DISEASE_GENE = {"DaG", "DuG", "DdG"}


def _ranks(values):
    """Average ranks (1..n), ties averaged — for rank-average ensembling."""
    order = np.argsort(values)
    r = np.empty(len(values), dtype=float)
    i = 0
    while i < len(values):
        j = i
        while j < len(values) and values[order[j]] == values[order[i]]:
            j += 1
        r[order[i:j]] = (i + 1 + j) / 2.0
        i = j
    return r


def run():
    log_lines = []

    def log(m=""):
        print(m, flush=True); log_lines.append(str(m))

    log("=" * 72)
    log("KG-EMBEDDING COMPLEMENTARITY  -  does it add to gene overlap?")
    log(f"Run: {datetime.datetime.now().isoformat(timespec='seconds')}")
    log("=" * 72)

    emb_path = ROOT / "data" / "kg_embeddings.npz"
    if not emb_path.exists():
        log("ERROR: data/kg_embeddings.npz missing — run validate_kg_model.py first.")
        sys.exit(1)

    triples, node2id, rel2id, id2node = load_hetionet_triples()
    model, n2i2, r2i2 = DistMultKGE.load(emb_path)
    tr = rel2id[TREATS]

    # Reconstruct the SAME held-out split as validate_kg_model.py (seed 7, 20%)
    rng = np.random.default_rng(7)
    ct_idx = np.where(triples[:, 1] == tr)[0]
    rng.shuffle(ct_idx)
    n_test = max(1, int(round(0.2 * len(ct_idx))))
    test = triples[ct_idx[:n_test]]
    known = set((int(h), int(t)) for h, t in triples[triples[:, 1] == tr][:, [0, 2]])
    diseases = np.unique(triples[triples[:, 1] == tr][:, 2])
    log(f"\nHeld-out test treats: {len(test)} | candidate diseases: {len(diseases)}")

    # Gene neighbour sets from biology edges
    cgene = defaultdict(set)
    dgene = defaultdict(set)
    cg_rels = {rel2id[m] for m in COMPOUND_GENE if m in rel2id}
    dg_rels = {rel2id[m] for m in DISEASE_GENE if m in rel2id}
    for h, r, t in triples:
        if r in cg_rels:
            cgene[int(h)].add(int(t))
        elif r in dg_rels:
            dgene[int(h)].add(int(t))

    def overlap(c, d):
        gc, gd = cgene.get(int(c)), dgene.get(int(d))
        if not gc or not gd:
            return 0.0
        return len(gc & gd) / np.sqrt(len(gc) * len(gd))   # cosine-like overlap

    # Build positive + negative pairs (20x negatives)
    rng2 = np.random.default_rng(11)
    pairs = []   # (compound, disease, label)
    for h, _, t in test:
        pairs.append((int(h), int(t), 1))
        got = 0
        while got < 20:
            d = int(diseases[rng2.integers(0, len(diseases))])
            if (int(h), d) not in known:
                pairs.append((int(h), d, 0)); got += 1

    H = np.array([p[0] for p in pairs]); D = np.array([p[1] for p in pairs])
    Y = [p[2] for p in pairs]
    kg = model.score(H, np.full(len(H), tr), D)
    ov = np.array([overlap(h, d) for h, d in zip(H, D)])
    ens = _ranks(kg) + _ranks(ov)          # rank-average ensemble

    au_ov = auroc(list(zip(ov.tolist(), Y)))
    au_kg = auroc(list(zip(kg.tolist(), Y)))
    au_ens = auroc(list(zip(ens.tolist(), Y)))

    log(f"\nAUROC on held-out treats (pos + 20x neg):")
    log(f"  gene-overlap only : {au_ov:.4f}")
    log(f"  kg-embedding only : {au_kg:.4f}")
    log(f"  ensemble (rank-avg): {au_ens:.4f}")
    best_single = max(au_ov, au_kg)
    lift = au_ens - best_single
    log(f"  ensemble lift over best single: {lift:+.4f}")

    sev = "INFO" if lift > 0.01 else ("WARN" if lift > -0.01 else "WARN")
    verdict = ("Ensemble beats both singles — the KG embedding carries orthogonal signal worth keeping."
               if lift > 0.01 else
               "Ensemble does not beat gene overlap — the KG embedding adds little here.")
    result = {
        "run_at": datetime.datetime.now().isoformat(timespec="seconds"),
        "held_out_treats": int(len(test)),
        "auroc": {"gene_overlap": round(au_ov, 4), "kg_embedding": round(au_kg, 4),
                  "ensemble": round(au_ens, 4)},
        "ensemble_lift_over_best_single": round(lift, 4),
        "verdict": verdict,
        "findings": [{
            "id": "KGE-03", "severity": sev,
            "title": "KG-embedding complementarity vs gene overlap",
            "detail": f"On held-out treat edges: gene-overlap AUROC {au_ov:.3f}, KG-embedding {au_kg:.3f}, "
                      f"rank-average ensemble {au_ens:.3f} (lift {lift:+.3f} over the best single). {verdict}",
            "capa": "Keep the KG embedding as an engine signal only if the ensemble lift is positive and "
                    "reproducible; otherwise prefer the simpler mechanistic score.",
        }],
    }
    (HERE / "kg_ensemble_results.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    (HERE / "kg_ensemble_run.log").write_text("\n".join(log_lines), encoding="utf-8")
    log(f"\n  {verdict}")
    log(f"  Wrote: {HERE / 'kg_ensemble_results.json'}")
    return result


if __name__ == "__main__":
    try:
        run()
    except Exception:
        import traceback; traceback.print_exc(); sys.exit(1)
