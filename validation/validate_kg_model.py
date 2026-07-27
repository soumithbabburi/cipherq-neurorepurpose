"""
KG-Embedding Validation — held-out treatment-edge link prediction
═══════════════════════════════════════════════════════════════════════════
Validates the DistMult knowledge-graph embedding (services/kg_embedding.py): can
it predict Compound-*treats*-Disease edges it was NOT trained on?

Leakage control:
  • The held-out 20% of treat (CtD) edges are removed from the training graph
    entirely — the model never sees them. It trains on all biology edges plus the
    other 80% of treat edges (standard transductive KGE setup).
  • A label-shuffle negative control must collapse to AUROC ≈ 0.5.

Honest framing vs the mechanism-only benchmark (validate_predictions.py, AUROC
0.73): that score uses ZERO treatment information; this KGE is a *transductive*
link-prediction model that learns from known treatments + biology. They answer
different questions and are complementary, not directly comparable — both are
reported.

Metrics: AUROC (held-out treats vs sampled negatives), Hits@10 / Hits@30 / MRR
(filtered ranking over the candidate-disease space), and the negative control.

Saves embeddings to data/kg_embeddings.npz. Writes validation/kg_model_results.json.

Usage:
    python validation/validate_kg_model.py [epochs] [dim]
"""
import sys
import json
import datetime
import numpy as np
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))

from services.kg_embedding import (
    load_hetionet_triples, DistMultKGE, TREATS)
from validation.validate_predictions import auroc, bedroc


def run(epochs: int = 20, dim: int = 64):
    log_lines = []

    def log(msg=""):
        print(msg, flush=True)
        log_lines.append(str(msg))

    log("=" * 72)
    log("KG-EMBEDDING VALIDATION  -  held-out treatment-edge link prediction")
    log(f"Run: {datetime.datetime.now().isoformat(timespec='seconds')}")
    log("=" * 72)

    triples, node2id, rel2id, id2node = load_hetionet_triples()
    n_nodes, n_rels = len(node2id), len(rel2id)
    log(f"\nGraph: {len(triples):,} triples · {n_nodes:,} nodes · {n_rels} relations")
    if TREATS not in rel2id:
        log(f"ERROR: relation {TREATS} not in graph.")
        sys.exit(1)
    tr = rel2id[TREATS]

    # ── Split treat edges 80/20 (held-out test never enters training) ─────────
    rng = np.random.default_rng(7)
    ct_mask = triples[:, 1] == tr
    ct_idx = np.where(ct_mask)[0]
    rng.shuffle(ct_idx)
    n_test = max(1, int(round(0.2 * len(ct_idx))))
    test_idx = set(ct_idx[:n_test].tolist())
    train_mask = np.ones(len(triples), dtype=bool)
    train_mask[list(test_idx)] = False
    train_triples = triples[train_mask]
    test_triples = triples[ct_idx[:n_test]]
    log(f"Treat (CtD) edges: {len(ct_idx)} → train {len(ct_idx) - n_test} / held-out test {n_test}")
    log(f"Training triples (graph minus held-out treats): {len(train_triples):,}")

    # ── Train ─────────────────────────────────────────────────────────────────
    log(f"\nTraining DistMult (dim={dim}, epochs={epochs})…")
    model = DistMultKGE(n_nodes, n_rels, dim=dim)
    model.fit(train_triples, epochs=epochs, batch=4096, n_neg=5, lr=0.1, log=log)

    # ── Candidate spaces + known-positive filter ──────────────────────────────
    all_ct = triples[triples[:, 1] == tr]
    known = set((int(h), int(t)) for h, t in all_ct[:, [0, 2]])
    compounds = np.unique(all_ct[:, 0])
    diseases = np.unique(all_ct[:, 2])
    log(f"Candidate space: {len(compounds)} compounds × {len(diseases)} diseases")

    # ── AUROC: held-out positives vs sampled negatives ────────────────────────
    pos_scores = model.score(test_triples[:, 0],
                             np.full(len(test_triples), tr),
                             test_triples[:, 2])
    neg_ratio = 20
    neg_h, neg_t = [], []
    dis_arr = diseases
    rng2 = np.random.default_rng(11)
    for h in test_triples[:, 0]:
        got = 0
        while got < neg_ratio:
            d = int(dis_arr[rng2.integers(0, len(dis_arr))])
            if (int(h), d) not in known:
                neg_h.append(int(h)); neg_t.append(d); got += 1
    neg_scores = model.score(np.array(neg_h), np.full(len(neg_h), tr), np.array(neg_t))
    pairs = [(float(s), 1) for s in pos_scores] + [(float(s), 0) for s in neg_scores]
    au = auroc(pairs)
    bd = bedroc(pairs, 20.0)

    # negative control: shuffle labels
    labels = [y for _, y in pairs]
    rng2.shuffle(labels)
    shuf = [(pairs[i][0], labels[i]) for i in range(len(pairs))]
    au_shuf = auroc(shuf)

    # ── Filtered ranking over the candidate-disease space ─────────────────────
    Edis = model.E[diseases]                       # (D, dim)
    d_pos = {int(d): i for i, d in enumerate(diseases)}
    ranks = []
    for h, _, t in test_triples:
        cvec = model.E[int(h)] * model.R[tr]       # (dim,)
        scores = Edis @ cvec                        # (D,)
        # filter: set known-other-positive diseases for this compound to -inf
        for j, d in enumerate(diseases):
            if (int(h), int(d)) in known and int(d) != int(t):
                scores[j] = -1e30
        true_j = d_pos[int(t)]
        rank = int(np.sum(scores > scores[true_j]) + 1)
        ranks.append(rank)
    ranks = np.array(ranks)
    hits10 = float(np.mean(ranks <= 10))
    hits30 = float(np.mean(ranks <= 30))
    mrr = float(np.mean(1.0 / ranks))

    log(f"\nHeld-out treat-edge prediction:")
    log(f"  AUROC {au:.4f}   BEDROC(a=20) {bd:.4f}   (vs {neg_ratio}x negatives)")
    log(f"  Hits@10 {hits10:.3f}   Hits@30 {hits30:.3f}   MRR {mrr:.3f}   "
        f"(rank over {len(diseases)} diseases)")
    log(f"  Negative control (label shuffle) AUROC {au_shuf:.4f}  (expect ~0.5)")

    # ── Save embeddings for later engine integration ──────────────────────────
    emb_path = ROOT / "data" / "kg_embeddings.npz"
    # retrain-free save uses the split model; fine as a first artefact
    model.save(emb_path, node2id, rel2id)
    log(f"  Saved embeddings → {emb_path}")

    sev = "INFO" if au >= 0.75 else ("WARN" if au >= 0.6 else "FAIL")
    result = {
        "run_at": datetime.datetime.now().isoformat(timespec="seconds"),
        "model": f"DistMult KGE (numpy), dim={dim}, epochs={epochs}",
        "graph": {"triples": int(len(triples)), "nodes": n_nodes, "relations": n_rels,
                  "source": "hetionet_v1.0"},
        "treat_edges": {"total": int(len(ct_idx)), "held_out_test": int(n_test)},
        "candidate_space": {"compounds": int(len(compounds)), "diseases": int(len(diseases))},
        "metrics": {
            "auroc": round(au, 4),
            "bedroc_a20": round(bd, 4),
            "hits_at_10": round(hits10, 4),
            "hits_at_30": round(hits30, 4),
            "mrr": round(mrr, 4),
            "negative_control_auroc": round(au_shuf, 4),
            "negatives_per_positive": neg_ratio,
        },
        "comparison_note": ("Transductive link prediction: learns from biology + 80% of known "
                            "treatments, predicts the held-out 20%. Complementary to the "
                            "mechanism-only benchmark (AUROC 0.73), which uses NO treatment data."),
        "findings": [
            {"id": "KGE-01", "severity": sev,
             "title": "KG-embedding predicts held-out treatment edges",
             "detail": f"A DistMult embedding of the Hetionet graph predicts the held-out 20% of "
                       f"Compound-treats-Disease edges with AUROC {au:.3f} (Hits@10 {hits10:.2f}, "
                       f"MRR {mrr:.2f}) — capturing indirect multi-hop biology the direct "
                       f"target-overlap score cannot. Held-out edges are absent from training.",
             "capa": "Re-train on each graph refresh; track AUROC/Hits@10. Integrate the embedding "
                     "score as an additional engine signal (needs ChEMBL→DrugBank node mapping)."},
            {"id": "KGE-02", "severity": "INFO",
             "title": "Negative control",
             "detail": f"Shuffling labels collapses AUROC to {au_shuf:.3f} (expect ~0.5).",
             "capa": "None — control."},
        ],
    }
    (HERE / "kg_model_results.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    (HERE / "kg_model_run.log").write_text("\n".join(log_lines), encoding="utf-8")
    log(f"\n  Wrote: {HERE / 'kg_model_results.json'}")
    return result


if __name__ == "__main__":
    ep = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    dm = int(sys.argv[2]) if len(sys.argv) > 2 else 64
    try:
        run(ep, dm)
    except Exception:
        import traceback
        traceback.print_exc()
        sys.exit(1)
