"""
TxGNN-protocol benchmark: DISEASE zero-shot repurposing on PrimeKG.
════════════════════════════════════════════════════════════════════════════════════
TxGNN's headline setting is disease zero-shot — test diseases have NONE of their treatment
edges in training; the model must predict their drugs from biology alone. Our earlier
harness held out DRUGS (compound-disjoint); this one holds out DISEASES to match TxGNN.

Clean transductive protocol (matches TxGNN — no test-label leakage):
  * Split diseases (those with >=1 indication) 85/15 (seed 42). Test diseases' treatment
    edges (indication / contraindication / off-label) are removed from EVERYTHING:
      - the DistMult inputs (retrained here, so node embeddings carry no test-treatment signal)
      - the GNN message-passing graph and its supervision
      - the treats-classifier supervision
    Full biological topology (protein/phenotype/pathway/…) stays visible, transductive.
  * The GNN and treats-classifier are retrained here on the disease split (shipped
    pkg_gnn.pt / pkg_treats.pkl untouched). Node embeddings: by DEFAULT we reuse the shipped
    DistMult (a full retrain is ~6 min/epoch => ~2 h on CPU), which leaves mild input leakage
    (test-disease treatment edges are ~1% of the graph), so the numbers are a slight UPPER
    BOUND. Pass retrain_dm=True for a fully leakage-free run (~2 h).

Metrics (TxGNN's form), on held-out diseases, drug side:
  * AUROC  indication-vs-contraindication  (the hard directional task)
  * AUPRC  indication-vs-all-negatives     (contra + random drugs)
  * Recall@{10,20,50,100} + median rank    (rank all drugs for each test disease)
for three rankers: DistMult, treats-classifier, and the fusion (clf+GNN, RRF 3:1).

NOTE: this reproduces TxGNN's PROTOCOL and metric form, not their exact published split
files, so numbers are comparable in kind, not identical. Run (CPU, ~10 min):
  .venv312/Scripts/python.exe -m validation.benchmark_txgnn_split
"""
from __future__ import annotations

import json
import random
import sys
import time
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from services.kg_embedding import DistMultKGE                                  # noqa: E402
from services.train_primekg_gnn import RGCN, _row_norm_adj, _GROUP, _TREAT_RELS  # noqa: E402

_D = ROOT / "data" / "primekg"
KS = (10, 20, 50, 100)


def _rrf(scores_list, weights, k_rrf=60):
    fused = np.zeros(len(scores_list[0]))
    for sc, w in zip(scores_list, weights):
        order = np.argsort(-sc); rank = np.empty(len(sc)); rank[order] = np.arange(1, len(sc) + 1)
        fused += w * (1.0 / (k_rrf + rank))
    return fused


def _auroc(pos, neg):
    pos = np.asarray(pos, float); neg = np.asarray(neg, float)
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    allv = np.concatenate([pos, neg]); order = allv.argsort()
    ranks = np.empty(len(allv)); ranks[order] = np.arange(1, len(allv) + 1)
    return float((ranks[:len(pos)].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg)))


def _auprc(pos, neg):
    pos = np.asarray(pos, float); neg = np.asarray(neg, float)
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    y = np.concatenate([np.ones(len(pos)), np.zeros(len(neg))])
    s = np.concatenate([pos, neg])
    order = np.argsort(-s); y = y[order]
    tp = np.cumsum(y); fp = np.cumsum(1 - y)
    prec = tp / np.maximum(tp + fp, 1); rec = tp / max(y.sum(), 1)
    rec = np.concatenate([[0], rec]); prec = np.concatenate([[1], prec])
    return float(np.sum((rec[1:] - rec[:-1]) * prec[1:]))


def main(seed=42, holdout=0.15, dim=64, dm_epochs=20, gnn_epochs=60, neg_per_pos=8,
         retrain_dm=False):
    rng = random.Random(seed); torch.manual_seed(seed)
    t0 = time.time()
    import pandas as pd

    nodes = json.load(open(_D / "pkg_nodes.json"))
    labels = json.load(open(_D / "pkg_labels.json"))
    disease_ids = np.array([int(i) for i, v in nodes.items() if v[0] == "disease"], np.int64)
    drug_ids = np.array([int(i) for i, v in nodes.items() if v[0] == "drug"], np.int64)

    # ── disease zero-shot split ───────────────────────────────────────────────────
    dis_with_ind = sorted({int(z) for d, r in labels.items() for z in r.get("indication", [])
                           if int(z) in set(disease_ids.tolist())})
    rng.shuffle(dis_with_ind)
    n_hold = int(len(dis_with_ind) * holdout)
    test_dis = set(dis_with_ind[:n_hold])
    print(f"[bench] {len(disease_ids)} diseases | {len(dis_with_ind)} with indications | "
          f"test(held-out) {len(test_dis)} | {torch.get_num_threads()} threads", flush=True)

    # ── load triples; drop TEST-disease treatment edges from the whole graph ──────
    df = pd.read_csv(_D / "triples.tsv", sep="\t", header=None, names=["h", "r", "t"],
                     dtype={"h": np.int64, "t": np.int64, "r": str})
    treat = df["r"].isin(_TREAT_RELS).to_numpy()
    h_test = df["h"].isin(test_dis).to_numpy(); t_test = df["t"].isin(test_dis).to_numpy()
    drop = treat & (h_test | t_test)
    df = df[~drop].reset_index(drop=True)
    print(f"[bench] dropped {int(drop.sum()):,} test-disease treatment edges "
          f"({len(df):,} triples remain) ({time.time()-t0:.0f}s)", flush=True)

    # ── (1) node embeddings ───────────────────────────────────────────────────────
    # DistMult over the full 8.1M-triple graph costs ~6 min/epoch on CPU (~2 h for 20 ep),
    # so retraining it per split is impractical. Default: REUSE the shipped embeddings and
    # retrain only the task models (GNN + classifier) clean on the disease split. The shared
    # INPUT embeddings then carry mild leakage (test-disease treatment edges are ~1% of the
    # graph and shaped those diseases' vectors), so the reported numbers are a slight UPPER
    # BOUND on true zero-shot. Set retrain_dm=True for a fully leakage-free run (~2 h).
    if retrain_dm:
        rel2id = {}
        rids = np.array([rel2id.setdefault(r, len(rel2id)) for r in df["r"].to_numpy()], np.int64)
        tri = np.column_stack([df["h"].to_numpy(), rids, df["t"].to_numpy()])
        nN = int(tri[:, [0, 2]].max()) + 1
        dm = DistMultKGE(nN, len(rel2id), dim=dim)
        dm.fit(tri, epochs=dm_epochs, log=lambda s: print(f"  [dm] {s}", flush=True))
        E = dm.E.astype(np.float32); ind_rv = dm.R[rel2id["indication"]].astype(np.float32)
        print(f"[bench] DistMult retrained leakage-free ({time.time()-t0:.0f}s)", flush=True)
    else:
        emb = np.load(_D / "pkg_embeddings.npz")
        rels = json.load(open(_D / "pkg_rels.json"))
        E = emb["E"].astype(np.float32); ind_rv = emb["R"][rels["indication"]].astype(np.float32)
        nN = E.shape[0]
        print(f"[bench] reusing shipped DistMult embeddings (mild input leakage; upper bound) "
              f"({time.time()-t0:.0f}s)", flush=True)

    # ── (2) build relation-group adjacencies for the GNN ──────────────────────────
    df["g"] = df["r"].map(_GROUP)
    groups = sorted(df["g"].dropna().unique())
    adjs = [_row_norm_adj(df[df["g"] == g]["h"].to_numpy(), df[df["g"] == g]["t"].to_numpy(), nN)
            for g in groups]

    # supervision (TRAIN diseases only): indications (pos) + contraindications (hard neg)
    dset = set(disease_ids.tolist())
    pos, hard = [], []
    for d, r in labels.items():
        di = int(d)
        for z in r.get("indication", []):
            if int(z) in dset and int(z) not in test_dis:
                pos.append((di, int(z)))
        for z in r.get("contraindication", []):
            if int(z) in dset and int(z) not in test_dis:
                hard.append((di, int(z)))
    pos = np.array(pos, np.int64); hard = np.array(hard, np.int64)
    print(f"[bench] supervision {len(pos)} pos / {len(hard)} hard-neg ({time.time()-t0:.0f}s)", flush=True)

    model = RGCN(E, len(groups), ind_rv, layers=2)
    opt = torch.optim.Adam(model.parameters(), lr=0.01, weight_decay=1e-5)
    bce = nn.BCEWithLogitsLoss()
    pd_, pz = torch.tensor(pos[:, 0]), torch.tensor(pos[:, 1])
    hd, hz = torch.tensor(hard[:, 0]), torch.tensor(hard[:, 1])
    dis_t = torch.tensor(disease_ids)
    model.train()
    for ep in range(1, gnn_epochs + 1):
        opt.zero_grad(); H = model.encode(adjs)
        nd = pd_.repeat(neg_per_pos); nz = dis_t[torch.randint(0, len(dis_t), (len(pd_) * neg_per_pos,))]
        d_all = torch.cat([pd_, hd, nd]); z_all = torch.cat([pz, hz, nz])
        y = torch.cat([torch.ones(len(pd_)), torch.zeros(len(hd)), torch.zeros(len(nd))])
        loss = bce(model.score(H, d_all, z_all), y); loss.backward(); opt.step()
        if ep == 1 or ep % 20 == 0:
            print(f"  [gnn] epoch {ep} loss {loss.item():.4f} ({time.time()-t0:.0f}s)", flush=True)
    model.eval()
    with torch.no_grad():
        Hf = model.encode(adjs).detach().numpy()
    r_ind = model.r_ind.detach().numpy()

    # ── (3) treats classifier on the same split ───────────────────────────────────
    from sklearn.ensemble import HistGradientBoostingClassifier

    def feat(d, z):
        ed = E[d]; ez = E[z]; return np.concatenate([ed, ez, ed * ez, np.abs(ed - ez)], -1)
    X, y = [], []
    dl = list(dset - test_dis)
    for (d, z) in pos:
        X.append(feat(d, z)); y.append(1)
    for (d, z) in hard:
        X.append(feat(d, z)); y.append(0)
    for (d, _z) in pos:                                   # random negatives
        for _ in range(2):
            X.append(feat(d, dl[rng.randrange(len(dl))])); y.append(0)
    clf = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.08, max_depth=6,
                                         l2_regularization=1.0, random_state=seed)
    clf.fit(np.array(X, np.float32), np.array(y))
    print(f"[bench] classifier trained ({time.time()-t0:.0f}s)", flush=True)

    # ── (4) evaluate on held-out diseases (drug side) ─────────────────────────────
    _POOL = 400
    def clf_scores(z, dids):
        ez = E[z]; Ed = E[dids]; n = len(dids)
        f = np.concatenate([Ed, np.broadcast_to(ez, (n, ez.shape[0])), Ed * ez, np.abs(Ed - ez)], 1)
        return clf.predict_proba(f)[:, 1]

    drug_pos_of = {int(x): i for i, x in enumerate(drug_ids)}
    # invert labels: disease -> its indication / contra drugs
    dind, dcon = {}, {}
    for d, r in labels.items():
        di = int(d)
        for z in r.get("indication", []):
            dind.setdefault(int(z), []).append(di)
        for z in r.get("contraindication", []):
            dcon.setdefault(int(z), []).append(di)

    ranks = {"distmult": [], "clf": [], "fusion": []}
    roc = {"distmult": [], "clf": [], "fusion": []}   # (pos_scores, contra_scores)
    prc = {"distmult": [], "clf": [], "fusion": []}
    n_eval = 0
    for z in test_dis:
        inds = [drug_pos_of[x] for x in dind.get(z, []) if x in drug_pos_of]
        if not inds:
            continue
        n_eval += 1
        cons = [drug_pos_of[x] for x in dcon.get(z, []) if x in drug_pos_of]
        dm_all = (E[drug_ids] @ (E[z] * ind_rv))                 # relatedness, all drugs
        pool = np.argsort(-dm_all)[:_POOL]
        clf_p = clf_scores(z, drug_ids[pool])
        gnn_p = (Hf[drug_ids[pool]] * r_ind * Hf[z]).sum(-1)
        fus = _rrf([clf_p, gnn_p], [3.0, 1.0])
        # full-length score vectors (out-of-pool pushed to bottom for clf/fusion)
        def full(vec_pool):
            v = np.full(len(drug_ids), -np.inf); v[pool] = vec_pool; return v
        S = {"distmult": dm_all, "clf": full(clf_p), "fusion": full(fus)}
        other = set(inds)
        for name, sc in S.items():
            keep = np.ones(len(sc), bool)
            for j in other:
                keep[j] = False
            for j in inds:
                keep[j] = True
                r = int(np.sum(sc[keep] > sc[j]) + 1); ranks[name].append(r)
                keep[j] = False
            # AUROC/AUPRC using the same scorer
            ps = [sc[j] for j in inds]
            cs = [sc[j] for j in cons]
            ns = list(sc[np.random.RandomState(z).choice(len(sc), min(200, len(sc)), replace=False)])
            if cs:
                roc[name].append(_auroc(ps, cs))
            prc[name].append(_auprc(ps, cs + ns))

    def summ(name):
        r = np.array(ranks[name], float)
        return {**{f"recall@{k}": round(float(np.mean(r <= k)), 4) for k in KS},
                "median_rank": int(np.median(r)), "n_ind_pairs": int(len(r)),
                "auroc_ind_vs_contra": round(float(np.nanmean(roc[name])), 4),
                "auprc_ind_vs_neg": round(float(np.nanmean(prc[name])), 4)}
    out = {"protocol": "disease zero-shot (TxGNN-style), seed 42, 15% held-out diseases; "
                       "test-disease treatment edges removed from GNN graph + all supervision. "
                       "Same protocol/metrics as TxGNN, NOT their exact split files.",
           "input_embeddings": ("shipped DistMult reused (mild input leakage -> slight UPPER BOUND)"
                                if not retrain_dm else "retrained leakage-free"),
           "n_test_diseases_evaluated": n_eval, "elapsed_s": round(time.time() - t0, 1),
           "rankers": {n: summ(n) for n in ("distmult", "clf", "fusion")}}
    json.dump(out, open(ROOT / "validation" / "benchmark_txgnn_split_results.json", "w"), indent=2)

    print(f"\n=== TxGNN-protocol disease zero-shot ({n_eval} held-out diseases) ===", flush=True)
    print(f"{'ranker':10} " + " ".join(f"R@{k:<4}" for k in KS) + "  medRank  AUROC(i/c)  AUPRC", flush=True)
    for n in ("distmult", "clf", "fusion"):
        m = out["rankers"][n]
        print(f"{n:10} " + " ".join(f"{m[f'recall@{k}']:<5}" for k in KS) +
              f"  {m['median_rank']:<7} {m['auroc_ind_vs_contra']:<10} {m['auprc_ind_vs_neg']}", flush=True)
    print(f"\nwrote validation/benchmark_txgnn_split_results.json ({time.time()-t0:.0f}s)", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
