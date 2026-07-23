"""
Signature-based connectivity engine  (CMap / LINCS family)
═══════════════════════════════════════════════════════════════════════════════
Backs the "Signature-based" approach in the platform's tagline with a real
connectivity computation instead of a marketing word.

The Connectivity Map idea (Lamb et al., Science 2006; Subramanian et al., Cell
2017 — L1000) is: a drug is a repurposing candidate for a disease when its
perturbation signature is ANTI-correlated with the disease signature — i.e. the
drug pushes genes in the OPPOSITE direction to the disease (it "reverses" the
disease state). A strongly negative connectivity score is the therapeutic hit.

This module implements that scoring as a signed, weighted cosine connectivity
over the genes shared by the two signatures — robust with the sparse signatures
we can build today, and identical in spirit to the L1000 score. The math is the
real algorithm; only the input signatures are pluggable:

  • Disease signature : disease-associated genes (gene_symbol + association
    weight) already resolved from Open Targets in the platform. Treated as
    "elevated in disease" (documented proxy for direction).
  • Drug signature    : the drug's mechanistic targets with action direction —
    INHIBITOR/ANTAGONIST ⇒ −1 (suppresses), AGONIST/ACTIVATOR ⇒ +1. These come
    straight from the mechanisms the platform already has.

Drop-in upgrade path: replace the two builders with real LINCS L1000 z-score
profiles and the exact same connectivity() function scores them — no other code
changes. That is the honest claim: the connectivity engine is real now and
sharpens as richer signatures are connected.
"""
from __future__ import annotations

import math
from typing import Dict, List, Optional

# Map ChEMBL/Open-Targets mechanism action types → signed perturbation direction.
_NEG = {"INHIBITOR", "ANTAGONIST", "BLOCKER", "NEGATIVE ALLOSTERIC MODULATOR",
        "NEGATIVE MODULATOR", "INVERSE AGONIST", "DISRUPTOR", "DEGRADER",
        "SUPPRESSOR", "ANTISENSE INHIBITOR"}
_POS = {"AGONIST", "PARTIAL AGONIST", "ACTIVATOR", "POSITIVE ALLOSTERIC MODULATOR",
        "POSITIVE MODULATOR", "OPENER", "STABILISER", "STABILIZER", "INDUCER"}


def _direction(action: Optional[str]) -> float:
    a = (action or "").strip().upper()
    if a in _NEG:
        return -1.0
    if a in _POS:
        return 1.0
    return 0.0          # MODULATOR / BINDING AGENT / unknown → no asserted direction


def drug_signature(targets: List[Dict]) -> Dict[str, float]:
    """
    Build a drug perturbation signature from mechanistic targets.
    targets: [{"gene": "PDE5A", "action": "INHIBITOR", "weight": 1.0?}, ...]
    Returns {gene: signed_weight}. Genes with no asserted direction are dropped.
    """
    sig: Dict[str, float] = {}
    for t in targets or []:
        gene = (t.get("gene") or t.get("gene_symbol") or "").strip().upper()
        if not gene:
            continue
        d = _direction(t.get("action") or t.get("action_type"))
        if d == 0.0:
            continue
        w = float(t.get("weight", 1.0) or 1.0)
        # if a gene appears twice, keep the strongest asserted direction
        if gene not in sig or abs(d * w) > abs(sig[gene]):
            sig[gene] = d * w
    return sig


def disease_signature(genes: List, default_weight: float = 1.0) -> Dict[str, float]:
    """
    Build a disease signature from associated genes (assumed elevated in disease).
    genes: ["APP", "BACE1", ...]  or  [{"gene": "APP", "weight": 0.8}, ...]
    Returns {gene: +weight}.
    """
    sig: Dict[str, float] = {}
    for g in genes or []:
        if isinstance(g, dict):
            gene = (g.get("gene") or g.get("gene_symbol") or "").strip().upper()
            w = float(g.get("weight", default_weight) or default_weight)
        else:
            gene, w = str(g).strip().upper(), default_weight
        if gene:
            sig[gene] = max(sig.get(gene, 0.0), abs(w))    # disease = elevated (+)
    return sig


def connectivity(disease_sig: Dict[str, float], drug_sig: Dict[str, float]) -> Dict:
    """
    Signed weighted cosine connectivity between a disease and a drug signature.

      connectivity ∈ [-1, 1]   (CMap convention)
        −1  drug perfectly REVERSES the disease signature  → best repurposing hit
         0  orthogonal / no shared signal
        +1  drug MIMICS the disease signature              → contraindicated

      reversal_score ∈ [0, 1] = max(0, −connectivity), the rankable "good" signal.
    """
    shared = [g for g in disease_sig if g in drug_sig]
    if not shared:
        return {"connectivity": 0.0, "reversal_score": 0.0, "n_shared": 0,
                "coverage": 0.0, "signature_score": 0.0, "n_disease_genes": len(disease_sig),
                "n_drug_genes": len(drug_sig), "top_reversed": [], "top_concordant": [],
                "interpretation": "No shared genes between disease and drug signatures."}

    num = sum(disease_sig[g] * drug_sig[g] for g in shared)
    den_d = math.sqrt(sum(disease_sig[g] ** 2 for g in shared))
    den_r = math.sqrt(sum(drug_sig[g] ** 2 for g in shared))
    cosine = num / (den_d * den_r) if den_d and den_r else 0.0
    cosine = max(-1.0, min(1.0, cosine))

    reversal = max(0.0, -cosine)
    coverage = len(shared) / max(1, len(disease_sig))
    # Confidence-weight the reversal by how much of the disease signature we cover
    # (a 1-gene overlap is weaker evidence than a 10-gene overlap). Saturating.
    conf = 1.0 - math.exp(-len(shared) / 3.0)
    signature_score = round(reversal * (0.4 + 0.6 * conf) * (0.6 + 0.4 * coverage), 4)

    contribs = sorted(
        ({"gene": g, "disease": round(disease_sig[g], 3), "drug": round(drug_sig[g], 3),
          "contribution": round(disease_sig[g] * drug_sig[g], 3)} for g in shared),
        key=lambda x: x["contribution"])
    top_reversed = [c for c in contribs if c["contribution"] < 0][:8]
    top_concordant = [c for c in reversed(contribs) if c["contribution"] > 0][:5]

    if cosine <= -0.5:
        interp = "Strong signature reversal — drug opposes the disease signature."
    elif cosine < -0.1:
        interp = "Partial signature reversal — therapeutically favourable direction."
    elif cosine <= 0.1:
        interp = "No net directional signal across shared genes."
    else:
        interp = "Concordant with disease signature — unfavourable (mimics disease)."

    return {
        "connectivity": round(cosine, 4),
        "reversal_score": round(reversal, 4),
        "signature_score": signature_score,
        "n_shared": len(shared),
        "coverage": round(coverage, 3),
        "n_disease_genes": len(disease_sig),
        "n_drug_genes": len(drug_sig),
        "top_reversed": top_reversed,
        "top_concordant": top_concordant,
        "interpretation": interp,
        "method": "Signed weighted-cosine connectivity (CMap/L1000 family)",
        "citations": ["Lamb et al., Science 2006", "Subramanian et al., Cell 2017"],
    }


def score_reversal(disease_genes: List, drug_targets: List[Dict]) -> Dict:
    """Convenience: build both signatures from platform data and score them."""
    return connectivity(disease_signature(disease_genes), drug_signature(drug_targets))


if __name__ == "__main__":
    import json
    # Disease "elevated" genes (e.g. a neuro-inflammatory signature)
    disease = ["TNF", "IL6", "BACE1", "APP", "NFKB1", "GSK3B"]
    print("De-novo connectivity demo (disease = elevated TNF/IL6/BACE1/APP/NFKB1/GSK3B)\n")

    reverser = [  # inhibits disease-elevated genes → should REVERSE
        {"gene": "TNF", "action": "INHIBITOR"},
        {"gene": "GSK3B", "action": "INHIBITOR"},
        {"gene": "BACE1", "action": "INHIBITOR"},
    ]
    mimic = [  # activates disease-elevated genes → should MIMIC (bad)
        {"gene": "TNF", "action": "AGONIST"},
        {"gene": "NFKB1", "action": "ACTIVATOR"},
    ]
    for label, tgts in [("Reverser (anti-inflammatory)", reverser), ("Mimic (pro-inflammatory)", mimic)]:
        r = score_reversal(disease, tgts)
        print(f"• {label}: connectivity={r['connectivity']}  reversal={r['reversal_score']}  "
              f"signature_score={r['signature_score']}  ({r['n_shared']}/{r['n_disease_genes']} genes)")
        print(f"    {r['interpretation']}")
    print("\nFull result for the reverser:")
    print(json.dumps(score_reversal(disease, reverser), indent=2))
