"""
Feature Extractor for ML Scoring
Computes real molecular and interaction features from drugs.json + drug_interactions.json + RDKit.
"""

import json
import logging
import os
from functools import lru_cache
from typing import Dict

logger = logging.getLogger(__name__)

_REPO_ROOT = os.path.dirname(os.path.abspath(__file__))

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, QED, rdMolDescriptors, AllChem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False


@lru_cache(maxsize=1)
def _load_drugs():
    path = os.path.join(_REPO_ROOT, "drugs.json")
    try:
        with open(path, "r", encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        return {}


@lru_cache(maxsize=1)
def _load_interactions():
    path = os.path.join(_REPO_ROOT, "drug_interactions.json")
    try:
        with open(path, "r", encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        return {}


def _rdkit_features(smiles: str) -> Dict:
    """Compute molecular features from SMILES using RDKit."""
    if not RDKIT_AVAILABLE or not smiles:
        return {}
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {}
        return {
            "molecular_weight": Descriptors.ExactMolWt(mol),
            "log_p": Descriptors.MolLogP(mol),
            "qed_score": QED.qed(mol),
            "tpsa": Descriptors.TPSA(mol),
            "rotatable_bonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
            "aromatic_rings": rdMolDescriptors.CalcNumAromaticRings(mol),
            "hba": rdMolDescriptors.CalcNumHBA(mol),
            "hbd": rdMolDescriptors.CalcNumHBD(mol),
            "heavy_atom_count": mol.GetNumHeavyAtoms(),
            "num_rings": rdMolDescriptors.CalcNumRings(mol),
            "num_heteroatoms": rdMolDescriptors.CalcNumHeteroatoms(mol),
            "fraction_sp3": rdMolDescriptors.CalcFractionCSP3(mol),
            "lipinski_violations": sum([
                Descriptors.ExactMolWt(mol) > 500,
                Descriptors.MolLogP(mol) > 5,
                rdMolDescriptors.CalcNumHBD(mol) > 5,
                rdMolDescriptors.CalcNumHBA(mol) > 10,
            ]),
        }
    except Exception as e:
        logger.debug(f"RDKit feature extraction failed: {e}")
        return {}


def _interaction_features(drug_name: str) -> Dict:
    """Compute target interaction features from drug_interactions.json."""
    interactions = _load_interactions()
    key = drug_name.lower()
    drug_data = interactions.get(key) or interactions.get(drug_name, {})

    targets = drug_data.get("targets", [])
    if not targets:
        return {}

    confidences = [t.get("confidence", 0.0) for t in targets if isinstance(t, dict)]
    if not confidences:
        return {}

    return {
        "num_targets": len(targets),
        "avg_confidence": sum(confidences) / len(confidences),
        "best_confidence": max(confidences),
        "worst_confidence": min(confidences),
        "primary_target_confidence": confidences[0] if confidences else 0.5,
        "target_diversity": min(1.0, len(set(t.get("gene_symbol", "") for t in targets if isinstance(t, dict))) / max(len(targets), 1)),
        "off_target_count": max(0, len(targets) - 1),
    }


def extract_all_features(drug_name: str, disease: str = None) -> Dict:
    """
    Extract ML features for a drug. Uses RDKit for molecular features,
    drug_interactions.json for target features, and drugs.json for metadata.
    Falls back to reasonable defaults for anything that can't be computed.
    """
    drugs = _load_drugs()
    drug_key = drug_name.lower()
    drug_data = drugs.get(drug_key) or drugs.get(drug_name, {})

    smiles = drug_data.get("smiles", "")

    # Start with RDKit molecular features
    features = _rdkit_features(smiles)

    # Overlay interaction features
    features.update(_interaction_features(drug_name))

    # Pull any stored properties from JSON
    for key in ["molecular_weight", "logp", "tpsa", "hba", "hbd", "rotatable_bonds"]:
        if key not in features and key in drug_data:
            mapped = {"logp": "log_p"}.get(key, key)
            features[mapped] = float(drug_data[key])

    # FDA approval / phase data
    fda_approved = drug_data.get("fda_approved", False)
    features["fda_approved"] = 1.0 if fda_approved else 0.0

    clinical_phase = drug_data.get("max_phase", 4.0 if fda_approved else 0.0)
    features["clinical_trial_phase"] = float(clinical_phase)
    features["years_approved"] = float(drug_data.get("years_on_market", 5.0 if fda_approved else 0.0))

    # Defaults for features that can't be computed from available data
    defaults = {
        "molecular_weight": 350.0,
        "log_p": 2.5,
        "qed_score": 0.7,
        "tpsa": 70.0,
        "rotatable_bonds": 4,
        "aromatic_rings": 2,
        "hba": 4,
        "hbd": 2,
        "lipinski_violations": 0,
        "fraction_sp3": 0.4,
        "heavy_atom_count": 25,
        "num_rings": 3,
        "num_heteroatoms": 5,
        "pains_alerts": 0,
        "synthetic_accessibility": 3.5,
        "num_targets": 2,
        "avg_confidence": 0.75,
        "best_confidence": 0.85,
        "worst_confidence": 0.6,
        "primary_target_confidence": 0.8,
        "target_diversity": 0.6,
        "off_target_count": 1,
        "avg_binding_affinity": -7.0,
        "best_binding_affinity": -8.0,
        "worst_binding_affinity": -6.0,
        "cyp450_interaction_count": 1,
        "receptor_binding_count": 1,
        "enzyme_interaction_count": 1,
        "pagerank_centrality": 0.05,
        "betweenness_centrality": 0.02,
        "closeness_centrality": 0.5,
        "degree_centrality": 0.1,
        "path_length_to_disease": 3.0,
        "pathway_connection_count": 2,
        "evidence_chain_strength": 0.7,
        "publication_count": 40,
        "network_clustering": 0.5,
        "fda_approved": 0.0,
        "years_approved": 0.0,
        "clinical_trial_phase": 0.0,
        "therapeutic_category_relevance": 0.7,
        "drug_drug_interaction_risk": 0.3,
        "adverse_event_score": 0.8,
        "safety_profile": 0.85,
        "clinical_trial_success_rate": 0.6,
        "market_presence_years": 5.0,
        "indication_similarity": 0.5,
        "disease_category_match": 0.7,
        "mechanism_disease_relevance": 0.7,
        "pathway_overlap_score": 0.5,
        "disease_gene_target_overlap": 0.5,
        "unmet_medical_need": 0.7,
        "competitive_landscape": 0.5,
        "repurposing_feasibility": 0.7,
    }

    for k, v in defaults.items():
        if k not in features:
            features[k] = v

    return features
