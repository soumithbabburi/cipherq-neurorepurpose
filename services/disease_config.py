"""
Disease-Specific Configuration Module
Defines optimization parameters, PBPK parameters, and target selection for each disease category
"""

import logging
from typing import Dict, List, Optional, Any
from dataclasses import dataclass

logger = logging.getLogger(__name__)

@dataclass
class DiseaseOptimizationProfile:
    """Disease-specific optimization profile"""
    name: str
    category: str
    optimization_weights: Dict[str, float]
    property_thresholds: Dict[str, Dict]
    target_tissues: List[str]
    preferred_targets: List[str]
    pbpk_adjustments: Dict[str, float]
    scoring_emphasis: str

DISEASE_PROFILES = {
    "Alzheimer's Disease": DiseaseOptimizationProfile(
        name="Alzheimer's Disease",
        category="neurological",
        optimization_weights={
            'bbb_penetration': 0.40,
            'drug_likeness': 0.20,
            'stability': 0.20,
            'selectivity': 0.20
        },
        property_thresholds={
            'molecular_weight': {'min': 150, 'max': 450, 'optimal': (200, 400)},
            'logp': {'min': 1.0, 'max': 4.0, 'optimal': (2.0, 3.5)},
            'tpsa': {'min': 20, 'max': 90, 'optimal': (40, 70)},
            'hbd': {'max': 3, 'optimal': (1, 2)},
            'hba': {'max': 7, 'optimal': (2, 5)},
            'rotatable_bonds': {'max': 8, 'optimal': (2, 6)}
        },
        target_tissues=['brain', 'cns'],
        preferred_targets=['AChE', 'BuChE', 'BACE1', 'NMDA', 'APP', 'MAPT', 'GSK3B', 'CDK5'],
        pbpk_adjustments={'kp_brain': 1.5, 'bbb_factor': 1.3},
        scoring_emphasis='cns_mpo'
    ),
    "Parkinson's Disease": DiseaseOptimizationProfile(
        name="Parkinson's Disease",
        category="neurological",
        optimization_weights={
            'bbb_penetration': 0.40,
            'drug_likeness': 0.20,
            'stability': 0.20,
            'selectivity': 0.20
        },
        property_thresholds={
            'molecular_weight': {'min': 150, 'max': 450, 'optimal': (200, 400)},
            'logp': {'min': 1.0, 'max': 4.0, 'optimal': (2.0, 3.5)},
            'tpsa': {'min': 20, 'max': 90, 'optimal': (40, 70)},
            'hbd': {'max': 3, 'optimal': (1, 2)},
            'hba': {'max': 7, 'optimal': (2, 5)},
            'rotatable_bonds': {'max': 8, 'optimal': (2, 6)}
        },
        target_tissues=['brain', 'substantia_nigra'],
        preferred_targets=['DDC', 'COMT', 'MAO-B', 'DRD2', 'DRD3', 'SNCA', 'LRRK2', 'PARK7'],
        pbpk_adjustments={'kp_brain': 1.5, 'bbb_factor': 1.3},
        scoring_emphasis='cns_mpo'
    ),
    "Multiple Sclerosis": DiseaseOptimizationProfile(
        name="Multiple Sclerosis",
        category="neurological",
        optimization_weights={
            'bbb_penetration': 0.35,
            'drug_likeness': 0.25,
            'stability': 0.20,
            'selectivity': 0.20
        },
        property_thresholds={
            'molecular_weight': {'min': 200, 'max': 500, 'optimal': (250, 450)},
            'logp': {'min': 0.5, 'max': 4.0, 'optimal': (1.5, 3.5)},
            'tpsa': {'min': 30, 'max': 100, 'optimal': (50, 80)},
            'hbd': {'max': 4, 'optimal': (1, 3)},
            'hba': {'max': 8, 'optimal': (3, 6)},
            'rotatable_bonds': {'max': 10, 'optimal': (3, 7)}
        },
        target_tissues=['brain', 'cns', 'immune'],
        preferred_targets=['S1PR1', 'CD20', 'ITK', 'BTK', 'VLA4', 'IL2RA'],
        pbpk_adjustments={'kp_brain': 1.3, 'bbb_factor': 1.2},
        scoring_emphasis='cns_mpo'
    ),
    "Type 2 Diabetes": DiseaseOptimizationProfile(
        name="Type 2 Diabetes",
        category="metabolic",
        optimization_weights={
            'bbb_penetration': 0.05,
            'drug_likeness': 0.35,
            'stability': 0.30,
            'selectivity': 0.30
        },
        property_thresholds={
            'molecular_weight': {'min': 200, 'max': 600, 'optimal': (300, 500)},
            'logp': {'min': -1.0, 'max': 5.0, 'optimal': (1.0, 4.0)},
            'tpsa': {'min': 40, 'max': 140, 'optimal': (60, 120)},
            'hbd': {'max': 5, 'optimal': (2, 4)},
            'hba': {'max': 10, 'optimal': (4, 8)},
            'rotatable_bonds': {'max': 12, 'optimal': (4, 8)}
        },
        target_tissues=['liver', 'pancreas', 'muscle', 'adipose'],
        preferred_targets=['INSR', 'GLP1R', 'DPP4', 'SGLT2', 'PPARG', 'AMPK', 'GCK', 'GCGR'],
        pbpk_adjustments={'kp_liver': 1.4, 'kp_pancreas': 1.3},
        scoring_emphasis='metabolic_stability'
    ),
    "Type 1 Diabetes": DiseaseOptimizationProfile(
        name="Type 1 Diabetes",
        category="metabolic",
        optimization_weights={
            'bbb_penetration': 0.05,
            'drug_likeness': 0.35,
            'stability': 0.30,
            'selectivity': 0.30
        },
        property_thresholds={
            'molecular_weight': {'min': 200, 'max': 600, 'optimal': (300, 500)},
            'logp': {'min': -1.0, 'max': 5.0, 'optimal': (1.0, 4.0)},
            'tpsa': {'min': 40, 'max': 140, 'optimal': (60, 120)},
            'hbd': {'max': 5, 'optimal': (2, 4)},
            'hba': {'max': 10, 'optimal': (4, 8)},
            'rotatable_bonds': {'max': 12, 'optimal': (4, 8)}
        },
        target_tissues=['pancreas', 'immune'],
        preferred_targets=['INSR', 'CD3', 'GAD65', 'IA2', 'ZNT8'],
        pbpk_adjustments={'kp_pancreas': 1.5},
        scoring_emphasis='metabolic_stability'
    ),
    "Hypertension": DiseaseOptimizationProfile(
        name="Hypertension",
        category="cardiovascular",
        optimization_weights={
            'bbb_penetration': 0.10,
            'drug_likeness': 0.35,
            'stability': 0.30,
            'selectivity': 0.25
        },
        property_thresholds={
            'molecular_weight': {'min': 200, 'max': 550, 'optimal': (250, 450)},
            'logp': {'min': 0.5, 'max': 5.0, 'optimal': (1.5, 4.0)},
            'tpsa': {'min': 40, 'max': 130, 'optimal': (60, 100)},
            'hbd': {'max': 4, 'optimal': (1, 3)},
            'hba': {'max': 9, 'optimal': (3, 7)},
            'rotatable_bonds': {'max': 10, 'optimal': (3, 8)}
        },
        target_tissues=['heart', 'kidney', 'vasculature'],
        preferred_targets=['ACE', 'ACE2', 'AT1R', 'ADRB1', 'ADRB2', 'CACNA1C', 'KCNH2', 'NOS3'],
        pbpk_adjustments={'kp_heart': 1.3, 'kp_kidney': 1.2},
        scoring_emphasis='cardiovascular_safety'
    ),
    "Heart Failure": DiseaseOptimizationProfile(
        name="Heart Failure",
        category="cardiovascular",
        optimization_weights={
            'bbb_penetration': 0.05,
            'drug_likeness': 0.35,
            'stability': 0.35,
            'selectivity': 0.25
        },
        property_thresholds={
            'molecular_weight': {'min': 200, 'max': 600, 'optimal': (300, 500)},
            'logp': {'min': 0.0, 'max': 5.0, 'optimal': (1.0, 4.0)},
            'tpsa': {'min': 50, 'max': 140, 'optimal': (70, 120)},
            'hbd': {'max': 5, 'optimal': (2, 4)},
            'hba': {'max': 10, 'optimal': (4, 8)},
            'rotatable_bonds': {'max': 12, 'optimal': (4, 9)}
        },
        target_tissues=['heart', 'kidney'],
        preferred_targets=['ACE', 'AT1R', 'ADRB1', 'ADRB2', 'NEP', 'SGLT2', 'MR', 'PDE3'],
        pbpk_adjustments={'kp_heart': 1.5, 'kp_kidney': 1.3},
        scoring_emphasis='cardiovascular_safety'
    ),
    "Lung Cancer": DiseaseOptimizationProfile(
        name="Lung Cancer",
        category="oncology",
        optimization_weights={
            'bbb_penetration': 0.05,
            'drug_likeness': 0.25,
            'stability': 0.30,
            'selectivity': 0.40
        },
        property_thresholds={
            'molecular_weight': {'min': 300, 'max': 700, 'optimal': (400, 600)},
            'logp': {'min': 1.0, 'max': 6.0, 'optimal': (2.0, 5.0)},
            'tpsa': {'min': 50, 'max': 150, 'optimal': (80, 130)},
            'hbd': {'max': 5, 'optimal': (2, 4)},
            'hba': {'max': 12, 'optimal': (5, 10)},
            'rotatable_bonds': {'max': 15, 'optimal': (5, 12)}
        },
        target_tissues=['lung', 'tumor'],
        preferred_targets=['EGFR', 'ALK', 'ROS1', 'KRAS', 'MET', 'BRAF', 'RET', 'NTRK', 'HER2'],
        pbpk_adjustments={'kp_lung': 1.5, 'tumor_penetration': 1.4},
        scoring_emphasis='tumor_selectivity'
    ),
    "Breast Cancer": DiseaseOptimizationProfile(
        name="Breast Cancer",
        category="oncology",
        optimization_weights={
            'bbb_penetration': 0.05,
            'drug_likeness': 0.25,
            'stability': 0.30,
            'selectivity': 0.40
        },
        property_thresholds={
            'molecular_weight': {'min': 300, 'max': 700, 'optimal': (400, 600)},
            'logp': {'min': 1.0, 'max': 6.0, 'optimal': (2.0, 5.0)},
            'tpsa': {'min': 50, 'max': 150, 'optimal': (80, 130)},
            'hbd': {'max': 5, 'optimal': (2, 4)},
            'hba': {'max': 12, 'optimal': (5, 10)},
            'rotatable_bonds': {'max': 15, 'optimal': (5, 12)}
        },
        target_tissues=['breast', 'tumor'],
        preferred_targets=['ESR1', 'HER2', 'ERBB2', 'CDK4', 'CDK6', 'PIK3CA', 'BRCA1', 'BRCA2', 'PARP'],
        pbpk_adjustments={'tumor_penetration': 1.4},
        scoring_emphasis='tumor_selectivity'
    ),
    "Rheumatoid Arthritis": DiseaseOptimizationProfile(
        name="Rheumatoid Arthritis",
        category="inflammatory",
        optimization_weights={
            'bbb_penetration': 0.05,
            'drug_likeness': 0.30,
            'stability': 0.30,
            'selectivity': 0.35
        },
        property_thresholds={
            'molecular_weight': {'min': 250, 'max': 650, 'optimal': (350, 550)},
            'logp': {'min': 0.5, 'max': 5.5, 'optimal': (1.5, 4.5)},
            'tpsa': {'min': 50, 'max': 140, 'optimal': (70, 120)},
            'hbd': {'max': 5, 'optimal': (2, 4)},
            'hba': {'max': 10, 'optimal': (4, 8)},
            'rotatable_bonds': {'max': 12, 'optimal': (4, 10)}
        },
        target_tissues=['joints', 'immune', 'synovium'],
        preferred_targets=['TNF', 'IL6R', 'JAK1', 'JAK3', 'BTK', 'SYK', 'CD20', 'CTLA4'],
        pbpk_adjustments={'kp_joint': 1.4, 'immune_distribution': 1.3},
        scoring_emphasis='anti_inflammatory'
    ),
    "Crohn's Disease": DiseaseOptimizationProfile(
        name="Crohn's Disease",
        category="inflammatory",
        optimization_weights={
            'bbb_penetration': 0.05,
            'drug_likeness': 0.30,
            'stability': 0.35,
            'selectivity': 0.30
        },
        property_thresholds={
            'molecular_weight': {'min': 250, 'max': 650, 'optimal': (350, 550)},
            'logp': {'min': 0.0, 'max': 5.0, 'optimal': (1.0, 4.0)},
            'tpsa': {'min': 60, 'max': 150, 'optimal': (80, 130)},
            'hbd': {'max': 5, 'optimal': (2, 4)},
            'hba': {'max': 11, 'optimal': (5, 9)},
            'rotatable_bonds': {'max': 12, 'optimal': (4, 10)}
        },
        target_tissues=['gut', 'colon', 'immune'],
        preferred_targets=['TNF', 'IL12B', 'IL23A', 'ITGA4', 'JAK1', 'S1PR1'],
        pbpk_adjustments={'kp_gut': 1.5, 'gut_absorption': 0.9},
        scoring_emphasis='gut_targeting'
    ),
    "HIV/AIDS": DiseaseOptimizationProfile(
        name="HIV/AIDS",
        category="infectious",
        optimization_weights={
            'bbb_penetration': 0.15,
            'drug_likeness': 0.30,
            'stability': 0.30,
            'selectivity': 0.25
        },
        property_thresholds={
            'molecular_weight': {'min': 300, 'max': 700, 'optimal': (400, 600)},
            'logp': {'min': 0.5, 'max': 5.5, 'optimal': (2.0, 4.5)},
            'tpsa': {'min': 50, 'max': 140, 'optimal': (70, 120)},
            'hbd': {'max': 5, 'optimal': (2, 4)},
            'hba': {'max': 12, 'optimal': (5, 10)},
            'rotatable_bonds': {'max': 15, 'optimal': (5, 12)}
        },
        target_tissues=['immune', 'lymphoid', 'cns'],
        preferred_targets=['HIV-RT', 'HIV-PR', 'HIV-IN', 'CCR5', 'CXCR4', 'CD4'],
        pbpk_adjustments={'immune_distribution': 1.4, 'cns_penetration': 1.2},
        scoring_emphasis='viral_inhibition'
    ),
    "COVID-19": DiseaseOptimizationProfile(
        name="COVID-19",
        category="infectious",
        optimization_weights={
            'bbb_penetration': 0.05,
            'drug_likeness': 0.35,
            'stability': 0.30,
            'selectivity': 0.30
        },
        property_thresholds={
            'molecular_weight': {'min': 250, 'max': 650, 'optimal': (350, 550)},
            'logp': {'min': 0.0, 'max': 5.0, 'optimal': (1.5, 4.0)},
            'tpsa': {'min': 50, 'max': 140, 'optimal': (70, 120)},
            'hbd': {'max': 5, 'optimal': (2, 4)},
            'hba': {'max': 11, 'optimal': (5, 9)},
            'rotatable_bonds': {'max': 12, 'optimal': (4, 10)}
        },
        target_tissues=['lung', 'respiratory'],
        preferred_targets=['3CLpro', 'PLpro', 'RdRp', 'ACE2', 'TMPRSS2', 'Spike'],
        pbpk_adjustments={'kp_lung': 1.5, 'respiratory_distribution': 1.4},
        scoring_emphasis='antiviral_potency'
    ),
    "Epilepsy": DiseaseOptimizationProfile(
        name="Epilepsy",
        category="neurological",
        optimization_weights={
            'bbb_penetration': 0.40,
            'drug_likeness': 0.25,
            'stability': 0.20,
            'selectivity': 0.15
        },
        property_thresholds={
            'molecular_weight': {'min': 150, 'max': 400, 'optimal': (200, 350)},
            'logp': {'min': 1.0, 'max': 4.0, 'optimal': (2.0, 3.5)},
            'tpsa': {'min': 20, 'max': 80, 'optimal': (35, 65)},
            'hbd': {'max': 2, 'optimal': (0, 2)},
            'hba': {'max': 6, 'optimal': (2, 4)},
            'rotatable_bonds': {'max': 6, 'optimal': (1, 4)}
        },
        target_tissues=['brain', 'cns'],
        preferred_targets=['GABA-A', 'SCN1A', 'SCN2A', 'CACNA1H', 'SV2A', 'AMPA', 'NMDA'],
        pbpk_adjustments={'kp_brain': 1.5, 'bbb_factor': 1.4},
        scoring_emphasis='cns_mpo'
    ),
    "Asthma": DiseaseOptimizationProfile(
        name="Asthma",
        category="respiratory",
        optimization_weights={
            'bbb_penetration': 0.05,
            'drug_likeness': 0.30,
            'stability': 0.35,
            'selectivity': 0.30
        },
        property_thresholds={
            'molecular_weight': {'min': 200, 'max': 550, 'optimal': (300, 450)},
            'logp': {'min': 0.5, 'max': 5.0, 'optimal': (2.0, 4.0)},
            'tpsa': {'min': 40, 'max': 120, 'optimal': (60, 100)},
            'hbd': {'max': 4, 'optimal': (1, 3)},
            'hba': {'max': 9, 'optimal': (4, 7)},
            'rotatable_bonds': {'max': 10, 'optimal': (3, 8)}
        },
        target_tissues=['lung', 'bronchi', 'airway'],
        preferred_targets=['ADRB2', 'GCR', 'ALOX5', 'LTC4S', 'PDE4', 'IL5', 'IL4RA', 'IgE'],
        pbpk_adjustments={'kp_lung': 1.5, 'inhalation_factor': 1.3},
        scoring_emphasis='respiratory_targeting'
    ),
    "COPD": DiseaseOptimizationProfile(
        name="COPD",
        category="respiratory",
        optimization_weights={
            'bbb_penetration': 0.05,
            'drug_likeness': 0.30,
            'stability': 0.35,
            'selectivity': 0.30
        },
        property_thresholds={
            'molecular_weight': {'min': 200, 'max': 550, 'optimal': (300, 450)},
            'logp': {'min': 0.5, 'max': 5.0, 'optimal': (2.0, 4.0)},
            'tpsa': {'min': 40, 'max': 120, 'optimal': (60, 100)},
            'hbd': {'max': 4, 'optimal': (1, 3)},
            'hba': {'max': 9, 'optimal': (4, 7)},
            'rotatable_bonds': {'max': 10, 'optimal': (3, 8)}
        },
        target_tissues=['lung', 'bronchi'],
        preferred_targets=['ADRB2', 'CHRM3', 'GCR', 'PDE4', 'NE', 'MMP9', 'MMP12'],
        pbpk_adjustments={'kp_lung': 1.5, 'inhalation_factor': 1.3},
        scoring_emphasis='respiratory_targeting'
    ),
    "Major Depression": DiseaseOptimizationProfile(
        name="Major Depression",
        category="psychiatric",
        optimization_weights={
            'bbb_penetration': 0.40,
            'drug_likeness': 0.25,
            'stability': 0.20,
            'selectivity': 0.15
        },
        property_thresholds={
            'molecular_weight': {'min': 150, 'max': 450, 'optimal': (200, 380)},
            'logp': {'min': 1.5, 'max': 4.5, 'optimal': (2.5, 4.0)},
            'tpsa': {'min': 20, 'max': 80, 'optimal': (30, 65)},
            'hbd': {'max': 2, 'optimal': (0, 2)},
            'hba': {'max': 6, 'optimal': (2, 5)},
            'rotatable_bonds': {'max': 7, 'optimal': (2, 5)}
        },
        target_tissues=['brain', 'cns'],
        preferred_targets=['SERT', 'NET', 'DAT', '5HT1A', '5HT2A', 'MAOA', 'MAOB', 'NK1R'],
        pbpk_adjustments={'kp_brain': 1.5, 'bbb_factor': 1.4},
        scoring_emphasis='cns_mpo'
    ),
    "Schizophrenia": DiseaseOptimizationProfile(
        name="Schizophrenia",
        category="psychiatric",
        optimization_weights={
            'bbb_penetration': 0.40,
            'drug_likeness': 0.25,
            'stability': 0.20,
            'selectivity': 0.15
        },
        property_thresholds={
            'molecular_weight': {'min': 200, 'max': 500, 'optimal': (250, 420)},
            'logp': {'min': 2.0, 'max': 5.0, 'optimal': (3.0, 4.5)},
            'tpsa': {'min': 20, 'max': 70, 'optimal': (25, 55)},
            'hbd': {'max': 2, 'optimal': (0, 1)},
            'hba': {'max': 5, 'optimal': (2, 4)},
            'rotatable_bonds': {'max': 6, 'optimal': (2, 5)}
        },
        target_tissues=['brain', 'cns'],
        preferred_targets=['DRD2', 'DRD3', '5HT2A', '5HT1A', 'H1', 'M1', 'ADRA1A'],
        pbpk_adjustments={'kp_brain': 1.5, 'bbb_factor': 1.4},
        scoring_emphasis='cns_mpo'
    ),
    "Chronic Kidney Disease": DiseaseOptimizationProfile(
        name="Chronic Kidney Disease",
        category="renal",
        optimization_weights={
            'bbb_penetration': 0.05,
            'drug_likeness': 0.35,
            'stability': 0.35,
            'selectivity': 0.25
        },
        property_thresholds={
            'molecular_weight': {'min': 200, 'max': 550, 'optimal': (250, 450)},
            'logp': {'min': -0.5, 'max': 4.0, 'optimal': (0.5, 3.0)},
            'tpsa': {'min': 50, 'max': 140, 'optimal': (70, 120)},
            'hbd': {'max': 5, 'optimal': (2, 4)},
            'hba': {'max': 10, 'optimal': (4, 8)},
            'rotatable_bonds': {'max': 10, 'optimal': (3, 8)}
        },
        target_tissues=['kidney'],
        preferred_targets=['RAAS', 'ACE', 'AT1R', 'SGLT2', 'MR', 'EGFR', 'VEGF'],
        pbpk_adjustments={'kp_kidney': 1.5, 'renal_clearance': 0.7},
        scoring_emphasis='renal_safety'
    ),
    "Osteoporosis": DiseaseOptimizationProfile(
        name="Osteoporosis",
        category="musculoskeletal",
        optimization_weights={
            'bbb_penetration': 0.05,
            'drug_likeness': 0.30,
            'stability': 0.35,
            'selectivity': 0.30
        },
        property_thresholds={
            'molecular_weight': {'min': 200, 'max': 600, 'optimal': (300, 500)},
            'logp': {'min': -1.0, 'max': 4.0, 'optimal': (0.0, 3.0)},
            'tpsa': {'min': 60, 'max': 160, 'optimal': (80, 140)},
            'hbd': {'max': 6, 'optimal': (2, 5)},
            'hba': {'max': 12, 'optimal': (5, 10)},
            'rotatable_bonds': {'max': 12, 'optimal': (4, 9)}
        },
        target_tissues=['bone'],
        preferred_targets=['RANKL', 'RANK', 'CTSK', 'SOST', 'PTH1R', 'VDR', 'ER'],
        pbpk_adjustments={'kp_bone': 1.5},
        scoring_emphasis='bone_targeting'
    ),
}

def get_disease_category(disease_name: str) -> str:
    """Map disease name to category for fallback handling"""
    disease_lower = disease_name.lower()
    
    neurological_keywords = ['alzheimer', 'parkinson', 'huntington', 'als', 'multiple sclerosis', 
                            'dementia', 'epilepsy', 'seizure', 'neuropathy', 'stroke']
    cardiovascular_keywords = ['heart', 'cardiac', 'hypertension', 'arrhythmia', 'atherosclerosis',
                              'coronary', 'myocardial', 'stroke', 'atrial']
    oncology_keywords = ['cancer', 'tumor', 'carcinoma', 'leukemia', 'lymphoma', 'melanoma',
                        'sarcoma', 'neoplasm', 'oncology']
    metabolic_keywords = ['diabetes', 'obesity', 'metabolic', 'thyroid', 'cholesterol']
    inflammatory_keywords = ['arthritis', 'crohn', 'colitis', 'lupus', 'inflammatory',
                            'psoriasis', 'autoimmune']
    infectious_keywords = ['hiv', 'aids', 'hepatitis', 'covid', 'tuberculosis', 'malaria',
                          'influenza', 'viral', 'bacterial', 'infection']
    respiratory_keywords = ['asthma', 'copd', 'pulmonary', 'lung', 'respiratory', 'bronchitis']
    psychiatric_keywords = ['depression', 'anxiety', 'schizophrenia', 'bipolar', 'psychiatric',
                           'ptsd', 'adhd', 'autism']
    renal_keywords = ['kidney', 'renal', 'nephro']
    
    if any(kw in disease_lower for kw in neurological_keywords):
        return 'neurological'
    elif any(kw in disease_lower for kw in cardiovascular_keywords):
        return 'cardiovascular'
    elif any(kw in disease_lower for kw in oncology_keywords):
        return 'oncology'
    elif any(kw in disease_lower for kw in metabolic_keywords):
        return 'metabolic'
    elif any(kw in disease_lower for kw in inflammatory_keywords):
        return 'inflammatory'
    elif any(kw in disease_lower for kw in infectious_keywords):
        return 'infectious'
    elif any(kw in disease_lower for kw in respiratory_keywords):
        return 'respiratory'
    elif any(kw in disease_lower for kw in psychiatric_keywords):
        return 'psychiatric'
    elif any(kw in disease_lower for kw in renal_keywords):
        return 'renal'
    else:
        return 'general'

def get_default_profile_for_category(category: str) -> DiseaseOptimizationProfile:
    """Get a default profile for a disease category"""
    category_defaults = {
        'neurological': DISEASE_PROFILES.get("Alzheimer's Disease"),
        'cardiovascular': DISEASE_PROFILES.get("Hypertension"),
        'oncology': DISEASE_PROFILES.get("Lung Cancer"),
        'metabolic': DISEASE_PROFILES.get("Type 2 Diabetes"),
        'inflammatory': DISEASE_PROFILES.get("Rheumatoid Arthritis"),
        'infectious': DISEASE_PROFILES.get("HIV/AIDS"),
        'respiratory': DISEASE_PROFILES.get("Asthma"),
        'psychiatric': DISEASE_PROFILES.get("Major Depression"),
        'renal': DISEASE_PROFILES.get("Chronic Kidney Disease"),
    }
    
    profile = category_defaults.get(category)
    if profile:
        return profile
    
    return DiseaseOptimizationProfile(
        name="General",
        category="general",
        optimization_weights={
            'bbb_penetration': 0.15,
            'drug_likeness': 0.35,
            'stability': 0.25,
            'selectivity': 0.25
        },
        property_thresholds={
            'molecular_weight': {'min': 200, 'max': 550, 'optimal': (300, 450)},
            'logp': {'min': 0.0, 'max': 5.0, 'optimal': (1.5, 4.0)},
            'tpsa': {'min': 40, 'max': 130, 'optimal': (60, 100)},
            'hbd': {'max': 5, 'optimal': (2, 4)},
            'hba': {'max': 10, 'optimal': (4, 8)},
            'rotatable_bonds': {'max': 10, 'optimal': (3, 8)}
        },
        target_tissues=['general'],
        preferred_targets=[],
        pbpk_adjustments={},
        scoring_emphasis='drug_likeness'
    )

def get_disease_profile(disease_name: str) -> DiseaseOptimizationProfile:
    """Get disease-specific optimization profile"""
    if disease_name in DISEASE_PROFILES:
        logger.info(f"Using specific profile for {disease_name}")
        return DISEASE_PROFILES[disease_name]
    
    category = get_disease_category(disease_name)
    logger.info(f"No specific profile for {disease_name}, using {category} category defaults")
    return get_default_profile_for_category(category)

def get_optimization_weights(disease_name: str) -> Dict[str, float]:
    """Get optimization weights for a disease"""
    profile = get_disease_profile(disease_name)
    return profile.optimization_weights

def get_property_thresholds(disease_name: str) -> Dict[str, Dict]:
    """Get property thresholds for a disease"""
    profile = get_disease_profile(disease_name)
    return profile.property_thresholds

def get_target_tissues(disease_name: str) -> List[str]:
    """Get target tissues for a disease"""
    profile = get_disease_profile(disease_name)
    return profile.target_tissues

def get_preferred_targets(disease_name: str) -> List[str]:
    """Get preferred protein targets for a disease"""
    profile = get_disease_profile(disease_name)
    return profile.preferred_targets

def get_pbpk_adjustments(disease_name: str) -> Dict[str, float]:
    """Get PBPK parameter adjustments for a disease"""
    profile = get_disease_profile(disease_name)
    return profile.pbpk_adjustments

def get_scoring_emphasis(disease_name: str) -> str:
    """Get scoring emphasis for a disease"""
    profile = get_disease_profile(disease_name)
    return profile.scoring_emphasis
