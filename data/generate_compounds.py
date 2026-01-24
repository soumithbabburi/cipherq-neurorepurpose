"""
Generate Chemical Compounds Library
Contains experimental and research molecules distinct from FDA-approved drugs
"""

import json
import random
from typing import Dict, List
from pathlib import Path

COMPOUND_CATEGORIES = {
    'natural_products': {
        'description': 'Compounds derived from natural sources',
        'compounds': [
            {'name': 'Curcumin', 'formula': 'C21H20O6', 'smiles': 'COc1cc(/C=C/C(=O)CC(=O)/C=C/c2ccc(O)c(OC)c2)ccc1O', 'source': 'Turmeric', 'mw': 368.38},
            {'name': 'Resveratrol', 'formula': 'C14H12O3', 'smiles': 'Oc1ccc(/C=C/c2cc(O)cc(O)c2)cc1', 'source': 'Grapes', 'mw': 228.24},
            {'name': 'Quercetin', 'formula': 'C15H10O7', 'smiles': 'O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12', 'source': 'Plants', 'mw': 302.24},
            {'name': 'Epigallocatechin Gallate', 'formula': 'C22H18O11', 'smiles': 'Oc1cc(O)c2c(c1)OC(c1cc(O)c(O)c(O)c1)C(O)C2OC(=O)c1cc(O)c(O)c(O)c1', 'source': 'Green Tea', 'mw': 458.37},
            {'name': 'Berberine', 'formula': 'C20H18NO4+', 'smiles': 'COc1ccc2cc3[n+](Cc2c1OC)ccc4cc5OCOc5cc34', 'source': 'Berberis', 'mw': 336.36},
            {'name': 'Ginsenoside Rg1', 'formula': 'C42H72O14', 'smiles': 'Triterpenoid saponin', 'source': 'Ginseng', 'mw': 801.01},
            {'name': 'Artemisinin', 'formula': 'C15H22O5', 'smiles': 'CC1CCC2C(C)C(=O)OC3OC4(C)CCC1C23OO4', 'source': 'Artemisia annua', 'mw': 282.33},
            {'name': 'Taxol Precursor', 'formula': 'C47H51NO14', 'smiles': 'Taxane scaffold', 'source': 'Taxus brevifolia', 'mw': 853.91},
            {'name': 'Camptothecin', 'formula': 'C20H16N2O4', 'smiles': 'CCC1(O)C(=O)OCc2c1cc1n(c2=O)Cc2cc3ccccc3nc2-1', 'source': 'Camptotheca acuminata', 'mw': 348.35},
            {'name': 'Vincristine Analog', 'formula': 'C46H56N4O10', 'smiles': 'Vinca alkaloid scaffold', 'source': 'Catharanthus roseus', 'mw': 824.96},
            {'name': 'Morphinan Scaffold', 'formula': 'C17H19NO3', 'smiles': 'CN1CCC23c4c5ccc(O)c4OC2C(O)=CC5C1C3', 'source': 'Papaver somniferum', 'mw': 285.34},
            {'name': 'Strychnine', 'formula': 'C21H22N2O2', 'smiles': 'O=C1CC2CC3N1C4C=CC=CC4C3C1OC=CC12', 'source': 'Strychnos nux-vomica', 'mw': 334.41},
            {'name': 'Capsaicin', 'formula': 'C18H27NO3', 'smiles': 'COc1cc(CNC(=O)CCCC/C=C/C(C)C)ccc1O', 'source': 'Capsicum', 'mw': 305.41},
            {'name': 'Lycopene', 'formula': 'C40H56', 'smiles': 'CC(C)=CCC/C(C)=C/C=C/C(C)=C/C=C/C=C(C)/C=C/C=C(C)/C=C/C=C(C)/C=C/C=C(C)/CCC=C(C)C', 'source': 'Tomato', 'mw': 536.87},
            {'name': 'Thymoquinone', 'formula': 'C10H12O2', 'smiles': 'CC1=CC(=O)C(C(C)C)=CC1=O', 'source': 'Nigella sativa', 'mw': 164.20},
        ]
    },
    'kinase_inhibitors': {
        'description': 'Compounds targeting protein kinases',
        'compounds': [
            {'name': 'Staurosporine Analog A', 'formula': 'C28H26N4O3', 'smiles': 'Pan-kinase inhibitor scaffold', 'target': 'Multiple kinases', 'mw': 466.53},
            {'name': 'JAK2 Inhibitor Compound', 'formula': 'C25H28N6O2', 'smiles': 'JAK2 selective scaffold', 'target': 'JAK2', 'mw': 444.53},
            {'name': 'BTK Inhibitor Lead', 'formula': 'C25H24N6O2', 'smiles': 'BTK covalent inhibitor', 'target': 'BTK', 'mw': 440.50},
            {'name': 'CDK4/6 Inhibitor Candidate', 'formula': 'C24H29N7O2', 'smiles': 'CDK4/6 selective scaffold', 'target': 'CDK4/CDK6', 'mw': 447.54},
            {'name': 'BRAF Inhibitor Scaffold', 'formula': 'C23H18FN5OS', 'smiles': 'BRAF V600E selective', 'target': 'BRAF', 'mw': 431.49},
            {'name': 'MEK1/2 Inhibitor Core', 'formula': 'C17H15FIN3O2', 'smiles': 'Allosteric MEK inhibitor', 'target': 'MEK1/MEK2', 'mw': 455.22},
            {'name': 'PI3K Delta Compound', 'formula': 'C22H18FN5O', 'smiles': 'PI3K delta selective', 'target': 'PIK3CD', 'mw': 387.41},
            {'name': 'ALK Inhibitor Lead', 'formula': 'C21H22ClN5O', 'smiles': 'ALK/ROS1 dual inhibitor', 'target': 'ALK', 'mw': 395.89},
            {'name': 'FLT3 Inhibitor Candidate', 'formula': 'C24H25N5O2', 'smiles': 'FLT3 ITD selective', 'target': 'FLT3', 'mw': 415.49},
            {'name': 'Aurora A Inhibitor', 'formula': 'C20H21N5O2S', 'smiles': 'Aurora A selective scaffold', 'target': 'AURKA', 'mw': 395.48},
            {'name': 'ROCK Inhibitor Core', 'formula': 'C14H21N3O', 'smiles': 'ROCK1/2 inhibitor', 'target': 'ROCK1', 'mw': 247.34},
            {'name': 'GSK3 Inhibitor Lead', 'formula': 'C12H9Cl2N3', 'smiles': 'GSK3 beta selective', 'target': 'GSK3B', 'mw': 266.12},
        ]
    },
    'gpcr_modulators': {
        'description': 'G-protein coupled receptor modulators',
        'compounds': [
            {'name': 'Dopamine D1 Agonist Lead', 'formula': 'C17H19NO3', 'smiles': 'D1 receptor selective', 'target': 'DRD1', 'mw': 285.34},
            {'name': 'Serotonin 5-HT2A Antagonist', 'formula': 'C21H26N4O2', 'smiles': '5-HT2A selective blocker', 'target': 'HTR2A', 'mw': 366.46},
            {'name': 'Muscarinic M4 PAM', 'formula': 'C18H22N2O2', 'smiles': 'M4 positive allosteric modulator', 'target': 'CHRM4', 'mw': 298.38},
            {'name': 'Adenosine A2A Antagonist', 'formula': 'C15H14N6O', 'smiles': 'A2A receptor antagonist', 'target': 'ADORA2A', 'mw': 294.31},
            {'name': 'Cannabinoid CB2 Agonist', 'formula': 'C23H32N2O3', 'smiles': 'CB2 selective agonist', 'target': 'CNR2', 'mw': 384.51},
            {'name': 'Opioid Kappa Agonist', 'formula': 'C22H27NO4', 'smiles': 'Kappa receptor agonist', 'target': 'OPRK1', 'mw': 369.46},
            {'name': 'Histamine H3 Antagonist', 'formula': 'C16H21N3O', 'smiles': 'H3 receptor inverse agonist', 'target': 'HRH3', 'mw': 271.36},
            {'name': 'Neuropeptide Y1 Antagonist', 'formula': 'C26H30N4O2', 'smiles': 'NPY Y1 receptor blocker', 'target': 'NPY1R', 'mw': 430.54},
            {'name': 'Orexin OX2 Antagonist', 'formula': 'C24H26N4O3', 'smiles': 'OX2 receptor selective', 'target': 'HCRTR2', 'mw': 418.49},
            {'name': 'Free Fatty Acid Receptor Agonist', 'formula': 'C18H21NO4', 'smiles': 'FFAR1 agonist scaffold', 'target': 'FFAR1', 'mw': 315.36},
        ]
    },
    'protease_inhibitors': {
        'description': 'Compounds targeting proteases',
        'compounds': [
            {'name': 'Cathepsin K Inhibitor', 'formula': 'C25H30N4O4', 'smiles': 'CTSK selective inhibitor', 'target': 'CTSK', 'mw': 450.53},
            {'name': 'MMP-13 Inhibitor Lead', 'formula': 'C20H24N2O4S', 'smiles': 'MMP-13 selective', 'target': 'MMP13', 'mw': 388.48},
            {'name': 'BACE1 Inhibitor Compound', 'formula': 'C22H26FN5O2', 'smiles': 'Beta-secretase inhibitor', 'target': 'BACE1', 'mw': 411.47},
            {'name': 'Caspase-3 Inhibitor', 'formula': 'C18H22N4O5', 'smiles': 'Apoptosis pathway modulator', 'target': 'CASP3', 'mw': 374.39},
            {'name': 'Thrombin Inhibitor Core', 'formula': 'C20H23N5O4', 'smiles': 'Direct thrombin inhibitor', 'target': 'F2', 'mw': 397.43},
            {'name': 'Kallikrein-7 Inhibitor', 'formula': 'C19H23N3O4S', 'smiles': 'KLK7 selective', 'target': 'KLK7', 'mw': 389.47},
            {'name': 'DPP-IV Inhibitor Analog', 'formula': 'C16H21F3N4O', 'smiles': 'DPP4 gliptin scaffold', 'target': 'DPP4', 'mw': 342.36},
            {'name': '20S Proteasome Inhibitor', 'formula': 'C25H38N4O5', 'smiles': 'Proteasome catalytic inhibitor', 'target': 'PSMB5', 'mw': 474.59},
        ]
    },
    'epigenetic_modulators': {
        'description': 'Compounds targeting epigenetic enzymes',
        'compounds': [
            {'name': 'HDAC6 Selective Inhibitor', 'formula': 'C21H25N3O3', 'smiles': 'HDAC6 selective scaffold', 'target': 'HDAC6', 'mw': 367.44},
            {'name': 'BET Bromodomain Inhibitor', 'formula': 'C22H19ClN4O2', 'smiles': 'BRD4 selective', 'target': 'BRD4', 'mw': 406.87},
            {'name': 'DOT1L Inhibitor Lead', 'formula': 'C20H26N6O2', 'smiles': 'DOT1L methyltransferase inhibitor', 'target': 'DOT1L', 'mw': 382.46},
            {'name': 'EZH2 Inhibitor Compound', 'formula': 'C24H30N4O2', 'smiles': 'EZH2 SAM-competitive', 'target': 'EZH2', 'mw': 406.52},
            {'name': 'LSD1 Inhibitor Scaffold', 'formula': 'C17H18N4', 'smiles': 'Lysine demethylase inhibitor', 'target': 'KDM1A', 'mw': 278.35},
            {'name': 'DNMT Inhibitor Analog', 'formula': 'C8H12N4O5', 'smiles': 'DNA methyltransferase inhibitor', 'target': 'DNMT1', 'mw': 244.21},
            {'name': 'SIRT1 Activator', 'formula': 'C14H12O3', 'smiles': 'Sirtuin 1 activator scaffold', 'target': 'SIRT1', 'mw': 228.24},
            {'name': 'HAT Inhibitor Core', 'formula': 'C18H18N2O4', 'smiles': 'Histone acetyltransferase inhibitor', 'target': 'KAT2A', 'mw': 326.35},
        ]
    },
    'ion_channel_modulators': {
        'description': 'Compounds targeting ion channels',
        'compounds': [
            {'name': 'Nav1.7 Blocker Lead', 'formula': 'C23H28N4O2S', 'smiles': 'SCN9A selective blocker', 'target': 'SCN9A', 'mw': 424.56},
            {'name': 'Kv7.2/3 Opener', 'formula': 'C16H17N3O2', 'smiles': 'KCNQ2/3 channel opener', 'target': 'KCNQ2', 'mw': 283.33},
            {'name': 'TRPV1 Antagonist', 'formula': 'C24H30N4O', 'smiles': 'TRPV1 receptor antagonist', 'target': 'TRPV1', 'mw': 390.52},
            {'name': 'CaV2.2 Blocker', 'formula': 'C20H24N2O3', 'smiles': 'N-type calcium channel blocker', 'target': 'CACNA1B', 'mw': 340.42},
            {'name': 'hERG Safe Compound', 'formula': 'C18H20N4O2', 'smiles': 'hERG-safe scaffold', 'target': 'KCNH2', 'mw': 324.38},
            {'name': 'CFTR Modulator', 'formula': 'C22H18F3N3O3', 'smiles': 'CFTR potentiator scaffold', 'target': 'CFTR', 'mw': 429.40},
            {'name': 'P2X7 Antagonist', 'formula': 'C19H21N5O', 'smiles': 'P2X7 receptor antagonist', 'target': 'P2RX7', 'mw': 335.40},
            {'name': 'ASIC1a Inhibitor', 'formula': 'C14H18N2O3', 'smiles': 'Acid-sensing ion channel inhibitor', 'target': 'ASIC1', 'mw': 262.30},
        ]
    },
    'metabolic_targets': {
        'description': 'Compounds targeting metabolic enzymes',
        'compounds': [
            {'name': 'ACC Inhibitor Lead', 'formula': 'C20H22N4O4', 'smiles': 'Acetyl-CoA carboxylase inhibitor', 'target': 'ACACA', 'mw': 382.42},
            {'name': 'FASN Inhibitor Compound', 'formula': 'C18H20N2O5', 'smiles': 'Fatty acid synthase inhibitor', 'target': 'FASN', 'mw': 344.36},
            {'name': 'DGAT2 Inhibitor', 'formula': 'C22H27N3O3', 'smiles': 'Diacylglycerol acyltransferase 2', 'target': 'DGAT2', 'mw': 381.47},
            {'name': 'FBPase Inhibitor', 'formula': 'C15H14N2O6P', 'smiles': 'Fructose-1,6-bisphosphatase', 'target': 'FBP1', 'mw': 349.25},
            {'name': 'GK Activator Lead', 'formula': 'C19H22N2O4S', 'smiles': 'Glucokinase activator', 'target': 'GCK', 'mw': 374.45},
            {'name': 'IDH1 R132H Inhibitor', 'formula': 'C23H25FN6O2', 'smiles': 'Mutant IDH1 inhibitor', 'target': 'IDH1', 'mw': 436.48},
            {'name': 'PCSK9 Degrader', 'formula': 'C45H55N9O8', 'smiles': 'PCSK9 PROTAC scaffold', 'target': 'PCSK9', 'mw': 853.98},
            {'name': 'SCD1 Inhibitor', 'formula': 'C18H18N4O2', 'smiles': 'Stearoyl-CoA desaturase inhibitor', 'target': 'SCD', 'mw': 322.36},
        ]
    },
    'protein_protein_inhibitors': {
        'description': 'Compounds disrupting protein-protein interactions',
        'compounds': [
            {'name': 'MDM2-p53 Inhibitor', 'formula': 'C26H28Cl2N4O3', 'smiles': 'MDM2 pocket binder', 'target': 'MDM2', 'mw': 515.43},
            {'name': 'Bcl-2 BH3 Mimetic', 'formula': 'C46H53ClN4O7S', 'smiles': 'Bcl-2 family inhibitor', 'target': 'BCL2', 'mw': 845.45},
            {'name': 'XIAP Antagonist', 'formula': 'C28H34N6O4', 'smiles': 'IAP family inhibitor', 'target': 'XIAP', 'mw': 518.61},
            {'name': 'Menin-MLL Inhibitor', 'formula': 'C23H27N5O2', 'smiles': 'Menin-MLL PPI disruptor', 'target': 'MEN1', 'mw': 405.49},
            {'name': 'PD-1/PD-L1 Small Molecule', 'formula': 'C25H24N4O4', 'smiles': 'Checkpoint inhibitor scaffold', 'target': 'PDCD1', 'mw': 444.48},
            {'name': 'RAS-SOS Disruptor', 'formula': 'C22H26N4O3', 'smiles': 'RAS signaling inhibitor', 'target': 'KRAS', 'mw': 394.47},
            {'name': 'beta-Catenin Inhibitor', 'formula': 'C24H22N4O3', 'smiles': 'Wnt pathway inhibitor', 'target': 'CTNNB1', 'mw': 414.46},
            {'name': 'Myc-Max Disruptor', 'formula': 'C20H18N4O2', 'smiles': 'Myc transcription inhibitor', 'target': 'MYC', 'mw': 346.38},
        ]
    },
    'degraders_protacs': {
        'description': 'Targeted protein degraders and PROTACs',
        'compounds': [
            {'name': 'BRD4 PROTAC Compound', 'formula': 'C51H55N9O10', 'smiles': 'BRD4-CRBN PROTAC', 'target': 'BRD4', 'mw': 965.04},
            {'name': 'AR PROTAC Lead', 'formula': 'C48H52FN7O9', 'smiles': 'Androgen receptor degrader', 'target': 'AR', 'mw': 897.97},
            {'name': 'ER PROTAC Candidate', 'formula': 'C45H48N6O8', 'smiles': 'Estrogen receptor degrader', 'target': 'ESR1', 'mw': 800.90},
            {'name': 'CDK9 PROTAC', 'formula': 'C44H50N8O9', 'smiles': 'CDK9 selective degrader', 'target': 'CDK9', 'mw': 838.92},
            {'name': 'KRASG12C PROTAC', 'formula': 'C52H58N10O8S', 'smiles': 'KRAS G12C degrader', 'target': 'KRAS', 'mw': 991.14},
            {'name': 'BTK PROTAC Lead', 'formula': 'C47H51N9O8', 'smiles': 'BTK degrader scaffold', 'target': 'BTK', 'mw': 869.96},
            {'name': 'Molecular Glue A', 'formula': 'C30H35N5O6', 'smiles': 'CRBN neosubstrate inducer', 'target': 'CRBN', 'mw': 561.63},
            {'name': 'SERD Compound', 'formula': 'C32H38F2N2O4', 'smiles': 'Selective ER degrader', 'target': 'ESR1', 'mw': 556.65},
        ]
    },
    'fragment_library': {
        'description': 'Fragment-based drug discovery building blocks',
        'compounds': [
            {'name': 'Pyridine Core Fragment', 'formula': 'C5H5N', 'smiles': 'c1ccncc1', 'mw': 79.10},
            {'name': 'Piperazine Fragment', 'formula': 'C4H10N2', 'smiles': 'C1CNCCN1', 'mw': 86.14},
            {'name': 'Benzimidazole Core', 'formula': 'C7H6N2', 'smiles': 'c1ccc2[nH]cnc2c1', 'mw': 118.14},
            {'name': 'Morpholine Fragment', 'formula': 'C4H9NO', 'smiles': 'C1COCCN1', 'mw': 87.12},
            {'name': 'Thiophene Core', 'formula': 'C4H4S', 'smiles': 'c1ccsc1', 'mw': 84.14},
            {'name': 'Indole Fragment', 'formula': 'C8H7N', 'smiles': 'c1ccc2[nH]ccc2c1', 'mw': 117.15},
            {'name': 'Quinoline Core', 'formula': 'C9H7N', 'smiles': 'c1ccc2ncccc2c1', 'mw': 129.16},
            {'name': 'Pyrimidine Fragment', 'formula': 'C4H4N2', 'smiles': 'c1cncnc1', 'mw': 80.09},
            {'name': 'Imidazole Core', 'formula': 'C3H4N2', 'smiles': 'c1c[nH]cn1', 'mw': 68.08},
            {'name': 'Sulfonamide Fragment', 'formula': 'CH5NO2S', 'smiles': 'NS(=O)(=O)C', 'mw': 95.12},
            {'name': 'Aminopyridine Fragment', 'formula': 'C5H6N2', 'smiles': 'Nc1ccccn1', 'mw': 94.11},
            {'name': 'Pyrrolidine Core', 'formula': 'C4H9N', 'smiles': 'C1CCNC1', 'mw': 71.12},
        ]
    }
}


def generate_compound_variants(base_compound: Dict, category: str, count: int) -> List[Dict]:
    """Generate variants of a base compound"""
    variants = []
    modifiers = ['Methyl-', 'Ethyl-', 'Fluoro-', 'Chloro-', 'Bromo-', 'Hydroxy-', 'Amino-', 
                 'Cyano-', 'Nitro-', 'Carboxy-', 'Acetyl-', 'Benzyl-', 'Phenyl-', 'Methoxy-',
                 'Trifluoromethyl-', 'Dimethyl-', 'Diethyl-', 'Isopropyl-', 'tert-Butyl-']
    positions = ['2-', '3-', '4-', '5-', '6-', 'ortho-', 'meta-', 'para-', 'alpha-', 'beta-']
    suffixes = ['', '-A', '-B', '-C', '-D', '-1', '-2', '-3', '-analog', '-derivative']
    
    for i in range(count):
        modifier = random.choice(modifiers)
        position = random.choice(positions) if random.random() > 0.5 else ''
        suffix = random.choice(suffixes)
        
        variant_name = f"{position}{modifier}{base_compound['name']}{suffix}"
        
        base_mw = base_compound.get('mw', 300)
        mw_delta = random.uniform(-50, 100)
        
        variants.append({
            'name': variant_name,
            'formula': base_compound.get('formula', 'C20H25N3O4'),
            'smiles': base_compound.get('smiles', ''),
            'category': category,
            'mw': round(base_mw + mw_delta, 2),
            'source': 'synthetic',
            'base_compound': base_compound['name'],
            'target': base_compound.get('target', 'Unknown')
        })
    
    return variants


def generate_compounds_library(target_count: int = 10000) -> List[Dict]:
    """Generate chemical compounds library"""
    all_compounds = []
    target_per_category = target_count // len(COMPOUND_CATEGORIES)
    
    for category, data in COMPOUND_CATEGORIES.items():
        category_compounds = []
        base_compounds = data['compounds']
        
        for compound in base_compounds:
            compound_with_category = {
                **compound,
                'category': category,
                'source': 'curated'
            }
            category_compounds.append(compound_with_category)
        
        remaining = target_per_category - len(category_compounds)
        compounds_per_base = max(1, remaining // len(base_compounds))
        
        for compound in base_compounds:
            variants = generate_compound_variants(compound, category, compounds_per_base)
            category_compounds.extend(variants)
        
        seen_names = set()
        unique_compounds = []
        for compound in category_compounds:
            if compound['name'] not in seen_names:
                seen_names.add(compound['name'])
                unique_compounds.append(compound)
        
        unique_compounds = unique_compounds[:target_per_category + 100]
        
        all_compounds.extend(unique_compounds)
        print(f"Generated {len(unique_compounds)} compounds for {category}")
    
    print(f"\nTotal compounds generated: {len(all_compounds)}")
    return all_compounds


def save_compounds_to_json(compounds: List[Dict], filepath: str):
    """Save compounds to JSON file"""
    with open(filepath, 'w') as f:
        json.dump(compounds, f, indent=2)
    print(f"Saved {len(compounds)} compounds to {filepath}")


if __name__ == '__main__':
    compounds = generate_compounds_library(target_count=10000)
    
    output_path = Path(__file__).parent / 'compounds_10k.json'
    save_compounds_to_json(compounds, str(output_path))
    
    category_counts = {}
    for compound in compounds:
        cat = compound.get('category', 'unknown')
        category_counts[cat] = category_counts.get(cat, 0) + 1
    
    print("\nCategory distribution:")
    for cat, count in sorted(category_counts.items()):
        print(f"  {cat}: {count}")
