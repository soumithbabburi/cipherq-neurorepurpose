"""
FDA Drug Database Service
Clean, structured service for loading and querying FDA-approved drugs
Integrates with ChEMBL API for real drug data
"""

import os
import logging
import requests
from typing import List, Dict, Optional, Tuple
from datetime import datetime

logger = logging.getLogger(__name__)

# ChEMBL API base URL
CHEMBL_API_BASE = "https://www.ebi.ac.uk/chembl/api/data"


class DrugDatabaseService:
    """
    Service for managing FDA-approved drug database
    Loads from ChEMBL API and stores in PostgreSQL
    Falls back to JSON storage when database is unavailable
    """
    
    # Hardcoded FDA drugs for instant results (40k simulation)
    HARDCODED_DRUGS = {
        'cardiovascular': [
            {'name': 'Lisinopril', 'chembl_id': 'CHEMBL1237', 'drug_class': 'ACE Inhibitor', 'smiles': 'CC(C)C[C@H](N[C@@H](CCc1ccccc1)C(=O)O)C(=O)N1CCC[C@H]1C(=O)O'},
            {'name': 'Atorvastatin', 'chembl_id': 'CHEMBL1487', 'drug_class': 'Statin', 'smiles': 'CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccc(F)cc2)c(-c2ccccc2)n1CC[C@H](O)C[C@@H](O)CC(=O)O'},
            {'name': 'Metoprolol', 'chembl_id': 'CHEMBL62', 'drug_class': 'Beta Blocker', 'smiles': 'COCCc1ccc(OCC(O)CNC(C)C)cc1'},
            {'name': 'Amlodipine', 'chembl_id': 'CHEMBL22', 'drug_class': 'Calcium Channel Blocker', 'smiles': 'CCOC(=O)C1=C(COCCN)NC(C)=C(C(=O)OC)C1c1ccccc1Cl'},
            {'name': 'Losartan', 'chembl_id': 'CHEMBL191', 'drug_class': 'ARB', 'smiles': 'CCCCc1nc(Cl)c(CO)n1Cc1ccc(-c2ccccc2-c2nn[nH]n2)cc1'},
            {'name': 'Hydrochlorothiazide', 'chembl_id': 'CHEMBL530', 'drug_class': 'Thiazide Diuretic', 'smiles': 'NS(=O)(=O)c1cc2c(cc1Cl)NC(Cl)NS2(=O)=O'},
            {'name': 'Warfarin', 'chembl_id': 'CHEMBL1464', 'drug_class': 'Anticoagulant', 'smiles': 'CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O'},
            {'name': 'Clopidogrel', 'chembl_id': 'CHEMBL1771', 'drug_class': 'Antiplatelet', 'smiles': 'COC(=O)[C@H](c1ccccc1Cl)N1CCc2sccc2C1'},
            {'name': 'Furosemide', 'chembl_id': 'CHEMBL25', 'drug_class': 'Loop Diuretic', 'smiles': 'NS(=O)(=O)c1cc(C(=O)O)c(NCc2ccco2)cc1Cl'},
            {'name': 'Digoxin', 'chembl_id': 'CHEMBL1751', 'drug_class': 'Cardiac Glycoside', 'smiles': 'C[C@H]1O[C@H](O[C@@H]2[C@@H](O)C[C@H](O[C@@H]3[C@@H](O)C[C@H](O[C@@H]4CC[C@@]5(C)[C@H](CC[C@@H]6[C@@H]5CC[C@@]5(C)[C@@H](C7=CC(=O)OC7)CC[C@@H]65)C4)O[C@H]3C)O[C@H]2C)C[C@@H](O)[C@@H]1O'},
        ],
        'diabetes': [
            {'name': 'Metformin', 'chembl_id': 'CHEMBL1431', 'drug_class': 'Biguanide', 'smiles': 'CN(C)C(=N)N=C(N)N'},
            {'name': 'Glyburide', 'chembl_id': 'CHEMBL434', 'drug_class': 'Sulfonylurea', 'smiles': 'COc1ccc(Cl)cc1C(=O)NCCc1ccc(S(=O)(=O)NC(=O)NC2CCCCC2)cc1'},
            {'name': 'Pioglitazone', 'chembl_id': 'CHEMBL595', 'drug_class': 'Thiazolidinedione', 'smiles': 'CCc1ccc(CCOc2ccc(CC3SC(=O)NC3=O)cc2)nc1'},
            {'name': 'Sitagliptin', 'chembl_id': 'CHEMBL1422', 'drug_class': 'DPP-4 Inhibitor', 'smiles': 'N[C@@H](CC(=O)N1CCn2c(nnc2C(F)(F)F)C1)Cc1cc(F)c(F)cc1F'},
            {'name': 'Empagliflozin', 'chembl_id': 'CHEMBL2107027', 'drug_class': 'SGLT2 Inhibitor', 'smiles': 'OC[C@H]1O[C@@H](c2ccc(Cl)c(Cc3ccc(O[C@H]4CCCC4)cc3)c2)[C@H](O)[C@@H](O)[C@@H]1O'},
            {'name': 'Glipizide', 'chembl_id': 'CHEMBL428', 'drug_class': 'Sulfonylurea', 'smiles': 'Cc1cnc(C(=O)NCCc2ccc(S(=O)(=O)NC(=O)NC3CCCCC3)cc2)cn1'},
            {'name': 'Rosiglitazone', 'chembl_id': 'CHEMBL121', 'drug_class': 'Thiazolidinedione', 'smiles': 'CN(CCOc1ccc(CC2SC(=O)NC2=O)cc1)c1ccccn1'},
            {'name': 'Liraglutide', 'chembl_id': 'CHEMBL1201534', 'drug_class': 'GLP-1 Agonist', 'smiles': 'CCCCCCCCCCCCCCCC(=O)N[C@@H](CCC(=O)O)C(=O)N[C@@H](CC1=CN=CN1)C(=O)O'},
            {'name': 'Canagliflozin', 'chembl_id': 'CHEMBL2048484', 'drug_class': 'SGLT2 Inhibitor', 'smiles': 'Cc1ccc(C2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)cc1Cc1ccc(-c2ccc(F)cc2)s1'},
            {'name': 'Repaglinide', 'chembl_id': 'CHEMBL1272', 'drug_class': 'Meglitinide', 'smiles': 'CCOc1cc(CC(=O)N[C@@H](CC(C)C)c2ccccc2N2CCCCC2)ccc1C(=O)O'},
        ],
        'anti_inflammatory': [
            {'name': 'Ibuprofen', 'chembl_id': 'CHEMBL521', 'drug_class': 'NSAID', 'smiles': 'CC(C)Cc1ccc(C(C)C(=O)O)cc1'},
            {'name': 'Naproxen', 'chembl_id': 'CHEMBL137', 'drug_class': 'NSAID', 'smiles': 'COc1ccc2cc(C(C)C(=O)O)ccc2c1'},
            {'name': 'Celecoxib', 'chembl_id': 'CHEMBL118', 'drug_class': 'COX-2 Inhibitor', 'smiles': 'Cc1ccc(-c2cc(C(F)(F)F)nn2-c2ccc(S(N)(=O)=O)cc2)cc1'},
            {'name': 'Prednisone', 'chembl_id': 'CHEMBL635', 'drug_class': 'Corticosteroid', 'smiles': 'CC12CCC(=O)C=C1CCC1C2C(O)CC2(C)C1CCC2(O)C(=O)CO'},
            {'name': 'Aspirin', 'chembl_id': 'CHEMBL25', 'drug_class': 'NSAID', 'smiles': 'CC(=O)Oc1ccccc1C(=O)O'},
            {'name': 'Diclofenac', 'chembl_id': 'CHEMBL139', 'drug_class': 'NSAID', 'smiles': 'OC(=O)Cc1ccccc1Nc1c(Cl)cccc1Cl'},
            {'name': 'Indomethacin', 'chembl_id': 'CHEMBL6', 'drug_class': 'NSAID', 'smiles': 'COc1ccc2c(c1)c(CC(=O)O)c(C)n2C(=O)c1ccc(Cl)cc1'},
            {'name': 'Meloxicam', 'chembl_id': 'CHEMBL599', 'drug_class': 'NSAID', 'smiles': 'Cc1cnc(NC(=O)C2=C(O)c3ccccc3S(=O)(=O)N2C)s1'},
            {'name': 'Dexamethasone', 'chembl_id': 'CHEMBL384467', 'drug_class': 'Corticosteroid', 'smiles': 'CC1CC2C3CCC4=CC(=O)C=CC4(C)C3(F)C(O)CC2(C)C1(O)C(=O)CO'},
            {'name': 'Piroxicam', 'chembl_id': 'CHEMBL527', 'drug_class': 'NSAID', 'smiles': 'CN1C(C(=O)Nc2ccccn2)=C(O)c2ccccc2S1(=O)=O'},
        ],
        'neurological': [
            {'name': 'Donepezil', 'chembl_id': 'CHEMBL502', 'drug_class': 'AChE Inhibitor', 'smiles': 'COc1cc2CC(CC3CCN(Cc4ccccc4)CC3)C(=O)c2cc1OC'},
            {'name': 'Memantine', 'chembl_id': 'CHEMBL807', 'drug_class': 'NMDA Antagonist', 'smiles': 'CC12CC3CC(C)(C1)CC(N)(C3)C2'},
            {'name': 'Levodopa', 'chembl_id': 'CHEMBL1028', 'drug_class': 'Dopamine Precursor', 'smiles': 'N[C@@H](Cc1ccc(O)c(O)c1)C(=O)O'},
            {'name': 'Carbidopa', 'chembl_id': 'CHEMBL1201228', 'drug_class': 'Decarboxylase Inhibitor', 'smiles': 'CC(Cc1ccc(O)c(O)c1)(NN)C(=O)O'},
            {'name': 'Rivastigmine', 'chembl_id': 'CHEMBL95', 'drug_class': 'AChE Inhibitor', 'smiles': 'CCN(C)C(=O)Oc1cccc(C(C)N(C)C)c1'},
            {'name': 'Galantamine', 'chembl_id': 'CHEMBL659', 'drug_class': 'AChE Inhibitor', 'smiles': 'COc1ccc2C3C=C[C@@H]4O[C@@H]5C=C(C2c1)N(C)[C@@H]5[C@@H]34'},
            {'name': 'Pramipexole', 'chembl_id': 'CHEMBL1178', 'drug_class': 'Dopamine Agonist', 'smiles': 'CCCN[C@H]1CCc2nc(N)sc2C1'},
            {'name': 'Ropinirole', 'chembl_id': 'CHEMBL829', 'drug_class': 'Dopamine Agonist', 'smiles': 'CCCN(CCC)CCc1cccc2NC(=O)Cc12'},
            {'name': 'Selegiline', 'chembl_id': 'CHEMBL972', 'drug_class': 'MAO-B Inhibitor', 'smiles': 'C#CCN(C)[C@H](C)Cc1ccccc1'},
            {'name': 'Rasagiline', 'chembl_id': 'CHEMBL933', 'drug_class': 'MAO-B Inhibitor', 'smiles': 'C#CCN[C@H]1CCc2ccccc21'},
        ],
        'psychiatric': [
            {'name': 'Sertraline', 'chembl_id': 'CHEMBL809', 'drug_class': 'SSRI', 'smiles': 'CN[C@H]1CC[C@@H](c2ccc(Cl)c(Cl)c2)c2ccccc21'},
            {'name': 'Fluoxetine', 'chembl_id': 'CHEMBL41', 'drug_class': 'SSRI', 'smiles': 'CNCCC(Oc1ccc(C(F)(F)F)cc1)c1ccccc1'},
            {'name': 'Escitalopram', 'chembl_id': 'CHEMBL1508', 'drug_class': 'SSRI', 'smiles': 'CN(C)CCCC1(c2ccc(F)cc2)OCc2cc(C#N)ccc21'},
            {'name': 'Venlafaxine', 'chembl_id': 'CHEMBL637', 'drug_class': 'SNRI', 'smiles': 'COc1ccc(C(CN(C)C)C2(O)CCCCC2)cc1'},
            {'name': 'Duloxetine', 'chembl_id': 'CHEMBL1175', 'drug_class': 'SNRI', 'smiles': 'CNCC[C@H](Oc1cccc2ccccc12)c1cccs1'},
            {'name': 'Bupropion', 'chembl_id': 'CHEMBL894', 'drug_class': 'NDRI', 'smiles': 'CC(NC(C)(C)C)C(=O)c1cccc(Cl)c1'},
            {'name': 'Quetiapine', 'chembl_id': 'CHEMBL716', 'drug_class': 'Atypical Antipsychotic', 'smiles': 'OCCOCCN1CCN(c2c3ccccc3Sc3ccccc23)CC1'},
            {'name': 'Risperidone', 'chembl_id': 'CHEMBL45', 'drug_class': 'Atypical Antipsychotic', 'smiles': 'Cc1nc2n(c1C(=O)CCN1CCC(c3noc4cc(F)ccc34)CC1)CCCC2'},
            {'name': 'Aripiprazole', 'chembl_id': 'CHEMBL1112', 'drug_class': 'Atypical Antipsychotic', 'smiles': 'Clc1cccc(N2CCN(CCCCOc3ccc4c(c3)CCC(=O)N4)CC2)c1Cl'},
            {'name': 'Olanzapine', 'chembl_id': 'CHEMBL715', 'drug_class': 'Atypical Antipsychotic', 'smiles': 'Cc1cc2c(s1)Nc1ccccc1N=C2N1CCN(C)CC1'},
        ],
        'antibiotic': [
            {'name': 'Amoxicillin', 'chembl_id': 'CHEMBL1082', 'drug_class': 'Penicillin', 'smiles': 'CC1(C)S[C@@H]2[C@H](NC(=O)[C@H](N)c3ccc(O)cc3)C(=O)N2[C@H]1C(=O)O'},
            {'name': 'Azithromycin', 'chembl_id': 'CHEMBL529', 'drug_class': 'Macrolide', 'smiles': 'CC[C@@H]1OC(=O)[C@H](C)[C@@H](O[C@H]2C[C@@](C)(OC)[C@@H](O)[C@H](C)O2)[C@H](C)[C@@H](O[C@@H]2O[C@H](C)C[C@@H](N(C)C)[C@H]2O)[C@](C)(O)C[C@@H](C)CN(C)[C@H](C)[C@@H](O)[C@]1(C)O'},
            {'name': 'Ciprofloxacin', 'chembl_id': 'CHEMBL8', 'drug_class': 'Fluoroquinolone', 'smiles': 'OC(=O)C1=CN(C2CC2)c2cc(N3CCNCC3)c(F)cc2C1=O'},
            {'name': 'Doxycycline', 'chembl_id': 'CHEMBL1560', 'drug_class': 'Tetracycline', 'smiles': 'C[C@@H]1C(O)=C(C(N)=O)C(=O)[C@@]2(O)C(O)=C3C(=O)c4c(O)cccc4[C@@](C)(O)[C@H]3[C@H](O)[C@@H]12'},
            {'name': 'Levofloxacin', 'chembl_id': 'CHEMBL33', 'drug_class': 'Fluoroquinolone', 'smiles': 'C[C@H]1COc2c(N3CCN(C)CC3)c(F)cc3c(=O)c(C(=O)O)cn1c23'},
            {'name': 'Cephalexin', 'chembl_id': 'CHEMBL1727', 'drug_class': 'Cephalosporin', 'smiles': 'C[C@@H]1[C@@H](N2C(=O)[C@@H](NC(=O)[C@H](N)c3ccccc3)C2S1(=O)=O)C(=O)O'},
            {'name': 'Metronidazole', 'chembl_id': 'CHEMBL137', 'drug_class': 'Nitroimidazole', 'smiles': 'Cc1ncc([N+](=O)[O-])n1CCO'},
            {'name': 'Trimethoprim', 'chembl_id': 'CHEMBL22', 'drug_class': 'Antifolate', 'smiles': 'COc1cc(Cc2cnc(N)nc2N)cc(OC)c1OC'},
            {'name': 'Clindamycin', 'chembl_id': 'CHEMBL398', 'drug_class': 'Lincosamide', 'smiles': 'CCC[C@@H]1CC(=O)N(C)[C@H](C(C)C)[C@@H](O)[C@@H](C)O[C@H]1[C@H](O)[C@H](O)[C@@H](NC(=O)C(Cl)Cl)C(C)C'},
            {'name': 'Vancomycin', 'chembl_id': 'CHEMBL262777', 'drug_class': 'Glycopeptide', 'smiles': 'CC(C)CC(NC(=O)C(CC(N)=O)NC(=O)C(NC(=O)C1CC(O)CN1C(=O)C(NC(=O)C(CC(C)C)NC(=O)C(Cc1ccc(O)cc1)NC=O)C(O)c1ccc(O)c(Cl)c1)C(O)c1ccc(O)c(c1)-c1cc(O)cc(O)c1C(O)C(NC(=O)C(CC(O)=O)NC=O)C(=O)O)C(=O)O'},
        ],
        'antiviral': [
            {'name': 'Acyclovir', 'chembl_id': 'CHEMBL184', 'drug_class': 'Nucleoside Analog', 'smiles': 'Nc1nc2c(ncn2COCCO)c(=O)[nH]1'},
            {'name': 'Oseltamivir', 'chembl_id': 'CHEMBL1229', 'drug_class': 'Neuraminidase Inhibitor', 'smiles': 'CCOC(=O)C1=C[C@@H](OC(CC)CC)[C@H](NC(C)=O)[C@@H](N)C1'},
            {'name': 'Remdesivir', 'chembl_id': 'CHEMBL4065906', 'drug_class': 'Nucleotide Analog', 'smiles': 'CCC(CC)COC(=O)C(C)NP(=O)(OCC1OC(C#N)(c2ccc3c(N)ncnc3c2)C(O)C1O)Oc1ccccc1'},
            {'name': 'Sofosbuvir', 'chembl_id': 'CHEMBL1259059', 'drug_class': 'NS5B Inhibitor', 'smiles': 'CC(C)OC(=O)C(C)NP(=O)(OCC1OC(n2ccc(=O)[nH]c2=O)C(F)C1O)Oc1ccccc1'},
            {'name': 'Ritonavir', 'chembl_id': 'CHEMBL163', 'drug_class': 'Protease Inhibitor', 'smiles': 'CC(C)[C@H](NC(=O)N(C)Cc1csc(n1)-c1ccc(C)cc1)C(=O)N[C@@H](C[C@H](O)[C@H](Cc1ccccc1)NC(=O)OCc1cncs1)Cc1ccccc1'},
            {'name': 'Lopinavir', 'chembl_id': 'CHEMBL729', 'drug_class': 'Protease Inhibitor', 'smiles': 'CC(C)[C@H](NC(=O)[C@H](Cc1ccccc1)NC(=O)c1cnc2ccccc2c1)C(=O)N[C@H](C[C@H](O)[C@H](Cc1ccccc1)NC(=O)COc1c(C)cccc1C)Cc1ccccc1'},
            {'name': 'Valacyclovir', 'chembl_id': 'CHEMBL1619', 'drug_class': 'Nucleoside Analog', 'smiles': 'CC(C)[C@H](N)C(=O)OCCOCn1cnc2c1nc(N)[nH]c2=O'},
            {'name': 'Famciclovir', 'chembl_id': 'CHEMBL882', 'drug_class': 'Nucleoside Analog', 'smiles': 'CC(=O)OCCn1cnc2c(N)nc(N)nc21'},
            {'name': 'Zanamivir', 'chembl_id': 'CHEMBL505', 'drug_class': 'Neuraminidase Inhibitor', 'smiles': 'CC(=O)N[C@@H]1[C@@H](N=C(N)N)C=C(C(=O)O)O[C@H]1[C@H](O)[C@H](O)CO'},
            {'name': 'Tenofovir', 'chembl_id': 'CHEMBL505', 'drug_class': 'Nucleotide Analog', 'smiles': 'CC(Cn1cnc2c(N)ncnc21)OCP(=O)(O)O'},
        ],
        'cancer': [
            {'name': 'Imatinib', 'chembl_id': 'CHEMBL941', 'drug_class': 'Tyrosine Kinase Inhibitor', 'smiles': 'Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1'},
            {'name': 'Erlotinib', 'chembl_id': 'CHEMBL553', 'drug_class': 'EGFR Inhibitor', 'smiles': 'COCCOc1cc2ncnc(Nc3cccc(C#C)c3)c2cc1OCCOC'},
            {'name': 'Tamoxifen', 'chembl_id': 'CHEMBL83', 'drug_class': 'SERM', 'smiles': 'CCC(=C(c1ccccc1)c1ccc(OCCN(C)C)cc1)c1ccccc1'},
            {'name': 'Methotrexate', 'chembl_id': 'CHEMBL34259', 'drug_class': 'Antimetabolite', 'smiles': 'CN(Cc1cnc2nc(N)nc(N)c2n1)c1ccc(C(=O)N[C@@H](CCC(=O)O)C(=O)O)cc1'},
            {'name': 'Cisplatin', 'chembl_id': 'CHEMBL11359', 'drug_class': 'Platinum Agent', 'smiles': '[NH3][Pt]([NH3])(Cl)Cl'},
            {'name': 'Paclitaxel', 'chembl_id': 'CHEMBL428647', 'drug_class': 'Taxane', 'smiles': 'CC(=O)OC1C(=O)C2(C)C(O)CC3OC3(C)C(OC(=O)c3ccccc3)C(O)(C(C)=C1C)C2(OC(C)=O)C(=O)C(O)C(NC(=O)c1ccccc1)c1ccccc1'},
            {'name': 'Doxorubicin', 'chembl_id': 'CHEMBL53463', 'drug_class': 'Anthracycline', 'smiles': 'COc1cccc2C(=O)c3c(O)c4C[C@](O)(C(=O)CO)C[C@H](O[C@@H]5C[C@H](N)[C@H](O)[C@H](C)O5)c4c(O)c3C(=O)c12'},
            {'name': 'Capecitabine', 'chembl_id': 'CHEMBL1773', 'drug_class': 'Antimetabolite', 'smiles': 'CCCCCC(=O)Oc1nc(N)n([C@@H]2O[C@H](C)[C@@H](O)[C@H]2O)c1=O'},
            {'name': 'Gefitinib', 'chembl_id': 'CHEMBL939', 'drug_class': 'EGFR Inhibitor', 'smiles': 'COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1'},
            {'name': 'Sorafenib', 'chembl_id': 'CHEMBL1336', 'drug_class': 'Multi-kinase Inhibitor', 'smiles': 'CNC(=O)c1cc(Oc2ccc(NC(=O)Nc3ccc(Cl)c(C(F)(F)F)c3)cc2)ccn1'},
        ],
        'pain': [
            {'name': 'Morphine', 'chembl_id': 'CHEMBL70', 'drug_class': 'Opioid', 'smiles': 'CN1CC[C@]23c4c5ccc(O)c4O[C@H]2C(O)=C[C@@H](O)[C@@H]3[C@H]1C5'},
            {'name': 'Oxycodone', 'chembl_id': 'CHEMBL656', 'drug_class': 'Opioid', 'smiles': 'CN1CC[C@]23c4c5ccc(O)c4O[C@H]2C(=O)CC[C@@H]3[C@H]1C5'},
            {'name': 'Tramadol', 'chembl_id': 'CHEMBL1102', 'drug_class': 'Opioid', 'smiles': 'COc1ccccc1C1(O)CCCCC1CN(C)C'},
            {'name': 'Gabapentin', 'chembl_id': 'CHEMBL940', 'drug_class': 'Anticonvulsant', 'smiles': 'NCC1(CC(=O)O)CCCCC1'},
            {'name': 'Pregabalin', 'chembl_id': 'CHEMBL1059', 'drug_class': 'Anticonvulsant', 'smiles': 'CC(C)C[C@H](CN)CC(=O)O'},
            {'name': 'Acetaminophen', 'chembl_id': 'CHEMBL112', 'drug_class': 'Analgesic', 'smiles': 'CC(=O)Nc1ccc(O)cc1'},
            {'name': 'Codeine', 'chembl_id': 'CHEMBL1613', 'drug_class': 'Opioid', 'smiles': 'COc1ccc2CC3N(C)CCC4=C[C@H](O)[C@@H]5Oc1c2[C@@]45C3'},
            {'name': 'Hydrocodone', 'chembl_id': 'CHEMBL1276', 'drug_class': 'Opioid', 'smiles': 'COc1ccc2CC3N(C)CCC4=CC(=O)[C@@H]5Oc1c2[C@@]45C3'},
            {'name': 'Fentanyl', 'chembl_id': 'CHEMBL596', 'drug_class': 'Opioid', 'smiles': 'CCC(=O)N(c1ccccc1)C1CCN(CCc2ccccc2)CC1'},
            {'name': 'Lidocaine', 'chembl_id': 'CHEMBL54', 'drug_class': 'Local Anesthetic', 'smiles': 'CCN(CC)CC(=O)Nc1c(C)cccc1C'},
        ],
    }
    
    # Therapeutic category mappings for ChEMBL ATC codes
    ATC_CATEGORY_MAP = {
        'A': 'gastrointestinal',
        'B': 'cardiovascular',  # Blood and blood forming organs
        'C': 'cardiovascular',
        'D': 'other',  # Dermatologicals
        'G': 'endocrine',  # Genito-urinary and sex hormones
        'H': 'endocrine',  # Systemic hormonal preparations
        'J': 'antibiotic',  # Anti-infectives
        'L': 'cancer',  # Antineoplastic and immunomodulating
        'M': 'anti_inflammatory',  # Musculo-skeletal
        'N': 'neurological',  # Nervous system
        'P': 'antiviral',  # Antiparasitic products
        'R': 'respiratory',
        'S': 'other',  # Sensory organs
        'V': 'other'  # Various
    }
    
    # Subcategory refinements
    ATC_SUBCATEGORY_MAP = {
        'A10': 'diabetes',  # Drugs used in diabetes
        'J01': 'antibiotic',  # Antibacterials
        'J05': 'antiviral',  # Antivirals
        'L01': 'cancer',  # Antineoplastic agents
        'L04': 'immunology',  # Immunosuppressants
        'M01': 'anti_inflammatory',  # Anti-inflammatory and antirheumatic
        'N02': 'pain',  # Analgesics
        'N03': 'neurological',  # Antiepileptics
        'N04': 'neurological',  # Anti-parkinson drugs
        'N05': 'psychiatric',  # Psycholeptics
        'N06': 'psychiatric',  # Psychoanaleptics
        'N07': 'neurological',  # Other nervous system drugs
    }
    
    def __init__(self):
        self.session = requests.Session()
        self.session.headers.update({
            'Accept': 'application/json',
            'User-Agent': 'CipherQ Drug Repurposing Platform'
        })
        self._db_session = None
        
    def _get_db_session(self):
        """Get database session lazily"""
        if self._db_session is None:
            try:
                from database.connection import db_session
                return db_session()
            except Exception as e:
                logger.warning(f"Database not available: {e}")
                return None
        return self._db_session
    
    def fetch_approved_drugs_from_chembl(
        self, 
        limit: int = 1000, 
        offset: int = 0
    ) -> List[Dict]:
        """
        Fetch FDA-approved drugs from ChEMBL API
        
        Args:
            limit: Number of drugs to fetch per request
            offset: Pagination offset
            
        Returns:
            List of drug dictionaries
        """
        try:
            url = f"{CHEMBL_API_BASE}/molecule.json"
            params = {
                'max_phase': 4,  # Phase 4 = Approved
                'limit': limit,
                'offset': offset,
                'molecule_type': 'Small molecule'
            }
            
            response = self.session.get(url, params=params, timeout=30)
            response.raise_for_status()
            
            data = response.json()
            molecules = data.get('molecules', [])
            
            drugs = []
            for mol in molecules:
                drug = self._parse_chembl_molecule(mol)
                if drug:
                    drugs.append(drug)
            
            logger.info(f"Fetched {len(drugs)} drugs from ChEMBL (offset={offset})")
            return drugs
            
        except Exception as e:
            logger.error(f"Error fetching drugs from ChEMBL: {e}")
            return []
    
    def _parse_chembl_molecule(self, mol: Dict) -> Optional[Dict]:
        """Parse ChEMBL molecule data into drug dictionary"""
        try:
            # Get structure info
            structures = mol.get('molecule_structures') or {}
            properties = mol.get('molecule_properties') or {}
            
            # Determine therapeutic category from ATC codes
            atc_classifications = mol.get('atc_classifications', [])
            category = self._determine_category(atc_classifications)
            
            drug = {
                'name': mol.get('pref_name') or mol.get('molecule_chembl_id'),
                'generic_name': mol.get('pref_name'),
                'chembl_id': mol.get('molecule_chembl_id'),
                'smiles': structures.get('canonical_smiles'),
                'inchi': structures.get('standard_inchi'),
                'inchi_key': structures.get('standard_inchi_key'),
                'molecular_formula': properties.get('full_molformula'),
                'molecular_weight': self._safe_float(properties.get('full_mwt')),
                'logp': self._safe_float(properties.get('alogp')),
                'polar_surface_area': self._safe_float(properties.get('psa')),
                'h_bond_donors': self._safe_int(properties.get('hbd')),
                'h_bond_acceptors': self._safe_int(properties.get('hba')),
                'rotatable_bonds': self._safe_int(properties.get('rtb')),
                'therapeutic_category': category,
                'drug_class': mol.get('molecule_type'),
                'max_phase': mol.get('max_phase', 4),
                'approval_status': 'fda_approved' if mol.get('max_phase') == 4 else 'investigational',
                'black_box_warning': mol.get('black_box_warning', False),
                'atc_codes': atc_classifications
            }
            
            return drug
            
        except Exception as e:
            logger.warning(f"Error parsing molecule: {e}")
            return None
    
    def _determine_category(self, atc_codes: List[str]) -> str:
        """Determine therapeutic category from ATC codes"""
        if not atc_codes:
            return 'other'
        
        # Check subcategories first (more specific)
        for atc in atc_codes:
            prefix3 = atc[:3] if len(atc) >= 3 else ''
            if prefix3 in self.ATC_SUBCATEGORY_MAP:
                return self.ATC_SUBCATEGORY_MAP[prefix3]
        
        # Fall back to main category
        for atc in atc_codes:
            first_letter = atc[0] if atc else ''
            if first_letter in self.ATC_CATEGORY_MAP:
                return self.ATC_CATEGORY_MAP[first_letter]
        
        return 'other'
    
    def _safe_float(self, value) -> Optional[float]:
        """Safely convert to float"""
        try:
            return float(value) if value is not None else None
        except (ValueError, TypeError):
            return None
    
    def _safe_int(self, value) -> Optional[int]:
        """Safely convert to int"""
        try:
            return int(value) if value is not None else None
        except (ValueError, TypeError):
            return None
    
    def load_all_approved_drugs(self, max_drugs: int = 40000) -> int:
        """
        Load all FDA-approved drugs from ChEMBL into database
        
        Args:
            max_drugs: Maximum number of drugs to load
            
        Returns:
            Number of drugs loaded
        """
        from database.connection import db_session
        from database.models import FDADrug
        
        total_loaded = 0
        batch_size = 1000
        offset = 0
        
        logger.info(f"Starting to load up to {max_drugs} FDA-approved drugs...")
        
        while total_loaded < max_drugs:
            drugs = self.fetch_approved_drugs_from_chembl(
                limit=batch_size, 
                offset=offset
            )
            
            if not drugs:
                logger.info("No more drugs to fetch from ChEMBL")
                break
            
            # Insert into database
            with db_session() as session:
                for drug_data in drugs:
                    # Check if drug already exists
                    existing = session.query(FDADrug).filter_by(
                        chembl_id=drug_data.get('chembl_id')
                    ).first()
                    
                    if not existing:
                        drug = FDADrug(**drug_data)
                        session.add(drug)
                        total_loaded += 1
            
            offset += batch_size
            
            if total_loaded >= max_drugs:
                break
            
            logger.info(f"Loaded {total_loaded} drugs so far...")
        
        logger.info(f"Finished loading {total_loaded} FDA-approved drugs")
        return total_loaded
    
    def search_drugs(
        self, 
        query: str = None,
        category: str = None,
        limit: int = 50
    ) -> List[Dict]:
        """
        Search drugs in database with fallback to hardcoded drugs
        
        Args:
            query: Search term (name, SMILES, etc.)
            category: Therapeutic category filter
            limit: Maximum results
            
        Returns:
            List of matching drugs
        """
        try:
            from database.connection import db_session
            from database.models import FDADrug
            
            with db_session() as session:
                q = session.query(FDADrug)
                
                if query:
                    search = f"%{query.lower()}%"
                    q = q.filter(
                        (FDADrug.name.ilike(search)) |
                        (FDADrug.generic_name.ilike(search)) |
                        (FDADrug.chembl_id.ilike(search))
                    )
                
                if category:
                    q = q.filter(FDADrug.therapeutic_category == category.lower())
                
                drugs = q.limit(limit).all()
                return [drug.to_dict() for drug in drugs]
                
        except Exception as e:
            logger.warning(f"Database unavailable, using hardcoded drugs: {e}")
            return self._search_hardcoded_drugs(query, category, limit)
    
    def _search_hardcoded_drugs(
        self, 
        query: str = None,
        category: str = None,
        limit: int = 50
    ) -> List[Dict]:
        """Search in hardcoded drugs when database is unavailable"""
        results = []
        
        categories_to_search = [category.lower().replace(' ', '_')] if category else list(self.HARDCODED_DRUGS.keys())
        
        for cat in categories_to_search:
            if cat not in self.HARDCODED_DRUGS:
                continue
                
            for drug in self.HARDCODED_DRUGS[cat]:
                if query:
                    query_lower = query.lower()
                    if not (query_lower in drug['name'].lower() or 
                            query_lower in drug.get('drug_class', '').lower() or
                            query_lower in drug.get('chembl_id', '').lower()):
                        continue
                
                results.append({
                    'name': drug['name'],
                    'generic_name': drug['name'],
                    'chembl_id': drug.get('chembl_id'),
                    'smiles': drug.get('smiles'),
                    'drug_class': drug.get('drug_class'),
                    'therapeutic_category': cat,
                    'approval_status': 'fda_approved',
                    'molecular_weight': self._estimate_mw(drug.get('smiles', '')),
                    'source': 'hardcoded'
                })
                
                if len(results) >= limit:
                    break
            
            if len(results) >= limit:
                break
        
        return results
    
    def _estimate_mw(self, smiles: str) -> Optional[float]:
        """Estimate molecular weight from SMILES"""
        if not smiles:
            return None
        try:
            from rdkit import Chem
            from rdkit.Chem import Descriptors
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                return round(Descriptors.MolWt(mol), 2)
        except:
            pass
        return None
    
    def get_drugs_by_category(self, category: str, limit: int = 100) -> List[Dict]:
        """Get all drugs in a therapeutic category"""
        return self.search_drugs(category=category, limit=limit)
    
    def get_drug_by_name(self, name: str) -> Optional[Dict]:
        """Get a specific drug by name"""
        results = self.search_drugs(query=name, limit=1)
        return results[0] if results else None
    
    def get_drug_count(self) -> int:
        """Get total number of drugs in database"""
        try:
            from database.connection import db_session
            from database.models import FDADrug
            
            with db_session() as session:
                return session.query(FDADrug).count()
        except Exception as e:
            logger.warning(f"Database unavailable for drug count, using hardcoded: {e}")
            return sum(len(drugs) for drugs in self.HARDCODED_DRUGS.values())
    
    def get_category_counts(self) -> Dict[str, int]:
        """Get drug counts by therapeutic category"""
        try:
            from database.connection import db_session
            from database.models import FDADrug
            from sqlalchemy import func
            
            with db_session() as session:
                results = session.query(
                    FDADrug.therapeutic_category,
                    func.count(FDADrug.id)
                ).group_by(FDADrug.therapeutic_category).all()
                
                return {cat: count for cat, count in results}
        except Exception as e:
            logger.warning(f"Database unavailable for category counts, using hardcoded: {e}")
            return {cat: len(drugs) for cat, drugs in self.HARDCODED_DRUGS.items()}
    
    def get_all_categories(self) -> List[str]:
        """Get all available therapeutic categories"""
        return list(self.HARDCODED_DRUGS.keys())


# Singleton instance
_drug_service = None

def get_drug_database_service() -> DrugDatabaseService:
    """Get singleton drug database service"""
    global _drug_service
    if _drug_service is None:
        _drug_service = DrugDatabaseService()
    return _drug_service
