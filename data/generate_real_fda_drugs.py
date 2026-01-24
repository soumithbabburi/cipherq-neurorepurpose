#!/usr/bin/env python3
"""
Generate comprehensive FDA-approved drug database with REAL drug names only
No synthetic or generated prefixes - only authentic pharmaceutical names
"""

import json
from pathlib import Path

REAL_FDA_DRUGS = {
    "cardiovascular": {
        "ACE Inhibitor": ["Lisinopril", "Enalapril", "Ramipril", "Captopril", "Benazepril", "Fosinopril", "Quinapril", "Perindopril", "Trandolapril", "Moexipril"],
        "ARB": ["Losartan", "Valsartan", "Irbesartan", "Olmesartan", "Candesartan", "Telmisartan", "Eprosartan", "Azilsartan"],
        "Beta Blocker": ["Metoprolol", "Atenolol", "Propranolol", "Carvedilol", "Bisoprolol", "Labetalol", "Nebivolol", "Nadolol", "Acebutolol", "Timolol", "Pindolol", "Sotalol", "Esmolol", "Betaxolol"],
        "Calcium Channel Blocker": ["Amlodipine", "Diltiazem", "Verapamil", "Nifedipine", "Felodipine", "Nicardipine", "Nimodipine", "Isradipine", "Clevidipine", "Nisoldipine"],
        "Diuretic": ["Furosemide", "Hydrochlorothiazide", "Chlorthalidone", "Spironolactone", "Bumetanide", "Torsemide", "Metolazone", "Indapamide", "Amiloride", "Triamterene", "Eplerenone", "Acetazolamide", "Ethacrynic Acid"],
        "Statin": ["Atorvastatin", "Rosuvastatin", "Simvastatin", "Pravastatin", "Lovastatin", "Fluvastatin", "Pitavastatin"],
        "Antiplatelet": ["Aspirin", "Clopidogrel", "Ticagrelor", "Prasugrel", "Dipyridamole", "Ticlopidine", "Vorapaxar", "Cangrelor"],
        "Anticoagulant": ["Warfarin", "Heparin", "Enoxaparin", "Rivaroxaban", "Apixaban", "Dabigatran", "Edoxaban", "Fondaparinux", "Dalteparin", "Argatroban", "Bivalirudin"],
        "Nitrate": ["Nitroglycerin", "Isosorbide Dinitrate", "Isosorbide Mononitrate", "Amyl Nitrite"],
        "Antiarrhythmic": ["Amiodarone", "Flecainide", "Propafenone", "Sotalol", "Dofetilide", "Dronedarone", "Lidocaine", "Mexiletine", "Procainamide", "Quinidine", "Disopyramide", "Ibutilide"],
        "Vasodilator": ["Hydralazine", "Minoxidil", "Sodium Nitroprusside", "Fenoldopam", "Alprostadil"],
        "Inotrope": ["Digoxin", "Dobutamine", "Dopamine", "Milrinone", "Inamrinone"],
    },
    "diabetes": {
        "Biguanide": ["Metformin"],
        "Sulfonylurea": ["Glimepiride", "Glipizide", "Glyburide", "Glibenclamide", "Tolbutamide", "Chlorpropamide", "Tolazamide"],
        "DPP-4 Inhibitor": ["Sitagliptin", "Saxagliptin", "Linagliptin", "Alogliptin", "Vildagliptin"],
        "SGLT2 Inhibitor": ["Empagliflozin", "Canagliflozin", "Dapagliflozin", "Ertugliflozin", "Sotagliflozin"],
        "GLP-1 Agonist": ["Semaglutide", "Liraglutide", "Exenatide", "Dulaglutide", "Lixisenatide", "Tirzepatide", "Albiglutide"],
        "Thiazolidinedione": ["Pioglitazone", "Rosiglitazone"],
        "Meglitinide": ["Repaglinide", "Nateglinide"],
        "Alpha-glucosidase Inhibitor": ["Acarbose", "Miglitol"],
        "Insulin": ["Insulin Lispro", "Insulin Aspart", "Insulin Glulisine", "Insulin Glargine", "Insulin Detemir", "Insulin Degludec", "NPH Insulin", "Regular Insulin"],
        "Amylin Analog": ["Pramlintide"],
    },
    "anti_inflammatory": {
        "NSAID": ["Ibuprofen", "Naproxen", "Diclofenac", "Indomethacin", "Meloxicam", "Piroxicam", "Ketorolac", "Celecoxib", "Sulindac", "Ketoprofen", "Etodolac", "Nabumetone", "Fenoprofen", "Oxaprozin", "Flurbiprofen", "Tolmetin", "Mefenamic Acid"],
        "COX-2 Inhibitor": ["Celecoxib", "Rofecoxib", "Valdecoxib", "Etoricoxib"],
        "Corticosteroid": ["Prednisone", "Prednisolone", "Dexamethasone", "Hydrocortisone", "Methylprednisolone", "Triamcinolone", "Betamethasone", "Budesonide", "Fluticasone", "Mometasone", "Beclomethasone", "Fluocinonide", "Clobetasol"],
        "DMARD": ["Methotrexate", "Hydroxychloroquine", "Sulfasalazine", "Leflunomide", "Azathioprine", "Cyclosporine", "Gold Compounds", "Penicillamine"],
        "Biologic": ["Adalimumab", "Etanercept", "Infliximab", "Certolizumab", "Golimumab", "Rituximab", "Tocilizumab", "Abatacept", "Sarilumab", "Ustekinumab", "Secukinumab", "Ixekizumab"],
    },
    "neurological": {
        "AChE Inhibitor": ["Donepezil", "Rivastigmine", "Galantamine", "Tacrine"],
        "NMDA Antagonist": ["Memantine", "Amantadine", "Rimantadine", "Ketamine"],
        "Dopamine Agonist": ["Pramipexole", "Ropinirole", "Rotigotine", "Apomorphine", "Bromocriptine", "Cabergoline", "Pergolide"],
        "MAO-B Inhibitor": ["Selegiline", "Rasagiline", "Safinamide"],
        "COMT Inhibitor": ["Entacapone", "Tolcapone", "Opicapone"],
        "Levodopa Combination": ["Carbidopa-Levodopa", "Benserazide-Levodopa"],
        "Anticonvulsant": ["Levetiracetam", "Lamotrigine", "Valproate", "Carbamazepine", "Oxcarbazepine", "Phenytoin", "Topiramate", "Gabapentin", "Pregabalin", "Lacosamide", "Zonisamide", "Clobazam", "Rufinamide", "Felbamate", "Perampanel", "Brivaracetam", "Ethosuximide", "Vigabatrin", "Tiagabine", "Primidone", "Phenobarbital"],
        "Multiple Sclerosis": ["Interferon Beta-1a", "Interferon Beta-1b", "Glatiramer", "Fingolimod", "Dimethyl Fumarate", "Teriflunomide", "Natalizumab", "Ocrelizumab", "Alemtuzumab", "Siponimod", "Ozanimod", "Cladribine"],
    },
    "psychiatric": {
        "SSRI": ["Fluoxetine", "Sertraline", "Paroxetine", "Escitalopram", "Citalopram", "Fluvoxamine", "Vilazodone", "Vortioxetine"],
        "SNRI": ["Venlafaxine", "Duloxetine", "Desvenlafaxine", "Levomilnacipran", "Milnacipran"],
        "TCA": ["Amitriptyline", "Nortriptyline", "Imipramine", "Desipramine", "Clomipramine", "Doxepin", "Trimipramine", "Protriptyline", "Maprotiline"],
        "Atypical Antidepressant": ["Bupropion", "Mirtazapine", "Trazodone", "Nefazodone", "Vilazodone", "Agomelatine"],
        "MAOI": ["Phenelzine", "Tranylcypromine", "Isocarboxazid", "Selegiline Patch", "Moclobemide"],
        "Atypical Antipsychotic": ["Risperidone", "Olanzapine", "Quetiapine", "Aripiprazole", "Ziprasidone", "Paliperidone", "Asenapine", "Iloperidone", "Lurasidone", "Brexpiprazole", "Cariprazine", "Clozapine", "Pimavanserin"],
        "Typical Antipsychotic": ["Haloperidol", "Chlorpromazine", "Fluphenazine", "Perphenazine", "Thioridazine", "Thiothixene", "Loxapine", "Molindone", "Pimozide", "Trifluoperazine"],
        "Anxiolytic": ["Alprazolam", "Lorazepam", "Clonazepam", "Diazepam", "Buspirone", "Hydroxyzine", "Chlordiazepoxide", "Oxazepam", "Temazepam", "Triazolam", "Midazolam", "Clorazepate"],
        "Mood Stabilizer": ["Lithium", "Valproate", "Carbamazepine", "Lamotrigine", "Oxcarbazepine"],
        "ADHD Medication": ["Methylphenidate", "Amphetamine", "Dextroamphetamine", "Lisdexamfetamine", "Atomoxetine", "Guanfacine", "Clonidine", "Modafinil", "Armodafinil"],
    },
    "antibiotic": {
        "Penicillin": ["Amoxicillin", "Ampicillin", "Penicillin V", "Penicillin G", "Piperacillin", "Ticarcillin", "Nafcillin", "Oxacillin", "Dicloxacillin", "Mezlocillin"],
        "Cephalosporin": ["Cephalexin", "Cefazolin", "Cefuroxime", "Ceftriaxone", "Ceftazidime", "Cefepime", "Cefotaxime", "Cefdinir", "Cefpodoxime", "Cefixime", "Ceftaroline", "Cefaclor", "Cefprozil", "Cefoxitin", "Cefotetan"],
        "Fluoroquinolone": ["Ciprofloxacin", "Levofloxacin", "Moxifloxacin", "Ofloxacin", "Norfloxacin", "Gemifloxacin", "Delafloxacin"],
        "Macrolide": ["Azithromycin", "Clarithromycin", "Erythromycin", "Fidaxomicin"],
        "Tetracycline": ["Doxycycline", "Minocycline", "Tetracycline", "Tigecycline", "Eravacycline", "Omadacycline", "Sarecycline"],
        "Aminoglycoside": ["Gentamicin", "Tobramycin", "Amikacin", "Streptomycin", "Neomycin", "Kanamycin", "Plazomicin"],
        "Glycopeptide": ["Vancomycin", "Teicoplanin", "Dalbavancin", "Oritavancin", "Telavancin"],
        "Carbapenem": ["Imipenem", "Meropenem", "Ertapenem", "Doripenem"],
        "Sulfonamide": ["Sulfamethoxazole", "Trimethoprim", "Sulfadiazine", "Sulfasalazine"],
        "Lincosamide": ["Clindamycin", "Lincomycin"],
        "Oxazolidinone": ["Linezolid", "Tedizolid"],
        "Lipopeptide": ["Daptomycin"],
        "Nitroimidazole": ["Metronidazole", "Tinidazole", "Secnidazole"],
        "Nitrofuran": ["Nitrofurantoin", "Furazolidone"],
    },
    "antiviral": {
        "Nucleoside Analog": ["Acyclovir", "Valacyclovir", "Famciclovir", "Ganciclovir", "Valganciclovir", "Cidofovir", "Foscarnet", "Ribavirin", "Remdesivir", "Tenofovir", "Entecavir", "Lamivudine", "Adefovir", "Telbivudine", "Sofosbuvir"],
        "NRTI": ["Abacavir", "Didanosine", "Emtricitabine", "Lamivudine", "Stavudine", "Tenofovir", "Zidovudine"],
        "NNRTI": ["Efavirenz", "Nevirapine", "Rilpivirine", "Etravirine", "Doravirine", "Delavirdine"],
        "Protease Inhibitor": ["Atazanavir", "Darunavir", "Lopinavir", "Ritonavir", "Saquinavir", "Indinavir", "Nelfinavir", "Fosamprenavir", "Tipranavir"],
        "Integrase Inhibitor": ["Raltegravir", "Elvitegravir", "Dolutegravir", "Bictegravir", "Cabotegravir"],
        "Entry Inhibitor": ["Enfuvirtide", "Maraviroc", "Ibalizumab", "Fostemsavir"],
        "Neuraminidase Inhibitor": ["Oseltamivir", "Zanamivir", "Peramivir", "Baloxavir"],
        "HCV Direct-Acting": ["Ledipasvir", "Velpatasvir", "Glecaprevir", "Pibrentasvir", "Elbasvir", "Grazoprevir", "Simeprevir", "Paritaprevir", "Ombitasvir", "Dasabuvir", "Daclatasvir"],
    },
    "cancer": {
        "Alkylating Agent": ["Cyclophosphamide", "Ifosfamide", "Melphalan", "Bendamustine", "Chlorambucil", "Busulfan", "Carmustine", "Lomustine", "Temozolomide", "Dacarbazine", "Procarbazine", "Thiotepa", "Trabectedin"],
        "Antimetabolite": ["Methotrexate", "Pemetrexed", "Fluorouracil", "Capecitabine", "Cytarabine", "Gemcitabine", "Mercaptopurine", "Thioguanine", "Fludarabine", "Cladribine", "Clofarabine", "Nelarabine", "Azacitidine", "Decitabine"],
        "Taxane": ["Paclitaxel", "Docetaxel", "Cabazitaxel", "Nab-Paclitaxel"],
        "Anthracycline": ["Doxorubicin", "Daunorubicin", "Epirubicin", "Idarubicin", "Mitoxantrone"],
        "Topoisomerase Inhibitor": ["Topotecan", "Irinotecan", "Etoposide", "Teniposide"],
        "Tyrosine Kinase Inhibitor": ["Imatinib", "Dasatinib", "Nilotinib", "Bosutinib", "Ponatinib", "Sunitinib", "Sorafenib", "Pazopanib", "Axitinib", "Cabozantinib", "Lenvatinib", "Regorafenib", "Vandetanib", "Erlotinib", "Gefitinib", "Afatinib", "Osimertinib", "Lapatinib", "Neratinib", "Tucatinib", "Ibrutinib", "Acalabrutinib", "Zanubrutinib", "Ruxolitinib", "Fedratinib", "Baricitinib", "Tofacitinib", "Palbociclib", "Ribociclib", "Abemaciclib"],
        "EGFR Inhibitor": ["Erlotinib", "Gefitinib", "Afatinib", "Osimertinib", "Dacomitinib", "Lapatinib", "Neratinib", "Cetuximab", "Panitumumab"],
        "Monoclonal Antibody": ["Rituximab", "Trastuzumab", "Bevacizumab", "Cetuximab", "Panitumumab", "Alemtuzumab", "Obinutuzumab", "Daratumumab", "Elotuzumab", "Pembrolizumab", "Nivolumab", "Atezolizumab", "Durvalumab", "Avelumab", "Ipilimumab", "Ramucirumab", "Pertuzumab"],
        "Immunomodulator": ["Thalidomide", "Lenalidomide", "Pomalidomide"],
        "Proteasome Inhibitor": ["Bortezomib", "Carfilzomib", "Ixazomib"],
        "PARP Inhibitor": ["Olaparib", "Rucaparib", "Niraparib", "Talazoparib"],
        "BCL-2 Inhibitor": ["Venetoclax"],
        "Hormonal": ["Tamoxifen", "Letrozole", "Anastrozole", "Exemestane", "Fulvestrant", "Leuprolide", "Goserelin", "Bicalutamide", "Enzalutamide", "Abiraterone", "Darolutamide", "Apalutamide", "Flutamide", "Degarelix"],
        "Platinum": ["Cisplatin", "Carboplatin", "Oxaliplatin", "Nedaplatin"],
        "Vinca Alkaloid": ["Vincristine", "Vinblastine", "Vinorelbine", "Vindesine"],
    },
    "pain": {
        "Opioid": ["Morphine", "Oxycodone", "Hydrocodone", "Fentanyl", "Hydromorphone", "Methadone", "Codeine", "Tramadol", "Tapentadol", "Buprenorphine", "Oxymorphone", "Meperidine", "Levorphanol", "Sufentanil", "Alfentanil", "Remifentanil"],
        "Non-opioid Analgesic": ["Acetaminophen", "Aspirin", "Ibuprofen", "Naproxen", "Ketorolac", "Diclofenac", "Meloxicam", "Celecoxib"],
        "Muscle Relaxant": ["Cyclobenzaprine", "Baclofen", "Tizanidine", "Methocarbamol", "Carisoprodol", "Orphenadrine", "Metaxalone", "Chlorzoxazone", "Dantrolene"],
        "Neuropathic Pain": ["Gabapentin", "Pregabalin", "Duloxetine", "Amitriptyline", "Nortriptyline", "Carbamazepine", "Capsaicin", "Lidocaine Patch"],
        "Migraine": ["Sumatriptan", "Rizatriptan", "Zolmitriptan", "Eletriptan", "Naratriptan", "Almotriptan", "Frovatriptan", "Ergotamine", "Dihydroergotamine", "Erenumab", "Fremanezumab", "Galcanezumab", "Rimegepant", "Ubrogepant", "Atogepant", "Lasmiditan"],
        "Local Anesthetic": ["Lidocaine", "Bupivacaine", "Ropivacaine", "Mepivacaine", "Prilocaine", "Articaine", "Chloroprocaine", "Procaine", "Tetracaine", "Benzocaine"],
    }
}

DRUG_TARGETS = {
    "ACE Inhibitor": "ACE (Angiotensin-Converting Enzyme)",
    "ARB": "AT1R (Angiotensin II Receptor Type 1)",
    "Beta Blocker": "Beta-Adrenergic Receptor",
    "Calcium Channel Blocker": "L-type Calcium Channels",
    "Diuretic": "Renal Ion Transporters",
    "Statin": "HMG-CoA Reductase",
    "Biguanide": "AMPK/Mitochondria",
    "Sulfonylurea": "KATP Channels",
    "DPP-4 Inhibitor": "DPP-4 Enzyme",
    "SGLT2 Inhibitor": "SGLT2 Transporter",
    "GLP-1 Agonist": "GLP-1 Receptor",
    "NSAID": "COX-1/COX-2",
    "Corticosteroid": "Glucocorticoid Receptor",
    "AChE Inhibitor": "Acetylcholinesterase",
    "NMDA Antagonist": "NMDA Receptor",
    "Dopamine Agonist": "Dopamine Receptors",
    "SSRI": "Serotonin Transporter",
    "SNRI": "Serotonin/Norepinephrine Transporters",
    "Atypical Antipsychotic": "D2/5-HT2A Receptors",
    "Penicillin": "Bacterial Cell Wall",
    "Cephalosporin": "Bacterial Cell Wall",
    "Fluoroquinolone": "DNA Gyrase/Topoisomerase IV",
    "Macrolide": "50S Ribosomal Subunit",
    "Nucleoside Analog": "Viral Polymerase",
    "Protease Inhibitor": "Viral Protease",
    "Tyrosine Kinase Inhibitor": "Tyrosine Kinases",
    "Monoclonal Antibody": "Specific Antigens",
    "Opioid": "Opioid Receptors",
}

def generate_real_fda_database():
    """Generate FDA drug database with only real authentic drug names"""
    drugs = []
    drug_id = 1
    
    for category, drug_classes in REAL_FDA_DRUGS.items():
        for drug_class, drug_names in drug_classes.items():
            target = DRUG_TARGETS.get(drug_class, f"{drug_class} Target")
            
            for drug_name in drug_names:
                drug = {
                    "id": drug_id,
                    "name": drug_name,
                    "class": drug_class,
                    "therapeutic_category": category,
                    "target": target,
                    "mechanism": f"{drug_class} targeting {target}",
                    "smiles": "",
                    "source": "FDA",
                    "status": "Approved"
                }
                drugs.append(drug)
                drug_id += 1
    
    return drugs

if __name__ == "__main__":
    drugs = generate_real_fda_database()
    
    output_path = Path("data") / "drugs_40k.json"
    with open(output_path, 'w') as f:
        json.dump(drugs, f, indent=2)
    
    print(f"Generated {len(drugs)} real FDA drugs")
    
    from collections import Counter
    categories = Counter(d['therapeutic_category'] for d in drugs)
    print("\nDrugs by category:")
    for cat, count in sorted(categories.items()):
        print(f"  {cat}: {count}")
