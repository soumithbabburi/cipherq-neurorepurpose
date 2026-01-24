#!/usr/bin/env python3
"""
Expand FDA drug database with brand names, generics, and additional therapeutic drugs
Target: 5,000+ authentic FDA-approved drugs with real names
"""

import json
from pathlib import Path

ADDITIONAL_FDA_DRUGS = {
    "cardiovascular": {
        "ACE Inhibitor": [
            "Lisinopril", "Enalapril", "Ramipril", "Captopril", "Benazepril", "Fosinopril", 
            "Quinapril", "Perindopril", "Trandolapril", "Moexipril",
            "Prinivil", "Zestril", "Vasotec", "Altace", "Capoten", "Lotensin", "Monopril",
            "Accupril", "Aceon", "Mavik", "Univasc"
        ],
        "ARB": [
            "Losartan", "Valsartan", "Irbesartan", "Olmesartan", "Candesartan", "Telmisartan",
            "Eprosartan", "Azilsartan", "Cozaar", "Diovan", "Avapro", "Benicar", "Atacand",
            "Micardis", "Teveten", "Edarbi"
        ],
        "Beta Blocker": [
            "Metoprolol", "Atenolol", "Propranolol", "Carvedilol", "Bisoprolol", "Labetalol",
            "Nebivolol", "Nadolol", "Acebutolol", "Timolol", "Pindolol", "Sotalol", "Esmolol",
            "Betaxolol", "Lopressor", "Toprol-XL", "Tenormin", "Inderal", "Coreg", "Zebeta",
            "Trandate", "Bystolic", "Corgard", "Sectral", "Blocadren", "Visken", "Betapace",
            "Brevibloc", "Kerlone"
        ],
        "Calcium Channel Blocker": [
            "Amlodipine", "Diltiazem", "Verapamil", "Nifedipine", "Felodipine", "Nicardipine",
            "Nimodipine", "Isradipine", "Clevidipine", "Nisoldipine", "Norvasc", "Cardizem",
            "Calan", "Procardia", "Plendil", "Cardene", "Nimotop", "DynaCirc", "Cleviprex", "Sular"
        ],
        "Diuretic": [
            "Furosemide", "Hydrochlorothiazide", "Chlorthalidone", "Spironolactone", "Bumetanide",
            "Torsemide", "Metolazone", "Indapamide", "Amiloride", "Triamterene", "Eplerenone",
            "Acetazolamide", "Ethacrynic Acid", "Lasix", "Microzide", "Thalitone", "Aldactone",
            "Bumex", "Demadex", "Zaroxolyn", "Lozol", "Midamor", "Dyrenium", "Inspra", "Diamox", "Edecrin"
        ],
        "Statin": [
            "Atorvastatin", "Rosuvastatin", "Simvastatin", "Pravastatin", "Lovastatin",
            "Fluvastatin", "Pitavastatin", "Lipitor", "Crestor", "Zocor", "Pravachol",
            "Mevacor", "Lescol", "Livalo"
        ],
        "Antiplatelet": [
            "Aspirin", "Clopidogrel", "Ticagrelor", "Prasugrel", "Dipyridamole", "Ticlopidine",
            "Vorapaxar", "Cangrelor", "Plavix", "Brilinta", "Effient", "Persantine", "Ticlid",
            "Zontivity", "Kengreal"
        ],
        "Anticoagulant": [
            "Warfarin", "Heparin", "Enoxaparin", "Rivaroxaban", "Apixaban", "Dabigatran",
            "Edoxaban", "Fondaparinux", "Dalteparin", "Argatroban", "Bivalirudin", "Coumadin",
            "Lovenox", "Xarelto", "Eliquis", "Pradaxa", "Savaysa", "Arixtra", "Fragmin", "Angiomax"
        ],
        "Nitrate": [
            "Nitroglycerin", "Isosorbide Dinitrate", "Isosorbide Mononitrate", "Amyl Nitrite",
            "Nitrostat", "Nitro-Bid", "Isordil", "Imdur", "Monoket"
        ],
        "Antiarrhythmic": [
            "Amiodarone", "Flecainide", "Propafenone", "Dofetilide", "Dronedarone", "Lidocaine",
            "Mexiletine", "Procainamide", "Quinidine", "Disopyramide", "Ibutilide", "Adenosine",
            "Cordarone", "Pacerone", "Tambocor", "Rythmol", "Tikosyn", "Multaq", "Xylocaine",
            "Mexitil", "Pronestyl", "Norpace", "Corvert", "Adenocard"
        ],
        "Vasodilator": [
            "Hydralazine", "Minoxidil", "Sodium Nitroprusside", "Fenoldopam", "Alprostadil",
            "Apresoline", "Loniten", "Nitropress", "Corlopam", "Caverject"
        ],
        "Inotrope": [
            "Digoxin", "Dobutamine", "Dopamine", "Milrinone", "Lanoxin", "Dobutrex", "Intropin", "Primacor"
        ],
    },
    "diabetes": {
        "Biguanide": ["Metformin", "Glucophage", "Fortamet", "Glumetza", "Riomet"],
        "Sulfonylurea": [
            "Glimepiride", "Glipizide", "Glyburide", "Chlorpropamide", "Tolbutamide", "Tolazamide",
            "Amaryl", "Glucotrol", "DiaBeta", "Micronase", "Diabinese", "Orinase", "Tolinase"
        ],
        "DPP-4 Inhibitor": [
            "Sitagliptin", "Saxagliptin", "Linagliptin", "Alogliptin",
            "Januvia", "Onglyza", "Tradjenta", "Nesina"
        ],
        "SGLT2 Inhibitor": [
            "Empagliflozin", "Canagliflozin", "Dapagliflozin", "Ertugliflozin",
            "Jardiance", "Invokana", "Farxiga", "Steglatro"
        ],
        "GLP-1 Agonist": [
            "Semaglutide", "Liraglutide", "Exenatide", "Dulaglutide", "Lixisenatide", "Tirzepatide",
            "Ozempic", "Wegovy", "Victoza", "Byetta", "Bydureon", "Trulicity", "Adlyxin", "Mounjaro"
        ],
        "Thiazolidinedione": ["Pioglitazone", "Rosiglitazone", "Actos", "Avandia"],
        "Meglitinide": ["Repaglinide", "Nateglinide", "Prandin", "Starlix"],
        "Alpha-glucosidase Inhibitor": ["Acarbose", "Miglitol", "Precose", "Glyset"],
        "Insulin": [
            "Insulin Lispro", "Insulin Aspart", "Insulin Glulisine", "Insulin Glargine",
            "Insulin Detemir", "Insulin Degludec", "NPH Insulin", "Regular Insulin",
            "Humalog", "NovoLog", "Apidra", "Lantus", "Basaglar", "Toujeo", "Levemir",
            "Tresiba", "Humulin N", "Novolin N", "Humulin R", "Novolin R", "Fiasp", "Admelog"
        ],
        "Amylin Analog": ["Pramlintide", "Symlin"],
    },
    "anti_inflammatory": {
        "NSAID": [
            "Ibuprofen", "Naproxen", "Diclofenac", "Indomethacin", "Meloxicam", "Piroxicam",
            "Ketorolac", "Sulindac", "Ketoprofen", "Etodolac", "Nabumetone", "Fenoprofen",
            "Oxaprozin", "Flurbiprofen", "Tolmetin", "Mefenamic Acid", "Diflunisal",
            "Advil", "Motrin", "Aleve", "Naprosyn", "Voltaren", "Indocin", "Mobic", "Feldene",
            "Toradol", "Clinoril", "Orudis", "Lodine", "Relafen", "Nalfon", "Daypro", "Ansaid",
            "Tolectin", "Ponstel", "Dolobid"
        ],
        "COX-2 Inhibitor": ["Celecoxib", "Celebrex"],
        "Corticosteroid": [
            "Prednisone", "Prednisolone", "Dexamethasone", "Hydrocortisone", "Methylprednisolone",
            "Triamcinolone", "Betamethasone", "Budesonide", "Fluticasone", "Mometasone",
            "Beclomethasone", "Fluocinonide", "Clobetasol", "Halobetasol", "Desoximetasone",
            "Deltasone", "Orapred", "Decadron", "Cortef", "Medrol", "Kenalog", "Celestone",
            "Entocort", "Flovent", "Nasonex", "Qvar", "Temovate", "Ultravate"
        ],
        "DMARD": [
            "Methotrexate", "Hydroxychloroquine", "Sulfasalazine", "Leflunomide", "Azathioprine",
            "Trexall", "Plaquenil", "Azulfidine", "Arava", "Imuran"
        ],
        "Biologic": [
            "Adalimumab", "Etanercept", "Infliximab", "Certolizumab", "Golimumab", "Rituximab",
            "Tocilizumab", "Abatacept", "Sarilumab", "Ustekinumab", "Secukinumab", "Ixekizumab",
            "Tofacitinib", "Baricitinib", "Upadacitinib", "Risankizumab", "Guselkumab",
            "Humira", "Enbrel", "Remicade", "Cimzia", "Simponi", "Rituxan", "Actemra",
            "Orencia", "Kevzara", "Stelara", "Cosentyx", "Taltz", "Xeljanz", "Olumiant",
            "Rinvoq", "Skyrizi", "Tremfya"
        ],
    },
    "neurological": {
        "AChE Inhibitor": [
            "Donepezil", "Rivastigmine", "Galantamine",
            "Aricept", "Exelon", "Razadyne"
        ],
        "NMDA Antagonist": [
            "Memantine", "Amantadine", "Ketamine",
            "Namenda", "Symmetrel", "Ketalar"
        ],
        "Dopamine Agonist": [
            "Pramipexole", "Ropinirole", "Rotigotine", "Apomorphine", "Bromocriptine", "Cabergoline",
            "Mirapex", "Requip", "Neupro", "Apokyn", "Parlodel", "Dostinex"
        ],
        "MAO-B Inhibitor": ["Selegiline", "Rasagiline", "Safinamide", "Eldepryl", "Zelapar", "Azilect", "Xadago"],
        "COMT Inhibitor": ["Entacapone", "Tolcapone", "Opicapone", "Comtan", "Tasmar", "Ongentys"],
        "Levodopa Combination": ["Carbidopa-Levodopa", "Sinemet", "Rytary", "Stalevo", "Duopa"],
        "Anticonvulsant": [
            "Levetiracetam", "Lamotrigine", "Valproate", "Carbamazepine", "Oxcarbazepine",
            "Phenytoin", "Topiramate", "Gabapentin", "Pregabalin", "Lacosamide", "Zonisamide",
            "Clobazam", "Rufinamide", "Felbamate", "Perampanel", "Brivaracetam", "Ethosuximide",
            "Vigabatrin", "Tiagabine", "Primidone", "Phenobarbital", "Clonazepam",
            "Keppra", "Lamictal", "Depakote", "Tegretol", "Trileptal", "Dilantin", "Topamax",
            "Neurontin", "Lyrica", "Vimpat", "Zonegran", "Onfi", "Banzel", "Felbatol", "Fycompa",
            "Briviact", "Zarontin", "Sabril", "Gabitril", "Mysoline", "Luminal", "Klonopin"
        ],
        "Multiple Sclerosis": [
            "Interferon Beta-1a", "Interferon Beta-1b", "Glatiramer", "Fingolimod", "Dimethyl Fumarate",
            "Teriflunomide", "Natalizumab", "Ocrelizumab", "Alemtuzumab", "Siponimod", "Ozanimod", "Cladribine",
            "Avonex", "Rebif", "Betaseron", "Copaxone", "Gilenya", "Tecfidera", "Aubagio",
            "Tysabri", "Ocrevus", "Lemtrada", "Mayzent", "Zeposia", "Mavenclad"
        ],
    },
    "psychiatric": {
        "SSRI": [
            "Fluoxetine", "Sertraline", "Paroxetine", "Escitalopram", "Citalopram", "Fluvoxamine",
            "Vilazodone", "Vortioxetine", "Prozac", "Zoloft", "Paxil", "Lexapro", "Celexa",
            "Luvox", "Viibryd", "Trintellix"
        ],
        "SNRI": [
            "Venlafaxine", "Duloxetine", "Desvenlafaxine", "Levomilnacipran", "Milnacipran",
            "Effexor", "Cymbalta", "Pristiq", "Fetzima", "Savella"
        ],
        "TCA": [
            "Amitriptyline", "Nortriptyline", "Imipramine", "Desipramine", "Clomipramine",
            "Doxepin", "Trimipramine", "Protriptyline", "Maprotiline",
            "Elavil", "Pamelor", "Tofranil", "Norpramin", "Anafranil", "Sinequan", "Surmontil", "Vivactil"
        ],
        "Atypical Antidepressant": [
            "Bupropion", "Mirtazapine", "Trazodone", "Nefazodone",
            "Wellbutrin", "Remeron", "Desyrel", "Serzone"
        ],
        "MAOI": [
            "Phenelzine", "Tranylcypromine", "Isocarboxazid", "Selegiline Patch",
            "Nardil", "Parnate", "Marplan", "Emsam"
        ],
        "Atypical Antipsychotic": [
            "Risperidone", "Olanzapine", "Quetiapine", "Aripiprazole", "Ziprasidone", "Paliperidone",
            "Asenapine", "Iloperidone", "Lurasidone", "Brexpiprazole", "Cariprazine", "Clozapine", "Pimavanserin",
            "Risperdal", "Zyprexa", "Seroquel", "Abilify", "Geodon", "Invega", "Saphris",
            "Fanapt", "Latuda", "Rexulti", "Vraylar", "Clozaril", "Nuplazid"
        ],
        "Typical Antipsychotic": [
            "Haloperidol", "Chlorpromazine", "Fluphenazine", "Perphenazine", "Thioridazine",
            "Thiothixene", "Loxapine", "Molindone", "Pimozide", "Trifluoperazine",
            "Haldol", "Thorazine", "Prolixin", "Trilafon", "Mellaril", "Navane",
            "Loxitane", "Moban", "Orap", "Stelazine"
        ],
        "Anxiolytic": [
            "Alprazolam", "Lorazepam", "Clonazepam", "Diazepam", "Buspirone", "Hydroxyzine",
            "Chlordiazepoxide", "Oxazepam", "Temazepam", "Triazolam", "Midazolam", "Clorazepate",
            "Xanax", "Ativan", "Klonopin", "Valium", "BuSpar", "Vistaril", "Atarax",
            "Librium", "Serax", "Restoril", "Halcion", "Versed", "Tranxene"
        ],
        "Mood Stabilizer": ["Lithium", "Lithobid", "Eskalith"],
        "ADHD Medication": [
            "Methylphenidate", "Amphetamine", "Dextroamphetamine", "Lisdexamfetamine", "Atomoxetine",
            "Guanfacine", "Clonidine", "Modafinil", "Armodafinil",
            "Ritalin", "Concerta", "Adderall", "Dexedrine", "Vyvanse", "Strattera",
            "Intuniv", "Kapvay", "Provigil", "Nuvigil"
        ],
        "Sleep Aid": [
            "Zolpidem", "Eszopiclone", "Zaleplon", "Ramelteon", "Suvorexant", "Lemborexant",
            "Ambien", "Lunesta", "Sonata", "Rozerem", "Belsomra", "Dayvigo"
        ],
    },
    "antibiotic": {
        "Penicillin": [
            "Amoxicillin", "Ampicillin", "Penicillin V", "Penicillin G", "Piperacillin", "Ticarcillin",
            "Nafcillin", "Oxacillin", "Dicloxacillin", "Amoxil", "Principen", "Pen-Vee K",
            "Zosyn", "Timentin", "Nallpen", "Bactocill"
        ],
        "Cephalosporin": [
            "Cephalexin", "Cefazolin", "Cefuroxime", "Ceftriaxone", "Ceftazidime", "Cefepime",
            "Cefotaxime", "Cefdinir", "Cefpodoxime", "Cefixime", "Ceftaroline", "Cefaclor",
            "Cefprozil", "Cefoxitin", "Cefotetan", "Ceftobiprole",
            "Keflex", "Ancef", "Ceftin", "Rocephin", "Fortaz", "Maxipime", "Claforan",
            "Omnicef", "Vantin", "Suprax", "Teflaro", "Ceclor", "Cefzil", "Mefoxin"
        ],
        "Fluoroquinolone": [
            "Ciprofloxacin", "Levofloxacin", "Moxifloxacin", "Ofloxacin", "Norfloxacin",
            "Gemifloxacin", "Delafloxacin",
            "Cipro", "Levaquin", "Avelox", "Floxin", "Noroxin", "Factive", "Baxdela"
        ],
        "Macrolide": [
            "Azithromycin", "Clarithromycin", "Erythromycin", "Fidaxomicin",
            "Zithromax", "Z-Pak", "Biaxin", "Ery-Tab", "Dificid"
        ],
        "Tetracycline": [
            "Doxycycline", "Minocycline", "Tetracycline", "Tigecycline", "Eravacycline",
            "Omadacycline", "Sarecycline", "Vibramycin", "Minocin", "Tygacil", "Xerava", "Nuzyra", "Seysara"
        ],
        "Aminoglycoside": [
            "Gentamicin", "Tobramycin", "Amikacin", "Streptomycin", "Neomycin", "Kanamycin", "Plazomicin",
            "Garamycin", "Tobrex", "Amikin", "Neo-Fradin", "Kantrex", "Zemdri"
        ],
        "Glycopeptide": [
            "Vancomycin", "Teicoplanin", "Dalbavancin", "Oritavancin", "Telavancin",
            "Vancocin", "Dalvance", "Orbactiv", "Vibativ"
        ],
        "Carbapenem": [
            "Imipenem", "Meropenem", "Ertapenem", "Doripenem",
            "Primaxin", "Merrem", "Invanz", "Doribax"
        ],
        "Sulfonamide": [
            "Sulfamethoxazole", "Trimethoprim", "Sulfadiazine", "Sulfasalazine",
            "Bactrim", "Septra"
        ],
        "Lincosamide": ["Clindamycin", "Lincomycin", "Cleocin", "Lincocin"],
        "Oxazolidinone": ["Linezolid", "Tedizolid", "Zyvox", "Sivextro"],
        "Lipopeptide": ["Daptomycin", "Cubicin"],
        "Nitroimidazole": ["Metronidazole", "Tinidazole", "Secnidazole", "Flagyl", "Tindamax", "Solosec"],
        "Nitrofuran": ["Nitrofurantoin", "Macrobid", "Macrodantin"],
    },
    "antiviral": {
        "Nucleoside Analog": [
            "Acyclovir", "Valacyclovir", "Famciclovir", "Ganciclovir", "Valganciclovir", "Cidofovir",
            "Foscarnet", "Ribavirin", "Remdesivir", "Tenofovir", "Entecavir", "Lamivudine", "Adefovir",
            "Telbivudine", "Sofosbuvir",
            "Zovirax", "Valtrex", "Famvir", "Cytovene", "Valcyte", "Vistide", "Foscavir",
            "Copegus", "Veklury", "Viread", "Baraclude", "Epivir", "Hepsera", "Tyzeka", "Sovaldi"
        ],
        "HIV Antiretroviral": [
            "Abacavir", "Didanosine", "Emtricitabine", "Stavudine", "Zidovudine",
            "Efavirenz", "Nevirapine", "Rilpivirine", "Etravirine", "Doravirine", "Delavirdine",
            "Atazanavir", "Darunavir", "Lopinavir", "Ritonavir", "Saquinavir", "Indinavir",
            "Nelfinavir", "Fosamprenavir", "Tipranavir",
            "Raltegravir", "Elvitegravir", "Dolutegravir", "Bictegravir", "Cabotegravir",
            "Enfuvirtide", "Maraviroc", "Ibalizumab", "Fostemsavir",
            "Ziagen", "Videx", "Emtriva", "Zerit", "Retrovir", "Sustiva", "Viramune",
            "Edurant", "Intelence", "Pifeltro", "Rescriptor", "Reyataz", "Prezista",
            "Kaletra", "Norvir", "Invirase", "Crixivan", "Viracept", "Lexiva", "Aptivus",
            "Isentress", "Vitekta", "Tivicay", "Biktarvy", "Cabenuva", "Fuzeon", "Selzentry",
            "Trogarzo", "Rukobia"
        ],
        "Neuraminidase Inhibitor": [
            "Oseltamivir", "Zanamivir", "Peramivir", "Baloxavir",
            "Tamiflu", "Relenza", "Rapivab", "Xofluza"
        ],
        "HCV Direct-Acting": [
            "Ledipasvir", "Velpatasvir", "Glecaprevir", "Pibrentasvir", "Elbasvir", "Grazoprevir",
            "Simeprevir", "Paritaprevir", "Ombitasvir", "Dasabuvir", "Daclatasvir",
            "Harvoni", "Epclusa", "Mavyret", "Zepatier", "Olysio", "Viekira Pak", "Daklinza"
        ],
    },
    "cancer": {
        "Alkylating Agent": [
            "Cyclophosphamide", "Ifosfamide", "Melphalan", "Bendamustine", "Chlorambucil",
            "Busulfan", "Carmustine", "Lomustine", "Temozolomide", "Dacarbazine", "Procarbazine",
            "Thiotepa", "Trabectedin",
            "Cytoxan", "Ifex", "Alkeran", "Treanda", "Leukeran", "Myleran", "BiCNU",
            "CeeNU", "Temodar", "DTIC-Dome", "Matulane", "Thioplex", "Yondelis"
        ],
        "Antimetabolite": [
            "Methotrexate", "Pemetrexed", "Fluorouracil", "Capecitabine", "Cytarabine",
            "Gemcitabine", "Mercaptopurine", "Thioguanine", "Fludarabine", "Cladribine",
            "Clofarabine", "Nelarabine", "Azacitidine", "Decitabine",
            "Trexall", "Alimta", "Adrucil", "Xeloda", "Cytosar-U", "Gemzar", "Purinethol",
            "Tabloid", "Fludara", "Leustatin", "Clolar", "Arranon", "Vidaza", "Dacogen"
        ],
        "Taxane": [
            "Paclitaxel", "Docetaxel", "Cabazitaxel", "Nab-Paclitaxel",
            "Taxol", "Taxotere", "Jevtana", "Abraxane"
        ],
        "Anthracycline": [
            "Doxorubicin", "Daunorubicin", "Epirubicin", "Idarubicin", "Mitoxantrone",
            "Adriamycin", "Doxil", "Cerubidine", "Ellence", "Idamycin", "Novantrone"
        ],
        "Topoisomerase Inhibitor": [
            "Topotecan", "Irinotecan", "Etoposide", "Teniposide",
            "Hycamtin", "Camptosar", "VePesid", "Vumon"
        ],
        "Tyrosine Kinase Inhibitor": [
            "Imatinib", "Dasatinib", "Nilotinib", "Bosutinib", "Ponatinib", "Sunitinib",
            "Sorafenib", "Pazopanib", "Axitinib", "Cabozantinib", "Lenvatinib", "Regorafenib",
            "Vandetanib", "Erlotinib", "Gefitinib", "Afatinib", "Osimertinib", "Lapatinib",
            "Neratinib", "Tucatinib", "Ibrutinib", "Acalabrutinib", "Zanubrutinib", "Ruxolitinib",
            "Fedratinib", "Palbociclib", "Ribociclib", "Abemaciclib", "Encorafenib", "Binimetinib",
            "Trametinib", "Cobimetinib", "Dabrafenib", "Vemurafenib", "Lorlatinib", "Crizotinib",
            "Ceritinib", "Alectinib", "Brigatinib",
            "Gleevec", "Sprycel", "Tasigna", "Bosulif", "Iclusig", "Sutent", "Nexavar",
            "Votrient", "Inlyta", "Cabometyx", "Lenvima", "Stivarga", "Caprelsa", "Tarceva",
            "Iressa", "Gilotrif", "Tagrisso", "Tykerb", "Nerlynx", "Tukysa", "Imbruvica",
            "Calquence", "Brukinsa", "Jakafi", "Inrebic", "Ibrance", "Kisqali", "Verzenio",
            "Braftovi", "Mektovi", "Mekinist", "Cotellic", "Tafinlar", "Zelboraf", "Lorbrena",
            "Xalkori", "Zykadia", "Alecensa", "Alunbrig"
        ],
        "Monoclonal Antibody": [
            "Rituximab", "Trastuzumab", "Bevacizumab", "Cetuximab", "Panitumumab", "Alemtuzumab",
            "Obinutuzumab", "Daratumumab", "Elotuzumab", "Pembrolizumab", "Nivolumab",
            "Atezolizumab", "Durvalumab", "Avelumab", "Ipilimumab", "Ramucirumab", "Pertuzumab",
            "Trastuzumab Deruxtecan", "Enfortumab", "Sacituzumab", "Brentuximab",
            "Rituxan", "Herceptin", "Avastin", "Erbitux", "Vectibix", "Campath", "Gazyva",
            "Darzalex", "Empliciti", "Keytruda", "Opdivo", "Tecentriq", "Imfinzi", "Bavencio",
            "Yervoy", "Cyramza", "Perjeta", "Enhertu", "Padcev", "Trodelvy", "Adcetris"
        ],
        "Immunomodulator": [
            "Thalidomide", "Lenalidomide", "Pomalidomide",
            "Thalomid", "Revlimid", "Pomalyst"
        ],
        "Proteasome Inhibitor": [
            "Bortezomib", "Carfilzomib", "Ixazomib",
            "Velcade", "Kyprolis", "Ninlaro"
        ],
        "PARP Inhibitor": [
            "Olaparib", "Rucaparib", "Niraparib", "Talazoparib",
            "Lynparza", "Rubraca", "Zejula", "Talzenna"
        ],
        "BCL-2 Inhibitor": ["Venetoclax", "Venclexta"],
        "Hormonal": [
            "Tamoxifen", "Letrozole", "Anastrozole", "Exemestane", "Fulvestrant", "Leuprolide",
            "Goserelin", "Bicalutamide", "Enzalutamide", "Abiraterone", "Darolutamide",
            "Apalutamide", "Flutamide", "Degarelix",
            "Nolvadex", "Femara", "Arimidex", "Aromasin", "Faslodex", "Lupron", "Eligard",
            "Zoladex", "Casodex", "Xtandi", "Zytiga", "Nubeqa", "Erleada", "Eulexin", "Firmagon"
        ],
        "Platinum": [
            "Cisplatin", "Carboplatin", "Oxaliplatin",
            "Platinol", "Paraplatin", "Eloxatin"
        ],
        "Vinca Alkaloid": [
            "Vincristine", "Vinblastine", "Vinorelbine",
            "Oncovin", "Velban", "Navelbine"
        ],
    },
    "pain": {
        "Opioid": [
            "Morphine", "Oxycodone", "Hydrocodone", "Fentanyl", "Hydromorphone", "Methadone",
            "Codeine", "Tramadol", "Tapentadol", "Buprenorphine", "Oxymorphone", "Meperidine",
            "Levorphanol", "Sufentanil", "Alfentanil", "Remifentanil",
            "MS Contin", "Kadian", "OxyContin", "Percocet", "Vicodin", "Norco", "Duragesic",
            "Actiq", "Dilaudid", "Dolophine", "Tylenol #3", "Ultram", "Nucynta", "Subutex",
            "Belbuca", "Opana", "Demerol", "Sublimaze", "Ultiva"
        ],
        "Non-opioid Analgesic": [
            "Acetaminophen", "Aspirin",
            "Tylenol", "Bayer", "Excedrin"
        ],
        "Muscle Relaxant": [
            "Cyclobenzaprine", "Baclofen", "Tizanidine", "Methocarbamol", "Carisoprodol",
            "Orphenadrine", "Metaxalone", "Chlorzoxazone", "Dantrolene",
            "Flexeril", "Lioresal", "Zanaflex", "Robaxin", "Soma", "Norflex", "Skelaxin", "Parafon Forte", "Dantrium"
        ],
        "Neuropathic Pain": [
            "Gabapentin", "Pregabalin", "Duloxetine", "Capsaicin", "Lidocaine Patch",
            "Neurontin", "Lyrica", "Cymbalta", "Zostrix", "Lidoderm"
        ],
        "Migraine": [
            "Sumatriptan", "Rizatriptan", "Zolmitriptan", "Eletriptan", "Naratriptan",
            "Almotriptan", "Frovatriptan", "Ergotamine", "Dihydroergotamine", "Erenumab",
            "Fremanezumab", "Galcanezumab", "Rimegepant", "Ubrogepant", "Atogepant", "Lasmiditan",
            "Imitrex", "Maxalt", "Zomig", "Relpax", "Amerge", "Axert", "Frova",
            "Cafergot", "Migranal", "Aimovig", "Ajovy", "Emgality", "Nurtec", "Ubrelvy", "Qulipta", "Reyvow"
        ],
        "Local Anesthetic": [
            "Lidocaine", "Bupivacaine", "Ropivacaine", "Mepivacaine", "Prilocaine",
            "Articaine", "Chloroprocaine", "Procaine", "Tetracaine", "Benzocaine",
            "Xylocaine", "Marcaine", "Naropin", "Carbocaine", "Citanest", "Septocaine",
            "Nesacaine", "Novocain", "Pontocaine", "Americaine"
        ],
    }
}

DRUG_TARGETS = {
    "ACE Inhibitor": "ACE",
    "ARB": "AT1R",
    "Beta Blocker": "Beta-Adrenergic Receptor",
    "Calcium Channel Blocker": "L-type Calcium Channels",
    "Diuretic": "Renal Ion Transporters",
    "Statin": "HMG-CoA Reductase",
    "Antiplatelet": "Platelet Aggregation",
    "Anticoagulant": "Coagulation Factors",
    "Nitrate": "Vascular Smooth Muscle",
    "Antiarrhythmic": "Cardiac Ion Channels",
    "Vasodilator": "Vascular Smooth Muscle",
    "Inotrope": "Cardiac Contractility",
    "Biguanide": "AMPK",
    "Sulfonylurea": "KATP Channels",
    "DPP-4 Inhibitor": "DPP-4",
    "SGLT2 Inhibitor": "SGLT2",
    "GLP-1 Agonist": "GLP-1R",
    "Thiazolidinedione": "PPAR-gamma",
    "Meglitinide": "KATP Channels",
    "Alpha-glucosidase Inhibitor": "Alpha-glucosidase",
    "Insulin": "Insulin Receptor",
    "Amylin Analog": "Amylin Receptor",
    "NSAID": "COX-1/COX-2",
    "COX-2 Inhibitor": "COX-2",
    "Corticosteroid": "Glucocorticoid Receptor",
    "DMARD": "Immune Cells",
    "Biologic": "Inflammatory Cytokines",
    "AChE Inhibitor": "Acetylcholinesterase",
    "NMDA Antagonist": "NMDA Receptor",
    "Dopamine Agonist": "Dopamine Receptors",
    "MAO-B Inhibitor": "MAO-B",
    "COMT Inhibitor": "COMT",
    "Levodopa Combination": "Dopamine Synthesis",
    "Anticonvulsant": "Neuronal Ion Channels",
    "Multiple Sclerosis": "Immune Modulation",
    "SSRI": "Serotonin Transporter",
    "SNRI": "Serotonin/Norepinephrine Transporters",
    "TCA": "Monoamine Transporters",
    "Atypical Antidepressant": "Various Receptors",
    "MAOI": "Monoamine Oxidase",
    "Atypical Antipsychotic": "D2/5-HT2A Receptors",
    "Typical Antipsychotic": "D2 Receptor",
    "Anxiolytic": "GABA-A Receptor",
    "Mood Stabilizer": "Ion Channels/Signaling",
    "ADHD Medication": "Dopamine/Norepinephrine",
    "Sleep Aid": "GABA/Orexin",
    "Penicillin": "Bacterial Cell Wall",
    "Cephalosporin": "Bacterial Cell Wall",
    "Fluoroquinolone": "DNA Gyrase",
    "Macrolide": "50S Ribosome",
    "Tetracycline": "30S Ribosome",
    "Aminoglycoside": "30S Ribosome",
    "Glycopeptide": "Bacterial Cell Wall",
    "Carbapenem": "Bacterial Cell Wall",
    "Sulfonamide": "Folate Synthesis",
    "Lincosamide": "50S Ribosome",
    "Oxazolidinone": "50S Ribosome",
    "Lipopeptide": "Cell Membrane",
    "Nitroimidazole": "DNA",
    "Nitrofuran": "DNA",
    "Nucleoside Analog": "Viral Polymerase",
    "HIV Antiretroviral": "HIV Lifecycle",
    "Neuraminidase Inhibitor": "Neuraminidase",
    "HCV Direct-Acting": "HCV Proteins",
    "Alkylating Agent": "DNA",
    "Antimetabolite": "DNA Synthesis",
    "Taxane": "Microtubules",
    "Anthracycline": "DNA/Topoisomerase",
    "Topoisomerase Inhibitor": "Topoisomerase",
    "Tyrosine Kinase Inhibitor": "Tyrosine Kinases",
    "Monoclonal Antibody": "Specific Antigens",
    "Immunomodulator": "Immune System",
    "Proteasome Inhibitor": "Proteasome",
    "PARP Inhibitor": "PARP",
    "BCL-2 Inhibitor": "BCL-2",
    "Hormonal": "Hormone Receptors",
    "Platinum": "DNA",
    "Vinca Alkaloid": "Microtubules",
    "Opioid": "Opioid Receptors",
    "Non-opioid Analgesic": "COX/Central",
    "Muscle Relaxant": "CNS/Muscle",
    "Neuropathic Pain": "Neuronal Channels",
    "Migraine": "Serotonin/CGRP",
    "Local Anesthetic": "Sodium Channels",
}

def generate_expanded_fda_database():
    """Generate expanded FDA drug database with real authentic drug names"""
    drugs = []
    seen_names = set()
    drug_id = 1
    
    for category, drug_classes in ADDITIONAL_FDA_DRUGS.items():
        for drug_class, drug_names in drug_classes.items():
            target = DRUG_TARGETS.get(drug_class, f"{drug_class} Target")
            
            for drug_name in drug_names:
                if drug_name.lower() in seen_names:
                    continue
                seen_names.add(drug_name.lower())
                
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
    drugs = generate_expanded_fda_database()
    
    output_path = Path("data") / "drugs_40k.json"
    with open(output_path, 'w') as f:
        json.dump(drugs, f, indent=2)
    
    print(f"Generated {len(drugs)} real FDA drugs")
    
    from collections import Counter
    categories = Counter(d['therapeutic_category'] for d in drugs)
    print("\nDrugs by category:")
    for cat, count in sorted(categories.items()):
        print(f"  {cat}: {count}")
