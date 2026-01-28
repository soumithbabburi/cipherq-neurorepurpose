import json

def load_data(drugs_path, interactions_path):
    """Loads drug and interaction data from JSON files."""
    with open(drugs_path, 'r') as f:
        drugs_data = json.load(f)
    
    with open(interactions_path, 'r') as f:
        interactions_data = json.load(f)
        
    return drugs_data, interactions_data

def get_comprehensive_drug_profile(drug_name, drugs_data, interactions_data):
    """
    Combines molecular properties and protein interaction data for a specific drug.
    """
    # Normalize drug name for lookup
    drug_key = drug_name.lower()
    
    # Retrieve base molecular properties from drugs.json
    # Examples in data: '1,2-benzodiazepine', 'abacavir', 'abemaciclib'
    drug_info = drugs_data.get(drug_key)
    
    if not drug_info:
        return f"Drug '{drug_name}' not found in database."

    # Retrieve interaction/target data from drug_interactions.json
    # Examples in data: 'neratinib', 'imatinib', 'cisplatin'
    interactions = interactions_data.get(drug_key, [])

    # Compile the final profile
    profile = {
        "name": drug_info.get("name"),
        "molecular_details": {
            "smiles": drug_info.get("smiles"),
            "molecular_weight": drug_info.get("molecular_weight"),
            "log_p": drug_info.get("log_p"),
            "hbd": drug_info.get("hbd"),
            "hba": drug_info.get("hba"),
            "rotatable_bonds": drug_info.get("rotatable_bonds"),
            "approved": drug_info.get("approved"),
            "source": drug_info.get("source")
        },
        "interactions": []
    }

    # Add interaction details such as gene symbol and binding affinity
    for inter in interactions:
        profile["interactions"].append({
            "target_gene": inter.get("gene_symbol"),
            "protein_name": inter.get("protein_name"),
            "binding_affinity": inter.get("binding_affinity"),
            "confidence": inter.get("confidence_score"),
            "type": inter.get("interaction_type")
        })

    return profile

# --- Example Usage ---
# Path to your uploaded files
DRUGS_FILE = 'drugs.json'
INTERACTIONS_FILE = 'drug_interactions.json'

try:
    drugs, interactions = load_data(DRUGS_FILE, INTERACTIONS_FILE)
    
    # Example: Look up "Neratinib"
    # Neratinib targets include KDR, ERBB2, and EGFR
    drug_to_search = "neratinib"
    result = get_comprehensive_drug_profile(drug_to_search, drugs, interactions)
    
    print(json.dumps(result, indent=2))

except FileNotFoundError as e:
    print(f"Error: Ensure {DRUGS_FILE} and {INTERACTIONS_FILE} are in the directory.")
