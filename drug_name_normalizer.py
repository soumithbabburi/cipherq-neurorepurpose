"""
Drug Name Normalizer
Removes salt/form suffixes to get base drug name
"""

def normalize_drug_name(drug_name: str) -> str:
    """
    Normalize drug name by removing salt/form suffixes
    
    Examples:
        "metformin hydrochloride" → "metformin"
        "atorvastatin calcium" → "atorvastatin"
        "lisinopril, sterile" → "lisinopril"
    """
    
    # Common salt/form suffixes to remove
    suffixes = [
        ' hydrochloride',
        ' sodium',
        ' sulfate',
        ' calcium',
        ' potassium',
        ' maleate',
        ' tartrate',
        ' phosphate',
        ' acetate',
        ' citrate',
        ' mesylate',
        ' fumarate',
        ' succinate',
        ', sterile',
        ' trihydrate',
        ' monohydrate',
        ' anhydrous',
        ' dihydrate'
    ]
    
    normalized = drug_name.lower().strip()
    
    for suffix in suffixes:
        if normalized.endswith(suffix):
            normalized = normalized[:-len(suffix)]
            break  # Only remove one suffix
    
    return normalized.strip()


def get_base_drug_name(drug_name: str) -> str:
    """
    Get base drug name for matching in databases
    Same as normalize but returns it
    """
    return normalize_drug_name(drug_name)


__all__ = ['normalize_drug_name', 'get_base_drug_name']
