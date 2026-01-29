"""
DATA-DRIVEN DESCRIPTIONS - NO LLM!
Simple templates with protein function lookup
100% reliable, works every time
"""

import logging

logger = logging.getLogger(__name__)

def generate_docking_description_with_llm(drug_name: str, target_protein: str, binding_affinity: float, disease_name: str) -> str:
    """
    Generate description using templates - NO LLM, NO API
    Just protein functions + binding interpretation
    """
    
    # Binding strength
    if binding_affinity < -8:
        strength = "strong"
    elif binding_affinity < -6:
        strength = "moderate"
    else:
        strength = "weak to moderate"
    
    # Protein functions
    functions = {
        'PPARG': 'a nuclear receptor regulating glucose metabolism and inflammation',
        'PPARA': 'a nuclear receptor controlling lipid metabolism',  
        'ABCC8': 'the SUR1 subunit of K-ATP channels regulating insulin secretion',
        'KCNJ11': 'the Kir6.2 subunit of K-ATP channels in pancreatic beta cells',
        'DPP4': 'an enzyme regulating incretin hormones and glucose homeostasis',
        'ACHE': 'an enzyme breaking down acetylcholine in synapses',
        'MAOB': 'an enzyme metabolizing dopamine',
        'DRD2': 'a dopamine receptor controlling motor function'
    }
    
    func = functions.get(target_protein, 'a therapeutic target')
    
    # Disease-specific mechanisms
    mechanisms = {
        'alzheimer': f"{target_protein} modulation affects neuroinflammation and metabolic dysfunction implicated in Alzheimer's pathology",
        'parkinson': f"{target_protein} activity influences dopaminergic neurotransmission central to Parkinson's motor symptoms",
        'diabetes': f"{target_protein} regulation directly impacts glucose homeostasis and insulin sensitivity"
    }
    
    disease_key = 'alzheimer' if 'alzheimer' in disease_name.lower() else (
                  'parkinson' if 'parkinson' in disease_name.lower() else
                  'diabetes' if 'diabetes' in disease_name.lower() else None)
    
    mechanism = mechanisms.get(disease_key, f"{target_protein} represents a potential therapeutic target for {disease_name}")
    
    description = f"{drug_name} exhibits {strength} binding to {target_protein} ({binding_affinity:.2f} kcal/mol). {target_protein} is {func}. {mechanism}."
    
    logger.info(f"✅ Generated description for {drug_name}")
    return description


__all__ = ['generate_docking_description_with_llm']
