"""
LLM-Powered Descriptions - GEMINI FLASH (FREE!)
Uses Google's Gemini 1.5 Flash for scientific descriptions
"""

import os
import logging

logger = logging.getLogger(__name__)

def generate_docking_description_with_llm(drug_name: str, target_protein: str, binding_affinity: float, disease_name: str) -> str:
    """
    Generate scientific description using Gemini 1.5 Flash (FREE!)
    
    Args:
        drug_name: Name of the drug
        target_protein: Target protein gene symbol
        binding_affinity: Binding affinity in kcal/mol
        disease_name: Target disease
    
    Returns:
        Scientific description string
    """
    
    try:
        import google.generativeai as genai
        
        # Configure Gemini API
        api_key = os.getenv('GEMINI_API_KEY')
        if not api_key:
            logger.warning("GEMINI_API_KEY not found in environment")
            raise Exception("No API key")
        
        genai.configure(api_key=api_key)
        
        # Use Gemini 1.5 Flash (free tier)
        model = genai.GenerativeModel('gemini-1.5-flash')
        
        # Create prompt
        prompt = f"""Analyze this molecular docking result for drug repurposing:

Drug: {drug_name}
Target Protein: {target_protein}
Binding Affinity: {binding_affinity} kcal/mol
Disease Context: {disease_name}

Provide a 2-3 sentence scientific explanation covering:
1. What the binding affinity indicates about the interaction strength
2. What {target_protein} does biologically
3. How this mechanism could be relevant to {disease_name}

Be concise, scientific, and focus on the therapeutic rationale."""
        
        # Generate response
        response = model.generate_content(prompt)
        
        if response and response.text:
            description = response.text.strip()
            logger.info(f"✅ Gemini description generated for {drug_name}")
            return description
        else:
            raise Exception("Empty response from Gemini")
        
    except ImportError:
        logger.error("google-generativeai not installed. Run: pip install google-generativeai")
        return _fallback_description(drug_name, target_protein, binding_affinity, disease_name)
    
    except Exception as e:
        logger.warning(f"Gemini description failed: {e}")
        return _fallback_description(drug_name, target_protein, binding_affinity, disease_name)


def _fallback_description(drug_name: str, target_protein: str, binding_affinity: float, disease_name: str) -> str:
    """Fallback description when Gemini unavailable"""
    
    # Determine binding strength
    if binding_affinity < -8:
        strength = "strong"
        interpretation = "indicates high binding affinity and stable complex formation"
    elif binding_affinity < -6:
        strength = "moderate"
        interpretation = "suggests favorable binding interactions"
    else:
        strength = "weak"
        interpretation = "indicates modest binding potential"
    
    # Basic protein function lookup
    protein_functions = {
        'PPARG': 'a nuclear receptor regulating glucose metabolism and inflammation',
        'PPARA': 'a nuclear receptor controlling lipid metabolism',
        'DPP4': 'an enzyme regulating incretin hormones and glucose homeostasis',
        'ACHE': 'an enzyme breaking down acetylcholine in synapses',
        'BCHE': 'an enzyme involved in neurotransmitter metabolism',
        'MAOB': 'an enzyme metabolizing dopamine and other monoamines',
        'DRD2': 'a dopamine receptor involved in motor control and reward',
        'DRD3': 'a dopamine receptor in the mesolimbic pathway',
        'HMGCR': 'the rate-limiting enzyme in cholesterol biosynthesis',
        'ACE': 'an enzyme regulating blood pressure via angiotensin conversion',
        'ABCC8': 'a potassium channel subunit regulating insulin secretion',
        'KCNJ11': 'a potassium channel involved in insulin release',
        'GRIN1': 'an NMDA receptor subunit involved in synaptic plasticity'
    }
    
    protein_func = protein_functions.get(target_protein, 'a therapeutic protein target')
    
    return f"{drug_name} exhibits {strength} binding to {target_protein} ({binding_affinity:.2f} kcal/mol), which {interpretation}. {target_protein} is {protein_func}, making this interaction potentially relevant to {disease_name} pathophysiology."


__all__ = ['generate_docking_description_with_llm']
