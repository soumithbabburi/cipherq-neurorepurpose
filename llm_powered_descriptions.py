"""
LLM-Powered Docking Descriptions using Groq API
Generates scientific descriptions of molecular docking results
"""
import os
import logging
from groq import Groq

logger = logging.getLogger(__name__)

def generate_docking_description_with_llm(drug_name, target_protein, binding_affinity, disease_name):
    """
    Generate scientific description of docking results using Groq API.
    
    Args:
        drug_name: Name of the drug
        target_protein: Target protein (e.g., PPARG, DPP4)
        binding_affinity: Binding affinity in kcal/mol
        disease_name: Target disease
    
    Returns:
        Scientific description string
    """
    
    try:
        # Get Groq API key
        api_key = os.getenv('GROQ_API_KEY')
        
        if not api_key:
            raise ValueError("GROQ_API_KEY not found in environment")
        
        # Initialize Groq client
        client = Groq(api_key=api_key)
        
        # Create prompt
        prompt = f"""Explain this molecular docking result in 2-3 scientific sentences:

Drug: {drug_name}
Target Protein: {target_protein}
Binding Affinity: {binding_affinity} kcal/mol
Disease Context: {disease_name}

Explain:
1. What this binding affinity means (more negative = stronger binding)
2. How {drug_name} binding to {target_protein} could be therapeutic for {disease_name}
3. The biological significance

Be factual, concise, and scientific. Use the actual protein name and binding value."""
        
        # Call Groq API
        logger.info(f"Calling Groq API for docking description...")
        response = client.chat.completions.create(
            model="mixtral-8x7b-32768",
            messages=[{
                "role": "user",
                "content": prompt
            }],
            max_tokens=400,
            temperature=0.7
        )
        
        description = response.choices[0].message.content.strip()
        logger.info(f"✅ Groq description generated ({len(description)} chars)")
        
        return description
        
    except Exception as e:
        logger.warning(f"Groq API failed: {e}")
        
        # Fallback to simple description
        strength = "strong" if binding_affinity < -8 else ("moderate" if binding_affinity < -6 else "weak")
        return f"**{drug_name} → {target_protein}**: Binding affinity of {binding_affinity} kcal/mol indicates {strength} molecular interaction. This {strength} binding suggests potential therapeutic activity for {disease_name}."
