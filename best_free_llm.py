"""
BEST FREE LLM Integration
1. Gemini 1.5 Flash (FREE, fast, smart)
2. Groq (FREE, blazing fast)
"""
import os
import requests
import logging

logger = logging.getLogger(__name__)

def generate_docking_description_BEST_LLM(drug_name, target_protein, binding_affinity, disease_name):
    """
    Generate description using BEST free LLMs:
    1. Gemini 1.5 Flash (best quality, free)
    2. Groq Llama 3.3 70B (fastest, free)
    """
    
    prompt = f"""Explain this molecular docking result in 2-3 detailed scientific sentences:

**Drug:** {drug_name}
**Target Protein:** {target_protein}
**Binding Affinity:** {binding_affinity} kcal/mol (more negative = stronger)
**Disease:** {disease_name}

Explain:
1. What this binding affinity means and why it's significant
2. The specific biological mechanism: HOW does {drug_name} interact with {target_protein}?
3. Why this interaction could be therapeutic for {disease_name} - what disease processes does it affect?

Be specific about molecular mechanisms, pathways, and therapeutic effects. Use scientific terminology."""

    # === METHOD 1: Gemini 1.5 Flash (BEST!) ===
    try:
        gemini_key = os.getenv('GEMINI_API_KEY')
        
        if gemini_key:
            # Try Gemini 1.5 Flash first (best free model!)
            response = requests.post(
                f"https://generativelanguage.googleapis.com/v1beta/models/gemini-1.5-flash:generateContent?key={gemini_key}",
                headers={"Content-Type": "application/json"},
                json={
                    "contents": [{"parts": [{"text": prompt}]}],
                    "generationConfig": {
                        "temperature": 0.7,
                        "maxOutputTokens": 1024,
                        "topP": 0.9
                    }
                },
                timeout=20
            )
            
            if response.status_code == 200:
                result = response.json()
                description = result['candidates'][0]['content']['parts'][0]['text']
                logger.info("✅ Description from Gemini 1.5 Flash (BEST!)")
                return description
            else:
                logger.warning(f"Gemini 1.5 Flash: {response.status_code}")
    except Exception as e:
        logger.warning(f"Gemini 1.5 Flash failed: {e}")
    
    # === METHOD 2: Groq (FASTEST!) ===
    try:
        groq_key = os.getenv('GROQ_API_KEY')
        
        if groq_key:
            response = requests.post(
                "https://api.groq.com/openai/v1/chat/completions",
                headers={
                    "Authorization": f"Bearer {groq_key}",
                    "Content-Type": "application/json"
                },
                json={
                    "model": "llama-3.3-70b-versatile",
                    "messages": [{"role": "user", "content": prompt}],
                    "temperature": 0.7,
                    "max_tokens": 1024
                },
                timeout=15
            )
            
            if response.status_code == 200:
                result = response.json()
                description = result['choices'][0]['message']['content']
                logger.info("✅ Description from Groq Llama 3.3 70B (FASTEST!)")
                return description
            else:
                logger.warning(f"Groq: {response.status_code}")
    except Exception as e:
        logger.warning(f"Groq failed: {e}")
    
    # === FALLBACK: Basic description ===
    strength = "strong" if binding_affinity < -8 else ("moderate" if binding_affinity < -6 else "weak")
    return f"""**{drug_name} → {target_protein}**

The binding affinity of **{binding_affinity} kcal/mol** indicates **{strength} molecular interaction**. This suggests that {drug_name} can effectively bind to {target_protein}, which may modulate its activity and potentially provide therapeutic benefit for {disease_name}. More negative values indicate stronger, more stable binding.

*Get detailed AI analysis by adding GEMINI_API_KEY (free at https://makersuite.google.com/app/apikey) or GROQ_API_KEY (free at https://console.groq.com)*"""


# Simple wrapper for compatibility
def generate_docking_description_with_llm(drug_name, target_protein, binding_affinity, disease_name):
    """Wrapper for backward compatibility"""
    return generate_docking_description_BEST_LLM(drug_name, target_protein, binding_affinity, disease_name)
