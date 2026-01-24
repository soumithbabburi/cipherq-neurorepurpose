def generate_docking_description_with_llm(drug_name, target_protein, binding_affinity, disease_name):
    """
    Generate scientific description of docking results using Gemini API.
    NO HARDCODING - uses actual docking data.
    """
    import os
    import logging
    import requests
    
    logger = logging.getLogger(__name__)
    
    try:
        api_key = os.getenv('GEMINI_API_KEY')
        
        if not api_key:
            raise ValueError("No Gemini API key")
        
        prompt = f"""Explain the molecular docking result in 2-3 scientific sentences:

Drug: {drug_name}
Target Protein: {target_protein}
Binding Affinity: {binding_affinity} kcal/mol
Disease Context: {disease_name}

Explain:
1. What this binding affinity means (more negative = stronger binding)
2. How {drug_name} binding to {target_protein} could be therapeutic for {disease_name}
3. The biological significance

Be factual, concise, and scientific. Use the actual protein name and binding value."""

        # Call Gemini API with CORRECT model name
        api_url = f"https://generativelanguage.googleapis.com/v1beta/models/gemini-pro:generateContent?key={api_key}"
        
        response = requests.post(
            api_url,
            headers={"Content-Type": "application/json"},
            json={
                "contents": [{
                    "parts": [{"text": prompt}]
                }],
                "generationConfig": {
                    "temperature": 0.7,
                    "maxOutputTokens": 800
                }
            },
            timeout=30
        )
        
        logger.info(f"Gemini API call: {response.status_code}")
        
        if response.status_code == 200:
            result = response.json()
            description = result['candidates'][0]['content']['parts'][0]['text']
            logger.info(f"✅ Gemini description generated ({len(description)} chars)")
            return description
        else:
            logger.error(f"Gemini error {response.status_code}: {response.text[:300]}")
        
    except Exception as e:
        logger.info(f"Gemini API unavailable: {e}")
    
    # Simple fallback with actual data
    strength = "strong" if binding_affinity < -8 else ("moderate" if binding_affinity < -6 else "weak")
    return f"**{drug_name} → {target_protein}**: Binding affinity of {binding_affinity} kcal/mol indicates {strength} molecular interaction."