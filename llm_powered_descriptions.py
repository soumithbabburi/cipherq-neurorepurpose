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

        # Call Gemini API with CORRECT model names (try multiple)
        models_to_try = ["gemini-1.5-flash", "gemini-1.5-pro", "gemini-pro"]
        
        for model in models_to_try:
            try:
                api_url = f"https://generativelanguage.googleapis.com/v1beta/models/{model}:generateContent?key={api_key}"
                
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
                
                if response.status_code == 200:
                    result = response.json()
                    description = result['candidates'][0]['content']['parts'][0]['text']
                    logger.info(f"✅ Gemini description generated using {model} ({len(description)} chars)")
                    return description
                else:
                    logger.warning(f"{model} returned {response.status_code}, trying next...")
                    continue
                    
            except Exception as model_error:
                logger.warning(f"{model} failed: {model_error}, trying next...")
                continue
        
    except Exception as e:
        logger.info(f"Gemini API unavailable: {e}")
    
    # Simple fallback with actual data
    strength = "strong" if binding_affinity < -8 else ("moderate" if binding_affinity < -6 else "weak")
    return f"**{drug_name} → {target_protein}**: Binding affinity of {binding_affinity} kcal/mol indicates {strength} molecular interaction."
