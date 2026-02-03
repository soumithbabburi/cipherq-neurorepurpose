"""
SMART TOP N DRUG SELECTOR
Selects the top N most relevant drugs from a scored list
"""

import logging

logger = logging.getLogger(__name__)

def select_top_n_drugs(scored_drugs, n=3, min_score=0.0):
    """
    Select top N drugs from a scored list
    
    Args:
        scored_drugs: List of drug dicts with 'connection_score' field
        n: Number of drugs to select (default 3)
        min_score: Minimum score required (default 0.0)
    
    Returns:
        List of top N drugs
    """
    if not scored_drugs:
        logger.warning("⚠️ Empty drug list provided")
        return []
    
    # Filter by min_score
    qualifying_drugs = [d for d in scored_drugs if d.get('connection_score', 0) >= min_score]
    
    if not qualifying_drugs:
        logger.warning(f"⚠️ No drugs meet minimum score threshold ({min_score})")
        # Return top N anyway, even if below threshold
        qualifying_drugs = scored_drugs
    
    # Already sorted by disease_connection_filter, but ensure it
    qualifying_drugs.sort(key=lambda x: x.get('connection_score', 0), reverse=True)
    
    # Select top N
    top_n = qualifying_drugs[:n]
    
    logger.info(f"✅ SELECTED TOP {len(top_n)} DRUGS:")
    for i, drug in enumerate(top_n, 1):
        score = drug.get('connection_score', 0)
        name = drug.get('drug_name', 'Unknown')
        pathways = drug.get('matched_pathways', [])
        logger.info(f"   {i}. {name} (score: {score:.1f}, pathways: {len(pathways)})")
    
    return top_n

def select_top_drugs_with_diversity(scored_drugs, n=3, diversity_factor=0.3):
    """
    Select top N drugs with mechanism diversity
    Ensures selected drugs work through different pathways
    
    Args:
        scored_drugs: List of scored drugs
        n: Number to select
        diversity_factor: Weight for diversity (0-1, higher = more diverse)
    
    Returns:
        List of top N diverse drugs
    """
    if not scored_drugs or len(scored_drugs) <= n:
        return scored_drugs[:n]
    
    selected = []
    remaining = scored_drugs.copy()
    
    # Always pick the top scorer first
    selected.append(remaining.pop(0))
    
    # For remaining slots, balance score and diversity
    while len(selected) < n and remaining:
        best_candidate = None
        best_combined_score = -1
        
        for i, candidate in enumerate(remaining):
            # Base score
            base_score = candidate.get('connection_score', 0)
            
            # Diversity score: how different from already selected drugs
            diversity_score = calculate_pathway_diversity(candidate, selected)
            
            # Combined score
            combined = (1 - diversity_factor) * base_score + diversity_factor * diversity_score * 100
            
            if combined > best_combined_score:
                best_combined_score = combined
                best_candidate = (i, candidate)
        
        if best_candidate:
            idx, drug = best_candidate
            selected.append(drug)
            remaining.pop(idx)
    
    logger.info(f"✅ SELECTED {len(selected)} DIVERSE DRUGS:")
    for i, drug in enumerate(selected, 1):
        score = drug.get('connection_score', 0)
        name = drug.get('drug_name', 'Unknown')
        pathways = drug.get('matched_pathways', [])
        logger.info(f"   {i}. {name} (score: {score:.1f}, pathways: {pathways})")
    
    return selected

def calculate_pathway_diversity(candidate, selected_drugs):
    """
    Calculate how different a candidate's pathways are from selected drugs
    
    Returns:
        Float between 0 (identical) and 1 (completely different)
    """
    if not selected_drugs:
        return 1.0
    
    candidate_pathways = set(candidate.get('matched_pathways', []))
    
    if not candidate_pathways:
        return 0.5  # Unknown, give medium diversity
    
    # Calculate average Jaccard distance to selected drugs
    distances = []
    for selected in selected_drugs:
        selected_pathways = set(selected.get('matched_pathways', []))
        
        if not selected_pathways:
            distances.append(0.5)
            continue
        
        # Jaccard distance = 1 - (intersection / union)
        intersection = len(candidate_pathways & selected_pathways)
        union = len(candidate_pathways | selected_pathways)
        
        if union == 0:
            distance = 0.5
        else:
            distance = 1 - (intersection / union)
        
        distances.append(distance)
    
    # Average distance
    return sum(distances) / len(distances)

def format_selection_summary(top_drugs, target_disease, source_category=None):
    """
    Create a human-readable summary of selected drugs
    
    Returns:
        String with formatted summary
    """
    if not top_drugs:
        return "No drugs selected."
    
    summary = []
    summary.append(f"🎯 TOP {len(top_drugs)} DRUGS")
    
    if source_category:
        summary.append(f"   Category: {source_category} → Target: {target_disease}")
    else:
        summary.append(f"   Target: {target_disease}")
    
    summary.append("")
    
    for i, drug in enumerate(top_drugs, 1):
        name = drug.get('drug_name', 'Unknown')
        score = drug.get('connection_score', 0)
        pathways = drug.get('matched_pathways', [])
        targets = drug.get('num_targets', 0)
        
        summary.append(f"{i}. {name}")
        summary.append(f"   Score: {score:.1f}/100")
        summary.append(f"   Targets: {targets} genes")
        
        if pathways:
            pathway_str = ', '.join(pathways[:3])
            if len(pathways) > 3:
                pathway_str += f" (+{len(pathways)-3} more)"
            summary.append(f"   Pathways: {pathway_str}")
        
        summary.append("")
    
    return "\n".join(summary)

# Convenience function for common use case
def get_top_3_drugs(all_drugs, target_disease, source_category=None, use_diversity=False):
    """
    One-line function to get top 3 drugs
    
    Args:
        all_drugs: List of drug dicts (will be scored)
        target_disease: Target disease name
        source_category: Source category (optional)
        use_diversity: Use diversity selection (default False)
    
    Returns:
        Tuple: (top_3_drugs, all_scored_drugs)
    """
    from disease_connection_filter import filter_drugs_by_disease_connection
    
    # Score all drugs (min_score=0 means return all)
    all_scored = filter_drugs_by_disease_connection(
        all_drugs,
        target_disease=target_disease,
        source_category=source_category,
        min_score=0.0,  # Return all with scores
        auto_enrich=True
    )
    
    # Select top 3
    if use_diversity:
        top_3 = select_top_drugs_with_diversity(all_scored, n=3)
    else:
        top_3 = select_top_n_drugs(all_scored, n=3)
    
    return top_3, all_scored

# Testing
if __name__ == "__main__":
    print("Testing top N selector...")
    
    # Mock data
    mock_drugs = [
        {'drug_name': 'DrugA', 'connection_score': 85, 'matched_pathways': ['metabolism', 'inflammation']},
        {'drug_name': 'DrugB', 'connection_score': 72, 'matched_pathways': ['metabolism', 'autophagy']},
        {'drug_name': 'DrugC', 'connection_score': 68, 'matched_pathways': ['synaptic', 'neurotransmitter']},
        {'drug_name': 'DrugD', 'connection_score': 65, 'matched_pathways': ['inflammation', 'oxidative_stress']},
        {'drug_name': 'DrugE', 'connection_score': 45, 'matched_pathways': ['metabolism']},
    ]
    
    print("\n=== STANDARD TOP 3 ===")
    top_3 = select_top_n_drugs(mock_drugs, n=3)
    print(format_selection_summary(top_3, 'Alzheimers', 'Cardiovascular'))
    
    print("\n=== DIVERSE TOP 3 ===")
    top_3_diverse = select_top_drugs_with_diversity(mock_drugs, n=3, diversity_factor=0.3)
    print(format_selection_summary(top_3_diverse, 'Alzheimers', 'Cardiovascular'))
