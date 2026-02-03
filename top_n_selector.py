"""
Top N Selector - Selects best drugs from scored list
"""

import logging

logger = logging.getLogger(__name__)


def select_top_n_drugs(scored_drugs: list, n: int = 3) -> list:
    """
    Select top N drugs from scored list
    
    Args:
        scored_drugs: List of drugs with 'repurposing_score'
        n: Number of top drugs to return
    
    Returns:
        List of top N drugs
    """
    if not scored_drugs:
        return []
    
    # Sort by score
    sorted_drugs = sorted(scored_drugs, key=lambda x: x.get('repurposing_score', 0), reverse=True)
    
    # Return top N
    top_n = sorted_drugs[:n]
    
    logger.info(f"Selected top {len(top_n)} drugs from {len(scored_drugs)} total")
    
    return top_n


def get_top_3_drugs(scored_drugs: list) -> list:
    """Convenience function to get top 3"""
    return select_top_n_drugs(scored_drugs, n=3)


__all__ = ['select_top_n_drugs', 'get_top_3_drugs']
