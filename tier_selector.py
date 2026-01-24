"""
Tier Selection UI Module
Standalone component for drug database tier selection
"""
import streamlit as st
from typing import Tuple, List

def render_tier_selector() -> Tuple[int, List[str]]:
    """
    Render tier selection UI in sidebar
    Returns: (selected_tier, selected_categories)
    """
    
    st.sidebar.markdown("---")
    st.sidebar.markdown("### 🗂️ Drug Database Tier")
    
    # Tier selection with descriptions
    tier_options = {
        "🟢 Tier 1: FDA-Approved (565 drugs)": 1,
        "🟡 Tier 2: Clinical Trials (20,000 drugs)": 2,
        "🔵 Tier 3: Research Compounds (100,000 drugs)": 3
    }
    
    selected_tier_text = st.sidebar.radio(
        "Select tier:",
        list(tier_options.keys()),
        index=0,  # Default to Tier 1
        label_visibility="collapsed",
        help="Choose drug database tier based on your research needs"
    )
    
    selected_tier = tier_options[selected_tier_text]
    
    # Show tier description
    tier_descriptions = {
        1: "✓ FDA-approved drugs only\n✓ Clinical-ready for repurposing\n✓ Full safety data available",
        2: "✓ FDA + clinical trials\n✓ Investigational compounds\n✓ Research & development",
        3: "✓ Complete research database\n✓ Tool compounds included\n✓ Maximum chemical space"
    }
    
    st.sidebar.caption(tier_descriptions[selected_tier])
    
    # Category filter
    st.sidebar.markdown("### 🎯 Filter by Category")
    
    all_categories = st.sidebar.checkbox("All Categories", value=True)
    
    selected_categories = []
    
    if not all_categories:
        # Show individual category checkboxes
        categories = [
            'Cardiovascular',
            'Neurological',
            'Diabetes',
            'Psychiatric',
            'Anti-inflammatory',
            'Cancer',
            'Pain',
            'Antibiotic',
            'Antiviral'
        ]
        
        for cat in categories:
            if st.sidebar.checkbox(cat, value=False):
                selected_categories.append(cat)
        
        if selected_categories:
            st.sidebar.success(f"✅ {len(selected_categories)} selected")
    else:
        selected_categories = []  # Empty list means all categories
    
    return selected_tier, selected_categories


def get_tier_filtered_drugs(tier: int, categories: List[str] = None):
    """
    Query drugs filtered by tier and optionally by categories
    """
    try:
        import psycopg2
        from psycopg2.extras import RealDictCursor
        
        conn = psycopg2.connect(
            host="localhost",
            database="cipherq_repurpose",
            user="babburisoumith",
            password=""
        )
        
        cur = conn.cursor(cursor_factory=RealDictCursor)
        
        if categories and len(categories) > 0:
            # Filter by tier AND categories
            cur.execute("""
                SELECT name, therapeutic_category, drug_class, original_indication, 
                       fda_status, database_tier
                FROM drugs
                WHERE database_tier <= %s
                  AND therapeutic_category = ANY(%s)
                ORDER BY name
            """, (tier, categories))
        else:
            # Filter by tier only (all categories)
            cur.execute("""
                SELECT name, therapeutic_category, drug_class, original_indication,
                       fda_status, database_tier
                FROM drugs
                WHERE database_tier <= %s
                ORDER BY name
            """, (tier,))
        
        results = cur.fetchall()
        cur.close()
        conn.close()
        
        return [dict(row) for row in results]
        
    except Exception as e:
        st.error(f"Query failed: {e}")
        return []


__all__ = ['render_tier_selector', 'get_tier_filtered_drugs']