"""
CipherQ NeuroRepurpose - Drug Repurposing Platform
Complete working version with fixed categorization and top N selection
"""

import streamlit as st
import pandas as pd
import json
import logging
from pathlib import Path

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s:%(name)s:%(message)s')
logger = logging.getLogger(__name__)

# ============================================================================
# IMPORTS - All modules with error handling
# ============================================================================

# Core imports
try:
    from services.drug_categorizer import get_drugs_by_category
    CATEGORIZER_AVAILABLE = True
    logger.info("✅ Drug categorizer loaded")
except Exception as e:
    logger.warning(f"⚠️ Drug categorizer not available: {e}")
    CATEGORIZER_AVAILABLE = False

# Top N selector (NEW)
try:
    from top_n_selector import get_top_3_drugs, format_selection_summary, select_top_n_drugs
    TOP_N_SELECTOR_AVAILABLE = True
    logger.info("✅ Top N selector loaded")
except Exception as e:
    logger.warning(f"⚠️ Top N selector not available: {e}")
    TOP_N_SELECTOR_AVAILABLE = False

# Disease connection filter (UPDATED)
try:
    from disease_connection_filter import filter_drugs_by_disease_connection
    FILTER_AVAILABLE = True
    logger.info("✅ Disease connection filter loaded")
except Exception as e:
    logger.warning(f"⚠️ Disease connection filter not available: {e}")
    FILTER_AVAILABLE = False

# Evidence graph builder
try:
    from evidence_graph_builder import build_evidence_graph
    EVIDENCE_GRAPH_AVAILABLE = True
    logger.info("✅ Evidence graph builder loaded")
except Exception as e:
    logger.warning(f"⚠️ Evidence graph builder not available: {e}")
    EVIDENCE_GRAPH_AVAILABLE = False

# ============================================================================
# PAGE CONFIG
# ============================================================================

st.set_page_config(
    page_title="CipherQ Drug Repurposing",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# ============================================================================
# CUSTOM CSS
# ============================================================================

st.markdown("""
<style>
    .main-header {
        font-size: 2.5rem;
        font-weight: 700;
        color: #1f77b4;
        margin-bottom: 0.5rem;
    }
    .sub-header {
        font-size: 1.2rem;
        color: #666;
        margin-bottom: 2rem;
    }
    .metric-card {
        background-color: #f0f2f6;
        padding: 1rem;
        border-radius: 0.5rem;
        border-left: 4px solid #1f77b4;
    }
    .success-box {
        background-color: #d4edda;
        border: 1px solid #c3e6cb;
        border-radius: 0.25rem;
        padding: 1rem;
        margin: 1rem 0;
    }
    .info-box {
        background-color: #d1ecf1;
        border: 1px solid #bee5eb;
        border-radius: 0.25rem;
        padding: 1rem;
        margin: 1rem 0;
    }
</style>
""", unsafe_allow_html=True)

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def check_dependencies():
    """Check if all required modules are available"""
    deps = {
        "Drug Categorizer": CATEGORIZER_AVAILABLE,
        "Top N Selector": TOP_N_SELECTOR_AVAILABLE,
        "Disease Filter": FILTER_AVAILABLE,
        "Evidence Graph": EVIDENCE_GRAPH_AVAILABLE
    }
    
    missing = [name for name, available in deps.items() if not available]
    
    if missing:
        st.warning(f"⚠️ Some modules are not available: {', '.join(missing)}")
        return False
    return True

def format_drug_card(drug, rank):
    """Format a drug as a card display"""
    name = drug.get('drug_name', 'Unknown')
    score = drug.get('connection_score', 0)
    targets = drug.get('num_targets', 0)
    pathways = drug.get('matched_pathways', [])
    
    pathway_str = ', '.join(pathways[:3])
    if len(pathways) > 3:
        pathway_str += f" (+{len(pathways)-3} more)"
    
    return f"""
    **{rank}. {name}**
    - Score: {score:.1f}/100
    - Targets: {targets} genes
    - Pathways: {pathway_str if pathway_str else 'None'}
    """

# ============================================================================
# SIDEBAR
# ============================================================================

with st.sidebar:
    st.image("https://via.placeholder.com/150x50/1f77b4/ffffff?text=CipherQ", use_container_width=True)
    
    st.markdown("---")
    
    st.header("📖 About")
    st.write("""
    CipherQ NeuroRepurpose identifies promising drug repurposing candidates 
    using pathway-based analysis and evidence networks.
    """)
    
    st.markdown("---")
    
    st.header("🔧 System Status")
    deps_status = {
        "Categorizer": CATEGORIZER_AVAILABLE,
        "Selector": TOP_N_SELECTOR_AVAILABLE,
        "Filter": FILTER_AVAILABLE,
        "Evidence Graph": EVIDENCE_GRAPH_AVAILABLE
    }
    
    for name, status in deps_status.items():
        icon = "✅" if status else "❌"
        st.write(f"{icon} {name}")
    
    st.markdown("---")
    
    st.header("📊 Quick Stats")
    try:
        # Load drug categories
        with open('drug_therapeutic_categories.json', 'r') as f:
            categories = json.load(f)
        
        category_counts = {}
        for drug, cats in categories.items():
            for cat in cats:
                category_counts[cat] = category_counts.get(cat, 0) + 1
        
        st.metric("Total Drugs", len(categories))
        st.metric("Categories", len(category_counts))
        
    except Exception as e:
        st.write("Stats unavailable")
    
    st.markdown("---")
    
    with st.expander("ℹ️ Help & Tips"):
        st.write("""
        **How to use:**
        1. Select a source drug category
        2. Choose a target disease
        3. Click 'Analyze Drugs'
        4. Review top candidates
        5. Explore evidence network
        
        **Tips:**
        - Enable diversity for varied mechanisms
        - Check 'View All Scores' for full ranking
        - Download results as CSV
        """)

# ============================================================================
# MAIN APP
# ============================================================================

def main():
    """Main application logic"""
    
    # Header
    st.markdown('<h1 class="main-header">🧬 CipherQ Drug Repurposing</h1>', unsafe_allow_html=True)
    st.markdown('<p class="sub-header">Discover novel therapeutic applications for existing drugs</p>', unsafe_allow_html=True)
    
    # Check dependencies
    if not check_dependencies():
        st.error("""
        ❌ **Missing Required Modules**
        
        Please ensure all required files are uploaded to your GitHub repository:
        - top_n_selector.py
        - disease_connection_filter.py
        - services/drug_categorizer.py
        - evidence_graph_builder.py
        """)
        return
    
    st.markdown("---")
    
    # ========================================================================
    # SELECTION INTERFACE
    # ========================================================================
    
    st.subheader("🎯 Drug Analysis Configuration")
    
    col1, col2 = st.columns(2)
    
    with col1:
        selected_category = st.selectbox(
            "📦 Source Drug Category",
            [
                "Diabetic",
                "Cardiovascular", 
                "Parkinsons",
                "Alzheimers",
                "Pain",
                "Psychiatric",
                "Antibiotics",
                "Antivirals",
                "Antifungals",
                "Cancer",
                "Respiratory",
                "Gastrointestinal",
                "Hormones",
                "Dermatological",
                "Ophthalmological"
            ],
            help="Select the category of drugs to analyze for repurposing"
        )
    
    with col2:
        target_disease = st.selectbox(
            "🎯 Target Disease",
            [
                "Alzheimers",
                "Parkinsons",
                "Cardiovascular Disease",
                "Type 2 Diabetes",
                "Cancer",
                "Depression",
                "Anxiety",
                "Chronic Pain",
                "Inflammatory Disorders"
            ],
            help="Select the disease you want to find drug candidates for"
        )
    
    # Advanced options
    with st.expander("⚙️ Advanced Options"):
        col_a, col_b, col_c = st.columns(3)
        
        with col_a:
            num_candidates = st.slider(
                "Number of top drugs",
                min_value=1,
                max_value=10,
                value=3,
                help="How many top candidates to select"
            )
        
        with col_b:
            use_diversity = st.checkbox(
                "Use pathway diversity",
                value=False,
                help="Select drugs with diverse mechanisms of action"
            )
        
        with col_c:
            show_all_scores = st.checkbox(
                "Show all scores",
                value=True,
                help="Display scoring for all analyzed drugs"
            )
    
    # Analyze button
    analyze_button = st.button(
        "🔍 Analyze Drugs",
        type="primary",
        use_container_width=True
    )
    
    st.markdown("---")
    
    # ========================================================================
    # ANALYSIS EXECUTION
    # ========================================================================
    
    if analyze_button:
        
        # Progress tracking
        progress_bar = st.progress(0)
        status_text = st.empty()
        
        try:
            # ================================================================
            # STEP 1: Retrieve drugs from category
            # ================================================================
            status_text.text(f"📦 Retrieving {selected_category} drugs...")
            progress_bar.progress(20)
            
            all_drugs = get_drugs_by_category(selected_category)
            
            if not all_drugs:
                st.error(f"❌ No drugs found in {selected_category} category")
                st.info("This could mean the category file is not loaded correctly.")
                return
            
            logger.info(f"Retrieved {len(all_drugs)} drugs from {selected_category}")
            
            # ================================================================
            # STEP 2: Score all drugs and select top N
            # ================================================================
            status_text.text(f"⚖️ Scoring {len(all_drugs)} drugs for {target_disease} relevance...")
            progress_bar.progress(40)
            
            # Use the top_n_selector
            if num_candidates == 3 and not use_diversity:
                # Optimized path for common case
                top_drugs, all_scored = get_top_3_drugs(
                    all_drugs=all_drugs,
                    target_disease=target_disease,
                    source_category=selected_category,
                    use_diversity=use_diversity
                )
            else:
                # Custom selection
                all_scored = filter_drugs_by_disease_connection(
                    all_drugs,
                    target_disease=target_disease,
                    source_category=selected_category,
                    min_score=0.0,  # Return all with scores
                    auto_enrich=True
                )
                
                if use_diversity:
                    from top_n_selector import select_top_drugs_with_diversity
                    top_drugs = select_top_drugs_with_diversity(all_scored, n=num_candidates)
                else:
                    top_drugs = select_top_n_drugs(all_scored, n=num_candidates)
            
            progress_bar.progress(60)
            
            # ================================================================
            # STEP 3: Validate results
            # ================================================================
            if not top_drugs:
                st.warning("⚠️ No suitable drug candidates found")
                
                st.info("""
                **Possible reasons:**
                - Drugs may not have target data available
                - No biological connection to target disease detected
                - Try a different category or disease combination
                """)
                
                # Show diagnostic info
                with st.expander("🔍 Diagnostic Information"):
                    drugs_with_targets = len([d for d in all_scored if d.get('num_targets', 0) > 0])
                    drugs_with_score = len([d for d in all_scored if d.get('connection_score', 0) > 0])
                    
                    st.write(f"- Total drugs analyzed: {len(all_drugs)}")
                    st.write(f"- Drugs with target data: {drugs_with_targets}")
                    st.write(f"- Drugs with relevance score > 0: {drugs_with_score}")
                    
                    if drugs_with_targets < len(all_drugs) * 0.5:
                        st.warning("⚠️ Many drugs are missing target data. Check if drug_interactions.json is properly loaded.")
                
                return
            
            status_text.text(f"✅ Analysis complete!")
            progress_bar.progress(100)
            
            # Clear progress indicators
            import time
            time.sleep(0.5)
            progress_bar.empty()
            status_text.empty()
            
            # ================================================================
            # DISPLAY RESULTS
            # ================================================================
            
            st.success(f"✅ Found {len(top_drugs)} top candidate(s)!")
            
            st.markdown("---")
            
            # ================================================================
            # Summary Cards
            # ================================================================
            st.subheader("🎯 Top Candidates Summary")
            
            summary = format_selection_summary(top_drugs, target_disease, selected_category)
            st.code(summary, language=None)
            
            # Statistics
            st.subheader("📊 Analysis Statistics")
            
            col_stat1, col_stat2, col_stat3, col_stat4 = st.columns(4)
            
            with col_stat1:
                st.metric("Total Analyzed", len(all_drugs))
            
            with col_stat2:
                drugs_with_targets = len([d for d in all_scored if d.get('num_targets', 0) > 0])
                st.metric("With Targets", drugs_with_targets)
            
            with col_stat3:
                drugs_with_score = len([d for d in all_scored if d.get('connection_score', 0) > 0])
                st.metric("Scored > 0", drugs_with_score)
            
            with col_stat4:
                st.metric("Top Selected", len(top_drugs))
            
            st.markdown("---")
            
            # ================================================================
            # Evidence Network Graph
            # ================================================================
            st.subheader("🔬 Evidence Network Graph")
            st.write(f"Drug-target-pathway connections for top {len(top_drugs)} candidate(s)")
            
            if EVIDENCE_GRAPH_AVAILABLE:
                try:
                    build_evidence_graph(top_drugs, target_disease)
                except Exception as e:
                    st.error(f"❌ Error building evidence graph: {e}")
                    logger.exception("Evidence graph error")
                    
                    st.write("**Debug Information:**")
                    for i, drug in enumerate(top_drugs, 1):
                        with st.expander(f"Drug {i} data"):
                            st.json(drug)
            else:
                st.warning("Evidence graph builder not available")
            
            st.markdown("---")
            
            # ================================================================
            # Detailed Drug Information
            # ================================================================
            st.subheader("📋 Detailed Analysis")
            
            for i, drug in enumerate(top_drugs, 1):
                with st.expander(f"**{i}. {drug.get('drug_name', 'Unknown')}** - Complete Analysis"):
                    
                    # Two-column layout
                    detail_col1, detail_col2 = st.columns(2)
                    
                    with detail_col1:
                        st.markdown("**🎯 Scoring Breakdown**")
                        st.write(f"- **Total Score:** {drug.get('connection_score', 0):.1f} / 100")
                        st.write(f"- Pathway Score: {drug.get('pathway_score', 0):.1f}")
                        st.write(f"- Cross-Disease Bonus: {drug.get('cross_disease_score', 0):.1f}")
                        
                        st.markdown("**🧬 Mechanism of Action**")
                        st.write(f"- **Targets:** {drug.get('num_targets', 0)} genes")
                        
                        pathways = drug.get('matched_pathways', [])
                        if pathways:
                            st.write(f"- **Pathways:** {', '.join(pathways)}")
                        else:
                            st.write("- **Pathways:** None identified")
                    
                    with detail_col2:
                        st.markdown("**🔍 Target Details**")
                        details = drug.get('target_details', [])
                        if details:
                            for detail in details:
                                st.write(f"• {detail}")
                        else:
                            st.write("_No detailed target information available_")
                        
                        # Cross-disease pathways
                        cross_pathways = drug.get('cross_disease_pathways', [])
                        if cross_pathways:
                            st.markdown("**🔗 Cross-Disease Connections**")
                            st.write(f"Connected via: {', '.join(cross_pathways)}")
            
            st.markdown("---")
            
            # ================================================================
            # All Scored Drugs
            # ================================================================
            if show_all_scores:
                with st.expander("📊 View All Scored Drugs", expanded=False):
                    st.write(f"Complete ranking of all {len(all_scored)} analyzed drugs:")
                    
                    # Create DataFrame
                    df_data = []
                    for i, drug in enumerate(all_scored, 1):
                        score = drug.get('connection_score', 0)
                        name = drug.get('drug_name', 'Unknown')
                        targets = drug.get('num_targets', 0)
                        pathways = drug.get('matched_pathways', [])
                        
                        pathway_str = ', '.join(pathways[:3])
                        if len(pathways) > 3:
                            pathway_str += f" (+{len(pathways)-3})"
                        
                        df_data.append({
                            'Rank': i,
                            'Drug Name': name,
                            'Score': f"{score:.1f}",
                            'Targets': targets,
                            'Pathways': pathway_str if pathway_str else 'None'
                        })
                    
                    df = pd.DataFrame(df_data)
                    
                    # Display with styling
                    st.dataframe(
                        df,
                        use_container_width=True,
                        hide_index=True,
                        height=400
                    )
                    
                    # Download button
                    csv = df.to_csv(index=False)
                    st.download_button(
                        label="📥 Download Complete Ranking (CSV)",
                        data=csv,
                        file_name=f"{selected_category}_to_{target_disease}_complete_scores.csv",
                        mime="text/csv",
                        use_container_width=True
                    )
        
        except Exception as e:
            st.error(f"❌ An error occurred during analysis: {str(e)}")
            logger.exception("Analysis error")
            
            with st.expander("🐛 Error Details"):
                st.code(str(e))
    
    else:
        # ====================================================================
        # Initial state - show instructions
        # ====================================================================
        
        st.info("""
        👈 **Get started:** Select criteria above and click **'Analyze Drugs'**
        
        **What this tool does:**
        1. Retrieves drugs from your selected category
        2. Scores each drug based on biological relevance to target disease
        3. Selects the top candidates using pathway analysis
        4. Builds an evidence network showing molecular connections
        5. Provides detailed analysis and downloadable results
        """)
        
        # Example use case
        with st.expander("💡 Example Use Case"):
            st.markdown("""
            **Scenario:** Finding cardiovascular drugs for Alzheimer's disease
            
            **Steps:**
            1. Select **"Cardiovascular"** as source category
            2. Select **"Alzheimers"** as target disease
            3. Click **"Analyze Drugs"**
            
            **Results:**
            - System analyzes ~88 cardiovascular drugs
            - Scores each based on Alzheimer's pathway connections
            - Shows top 3 most promising candidates
            - Displays evidence network with drug-target-pathway connections
            
            **Why it works:**
            Many cardiovascular drugs affect pathways relevant to neurodegeneration:
            - Calcium channel blockers → neuroprotection
            - Statins → inflammation reduction
            - ACE inhibitors → oxidative stress reduction
            """)
        
        # Quick start guide
        with st.expander("🚀 Quick Start Guide"):
            st.markdown("""
            **For best results:**
            
            1. **Choose related categories**
               - Diabetic → Alzheimers (metabolism connection)
               - Cardiovascular → Parkinsons (circulation/inflammation)
               - Psychiatric → Pain (neurotransmitter pathways)
            
            2. **Adjust settings**
               - Use 3-5 candidates for focused analysis
               - Enable diversity for varied mechanisms
               - Check 'Show all scores' to see full ranking
            
            3. **Interpret results**
               - Higher scores = stronger biological connection
               - Review matched pathways for mechanism insights
               - Examine evidence network for validation
            
            4. **Export data**
               - Download complete rankings as CSV
               - Share top candidates for further review
            """)

# ============================================================================
# RUN APP
# ============================================================================

if __name__ == "__main__":
    main()
