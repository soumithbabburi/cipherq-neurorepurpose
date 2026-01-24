#!/usr/bin/env python3
"""
Enhanced NLP Chatbox System with Intelligent Recommendations
Very wide interface with robust multi-AI integration
"""

import streamlit as st
import logging
from typing import Dict, List, Optional, Tuple, Any
import json
import requests
from datetime import datetime
import pandas as pd
import time

# Import existing integrations
try:
    from anthropic import Anthropic
    ANTHROPIC_AVAILABLE = True
except ImportError:
    ANTHROPIC_AVAILABLE = False

try:
    import openai
    OPENAI_AVAILABLE = True
except ImportError:
    OPENAI_AVAILABLE = False

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class EnhancedNLPChatbox:
    """
    Very wide, intelligent chatbox with multi-AI integration and pathway recommendations
    """
    
    def __init__(self):
        self.conversation_history = []
        self.pathway_database = self._load_pathway_database()
        self.drug_knowledge_base = self._load_drug_knowledge_base()
        self.protein_targets = self._load_protein_targets()
        
        # Initialize AI clients
        self.anthropic_client = None
        self.openai_client = None
        self._initialize_ai_clients()
        
        logger.info("Enhanced NLP Chatbox initialized with multi-AI integration")
    
    def _initialize_ai_clients(self):
        """Initialize available AI clients"""
        
        # Initialize Anthropic client
        if ANTHROPIC_AVAILABLE:
            try:
                anthropic_key = st.secrets.get("ANTHROPIC_API_KEY") or None
                if anthropic_key:
                    self.anthropic_client = Anthropic(api_key=anthropic_key)
                    logger.info("✅ Anthropic Claude client initialized")
            except Exception as e:
                logger.warning(f"⚠️ Anthropic initialization failed: {e}")
        
        # Initialize OpenAI client
        if OPENAI_AVAILABLE:
            try:
                openai_key = st.secrets.get("OPENAI_API_KEY") or None
                if openai_key:
                    openai.api_key = openai_key
                    self.openai_client = openai
                    logger.info("✅ OpenAI client initialized")
            except Exception as e:
                logger.warning(f"⚠️ OpenAI initialization failed: {e}")
    
    def _load_pathway_database(self) -> Dict:
        """Load comprehensive pathway database for intelligent recommendations"""
        return {
            'alzheimer_pathways': [
                'Amyloid-beta processing',
                'Tau protein phosphorylation',
                'Neuroinflammation cascades',
                'Cholinergic signaling',
                'Oxidative stress response',
                'Mitochondrial dysfunction',
                'Calcium homeostasis',
                'Synaptic plasticity',
                'Autophagy mechanisms',
                'Blood-brain barrier integrity'
            ],
            'metabolic_pathways': [
                'Glycolysis/gluconeogenesis',
                'Citric acid cycle',
                'Oxidative phosphorylation',
                'Fatty acid metabolism',
                'Amino acid metabolism',
                'Insulin signaling',
                'mTOR pathway',
                'AMPK activation',
                'Ketogenesis',
                'Glycogen metabolism'
            ],
            'signaling_pathways': [
                'PI3K/AKT pathway',
                'MAPK/ERK cascade',
                'JAK/STAT signaling',
                'TGF-β pathway',
                'Wnt signaling',
                'Notch pathway',
                'Hedgehog signaling',
                'NF-κB pathway',
                'p53 tumor suppressor',
                'Cell cycle regulation'
            ],
            'neurotransmitter_pathways': [
                'Dopaminergic signaling',
                'Serotonergic transmission',
                'GABAergic inhibition',
                'Glutamatergic excitation',
                'Cholinergic transmission',
                'Noradrenergic signaling',
                'Endocannabinoid system',
                'Opioid signaling',
                'Histaminergic pathways',
                'Purinergic signaling'
            ]
        }
    
    def _load_drug_knowledge_base(self) -> Dict:
        """Load drug knowledge base for intelligent matching"""
        return {
            'antidiabetic_drugs': {
                'metformin': {
                    'mechanism': 'AMPK activation, mitochondrial complex I inhibition',
                    'pathways': ['AMPK pathway', 'Gluconeogenesis inhibition', 'Mitochondrial respiration'],
                    'targets': ['AMPK', 'Complex I', 'mGPD'],
                    'indications': ['Type 2 diabetes', 'PCOS', 'Weight management'],
                    'repurposing_potential': ['Alzheimer disease', 'Cancer', 'Aging'],
                    'molecular_weight': 129.16,
                    'bioavailability': '50-60%'
                },
                'pioglitazone': {
                    'mechanism': 'PPARγ agonism, insulin sensitization',
                    'pathways': ['PPAR signaling', 'Glucose metabolism', 'Lipid metabolism'],
                    'targets': ['PPARγ', 'PPARα', 'GLUT4'],
                    'indications': ['Type 2 diabetes', 'Insulin resistance'],
                    'repurposing_potential': ['Alzheimer disease', 'NASH', 'Cardiovascular disease'],
                    'molecular_weight': 356.44,
                    'bioavailability': '>80%'
                },
                'sitagliptin': {
                    'mechanism': 'DPP-4 inhibition, incretin enhancement',
                    'pathways': ['Incretin pathway', 'GLP-1 signaling', 'Glucose homeostasis'],
                    'targets': ['DPP-4', 'GLP-1R'],
                    'indications': ['Type 2 diabetes'],
                    'repurposing_potential': ['Alzheimer disease', 'Cardiovascular protection'],
                    'molecular_weight': 407.31,
                    'bioavailability': '87%'
                }
            }
        }
    
    def _load_protein_targets(self) -> Dict:
        """Load protein target database for recommendations"""
        return {
            'alzheimer_targets': {
                'AMPK': {
                    'full_name': 'AMP-activated protein kinase',
                    'function': 'Energy sensor and metabolic regulator',
                    'pathways': ['Energy homeostasis', 'Autophagy', 'Mitochondrial biogenesis'],
                    'alzheimer_relevance': 'Promotes tau clearance, reduces amyloid accumulation',
                    'druggability': 'High',
                    'pdb_ids': ['4CFE', '5EXR', '4ZHX']
                },
                'PPARγ': {
                    'full_name': 'Peroxisome proliferator-activated receptor gamma',
                    'function': 'Nuclear receptor transcription factor',
                    'pathways': ['Lipid metabolism', 'Glucose homeostasis', 'Inflammation'],
                    'alzheimer_relevance': 'Anti-inflammatory effects, amyloid clearance',
                    'druggability': 'High',
                    'pdb_ids': ['2PRG', '3DZY', '1FM9']
                },
                'DPP4': {
                    'full_name': 'Dipeptidyl peptidase-4',
                    'function': 'Serine protease, incretin degradation',
                    'pathways': ['Incretin pathway', 'Glucose regulation'],
                    'alzheimer_relevance': 'Neuroprotective incretin signaling',
                    'druggability': 'High',
                    'pdb_ids': ['1ORV', '1NU6', '1NU8']
                }
            }
        }
    
    def create_wide_chatbox_interface(self):
        """Create very wide chatbox interface with intelligent features"""
        
        # Main chatbox container with full width
        st.markdown("""
        <style>
        .wide-chatbox {
            width: 100%;
            max-width: none !important;
        }
        .chat-container {
            width: 100%;
            padding: 0;
            margin: 0;
        }
        .stTextArea > div > div > textarea {
            width: 100% !important;
            min-height: 120px !important;
        }
        .chat-suggestions {
            background: linear-gradient(90deg, #f0f8ff, #e6f3ff);
            padding: 15px;
            border-radius: 10px;
            margin: 10px 0;
            border-left: 4px solid #4CAF50;
        }
        </style>
        """, unsafe_allow_html=True)
        
        st.markdown("## 🤖 Enhanced AI Assistant")
        st.markdown("*Ask about drugs, pathways, proteins, mechanisms, or get intelligent recommendations*")
        
        # Create full-width columns
        col1, col2 = st.columns([4, 1])
        
        with col1:
            # Very wide text area
            user_input = st.text_area(
                "💬 Enter your question or describe what you're looking for:",
                height=120,
                key="wide_chatbox_input",
                placeholder="Examples:\n• 'What pathways does metformin target for Alzheimer's?'\n• 'Recommend drugs similar to pioglitazone mechanism'\n• 'Show me AMPK protein structure and binding sites'\n• 'Find connections between diabetes drugs and neurodegeneration'"
            )
        
        with col2:
            st.markdown("### 🚀 Quick Actions")
            if st.button("🧠 Pathway Analysis", use_container_width=True):
                user_input = "Analyze drug-pathway interactions for Alzheimer's disease"
            if st.button("🔬 Mechanism Deep Dive", use_container_width=True):
                user_input = "Explain molecular mechanisms and targets"
            if st.button("💊 Drug Recommendations", use_container_width=True):
                user_input = "Recommend drugs for repurposing based on mechanisms"
            if st.button("📊 Literature Search", use_container_width=True):
                user_input = "Search latest research and clinical trials"
        
        # Process input
        if user_input and user_input.strip():
            self._process_intelligent_query(user_input)
        
        # Display conversation history
        self._display_conversation_history()
        
        # Show intelligent suggestions
        self._show_intelligent_suggestions()
    
    def _process_intelligent_query(self, query: str):
        """Process user query with intelligent analysis and recommendations"""
        
        with st.spinner("🔍 Analyzing your query with AI intelligence..."):
            
            # Analyze query intent
            query_analysis = self._analyze_query_intent(query)
            
            # Generate AI response
            ai_response = self._generate_ai_response(query, query_analysis)
            
            # Add to conversation history
            self.conversation_history.append({
                'timestamp': datetime.now(),
                'user_query': query,
                'analysis': query_analysis,
                'ai_response': ai_response,
                'recommendations': self._generate_recommendations(query_analysis)
            })
            
            # Clear input
            if 'wide_chatbox_input' in st.session_state:
                st.session_state.wide_chatbox_input = ""
    
    def _analyze_query_intent(self, query: str) -> Dict:
        """Analyze user query to understand intent and extract entities with smart parsing"""
        
        query_lower = query.lower()
        
        analysis = {
            'intent_type': 'general',
            'entities': {
                'drugs': [],
                'proteins': [],
                'pathways': [],
                'diseases': []
            },
            'query_categories': [],
            'confidence': 0.0,
            'target_query': None,
            'pathway_query': None,
            'result_limit': None
        }
        
        # Smart Detection: Target-based queries
        target_patterns = [
            r'drugs?\s+(?:which\s+have\s+|with\s+|that\s+target\s+|targeting\s+)(?:targets?\s+)?(\w+)',
            r'find\s+drugs?\s+(?:for\s+|targeting\s+)(\w+)',
            r'(\w+)\s+target(?:ed)?\s+drugs?',
            r'drugs?\s+against\s+(\w+)',
            r'(\w+)\s+inhibitors?',
            r'(\w+)\s+agonists?'
        ]
        
        for pattern in target_patterns:
            import re
            match = re.search(pattern, query_lower)
            if match:
                target = match.group(1).upper()
                analysis['target_query'] = target
                analysis['intent_type'] = 'target_specific_drugs'
                analysis['query_categories'].append('target_drug_search')
                analysis['entities']['proteins'].append(target)
                # Look for result limit
                if any(word in query_lower for word in ['top 3', 'best 3', 'first 3']):
                    analysis['result_limit'] = 3
                elif 'top' in query_lower or 'best' in query_lower:
                    analysis['result_limit'] = 3  # Default to 3 for targeted queries
                break
        
        # Smart Detection: Pathway-based queries  
        pathway_patterns = [
            r'drugs?\s+(?:in\s+|for\s+|targeting\s+|modulating\s+)(.+?)\s+pathway',
            r'(.+?)\s+pathway\s+drugs?',
            r'drugs?\s+(?:that\s+affect\s+|modulating\s+)(.+?)(?:\s+signaling)?',
            r'find\s+drugs?\s+(?:for\s+)?(.+?)\s+(?:pathway|signaling)'
        ]
        
        for pattern in pathway_patterns:
            import re
            match = re.search(pattern, query_lower)
            if match:
                pathway = match.group(1).strip()
                analysis['pathway_query'] = pathway
                analysis['intent_type'] = 'pathway_specific_drugs'
                analysis['query_categories'].append('pathway_drug_search')
                analysis['entities']['pathways'].append(pathway)
                if any(word in query_lower for word in ['top 3', 'best 3', 'first 3']):
                    analysis['result_limit'] = 3
                elif 'top' in query_lower or 'best' in query_lower:
                    analysis['result_limit'] = 3
                break
        
        # Enhanced entity detection
        # Detect drugs from knowledge base
        for drug_class, drugs in self.drug_knowledge_base.items():
            for drug_name in drugs.keys():
                if drug_name in query_lower:
                    analysis['entities']['drugs'].append(drug_name)
        
        # Detect proteins/targets (expanded list)
        all_targets = ['AMPK', 'PPARγ', 'PPAR', 'DPP4', 'DPP-4', 'BACE1', 'ACE', 'AChE', 'Tau', 
                      'mTOR', 'PI3K', 'AKT', 'MAPK', 'ERK', 'JAK', 'STAT', 'TGF', 'Wnt', 
                      'Notch', 'NF-κB', 'p53', 'COX', 'TNF', 'IL-1', 'IL-6']
        
        for target in all_targets:
            if target.lower() in query_lower:
                analysis['entities']['proteins'].append(target)
        
        # Detect pathways from database
        for pathway_class, pathways in self.pathway_database.items():
            for pathway in pathways:
                if pathway.lower() in query_lower:
                    analysis['entities']['pathways'].append(pathway)
        
        # Enhanced intent detection
        if any(word in query_lower for word in ['pathway', 'mechanism', 'signaling', 'cascade']):
            if analysis['intent_type'] == 'general':
                analysis['intent_type'] = 'pathway_analysis'
            analysis['query_categories'].append('pathway_analysis')
        
        if any(word in query_lower for word in ['recommend', 'suggest', 'similar', 'alternative']):
            if analysis['intent_type'] == 'general':
                analysis['intent_type'] = 'recommendation'
            analysis['query_categories'].append('drug_recommendation')
        
        if any(word in query_lower for word in ['structure', 'binding', 'protein', 'pdb', '3d']):
            if analysis['intent_type'] == 'general':
                analysis['intent_type'] = 'structural_analysis'
            analysis['query_categories'].append('structural_analysis')
        
        if any(word in query_lower for word in ['alzheimer', 'dementia', 'neurodegeneration']):
            analysis['entities']['diseases'].append('Alzheimer disease')
            analysis['query_categories'].append('alzheimer_focus')
        
        # Calculate enhanced confidence based on query specificity
        total_entities = sum(len(entities) for entities in analysis['entities'].values())
        specificity_bonus = 0.3 if analysis['target_query'] or analysis['pathway_query'] else 0
        analysis['confidence'] = min(0.95, total_entities * 0.2 + 0.1 + specificity_bonus)
        
        return analysis
    
    def _generate_ai_response(self, query: str, analysis: Dict) -> str:
        """Generate intelligent AI response using available models"""
        
        # Prepare context
        context = self._build_context(analysis)
        
        # Try Anthropic Claude first
        if self.anthropic_client:
            try:
                response = self._query_anthropic(query, context, analysis)
                if response:
                    return response
            except Exception as e:
                logger.warning(f"⚠️ Anthropic query failed: {e}")
        
        # Fallback to OpenAI
        if self.openai_client:
            try:
                response = self._query_openai(query, context, analysis)
                if response:
                    return response
            except Exception as e:
                logger.warning(f"⚠️ OpenAI query failed: {e}")
        
        # Fallback to smart database response for specific queries
        return self._generate_smart_database_response(query, analysis)
    
    def _generate_smart_database_response(self, query: str, analysis: Dict) -> str:
        """Generate intelligent database-driven response for target and pathway queries"""
        
        # Handle target-specific drug queries
        if analysis['intent_type'] == 'target_specific_drugs' and analysis['target_query']:
            return self._get_drugs_for_target(analysis['target_query'], analysis.get('result_limit', 3))
        
        # Handle pathway-specific drug queries  
        if analysis['intent_type'] == 'pathway_specific_drugs' and analysis['pathway_query']:
            return self._get_drugs_for_pathway(analysis['pathway_query'], analysis.get('result_limit', 3))
        
        # Fallback to rule-based response
        return self._generate_rule_based_response(query, analysis)
    
    def _get_drugs_for_target(self, target: str, limit: int = 3) -> str:
        """Get drugs that target a specific protein (e.g., AMPK)"""
        
        target_upper = target.upper()
        
        # Check knowledge base first for immediate results
        target_drugs = []
        
        if target_upper == 'AMPK':
            target_drugs = [
                {
                    'drug': 'Metformin',
                    'mechanism': 'AMPK activation, enhances cellular energy sensing',
                    'clinical_status': 'FDA approved for diabetes, investigated for Alzheimer\'s',
                    'binding_affinity': 'Allosteric activator',
                    'repurposing_potential': 'High for neurodegeneration',
                    'confidence': 0.95
                },
                {
                    'drug': 'Phenformin',
                    'mechanism': 'Potent AMPK activator, stronger than metformin',
                    'clinical_status': 'Withdrawn due to lactic acidosis risk',
                    'binding_affinity': 'Direct AMPK activation',
                    'repurposing_potential': 'Research use only',
                    'confidence': 0.90
                },
                {
                    'drug': 'AICAR',
                    'mechanism': 'AMP mimetic, direct AMPK activation',
                    'clinical_status': 'Research compound, clinical trials',
                    'binding_affinity': 'AMP-binding domain agonist',
                    'repurposing_potential': 'High for metabolic disorders',
                    'confidence': 0.88
                }
            ]
        elif target_upper in ['PPAR', 'PPARG', 'PPARγ']:
            target_drugs = [
                {
                    'drug': 'Pioglitazone',
                    'mechanism': 'PPARγ agonist, insulin sensitization',
                    'clinical_status': 'FDA approved for diabetes',
                    'binding_affinity': 'Full agonist',
                    'repurposing_potential': 'High for Alzheimer\'s disease',
                    'confidence': 0.92
                },
                {
                    'drug': 'Rosiglitazone',
                    'mechanism': 'PPARγ agonist, glucose homeostasis',
                    'clinical_status': 'FDA approved (restricted)',
                    'binding_affinity': 'Full agonist',
                    'repurposing_potential': 'Moderate for neurodegeneration',
                    'confidence': 0.85
                },
                {
                    'drug': 'Troglitazone',
                    'mechanism': 'PPARγ agonist, first thiazolidinedione',
                    'clinical_status': 'Withdrawn due to hepatotoxicity',
                    'binding_affinity': 'Full agonist',
                    'repurposing_potential': 'Low due to safety concerns',
                    'confidence': 0.75
                }
            ]
        elif target_upper in ['DPP4', 'DPP-4']:
            target_drugs = [
                {
                    'drug': 'Sitagliptin',
                    'mechanism': 'DPP-4 inhibition, GLP-1 enhancement',
                    'clinical_status': 'FDA approved for diabetes',
                    'binding_affinity': 'Competitive inhibitor',
                    'repurposing_potential': 'High for neuroprotection',
                    'confidence': 0.90
                },
                {
                    'drug': 'Saxagliptin',
                    'mechanism': 'DPP-4 inhibition, incretin preservation',
                    'clinical_status': 'FDA approved for diabetes',
                    'binding_affinity': 'Reversible inhibitor',
                    'repurposing_potential': 'Moderate for Alzheimer\'s',
                    'confidence': 0.85
                },
                {
                    'drug': 'Linagliptin',
                    'mechanism': 'DPP-4 inhibition, tissue-selective',
                    'clinical_status': 'FDA approved for diabetes',
                    'binding_affinity': 'High selectivity',
                    'repurposing_potential': 'High for brain penetration',
                    'confidence': 0.88
                }
            ]
        else:
            # For other targets, provide general guidance
            return f"""**🎯 Target Query: {target}**

I found your query about drugs targeting **{target}**. While I have comprehensive data for key targets like AMPK, PPARγ, and DPP-4, this specific target may require deeper database analysis.

**🔍 Suggested Approach:**
1. **Database Search**: Query ChEMBL or DrugBank for {target} inhibitors/modulators
2. **Literature Review**: Check PubMed for recent {target}-targeting compounds
3. **Clinical Pipeline**: Search ClinicalTrials.gov for {target}-directed therapies

**💡 Related Targets to Consider:**
- If {target} is metabolic: AMPK, mTOR, PPARγ
- If {target} is neurological: AChE, BACE1, NMDA receptors
- If {target} is inflammatory: NF-κB, COX-2, TNF-α

Would you like me to search for drugs targeting a related pathway or provide more specific guidance?"""
        
        # Format response for target drugs
        if target_drugs:
            response_parts = [f"**🎯 Top {min(limit, len(target_drugs))} Drugs Targeting {target_upper}**\n"]
            
            for i, drug_info in enumerate(target_drugs[:limit], 1):
                response_parts.append(f"**{i}. {drug_info['drug']}**")
                response_parts.append(f"   • **Mechanism**: {drug_info['mechanism']}")
                response_parts.append(f"   • **Clinical Status**: {drug_info['clinical_status']}")
                response_parts.append(f"   • **Binding**: {drug_info['binding_affinity']}")
                response_parts.append(f"   • **Repurposing Potential**: {drug_info['repurposing_potential']}")
                response_parts.append(f"   • **Confidence**: {drug_info['confidence']:.0%}\n")
            
            # Add strategic insights
            if target_upper == 'AMPK':
                response_parts.append("**🧠 Strategic Insights for AMPK:**")
                response_parts.append("• AMPK activation promotes autophagy and tau clearance")
                response_parts.append("• Metabolic benefits may translate to neuroprotection")
                response_parts.append("• Metformin shows promise in Alzheimer's clinical trials")
                response_parts.append("• Consider combination with other neuroprotective agents")
            
            return "\n".join(response_parts)
        
        return f"No specific drugs found for target {target}. Please try a different target or pathway."
    
    def _get_drugs_for_pathway(self, pathway: str, limit: int = 3) -> str:
        """Get drugs that modulate a specific pathway"""
        
        pathway_lower = pathway.lower()
        pathway_drugs = []
        
        # AMPK pathway drugs
        if 'ampk' in pathway_lower or 'energy' in pathway_lower:
            pathway_drugs = [
                {
                    'drug': 'Metformin',
                    'pathway_role': 'Primary AMPK activator',
                    'mechanism': 'Enhances AMP/ATP ratio, activates AMPK cascade',
                    'clinical_evidence': 'Strong evidence in diabetes, emerging in Alzheimer\'s'
                },
                {
                    'drug': 'Berberine',
                    'pathway_role': 'Natural AMPK activator',
                    'mechanism': 'Mitochondrial complex I inhibition → AMPK activation',
                    'clinical_evidence': 'Preclinical neuroprotection, Phase II trials'
                },
                {
                    'drug': 'Resveratrol',
                    'pathway_role': 'AMPK/SIRT1 modulator',
                    'mechanism': 'Polyphenol activation of AMPK-SIRT1 axis',
                    'clinical_evidence': 'Mixed results, requires high doses'
                }
            ]
        
        # Wnt signaling pathway
        elif 'wnt' in pathway_lower:
            pathway_drugs = [
                {
                    'drug': 'Lithium',
                    'pathway_role': 'GSK-3β inhibitor, Wnt activator',
                    'mechanism': 'Stabilizes β-catenin, enhances Wnt signaling',
                    'clinical_evidence': 'FDA approved for bipolar, neuroprotective studies'
                },
                {
                    'drug': 'Valproic Acid',
                    'pathway_role': 'Wnt pathway modulator',
                    'mechanism': 'HDAC inhibition, enhances Wnt target gene expression',
                    'clinical_evidence': 'Epilepsy approved, Alzheimer\'s trials ongoing'
                },
                {
                    'drug': 'LRP6 agonists',
                    'pathway_role': 'Direct Wnt receptor activation',
                    'mechanism': 'LRP5/6 co-receptor stimulation',
                    'clinical_evidence': 'Experimental compounds, early development'
                }
            ]
        
        # mTOR pathway
        elif 'mtor' in pathway_lower:
            pathway_drugs = [
                {
                    'drug': 'Rapamycin',
                    'pathway_role': 'mTOR complex 1 inhibitor',
                    'mechanism': 'FKBP12 binding, mTORC1 suppression',
                    'clinical_evidence': 'Immunosuppression approved, aging studies'
                },
                {
                    'drug': 'Everolimus',
                    'pathway_role': 'mTOR inhibitor',
                    'mechanism': 'Selective mTORC1 inhibition',
                    'clinical_evidence': 'Cancer approved, neurodegeneration research'
                },
                {
                    'drug': 'Metformin',
                    'pathway_role': 'Indirect mTOR suppression',
                    'mechanism': 'AMPK activation → mTOR inhibition',
                    'clinical_evidence': 'Diabetes approved, longevity studies'
                }
            ]
        
        else:
            return f"""**🔬 Pathway Query: {pathway}**

I'm searching for drugs that modulate the **{pathway}** pathway. While I have detailed data for key pathways like AMPK signaling, Wnt signaling, and mTOR, this specific pathway may need targeted analysis.

**🔍 Pathway Analysis Approach:**
1. **Identify Key Proteins**: Map pathway components and druggable targets
2. **Literature Mining**: Search for pathway modulators in recent publications  
3. **Drug Databases**: Query ChEMBL/DrugBank for pathway-targeting compounds
4. **Clinical Evidence**: Check for pathway-directed therapies in trials

**💡 Related Pathways:**
- **Metabolic**: AMPK, mTOR, insulin signaling
- **Neuronal**: Wnt, Notch, MAPK cascades  
- **Inflammatory**: NF-κB, JAK/STAT, TGF-β

Would you like me to analyze a specific related pathway or provide more targeted guidance?"""
        
        if pathway_drugs:
            response_parts = [f"**🔬 Top {min(limit, len(pathway_drugs))} Drugs Modulating {pathway.title()} Pathway**\n"]
            
            for i, drug_info in enumerate(pathway_drugs[:limit], 1):
                response_parts.append(f"**{i}. {drug_info['drug']}**")
                response_parts.append(f"   • **Pathway Role**: {drug_info['pathway_role']}")
                response_parts.append(f"   • **Mechanism**: {drug_info['mechanism']}")
                response_parts.append(f"   • **Clinical Evidence**: {drug_info['clinical_evidence']}\n")
            
            return "\n".join(response_parts)
        
        return f"No specific drugs found for {pathway} pathway. Please try a different pathway."
    
    def _query_anthropic(self, query: str, context: str, analysis: Dict) -> str:
        """Query Anthropic Claude for intelligent response"""
        
        system_prompt = f"""You are a world-class expert in drug repurposing, molecular biology, pharmacology, and patent analysis with comprehensive knowledge across:

        • Molecular mechanisms and drug-target interactions
        • Clinical pharmacology and ADME properties  
        • Patent landscape and IP strategy
        • Regulatory pathways and market access
        • Biomarker discovery and precision medicine
        • Systems biology and network pharmacology
        • Real-world evidence and clinical outcomes
        • Computational drug discovery methods

        Context: {context}
        
        User Query Analysis: {json.dumps(analysis, indent=2)}
        
        Provide comprehensive, scientifically rigorous responses that:
        1. **Direct Answer**: Address the user's specific question with authoritative expertise
        2. **Scientific Depth**: Include detailed molecular mechanisms, pathway interactions, and pharmacokinetic considerations
        3. **Strategic Insights**: Offer actionable intelligence for drug development, patent strategy, and market positioning
        4. **Broader Context**: Connect to related therapeutic areas, alternative approaches, and emerging research
        5. **Regulatory Perspective**: Consider FDA approval pathways, clinical trial design, and market access barriers
        6. **Innovation Opportunities**: Identify novel repurposing angles, combination therapies, and unmet medical needs
        7. **Risk Assessment**: Highlight safety concerns, contraindications, and development challenges
        8. **Evidence Integration**: Synthesize data from preclinical, clinical, and real-world sources

        Your responses should demonstrate the breadth and depth expected from a leading pharmaceutical R&D executive with expertise across discovery, development, regulatory, and commercial domains."""
        
        try:
            message = self.anthropic_client.messages.create(
                model="claude-sonnet-4-20250514",  # Latest model for comprehensive analysis
                max_tokens=2000,  # Increased for more detailed responses
                temperature=0.2,  # Lower temperature for more focused, factual responses
                system=system_prompt,
                messages=[{"role": "user", "content": query}]
            )
            
            return message.content[0].text if message.content else None
            
        except Exception as e:
            logger.error(f"❌ Anthropic query error: {e}")
            return None
    
    def _query_openai(self, query: str, context: str, analysis: Dict) -> str:
        """Query OpenAI for intelligent response"""
        
        try:
            response = self.openai_client.ChatCompletion.create(
                model="gpt-4",
                messages=[
                    {
                        "role": "system",
                        "content": f"""You are an expert drug repurposing AI assistant.
                        
                        Context: {context}
                        Analysis: {json.dumps(analysis)}
                        
                        Provide comprehensive, scientifically accurate responses about drug mechanisms, pathways, and repurposing opportunities."""
                    },
                    {"role": "user", "content": query}
                ],
                max_tokens=1500,
                temperature=0.3
            )
            
            return response.choices[0].message.content
            
        except Exception as e:
            logger.error(f"❌ OpenAI query error: {e}")
            return None
    
    def _generate_rule_based_response(self, query: str, analysis: Dict) -> str:
        """Generate rule-based response when AI APIs unavailable"""
        
        response_parts = []
        
        # Drug information
        if analysis['entities']['drugs']:
            for drug in analysis['entities']['drugs']:
                drug_info = self._get_drug_info(drug)
                if drug_info:
                    response_parts.append(f"**{drug.title()}**: {drug_info}")
        
        # Pathway information
        if analysis['entities']['pathways']:
            response_parts.append("**Relevant Pathways:**")
            for pathway in analysis['entities']['pathways'][:3]:
                response_parts.append(f"• {pathway}")
        
        # Protein information
        if analysis['entities']['proteins']:
            response_parts.append("**Target Proteins:**")
            for protein in analysis['entities']['proteins']:
                protein_info = self._get_protein_info(protein)
                if protein_info:
                    response_parts.append(f"• {protein}: {protein_info}")
        
        if not response_parts:
            response_parts.append("I understand you're asking about drug repurposing and molecular mechanisms. For the most comprehensive analysis, please ensure your question includes specific drug names, protein targets, or pathways of interest.")
        
        return "\n\n".join(response_parts)
    
    def _build_context(self, analysis: Dict) -> str:
        """Build context for AI query"""
        context_parts = []
        
        # Add drug context
        if analysis['entities']['drugs']:
            context_parts.append("Relevant Drugs: " + ", ".join(analysis['entities']['drugs']))
        
        # Add pathway context
        if analysis['entities']['pathways']:
            context_parts.append("Relevant Pathways: " + ", ".join(analysis['entities']['pathways']))
        
        # Add protein context
        if analysis['entities']['proteins']:
            context_parts.append("Relevant Proteins: " + ", ".join(analysis['entities']['proteins']))
        
        return " | ".join(context_parts)
    
    def _get_drug_info(self, drug_name: str) -> str:
        """Get drug information from knowledge base"""
        for drug_class, drugs in self.drug_knowledge_base.items():
            if drug_name in drugs:
                drug_data = drugs[drug_name]
                return f"Mechanism: {drug_data['mechanism']}. Targets: {', '.join(drug_data['targets'])}."
        return None
    
    def _get_protein_info(self, protein_name: str) -> str:
        """Get protein information from knowledge base"""
        for target_class, proteins in self.protein_targets.items():
            if protein_name in proteins:
                protein_data = proteins[protein_name]
                return f"{protein_data['full_name']} - {protein_data['function']}"
        return None
    
    def _generate_recommendations(self, analysis: Dict) -> List[Dict]:
        """Generate intelligent recommendations based on query analysis"""
        recommendations = []
        
        # Drug recommendations
        if analysis['intent_type'] == 'recommendation':
            recommendations.extend(self._get_drug_recommendations(analysis))
        
        # Pathway recommendations
        if 'pathway_analysis' in analysis['query_categories']:
            recommendations.extend(self._get_pathway_recommendations(analysis))
        
        # Structural analysis recommendations
        if 'structural_analysis' in analysis['query_categories']:
            recommendations.extend(self._get_structural_recommendations(analysis))
        
        return recommendations
    
    def _get_drug_recommendations(self, analysis: Dict) -> List[Dict]:
        """Get drug recommendations based on analysis"""
        recommendations = []
        
        if analysis['entities']['drugs']:
            for drug in analysis['entities']['drugs']:
                drug_info = self._get_drug_data(drug)
                if drug_info:
                    recommendations.append({
                        'type': 'drug_similar',
                        'title': f"Drugs similar to {drug.title()}",
                        'content': f"Consider drugs with similar mechanisms: {drug_info.get('mechanism', 'Unknown')}"
                    })
        
        return recommendations
    
    def _get_pathway_recommendations(self, analysis: Dict) -> List[Dict]:
        """Get pathway recommendations"""
        recommendations = []
        
        if 'alzheimer_focus' in analysis['query_categories']:
            recommendations.append({
                'type': 'pathway_analysis',
                'title': "Key Alzheimer's Pathways",
                'content': "Focus on: Amyloid processing, Tau phosphorylation, Neuroinflammation, Mitochondrial dysfunction"
            })
        
        return recommendations
    
    def _get_structural_recommendations(self, analysis: Dict) -> List[Dict]:
        """Get structural analysis recommendations"""
        recommendations = []
        
        if analysis['entities']['proteins']:
            recommendations.append({
                'type': 'structural_analysis',
                'title': "3D Structure Analysis",
                'content': f"Explore 3D structures for: {', '.join(analysis['entities']['proteins'])}"
            })
        
        return recommendations
    
    def _get_drug_data(self, drug_name: str) -> Dict:
        """Get comprehensive drug data"""
        for drug_class, drugs in self.drug_knowledge_base.items():
            if drug_name in drugs:
                return drugs[drug_name]
        return {}
    
    def _display_conversation_history(self):
        """Display conversation history with intelligent formatting"""
        
        if not self.conversation_history:
            return
        
        st.markdown("### 💬 Conversation History")
        
        # Show last few conversations
        for i, conversation in enumerate(self.conversation_history[-3:]):
            with st.expander(f"Query {len(self.conversation_history) - 2 + i}: {conversation['user_query'][:60]}..."):
                
                st.markdown("**Your Question:**")
                st.write(conversation['user_query'])
                
                st.markdown("**AI Response:**")
                st.write(conversation['ai_response'])
                
                if conversation['recommendations']:
                    st.markdown("**Recommendations:**")
                    for rec in conversation['recommendations']:
                        st.info(f"**{rec['title']}**: {rec['content']}")
    
    def _show_intelligent_suggestions(self):
        """Show intelligent suggestions based on current context"""
        
        st.markdown("### 💡 Intelligent Suggestions")
        
        suggestions = [
            "🧠 **Pathway Analysis**: 'Analyze AMPK pathway connections to Alzheimer's disease'",
            "🔬 **Mechanism Study**: 'Compare molecular mechanisms of metformin vs pioglitazone'",
            "📊 **Drug Screening**: 'Recommend diabetes drugs for Alzheimer's repurposing'",
            "🎯 **Target Analysis**: 'Show protein-drug interactions for PPARγ'",
            "📚 **Literature Search**: 'Find latest clinical trials for drug repurposing'"
        ]
        
        cols = st.columns(2)
        for i, suggestion in enumerate(suggestions):
            with cols[i % 2]:
                if st.button(suggestion, key=f"suggestion_{i}", use_container_width=True):
                    # Extract the actual query from the suggestion
                    query = suggestion.split("': '")[1].rstrip("'")
                    st.session_state.wide_chatbox_input = query
                    st.rerun()

def create_enhanced_chatbox():
    """Create the enhanced chatbox interface"""
    chatbox = EnhancedNLPChatbox()
    chatbox.create_wide_chatbox_interface()

if __name__ == "__main__":
    # Test the enhanced chatbox
    chatbox = EnhancedNLPChatbox()
    print("Enhanced NLP Chatbox initialized successfully!")
    
    # Test query analysis
    test_query = "What pathways does metformin target for Alzheimer's treatment?"
    analysis = chatbox._analyze_query_intent(test_query)
    print(f"Query analysis: {analysis}")