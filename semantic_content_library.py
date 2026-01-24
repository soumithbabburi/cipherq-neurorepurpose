"""
Semantic Content Library
Pre-built knowledge base for intelligent drug repurposing conversations
"""
import psycopg2
import os
import logging
from typing import Dict, List, Optional
import json

logger = logging.getLogger(__name__)

class SemanticContentLibrary:
    """
    Comprehensive knowledge base for drug repurposing.
    Uses 200K pathways + 100K genes + all database relationships.
    """
    
    def __init__(self):
        self.conn = psycopg2.connect(
            host=os.getenv("DB_HOST", "localhost"),
            database=os.getenv("DB_NAME", "cipherq_repurpose"),
            user=os.getenv("DB_USER", "babburisoumith"),
            password=os.getenv("DB_PASSWORD", "")
        )
        self.cache = {}
        logger.info("✅ Semantic Content Library initialized")
    
    def get_drug_profile(self, drug_name: str) -> Dict:
        """Get complete drug profile with targets, pathways, diseases"""
        
        cursor = self.conn.cursor()
        
        profile = {
            'drug_name': drug_name,
            'targets': [],
            'pathways': [],
            'known_diseases': [],
            'repurposing_candidates': [],
            'mechanisms': []
        }
        
        # Get drug targets
        cursor.execute("""
            SELECT DISTINCT g.gene_symbol, g.gene_name, dpi.interaction_type
            FROM drug_protein_interactions dpi
            JOIN genes g ON dpi.gene_symbol = g.gene_symbol
            WHERE dpi.drug_name = %s
            LIMIT 20
        """, (drug_name,))
        
        targets = cursor.fetchall()
        profile['targets'] = [{'gene': gene, 'name': name, 'interaction': itype} for gene, name, itype in targets]
        
        # Get pathways (from 200K!)
        if targets:
            gene_symbols = [t[0] for t in targets]
            cursor.execute("""
                SELECT DISTINCT p.pathway_name, p.pathway_id, p.source
                FROM protein_pathway_members ppm
                JOIN pathways p ON ppm.pathway_id = p.pathway_id
                WHERE ppm.gene_symbol = ANY(%s)
                LIMIT 30
            """, (gene_symbols,))
            
            pathways = cursor.fetchall()
            profile['pathways'] = [{'name': name, 'id': pid, 'source': src} for name, pid, src in pathways]
        
        # Get known indications
        cursor.execute("""
            SELECT therapeutic_category, mechanism_of_action, fda_approved
            FROM drugs
            WHERE name = %s
        """, (drug_name,))
        
        drug_info = cursor.fetchone()
        if drug_info:
            profile['category'] = drug_info[0]
            profile['mechanism'] = drug_info[1]
            profile['fda_approved'] = drug_info[2]
        
        # Get repurposing opportunities
        cursor.execute("""
            SELECT disease_name, confidence_score, evidence_type
            FROM drug_disease_repurposing
            WHERE drug_name = %s
            ORDER BY confidence_score DESC
            LIMIT 10
        """, (drug_name,))
        
        repurposing = cursor.fetchall()
        profile['repurposing_candidates'] = [{'disease': dis, 'confidence': conf, 'evidence': ev} for dis, conf, ev in repurposing]
        
        cursor.close()
        return profile
    
    def get_disease_context(self, disease_name: str) -> Dict:
        """Get comprehensive disease context"""
        
        cursor = self.conn.cursor()
        
        context = {
            'disease_name': disease_name,
            'associated_genes': [],
            'pathways': [],
            'candidate_drugs': [],
            'clinical_trials': []
        }
        
        # Get associated genes (from 100K!)
        cursor.execute("""
            SELECT g.gene_symbol, g.gene_name
            FROM gene_disease_associations gda
            JOIN genes g ON gda.gene_symbol = g.gene_symbol
            WHERE gda.disease_name ILIKE %s
            LIMIT 50
        """, (f'%{disease_name}%',))
        
        genes = cursor.fetchall()
        context['associated_genes'] = [{'symbol': symbol, 'name': name} for symbol, name in genes]
        
        # Get disease pathways (from 200K + 15 associations!)
        cursor.execute("""
            SELECT p.pathway_name, p.pathway_id, pda.evidence_score
            FROM pathway_disease_associations pda
            JOIN pathways p ON pda.pathway_id = p.pathway_id
            WHERE pda.disease_name ILIKE %s
            ORDER BY pda.evidence_score DESC
            LIMIT 20
        """, (f'%{disease_name}%',))
        
        pathways = cursor.fetchall()
        context['pathways'] = [{'name': name, 'id': pid, 'score': score} for name, pid, score in pathways]
        
        # Get drug candidates
        cursor.execute("""
            SELECT drug_name, confidence_score, evidence_type
            FROM drug_disease_repurposing
            WHERE disease_name ILIKE %s
            ORDER BY confidence_score DESC
            LIMIT 20
        """, (f'%{disease_name}%',))
        
        drugs = cursor.fetchall()
        context['candidate_drugs'] = [{'drug': drug, 'confidence': conf, 'evidence': ev} for drug, conf, ev in drugs]
        
        # Get clinical trials
        cursor.execute("""
            SELECT nct_id, title, status, phase
            FROM clinical_trials
            WHERE disease ILIKE %s
            LIMIT 10
        """, (f'%{disease_name}%',))
        
        trials = cursor.fetchall()
        context['clinical_trials'] = [{'nct': nct, 'title': title, 'status': status, 'phase': phase} for nct, title, status, phase in trials]
        
        cursor.close()
        return context
    
    def get_pathway_information(self, pathway_name: str) -> Dict:
        """Get detailed pathway information"""
        
        cursor = self.conn.cursor()
        
        # Get pathway details
        cursor.execute("""
            SELECT pathway_id, pathway_name, source, description
            FROM pathways
            WHERE pathway_name ILIKE %s
            LIMIT 5
        """, (f'%{pathway_name}%',))
        
        pathway_info = cursor.fetchone()
        
        if not pathway_info:
            cursor.close()
            return {}
        
        pathway_id = pathway_info[0]
        
        info = {
            'pathway_id': pathway_id,
            'pathway_name': pathway_info[1],
            'source': pathway_info[2],
            'description': pathway_info[3],
            'proteins': [],
            'associated_diseases': [],
            'drugs_targeting': []
        }
        
        # Get proteins in pathway
        cursor.execute("""
            SELECT g.gene_symbol, g.gene_name
            FROM protein_pathway_members ppm
            JOIN genes g ON ppm.gene_symbol = g.gene_symbol
            WHERE ppm.pathway_id = %s
            LIMIT 50
        """, (pathway_id,))
        
        proteins = cursor.fetchall()
        info['proteins'] = [{'symbol': symbol, 'name': name} for symbol, name in proteins]
        
        # Get associated diseases
        cursor.execute("""
            SELECT disease_name, evidence_score
            FROM pathway_disease_associations
            WHERE pathway_id = %s
            ORDER BY evidence_score DESC
        """, (pathway_id,))
        
        diseases = cursor.fetchall()
        info['associated_diseases'] = [{'disease': dis, 'score': score} for dis, score in diseases]
        
        cursor.close()
        return info
    
    def search_knowledge_base(self, query: str) -> Dict:
        """
        Comprehensive search across all database tables.
        Returns relevant drugs, proteins, pathways, diseases.
        """
        
        cursor = self.conn.cursor()
        
        results = {
            'drugs': [],
            'pathways': [],
            'genes': [],
            'diseases': [],
            'interactions': []
        }
        
        query_pattern = f'%{query}%'
        
        # Search drugs
        cursor.execute("""
            SELECT name, therapeutic_category, mechanism_of_action
            FROM drugs
            WHERE name ILIKE %s OR therapeutic_category ILIKE %s OR mechanism_of_action ILIKE %s
            LIMIT 20
        """, (query_pattern, query_pattern, query_pattern))
        
        results['drugs'] = [{'name': name, 'category': cat, 'mechanism': mech} for name, cat, mech in cursor.fetchall()]
        
        # Search pathways (200K!)
        cursor.execute("""
            SELECT pathway_name, pathway_id, source
            FROM pathways
            WHERE pathway_name ILIKE %s
            LIMIT 20
        """, (query_pattern,))
        
        results['pathways'] = [{'name': name, 'id': pid, 'source': src} for name, pid, src in cursor.fetchall()]
        
        # Search genes (100K!)
        cursor.execute("""
            SELECT gene_symbol, gene_name
            FROM genes
            WHERE gene_symbol ILIKE %s OR gene_name ILIKE %s
            LIMIT 20
        """, (query_pattern, query_pattern))
        
        results['genes'] = [{'symbol': symbol, 'name': name} for symbol, name in cursor.fetchall()]
        
        # Search diseases (100K!)
        cursor.execute("""
            SELECT disease_name, disease_type
            FROM diseases
            WHERE disease_name ILIKE %s
            LIMIT 20
        """, (query_pattern,))
        
        results['diseases'] = [{'name': name, 'type': dtype} for name, dtype in cursor.fetchall()]
        
        cursor.close()
        return results
    
    def __del__(self):
        if self.conn:
            self.conn.close()


# Simple API for chatbox
def query_semantic_library(question: str) -> str:
    """
    Query the semantic library and return formatted context.
    Use this in your chatbox before sending to Gemini.
    """
    library = SemanticContentLibrary()
    results = library.search_knowledge_base(question)
    
    context = []
    
    if results['drugs']:
        context.append(f"**Drugs ({len(results['drugs'])}):**")
        for drug in results['drugs'][:5]:
            context.append(f"  • {drug['name']}: {drug['mechanism']}")
    
    if results['pathways']:
        context.append(f"\n**Pathways ({len(results['pathways'])} from 200K):**")
        for pw in results['pathways'][:5]:
            context.append(f"  • {pw['name']} ({pw['source']})")
    
    if results['genes']:
        context.append(f"\n**Genes ({len(results['genes'])} from 100K):**")
        for gene in results['genes'][:5]:
            context.append(f"  • {gene['symbol']}: {gene['name']}")
    
    if results['diseases']:
        context.append(f"\n**Diseases ({len(results['diseases'])} from 100K):**")
        for dis in results['diseases'][:5]:
            context.append(f"  • {dis['name']} ({dis['type']})")
    
    return "\n".join(context)