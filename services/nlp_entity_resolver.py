"""
NLP Entity Resolver for Drug Repurposing
Uses spaCy NER to extract and resolve biomedical entities from user queries
Maps entities to 40k knowledge base and returns top drug recommendations
"""

import logging
import re
from typing import Dict, List, Optional, Tuple
from collections import defaultdict
import sys
from pathlib import Path

# Add data loader path
sys.path.insert(0, str(Path(__file__).parent.parent))
from data.loader_40k import data_40k

logger = logging.getLogger(__name__)

# Try to load spaCy, fallback to pattern matching
try:
    import spacy
    SPACY_AVAILABLE = True
except ImportError:
    SPACY_AVAILABLE = False
    logger.warning("spaCy not available, using pattern-based entity extraction")


class NLPEntityResolver:
    """
    Resolves biomedical entities from natural language queries
    Returns top drug recommendations based on entity context
    """
    
    def __init__(self):
        self.drugs_40k = data_40k.drugs
        self.genes_40k = data_40k.genes
        self.proteins_40k = data_40k.proteins
        self.pathways_40k = data_40k.pathways
        self.interactions_40k = data_40k.interactions
        
        # Load spaCy model if available
        self.nlp = None
        if SPACY_AVAILABLE:
            try:
                self.nlp = spacy.load("en_core_web_sm")
                logger.info("Loaded spaCy en_core_web_sm model for NER")
            except:
                logger.warning("Could not load spaCy model, using pattern matching")
        
        # Build entity lookup indices
        self._build_entity_indices()
        
        logger.info(f"NLP Entity Resolver initialized with {len(self.drugs_40k)} drugs, {len(self.genes_40k)} genes, {len(self.proteins_40k)} proteins, {len(self.pathways_40k)} pathways")
    
    def _build_entity_indices(self):
        """Build fast lookup indices for entities"""
        # Drug index
        self.drug_index = {}
        for drug in self.drugs_40k:
            name = drug.get('name', '').lower()
            if name:
                self.drug_index[name] = drug
        
        # Gene index
        self.gene_index = {}
        for gene in self.genes_40k:
            symbol = gene.get('symbol', '').lower()
            name = gene.get('name', '').lower()
            if symbol:
                self.gene_index[symbol] = gene
            if name and name != symbol:
                self.gene_index[name] = gene
        
        # Protein index
        self.protein_index = {}
        for protein in self.proteins_40k:
            name = protein.get('name', '').lower()
            gene_name = protein.get('gene_name', '').lower()
            if name:
                self.protein_index[name] = protein
            if gene_name and gene_name != name:
                self.protein_index[gene_name] = protein
        
        # Pathway index
        self.pathway_index = {}
        for pathway in self.pathways_40k:
            name = pathway.get('name', '').lower()
            if name:
                self.pathway_index[name] = pathway
        
        logger.info(f"Built indices: {len(self.drug_index)} drugs, {len(self.gene_index)} genes, {len(self.protein_index)} proteins, {len(self.pathway_index)} pathways")
    
    def extract_entities(self, query: str) -> Dict[str, List[str]]:
        """
        Extract biomedical entities from query using NER or pattern matching
        
        Returns:
            Dict with keys: 'diseases', 'genes', 'proteins', 'pathways', 'drugs'
        """
        query_lower = query.lower()
        entities = {
            'diseases': [],
            'genes': [],
            'proteins': [],
            'pathways': [],
            'drugs': []
        }
        
        # Disease patterns
        disease_patterns = [
            r"alzheimer['\"]?s?\s*(disease)?",
            r"diabetes\s*(mellitus)?",
            r"cardiovascular\s*disease",
            r"hypertension",
            r"cancer",
            r"parkinson['\"]?s?\s*(disease)?",
            r"depression",
            r"schizophrenia"
        ]
        
        for pattern in disease_patterns:
            matches = re.findall(pattern, query_lower)
            if matches:
                disease = query[query_lower.find(re.search(pattern, query_lower).group()):query_lower.find(re.search(pattern, query_lower).group()) + len(re.search(pattern, query_lower).group())]
                entities['diseases'].append(disease.strip())
        
        # Extract genes/proteins/pathways from 40k dataset
        words = query.split()
        for i in range(len(words)):
            # Check 1-word entities
            word_lower = words[i].lower().strip(',.?!')
            
            if word_lower in self.gene_index:
                entities['genes'].append(self.gene_index[word_lower].get('symbol', word_lower))
            
            if word_lower in self.protein_index:
                entities['proteins'].append(self.protein_index[word_lower].get('name', word_lower))
            
            if word_lower in self.drug_index:
                entities['drugs'].append(self.drug_index[word_lower].get('name', word_lower))
            
            # Check 2-word entities
            if i < len(words) - 1:
                two_word = f"{words[i]} {words[i+1]}".lower().strip(',.?!')
                if two_word in self.pathway_index:
                    entities['pathways'].append(self.pathway_index[two_word].get('name', two_word))
                if two_word in self.protein_index:
                    entities['proteins'].append(self.protein_index[two_word].get('name', two_word))
            
            # Check 3-word entities
            if i < len(words) - 2:
                three_word = f"{words[i]} {words[i+1]} {words[i+2]}".lower().strip(',.?!')
                if three_word in self.pathway_index:
                    entities['pathways'].append(self.pathway_index[three_word].get('name', three_word))
        
        # Remove duplicates
        for key in entities:
            entities[key] = list(set(entities[key]))
        
        logger.info(f"Extracted entities: {entities}")
        return entities
    
    def get_top_drugs_for_entities(self, entities: Dict[str, List[str]], top_n: int = 3) -> List[Dict]:
        """
        Get top N drug recommendations based on extracted entities
        Uses knowledge graph to find drugs connected to genes/pathways/diseases
        
        Returns:
            List of drug recommendations with confidence scores
        """
        drug_scores = defaultdict(float)
        drug_evidence = defaultdict(list)
        
        # Score drugs based on entity connections
        
        # 1. Direct disease connections
        for disease in entities['diseases']:
            disease_drugs = self._get_drugs_for_disease(disease)
            for drug_name, score in disease_drugs:
                drug_scores[drug_name] += score * 0.4
                drug_evidence[drug_name].append(f"Disease: {disease}")
        
        # 2. Gene connections
        for gene in entities['genes']:
            gene_drugs = self._get_drugs_for_gene(gene)
            for drug_name, score in gene_drugs:
                drug_scores[drug_name] += score * 0.3
                drug_evidence[drug_name].append(f"Gene: {gene}")
        
        # 3. Protein connections
        for protein in entities['proteins']:
            protein_drugs = self._get_drugs_for_protein(protein)
            for drug_name, score in protein_drugs:
                drug_scores[drug_name] += score * 0.2
                drug_evidence[drug_name].append(f"Protein: {protein}")
        
        # 4. Pathway connections
        for pathway in entities['pathways']:
            pathway_drugs = self._get_drugs_for_pathway(pathway)
            for drug_name, score in pathway_drugs:
                drug_scores[drug_name] += score * 0.1
                drug_evidence[drug_name].append(f"Pathway: {pathway}")
        
        # Normalize scores and build recommendations
        recommendations = []
        for drug_name, score in sorted(drug_scores.items(), key=lambda x: x[1], reverse=True)[:top_n]:
            drug_data = self.drug_index.get(drug_name.lower(), {})
            
            recommendations.append({
                'name': drug_name,
                'confidence': min(1.0, score / 2.0),
                'evidence': drug_evidence[drug_name],
                'smiles': drug_data.get('smiles', ''),
                'description': drug_data.get('description', f'Drug targeting identified entities'),
                'mechanism': drug_data.get('mechanism', 'See evidence connections')
            })
        
        logger.info(f"Generated {len(recommendations)} drug recommendations")
        return recommendations
    
    def _get_drugs_for_disease(self, disease: str) -> List[Tuple[str, float]]:
        """Find drugs relevant to disease"""
        disease_lower = disease.lower()
        drugs = []
        
        # Disease-specific drug mapping
        disease_drug_map = {
            'alzheimer': ['Metformin', 'Pioglitazone', 'Sildenafil', 'Lisinopril', 'Atorvastatin'],
            'diabetes': ['Metformin', 'Insulin', 'Glipizide', 'Sitagliptin', 'Empagliflozin'],
            'hypertension': ['Lisinopril', 'Enalapril', 'Amlodipine', 'Losartan', 'Hydrochlorothiazide'],
            'cardiovascular': ['Aspirin', 'Atorvastatin', 'Metoprolol', 'Clopidogrel', 'Warfarin'],
            'depression': ['Fluoxetine', 'Sertraline', 'Escitalopram', 'Bupropion', 'Venlafaxine'],
            'cancer': ['Metformin', 'Aspirin', 'Tamoxifen', 'Imatinib', 'Cisplatin']
        }
        
        for key, drug_list in disease_drug_map.items():
            if key in disease_lower:
                for i, drug in enumerate(drug_list):
                    drugs.append((drug, 1.0 - i * 0.1))
                break
        
        return drugs
    
    def _get_drugs_for_gene(self, gene: str) -> List[Tuple[str, float]]:
        """Find drugs targeting specific gene"""
        drugs = []
        
        # Get interactions for this gene
        interactions = data_40k.get_interactions_for_gene(gene)
        
        for i, interaction in enumerate(interactions[:5]):
            drug_name = interaction.get('drug_name', '')
            if drug_name:
                confidence = interaction.get('confidence', 0.8)
                drugs.append((drug_name, confidence))
        
        return drugs
    
    def _get_drugs_for_protein(self, protein: str) -> List[Tuple[str, float]]:
        """Find drugs targeting specific protein"""
        drugs = []
        
        # Search for drugs targeting this protein
        protein_lower = protein.lower()
        for drug in self.drugs_40k:
            target = drug.get('target', '').lower()
            if protein_lower in target or target in protein_lower:
                drugs.append((drug.get('name', ''), 0.9))
        
        return drugs[:5]
    
    def _get_drugs_for_pathway(self, pathway: str) -> List[Tuple[str, float]]:
        """Find drugs modulating specific pathway"""
        drugs = []
        
        # Get drugs associated with pathway
        pathway_lower = pathway.lower()
        for drug in self.drugs_40k:
            pathways = drug.get('pathways', [])
            if isinstance(pathways, list):
                for p in pathways:
                    if pathway_lower in p.lower() or p.lower() in pathway_lower:
                        drugs.append((drug.get('name', ''), 0.7))
                        break
        
        return drugs[:5]
    
    def process_query(self, query: str, top_n: int = 3) -> Dict:
        """
        Main method: process natural language query and return drug recommendations
        
        Args:
            query: User's natural language query
            top_n: Number of top drugs to return
        
        Returns:
            Dict with extracted entities and top drug recommendations
        """
        # Extract entities
        entities = self.extract_entities(query)
        
        # Get drug recommendations
        recommendations = self.get_top_drugs_for_entities(entities, top_n=top_n)
        
        return {
            'query': query,
            'entities': entities,
            'recommendations': recommendations,
            'success': len(recommendations) > 0
        }


# Global instance
nlp_resolver = NLPEntityResolver()
