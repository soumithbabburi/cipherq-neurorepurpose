"""
Dynamic Target Resolver - Automatically discovers drug-target-pathway mappings
Queries ChEMBL, DrugBank, UniProt, and KEGG APIs for real-time data
"""

import requests
import json
import logging
from typing import Dict, List, Optional, Tuple
from pathlib import Path
import time

logger = logging.getLogger(__name__)

class DynamicTargetResolver:
    """
    Dynamically resolve drug targets using real biomedical APIs
    No hardcoded mappings - all data fetched from authoritative sources
    """
    
    def __init__(self, cache_dir: str = "data"):
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(exist_ok=True)
        
        # API endpoints
        self.chembl_base = "https://www.ebi.ac.uk/chembl/api/data"
        self.pubchem_base = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
        self.uniprot_base = "https://rest.uniprot.org/uniprotkb"
        self.pdb_base = "https://data.rcsb.org/rest/v1/core/entry"
        
        # Load cached data if available
        self.drugs_cache = self._load_cache("drugs_cache.json")
        self.targets_cache = self._load_cache("targets_cache.json")
        
    def _load_cache(self, filename: str) -> dict:
        """Load cached data from JSON file"""
        cache_file = self.cache_dir / filename
        if cache_file.exists():
            try:
                with open(cache_file, 'r') as f:
                    return json.load(f)
            except Exception as e:
                logger.warning(f"Could not load cache {filename}: {e}")
        return {}
    
    def _save_cache(self, data: dict, filename: str):
        """Save data to cache file"""
        cache_file = self.cache_dir / filename
        try:
            with open(cache_file, 'w') as f:
                json.dump(data, f, indent=2)
        except Exception as e:
            logger.warning(f"Could not save cache {filename}: {e}")
    
    def get_drug_targets_from_chembl(self, drug_name: str) -> List[Dict]:
        """Get drug targets from ChEMBL API and resolve to UniProt IDs"""
        try:
            # Search for drug by name
            search_url = f"{self.chembl_base}/molecule/search.json"
            params = {"q": drug_name, "limit": 5}
            
            response = requests.get(search_url, params=params, timeout=10)
            if response.status_code != 200:
                logger.warning(f"ChEMBL search failed for {drug_name}")
                return []
            
            molecules = response.json().get('molecules', [])
            if not molecules:
                return []
            
            # Get the first match (most relevant)
            molecule_id = molecules[0].get('molecule_chembl_id')
            
            # Get mechanisms of action (drug-target interactions)
            mech_url = f"{self.chembl_base}/mechanism.json"
            params = {"molecule_chembl_id": molecule_id}
            
            response = requests.get(mech_url, params=params, timeout=10)
            if response.status_code != 200:
                return []
            
            mechanisms = response.json().get('mechanisms', [])
            
            targets = []
            for mech in mechanisms:
                chembl_target_id = mech.get('target_chembl_id', '')
                
                # Resolve ChEMBL target ID to protein details
                target_details = self._resolve_chembl_target_to_protein(chembl_target_id)
                
                if target_details:
                    target_info = {
                        'target_name': target_details.get('gene_name', target_details.get('name', '')),
                        'protein_name': target_details.get('name', ''),
                        'uniprot_id': target_details.get('uniprot_id', ''),
                        'mechanism': mech.get('mechanism_of_action', ''),
                        'action_type': mech.get('action_type', ''),
                        'target_type': mech.get('target_type', ''),
                        'source': 'chembl'
                    }
                    targets.append(target_info)
            
            logger.info(f"Found {len(targets)} resolved targets for {drug_name} from ChEMBL")
            return targets
            
        except Exception as e:
            logger.error(f"ChEMBL API error for {drug_name}: {e}")
            return []
    
    def _resolve_chembl_target_to_protein(self, chembl_target_id: str) -> Optional[Dict]:
        """Resolve ChEMBL target ID to protein information with UniProt ID"""
        try:
            # Get target details from ChEMBL
            target_url = f"{self.chembl_base}/target/{chembl_target_id}.json"
            response = requests.get(target_url, timeout=10)
            
            if response.status_code != 200:
                return None
            
            target_data = response.json()
            
            # Extract protein components
            target_components = target_data.get('target_components', [])
            if not target_components:
                return None
            
            # Get first component (primary target)
            component = target_components[0]
            
            # Extract UniProt ID and gene name
            accessions = component.get('target_component_xrefs', [])
            uniprot_id = None
            gene_name = None
            
            for xref in accessions:
                if xref.get('xref_src_db') == 'UniProt':
                    uniprot_id = xref.get('xref_id', '')
                    break
            
            # Get gene names
            gene_names = component.get('target_component_synonyms', [])
            if gene_names:
                gene_name = gene_names[0].get('component_synonym', '')
            
            protein_name = component.get('component_description', '')
            
            return {
                'name': protein_name,
                'gene_name': gene_name,
                'uniprot_id': uniprot_id,
                'source': 'chembl'
            }
            
        except Exception as e:
            logger.debug(f"Failed to resolve ChEMBL target {chembl_target_id}: {e}")
            return None
    
    def get_target_protein_info(self, target_name: str) -> Optional[Dict]:
        """Get protein information from UniProt"""
        try:
            # Search UniProt for the target
            search_url = f"{self.uniprot_base}/search"
            params = {
                "query": target_name,
                "format": "json",
                "size": 1
            }
            
            response = requests.get(search_url, params=params, timeout=10)
            if response.status_code != 200:
                return None
            
            results = response.json().get('results', [])
            if not results:
                return None
            
            protein = results[0]
            
            # Extract relevant information
            protein_info = {
                'uniprot_id': protein.get('primaryAccession', ''),
                'name': protein.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', target_name),
                'organism': protein.get('organism', {}).get('scientificName', ''),
                'gene': protein.get('genes', [{}])[0].get('geneName', {}).get('value', '') if protein.get('genes') else '',
                'source': 'uniprot'
            }
            
            return protein_info
            
        except Exception as e:
            logger.error(f"UniProt API error for {target_name}: {e}")
            return None
    
    def get_pdb_structures_for_target(self, uniprot_id: str) -> List[str]:
        """Get PDB structures associated with UniProt ID"""
        try:
            # Query RCSB PDB for structures
            search_url = "https://search.rcsb.org/rcsbsearch/v2/query"
            
            query = {
                "query": {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                        "operator": "exact_match",
                        "value": uniprot_id
                    }
                },
                "return_type": "entry",
                "request_options": {
                    "return_all_hits": False,
                    "results_content_type": ["experimental"],
                    "sort": [{"sort_by": "score", "direction": "desc"}]
                }
            }
            
            response = requests.post(search_url, json=query, timeout=15)
            if response.status_code != 200:
                return []
            
            results = response.json().get('result_set', [])
            pdb_ids = [result.get('identifier', '') for result in results[:5]]  # Top 5 structures
            
            logger.info(f"Found {len(pdb_ids)} PDB structures for UniProt {uniprot_id}")
            return pdb_ids
            
        except Exception as e:
            logger.error(f"PDB search error for {uniprot_id}: {e}")
            return []
    
    def resolve_drug_target_dynamically(self, drug_name: str) -> Optional[Dict]:
        """
        Dynamically resolve drug to target protein with PDB structure
        Returns complete mapping: drug -> target -> protein -> PDB
        Fallback to 500-drug dataset if ChEMBL fails
        """
        
        # Check cache first
        drug_key = drug_name.upper().replace(' ', '_')
        if drug_key in self.drugs_cache:
            logger.info(f"Using cached target for {drug_name}")
            return self.drugs_cache[drug_key]
        
        try:
            # Step 1: Try ChEMBL API for drug targets
            targets = self.get_drug_targets_from_chembl(drug_name)
            
            if not targets:
                logger.info(f"ChEMBL found no targets for {drug_name}, trying dataset fallback...")
                # Fallback: Load from 500-drug dataset
                return self._resolve_from_dataset(drug_name)
            
            # Step 2: Get protein info from ChEMBL target (already has UniProt ID)
            primary_target = targets[0]
            uniprot_id = primary_target.get('uniprot_id', '')
            protein_name = primary_target.get('protein_name', '')
            gene_name = primary_target.get('target_name', '')
            
            if not uniprot_id:
                logger.warning(f"No UniProt ID for {drug_name}, trying dataset fallback...")
                return self._resolve_from_dataset(drug_name)
            
            # Step 3: Get PDB structures for the protein
            pdb_structures = self.get_pdb_structures_for_target(uniprot_id)
            
            if not pdb_structures:
                logger.warning(f"No PDB structures found for UniProt {uniprot_id}")
            
            # Compile complete mapping
            result = {
                'drug_name': drug_name,
                'targets': targets,
                'primary_target': {
                    'name': gene_name or protein_name,
                    'mechanism': primary_target.get('mechanism', ''),
                    'action_type': primary_target.get('action_type', '')
                },
                'protein': {
                    'name': protein_name,
                    'gene': gene_name,
                    'uniprot_id': uniprot_id,
                    'source': 'chembl'
                },
                'pdb_structures': pdb_structures,
                'preferred_pdb': pdb_structures[0] if pdb_structures else None
            }
            
            # Cache the result
            self.drugs_cache[drug_key] = result
            self._save_cache(self.drugs_cache, "drugs_cache.json")
            
            logger.info(f"✅ Dynamic resolution: {drug_name} -> {protein_name} -> PDB {result.get('preferred_pdb')}")
            return result
            
        except Exception as e:
            logger.error(f"Dynamic target resolution failed for {drug_name}: {e}")
            # Final fallback to dataset
            return self._resolve_from_dataset(drug_name)
    
    def _resolve_from_dataset(self, drug_name: str) -> Optional[Dict]:
        """Fallback: Resolve drug target from 500-drug dataset"""
        try:
            import json
            from pathlib import Path
            
            drugs_file = Path("data") / "drugs_500.json"
            # Use drug_target_proteins.json instead of generic proteins_500.json
            proteins_file = Path("data") / "drug_target_proteins.json"
            
            # Load datasets
            with open(drugs_file, 'r') as f:
                drugs = json.load(f)
            
            # Load drug target proteins
            if proteins_file.exists():
                with open(proteins_file, 'r') as f:
                    proteins = json.load(f)
            else:
                logger.warning("drug_target_proteins.json not found, using generic proteins")
                proteins_file = Path("data") / "proteins_500.json"
                with open(proteins_file, 'r') as f:
                    proteins = json.load(f)
            
            # Find drug in dataset
            drug_data = None
            for drug in drugs:
                if drug.get('name', '').lower() == drug_name.lower():
                    drug_data = drug
                    break
            
            if not drug_data:
                logger.warning(f"Drug {drug_name} not found in dataset")
                return None
            
            target_gene = drug_data.get('target', '')
            if not target_gene:
                logger.warning(f"No target for {drug_name} in dataset")
                return None
            
            # Find protein with matching gene name
            protein_data = None
            for protein in proteins:
                if protein.get('gene_name', '').upper() == target_gene.upper():
                    protein_data = protein
                    break
            
            if not protein_data:
                logger.info(f"Protein data not found for target {target_gene}, using target name only")
                # Return minimal data with just the target
                return {
                    'drug_name': drug_name,
                    'primary_target': {'name': target_gene},
                    'protein': {'gene': target_gene, 'name': target_gene},
                    'pdb_structures': [],
                    'preferred_pdb': None,
                    'source': 'dataset_minimal'
                }
            
            # Return complete protein data from dataset
            result = {
                'drug_name': drug_name,
                'primary_target': {'name': target_gene},
                'protein': {
                    'name': protein_data.get('name', ''),
                    'gene': protein_data.get('gene_name', ''),
                    'uniprot_id': protein_data.get('uniprot_id', ''),
                    'source': 'dataset'
                },
                'pdb_structures': protein_data.get('pdb_ids', []),
                'preferred_pdb': protein_data.get('pdb_ids', [None])[0],
                'source': 'dataset'
            }
            
            logger.info(f"✅ Dataset fallback: {drug_name} -> {target_gene} -> PDB {result.get('preferred_pdb')}")
            return result
            
        except Exception as e:
            logger.error(f"Dataset fallback failed for {drug_name}: {e}")
            return None
    
    def batch_resolve_drugs(self, drug_names: List[str], delay: float = 0.5) -> Dict[str, Dict]:
        """
        Resolve multiple drugs to targets (batch processing)
        """
        results = {}
        
        for i, drug_name in enumerate(drug_names):
            logger.info(f"Resolving {i+1}/{len(drug_names)}: {drug_name}")
            
            result = self.resolve_drug_target_dynamically(drug_name)
            if result:
                results[drug_name] = result
            
            # Rate limiting
            if delay > 0:
                time.sleep(delay)
        
        logger.info(f"✅ Batch resolved {len(results)}/{len(drug_names)} drugs")
        return results
