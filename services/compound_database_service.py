"""
Chemical Compound Database Service
Clean, structured service for loading and querying chemical compounds
Integrates with PubChem API for compound data
"""

import os
import logging
import requests
from typing import List, Dict, Optional
from datetime import datetime

logger = logging.getLogger(__name__)

# PubChem API base URLs
PUBCHEM_API_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
PUBCHEM_COMPOUND_URL = f"{PUBCHEM_API_BASE}/compound"


class CompoundDatabaseService:
    """
    Service for managing chemical compounds database
    Loads from PubChem API and stores in PostgreSQL
    """
    
    def __init__(self):
        self.session = requests.Session()
        self.session.headers.update({
            'Accept': 'application/json',
            'User-Agent': 'CipherQ Drug Discovery Platform'
        })
        self._db_session = None
    
    def fetch_compound_by_name(self, name: str) -> Optional[Dict]:
        """
        Fetch compound from PubChem by name
        
        Args:
            name: Compound name to search
            
        Returns:
            Compound dictionary or None
        """
        try:
            # First get CID from name
            url = f"{PUBCHEM_COMPOUND_URL}/name/{name}/cids/JSON"
            response = self.session.get(url, timeout=15)
            
            if response.status_code != 200:
                return None
            
            data = response.json()
            cids = data.get('IdentifierList', {}).get('CID', [])
            
            if not cids:
                return None
            
            cid = cids[0]
            return self.fetch_compound_by_cid(cid)
            
        except Exception as e:
            logger.warning(f"Error fetching compound {name}: {e}")
            return None
    
    def fetch_compound_by_cid(self, cid: int) -> Optional[Dict]:
        """
        Fetch compound details from PubChem by CID
        
        Args:
            cid: PubChem Compound ID
            
        Returns:
            Compound dictionary or None
        """
        try:
            # Fetch properties
            properties = [
                'MolecularFormula', 'MolecularWeight', 'CanonicalSMILES',
                'IsomericSMILES', 'InChI', 'InChIKey', 'XLogP',
                'ExactMass', 'MonoisotopicMass', 'TPSA', 'Complexity',
                'HBondDonorCount', 'HBondAcceptorCount', 'RotatableBondCount',
                'HeavyAtomCount', 'AtomStereoCount', 'DefinedAtomStereoCount',
                'IUPACName'
            ]
            
            url = f"{PUBCHEM_COMPOUND_URL}/cid/{cid}/property/{','.join(properties)}/JSON"
            response = self.session.get(url, timeout=15)
            
            if response.status_code != 200:
                return None
            
            data = response.json()
            props = data.get('PropertyTable', {}).get('Properties', [{}])[0]
            
            compound = {
                'name': props.get('IUPACName') or f"CID_{cid}",
                'pubchem_cid': cid,
                'smiles': props.get('CanonicalSMILES'),
                'inchi': props.get('InChI'),
                'inchi_key': props.get('InChIKey'),
                'molecular_formula': props.get('MolecularFormula'),
                'molecular_weight': self._safe_float(props.get('MolecularWeight')),
                'logp': self._safe_float(props.get('XLogP')),
                'polar_surface_area': self._safe_float(props.get('TPSA')),
                'h_bond_donors': self._safe_int(props.get('HBondDonorCount')),
                'h_bond_acceptors': self._safe_int(props.get('HBondAcceptorCount')),
                'rotatable_bonds': self._safe_int(props.get('RotatableBondCount')),
                'heavy_atom_count': self._safe_int(props.get('HeavyAtomCount')),
                'source': 'PubChem',
                'source_id': str(cid)
            }
            
            # Calculate drug-likeness
            compound['lipinski_violations'] = self._calculate_lipinski_violations(compound)
            compound['druglikeness_score'] = self._calculate_druglikeness(compound)
            
            return compound
            
        except Exception as e:
            logger.warning(f"Error fetching compound CID {cid}: {e}")
            return None
    
    def fetch_compound_by_smiles(self, smiles: str) -> Optional[Dict]:
        """
        Fetch compound from PubChem by SMILES
        
        Args:
            smiles: SMILES string
            
        Returns:
            Compound dictionary or None
        """
        try:
            # Get CID from SMILES
            url = f"{PUBCHEM_COMPOUND_URL}/smiles/{smiles}/cids/JSON"
            response = self.session.get(url, timeout=15)
            
            if response.status_code != 200:
                # Calculate properties locally if not in PubChem
                return self._calculate_properties_locally(smiles)
            
            data = response.json()
            cids = data.get('IdentifierList', {}).get('CID', [])
            
            if not cids:
                return self._calculate_properties_locally(smiles)
            
            return self.fetch_compound_by_cid(cids[0])
            
        except Exception as e:
            logger.warning(f"Error fetching compound by SMILES: {e}")
            return self._calculate_properties_locally(smiles)
    
    def _calculate_properties_locally(self, smiles: str) -> Optional[Dict]:
        """Calculate compound properties locally using RDKit"""
        try:
            from rdkit import Chem
            from rdkit.Chem import Descriptors, Crippen, Lipinski
            
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None
            
            compound = {
                'name': f"Custom_{hash(smiles) % 10000}",
                'smiles': smiles,
                'molecular_weight': Descriptors.ExactMolWt(mol),
                'logp': Crippen.MolLogP(mol),
                'polar_surface_area': Descriptors.TPSA(mol),
                'h_bond_donors': Lipinski.NumHDonors(mol),
                'h_bond_acceptors': Lipinski.NumHAcceptors(mol),
                'rotatable_bonds': Descriptors.NumRotatableBonds(mol),
                'heavy_atom_count': mol.GetNumHeavyAtoms(),
                'aromatic_rings': Descriptors.NumAromaticRings(mol),
                'source': 'Local',
                'source_id': None
            }
            
            compound['lipinski_violations'] = self._calculate_lipinski_violations(compound)
            compound['druglikeness_score'] = self._calculate_druglikeness(compound)
            
            return compound
            
        except Exception as e:
            logger.error(f"Error calculating properties locally: {e}")
            return None
    
    def _calculate_lipinski_violations(self, compound: Dict) -> int:
        """Calculate Lipinski Rule of 5 violations"""
        violations = 0
        
        mw = compound.get('molecular_weight') or 0
        logp = compound.get('logp') or 0
        hbd = compound.get('h_bond_donors') or 0
        hba = compound.get('h_bond_acceptors') or 0
        
        if mw > 500:
            violations += 1
        if logp > 5:
            violations += 1
        if hbd > 5:
            violations += 1
        if hba > 10:
            violations += 1
        
        return violations
    
    def _calculate_druglikeness(self, compound: Dict) -> float:
        """
        Calculate drug-likeness score (0-100)
        Based on Lipinski, Veber, and other rules
        """
        score = 100.0
        
        mw = compound.get('molecular_weight') or 0
        logp = compound.get('logp') or 0
        psa = compound.get('polar_surface_area') or 0
        hbd = compound.get('h_bond_donors') or 0
        hba = compound.get('h_bond_acceptors') or 0
        rb = compound.get('rotatable_bonds') or 0
        
        # Lipinski penalties
        if mw > 500:
            score -= 15
        elif mw > 600:
            score -= 30
        
        if logp > 5:
            score -= 15
        elif logp < -1:
            score -= 10
        
        if hbd > 5:
            score -= 10
        
        if hba > 10:
            score -= 10
        
        # Veber penalties
        if psa > 140:
            score -= 15
        
        if rb > 10:
            score -= 10
        
        return max(0, min(100, score))
    
    def _safe_float(self, value) -> Optional[float]:
        """Safely convert to float"""
        try:
            return float(value) if value is not None else None
        except (ValueError, TypeError):
            return None
    
    def _safe_int(self, value) -> Optional[int]:
        """Safely convert to int"""
        try:
            return int(value) if value is not None else None
        except (ValueError, TypeError):
            return None
    
    def save_compound(self, compound: Dict) -> bool:
        """Save compound to database"""
        from database.connection import db_session
        from database.models import ChemicalCompound
        
        try:
            with db_session() as session:
                # Check if compound already exists
                existing = None
                if compound.get('pubchem_cid'):
                    existing = session.query(ChemicalCompound).filter_by(
                        pubchem_cid=compound['pubchem_cid']
                    ).first()
                
                if not existing:
                    db_compound = ChemicalCompound(**compound)
                    session.add(db_compound)
                    return True
                return False
                
        except Exception as e:
            logger.error(f"Error saving compound: {e}")
            return False
    
    def search_compounds(
        self,
        query: str = None,
        min_mw: float = None,
        max_mw: float = None,
        min_druglikeness: float = None,
        limit: int = 50
    ) -> List[Dict]:
        """
        Search compounds in database
        
        Args:
            query: Search term (name, SMILES)
            min_mw: Minimum molecular weight
            max_mw: Maximum molecular weight
            min_druglikeness: Minimum drug-likeness score
            limit: Maximum results
            
        Returns:
            List of matching compounds
        """
        from database.connection import db_session
        from database.models import ChemicalCompound
        
        try:
            with db_session() as session:
                q = session.query(ChemicalCompound)
                
                if query:
                    search = f"%{query.lower()}%"
                    q = q.filter(
                        (ChemicalCompound.name.ilike(search)) |
                        (ChemicalCompound.smiles.ilike(search))
                    )
                
                if min_mw is not None:
                    q = q.filter(ChemicalCompound.molecular_weight >= min_mw)
                
                if max_mw is not None:
                    q = q.filter(ChemicalCompound.molecular_weight <= max_mw)
                
                if min_druglikeness is not None:
                    q = q.filter(ChemicalCompound.druglikeness_score >= min_druglikeness)
                
                compounds = q.limit(limit).all()
                return [c.to_dict() for c in compounds]
                
        except Exception as e:
            logger.error(f"Error searching compounds: {e}")
            return []
    
    def get_compound_count(self) -> int:
        """Get total number of compounds in database"""
        from database.connection import db_session
        from database.models import ChemicalCompound
        
        try:
            with db_session() as session:
                return session.query(ChemicalCompound).count()
        except Exception as e:
            logger.error(f"Error getting compound count: {e}")
            return 0
    
    def load_druglike_compounds_from_pubchem(
        self, 
        max_compounds: int = 10000,
        min_druglikeness: float = 70
    ) -> int:
        """
        Load drug-like compounds from PubChem
        Uses PubChem's property filters
        
        Args:
            max_compounds: Maximum compounds to load
            min_druglikeness: Minimum drug-likeness score
            
        Returns:
            Number of compounds loaded
        """
        logger.info(f"Loading up to {max_compounds} drug-like compounds from PubChem...")
        
        # Note: This is a simplified version. For production, you'd use
        # PubChem's bulk download or power user gateway
        total_loaded = 0
        
        # Example CIDs of drug-like compounds (in production, use PUG REST queries)
        # For now, we'll load a sample set
        sample_cids = list(range(1, min(max_compounds, 1000) + 1))
        
        for cid in sample_cids:
            compound = self.fetch_compound_by_cid(cid)
            if compound and compound.get('druglikeness_score', 0) >= min_druglikeness:
                if self.save_compound(compound):
                    total_loaded += 1
            
            if total_loaded >= max_compounds:
                break
        
        logger.info(f"Loaded {total_loaded} drug-like compounds from PubChem")
        return total_loaded


# Singleton instance
_compound_service = None

def get_compound_database_service() -> CompoundDatabaseService:
    """Get singleton compound database service"""
    global _compound_service
    if _compound_service is None:
        _compound_service = CompoundDatabaseService()
    return _compound_service
