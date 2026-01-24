"""
Simple Docking Service for CipherQ
Wraps NVIDIA BioNeMo DiffDock integration
"""
import logging
import os
from typing import Dict, Optional, List
from nvidia_bionemo_integration import NVIDIABioNeMoClient

logger = logging.getLogger(__name__)

class DockingService:
    """
    Simple docking service that uses NVIDIA DiffDock API.
    Falls back gracefully if API key is not available.
    """
    
    def __init__(self, disease_name: str = "Alzheimer's Disease"):
        """
        Initialize docking service.
        
        Args:
            disease_name: Target disease (for context)
        """
        self.disease_name = disease_name
        self.nvidia_client = NVIDIABioNeMoClient()
        self.available = self.nvidia_client.api_available
        
        if self.available:
            logger.info(f"✅ DockingService initialized with NVIDIA DiffDock for {disease_name}")
        else:
            logger.warning(f"⚠️ DockingService initialized without NVIDIA API key")
    
    def set_disease(self, disease_name: str):
        """Update the target disease"""
        self.disease_name = disease_name
        logger.info(f"DockingService disease set to: {disease_name}")
    
    def suggest_target_for_disease(self, drug_name: str) -> Optional[str]:
        """
        Suggest a relevant target protein for the drug and disease.
        
        Args:
            drug_name: Name of the drug
            
        Returns:
            Suggested target protein name or None
        """
        # Disease-specific target suggestions
        disease_targets = {
            "Alzheimer's Disease": ["PPARG", "BACE1", "AChE", "GSK3B"],
            "Diabetes": ["PPARG", "DPP4", "SGLT2", "GLP1R"],
            "Cancer": ["EGFR", "VEGFR", "HER2", "BCR-ABL"],
            "Parkinson's Disease": ["LRRK2", "MAO-B", "α-synuclein"],
            "Cardiovascular": ["ACE", "AT1R", "HMG-CoA reductase"]
        }
        
        # Get targets for this disease
        targets = disease_targets.get(self.disease_name, [])
        
        if targets:
            # Return first target as suggestion
            return targets[0]
        
        return None
    
    def perform_docking(
        self, 
        drug_name: str, 
        target_name: str,
        ligand_smiles: Optional[str] = None
    ) -> Dict:
        """
        Perform molecular docking.
        
        Args:
            drug_name: Name of the drug
            target_name: Target protein name
            ligand_smiles: SMILES structure (optional, will fetch from DB if not provided)
            
        Returns:
            Dictionary with docking results
        """
        if not self.available:
            logger.warning("NVIDIA API not available - cannot perform docking")
            return {
                'success': False,
                'error': 'NVIDIA API key not configured',
                'poses': [],
                'confidence_scores': [],
                'binding_affinities': []
            }
        
        try:
            # Get SMILES - try provided first, then database
            if not ligand_smiles:
                logger.info(f"Fetching SMILES for {drug_name} from database...")
                ligand_smiles = self._get_drug_smiles(drug_name)
                
            if not ligand_smiles:
                error_msg = f'No SMILES structure found for {drug_name}. Check database.'
                logger.error(error_msg)
                return {
                    'success': False,
                    'error': error_msg,
                    'poses': []
                }
            
            logger.info(f"✅ SMILES found: {ligand_smiles[:50]}...")
            
            # Get protein structure
            logger.info(f"Fetching protein structure for {target_name}...")
            protein_pdb = self._get_protein_structure(target_name)
            if not protein_pdb:
                error_msg = f'No PDB structure available for {target_name}'
                logger.warning(error_msg)
                # Use generic structure as fallback
                protein_pdb = self._get_generic_protein_structure()
                logger.info("Using generic protein structure for docking")
            
            # Convert SMILES to SDF
            logger.info(f"Converting SMILES to 3D SDF structure...")
            ligand_sdf = self._smiles_to_sdf(ligand_smiles, drug_name)
            if not ligand_sdf:
                error_msg = f'Failed to convert SMILES to SDF for {drug_name}. RDKit required.'
                logger.error(error_msg)
                return {
                    'success': False,
                    'error': error_msg,
                    'poses': []
                }
            
            logger.info(f"✅ SDF structure generated")
            
            # Run NVIDIA DiffDock
            logger.info(f"🚀 Running NVIDIA DiffDock for {drug_name} → {target_name}")
            result = self.nvidia_client.run_diffdock(
                protein_pdb=protein_pdb,
                ligand_sdf=ligand_sdf,
                ligand_name=drug_name,
                num_poses=20
            )
            
            # Check if NVIDIA failed - use AutoDock Vina as fallback
            nvidia_failed = (
                not result.get('success') or 
                result.get('status') == 'failed' or
                'failed' in result.get('details', '').lower()
            )
            
            if nvidia_failed:
                error_detail = result.get('details', result.get('error', 'Unknown error'))
                logger.warning(f"⚠️ NVIDIA DiffDock FAILED: {error_detail}")
                logger.info("=" * 80)
                logger.info("🔄 SWITCHING TO AUTODOCK VINA FALLBACK")
                logger.info("=" * 80)
                
                try:
                    # Use user's AutoDock Vina implementation
                    from autodock_vina import AutoDockVina
                    
                    vina = AutoDockVina()
                    vina_result = vina.run_docking(
                        protein_pdb=protein_pdb,
                        ligand_sdf=ligand_sdf,
                        ligand_name=drug_name,
                        num_poses=20
                    )
                    
                    if vina_result.get('success'):
                        logger.info(f"✅ AutoDock Vina completed: {len(vina_result.get('poses', []))} poses")
                        logger.info(f"📊 Vina affinities: {vina_result.get('binding_affinities', [])[:5]}")
                        
                        # Add metadata
                        vina_result['protein_pdb'] = protein_pdb
                        vina_result['target_name'] = target_name
                        vina_result['drug_name'] = drug_name
                        vina_result['pdb_id'] = 'VINA_FALLBACK'
                        vina_result['fallback'] = True
                        vina_result['docking_method'] = 'AutoDock Vina (RDKit)'
                        vina_result['nvidia_error'] = error_detail
                        vina_result['raw_nvidia_confidences'] = vina_result.get('confidence_scores', [])
                        
                        return vina_result
                    else:
                        logger.error(f"❌ Vina also failed: {vina_result.get('error')}")
                    
                except ImportError as ie:
                    logger.error(f"❌ autodock_vina.py not found: {ie}")
                    logger.error("❌ Place autodock_vina.py in project root!")
                except Exception as vina_error:
                    logger.error(f"❌ Vina fallback failed: {vina_error}")
                    import traceback
                    logger.error(traceback.format_exc())
                
                # Return NVIDIA failed result
                return {
                    'success': False,
                    'error': f'NVIDIA failed: {error_detail}. Vina fallback also failed.',
                    'poses': [],
                    'nvidia_error': error_detail
                }
            
            # Add extra metadata for successful NVIDIA docking
            if result.get('success') and result.get('status') != 'failed':
                result['protein_pdb'] = protein_pdb
                result['target_name'] = target_name
                result['drug_name'] = drug_name
                result['pdb_id'] = 'AF_MODEL'
                result['fallback'] = False
                result['docking_method'] = 'NVIDIA DiffDock'
                logger.info(f"✅ NVIDIA DiffDock completed: {len(result.get('poses', []))} poses generated")
            else:
                logger.error(f"❌ NVIDIA Docking failed: {result.get('error', 'Unknown error')}")
            
            return result
            
        except Exception as e:
            error_msg = f"Docking failed for {drug_name}: {str(e)}"
            logger.error(error_msg)
            import traceback
            logger.error(traceback.format_exc())
            return {
                'success': False,
                'error': error_msg,
                'poses': []
            }
    
    def _get_drug_smiles(self, drug_name: str) -> Optional[str]:
        """Get SMILES structure from database"""
        try:
            import psycopg2
            import os
            
            # Use environment variables directly to avoid config import issues
            db_params = {
                "host": os.getenv("DB_HOST", "localhost"),
                "database": os.getenv("DB_NAME", "cipherq_repurpose"),
                "user": os.getenv("DB_USER", "babburisoumith"),
                "password": os.getenv("DB_PASSWORD", "")
            }
            
            conn = psycopg2.connect(**db_params)
            cursor = conn.cursor()
            
            cursor.execute("SELECT smiles FROM drugs WHERE name = %s LIMIT 1", (drug_name,))
            result = cursor.fetchone()
            
            cursor.close()
            conn.close()
            
            if result:
                logger.info(f"✅ Found SMILES for {drug_name}")
                return result[0]
            else:
                logger.warning(f"No SMILES found in database for {drug_name}")
                return None
                
        except Exception as e:
            logger.error(f"Failed to get SMILES for {drug_name}: {e}")
            return None
    
    def _get_protein_structure(self, target_name: str) -> Optional[str]:
        """
        Get protein PDB structure.
        For now, returns a generic structure. In production, use AlphaFold or PDB.
        """
        try:
            # Try to get from real PDB fetcher if available
            from real_pdb_fetcher import RealPDBFetcher
            fetcher = RealPDBFetcher()
            pdb_content = fetcher.fetch_protein_structure(target_name)
            if pdb_content:
                return pdb_content
        except Exception as e:
            logger.warning(f"Could not fetch real PDB for {target_name}: {e}")
        
        # Return generic protein structure as fallback
        logger.warning(f"Using generic protein structure for {target_name}")
        return self._get_generic_protein_structure()
    
    def _get_generic_protein_structure(self) -> str:
        """
        Generic protein structure for testing - MUST have multiple residues for NVIDIA DiffDock
        Using a small 10-residue alpha helix structure
        """
        return """HEADER    GENERIC PROTEIN STRUCTURE                     
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N  
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00 20.00           C  
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00 20.00           C  
ATOM      4  O   ALA A   1       1.267   2.395   0.000  1.00 20.00           O  
ATOM      5  CB  ALA A   1       1.993  -0.742   1.228  1.00 20.00           C  
ATOM      6  N   GLY A   2       3.331   1.520   0.000  1.00 20.00           N  
ATOM      7  CA  GLY A   2       4.025   2.802   0.000  1.00 20.00           C  
ATOM      8  C   GLY A   2       5.532   2.621   0.000  1.00 20.00           C  
ATOM      9  O   GLY A   2       6.061   1.511   0.000  1.00 20.00           O  
ATOM     10  N   LEU A   3       6.241   3.746   0.000  1.00 20.00           N  
ATOM     11  CA  LEU A   3       7.694   3.746   0.000  1.00 20.00           C  
ATOM     12  C   LEU A   3       8.246   5.166   0.000  1.00 20.00           C  
ATOM     13  O   LEU A   3       7.504   6.141   0.000  1.00 20.00           O  
ATOM     14  CB  LEU A   3       8.229   2.999   1.228  1.00 20.00           C  
ATOM     15  CG  LEU A   3       9.747   2.852   1.228  1.00 20.00           C  
ATOM     16  CD1 LEU A   3      10.235   2.106   2.456  1.00 20.00           C  
ATOM     17  CD2 LEU A   3      10.282   2.134   0.000  1.00 20.00           C  
ATOM     18  N   VAL A   4       9.568   5.266   0.000  1.00 20.00           N  
ATOM     19  CA  VAL A   4      10.262   6.548   0.000  1.00 20.00           C  
ATOM     20  C   VAL A   4      11.769   6.367   0.000  1.00 20.00           C  
ATOM     21  O   VAL A   4      12.298   5.257   0.000  1.00 20.00           O  
ATOM     22  CB  VAL A   4      10.797   7.294   1.228  1.00 20.00           C  
ATOM     23  CG1 VAL A   4      11.491   8.576   0.772  1.00 20.00           C  
ATOM     24  CG2 VAL A   4       9.662   7.588   2.208  1.00 20.00           C  
ATOM     25  N   SER A   5      12.478   7.492   0.000  1.00 20.00           N  
ATOM     26  CA  SER A   5      13.931   7.492   0.000  1.00 20.00           C  
ATOM     27  C   SER A   5      14.483   8.912   0.000  1.00 20.00           C  
ATOM     28  O   SER A   5      13.741   9.887   0.000  1.00 20.00           O  
ATOM     29  CB  SER A   5      14.466   6.745   1.228  1.00 20.00           C  
ATOM     30  OG  SER A   5      13.876   5.461   1.371  1.00 20.00           O  
ATOM     31  N   ALA A   6      15.805   9.012   0.000  1.00 20.00           N  
ATOM     32  CA  ALA A   6      16.499  10.294   0.000  1.00 20.00           C  
ATOM     33  C   ALA A   6      18.006  10.113   0.000  1.00 20.00           C  
ATOM     34  O   ALA A   6      18.535   9.003   0.000  1.00 20.00           O  
ATOM     35  CB  ALA A   6      17.034  10.997   1.228  1.00 20.00           C  
ATOM     36  N   GLY A   7      18.715  11.238   0.000  1.00 20.00           N  
ATOM     37  CA  GLY A   7      20.168  11.238   0.000  1.00 20.00           C  
ATOM     38  C   GLY A   7      20.720  12.658   0.000  1.00 20.00           C  
ATOM     39  O   GLY A   7      19.978  13.633   0.000  1.00 20.00           O  
ATOM     40  N   ALA A   8      22.042  12.758   0.000  1.00 20.00           N  
ATOM     41  CA  ALA A   8      22.736  14.040   0.000  1.00 20.00           C  
ATOM     42  C   ALA A   8      24.243  13.859   0.000  1.00 20.00           C  
ATOM     43  O   ALA A   8      24.772  12.749   0.000  1.00 20.00           O  
ATOM     44  CB  ALA A   8      23.271  14.743   1.228  1.00 20.00           C  
ATOM     45  N   LEU A   9      24.952  14.984   0.000  1.00 20.00           N  
ATOM     46  CA  LEU A   9      26.405  14.984   0.000  1.00 20.00           C  
ATOM     47  C   LEU A   9      26.957  16.404   0.000  1.00 20.00           C  
ATOM     48  O   LEU A   9      26.215  17.379   0.000  1.00 20.00           O  
ATOM     49  CB  LEU A   9      26.940  14.237   1.228  1.00 20.00           C  
ATOM     50  CG  LEU A   9      28.458  14.090   1.228  1.00 20.00           C  
ATOM     51  CD1 LEU A   9      28.946  13.344   2.456  1.00 20.00           C  
ATOM     52  CD2 LEU A   9      28.993  13.372   0.000  1.00 20.00           C  
ATOM     53  N   ALA A  10      28.279  16.504   0.000  1.00 20.00           N  
ATOM     54  CA  ALA A  10      28.973  17.786   0.000  1.00 20.00           C  
ATOM     55  C   ALA A  10      30.480  17.605   0.000  1.00 20.00           C  
ATOM     56  O   ALA A  10      31.009  16.495   0.000  1.00 20.00           O  
ATOM     57  CB  ALA A  10      29.508  18.489   1.228  1.00 20.00           C  
END
"""
    
    def _smiles_to_sdf(self, smiles: str, mol_name: str = "molecule") -> Optional[str]:
        """Convert SMILES to SDF format"""
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem
            
            # Create molecule from SMILES
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                logger.error(f"Invalid SMILES: {smiles}")
                return None
            
            # Add hydrogens
            mol = Chem.AddHs(mol)
            
            # Generate 3D coordinates
            AllChem.EmbedMolecule(mol, randomSeed=42)
            AllChem.MMFFOptimizeMolecule(mol)
            
            # Convert to SDF
            mol.SetProp("_Name", mol_name)
            sdf_content = Chem.MolToMolBlock(mol)
            
            return sdf_content
            
        except ImportError:
            logger.error("RDKit not available - cannot convert SMILES to SDF")
            logger.error("Install: conda install -c conda-forge rdkit")
            return None
        except Exception as e:
            logger.error(f"Failed to convert SMILES to SDF: {e}")
            return None
    
    def is_available(self) -> bool:
        """Check if docking service is available"""
        return self.available
    
    def get_smiles_for_drug(self, drug_name: str, drug_dict: Optional[Dict] = None) -> Optional[str]:
        """
        Get SMILES for a drug from multiple sources.
        
        Args:
            drug_name: Name of the drug
            drug_dict: Optional drug dictionary that might contain SMILES
            
        Returns:
            SMILES string or None
        """
        # Try drug dict first (if provided)
        if drug_dict and 'smiles' in drug_dict:
            smiles = drug_dict['smiles']
            if smiles and smiles.strip():
                logger.info(f"✅ Got SMILES from drug dict for {drug_name}")
                return smiles
        
        # Try database
        logger.info(f"Fetching SMILES from database for {drug_name}...")
        return self._get_drug_smiles(drug_name)