"""
Database Models for CipherQ
Clean SQLAlchemy ORM models for FDA drugs and chemical compounds
"""

from datetime import datetime
from typing import Optional, List
from sqlalchemy import (
    Column, Integer, String, Float, Boolean, Text, DateTime,
    ForeignKey, Index, JSON, Enum as SQLEnum
)
from sqlalchemy.orm import declarative_base, relationship
import enum

Base = declarative_base()


class TherapeuticCategory(enum.Enum):
    """Therapeutic drug categories"""
    CARDIOVASCULAR = "cardiovascular"
    DIABETES = "diabetes"
    ANTI_INFLAMMATORY = "anti_inflammatory"
    NEUROLOGICAL = "neurological"
    PSYCHIATRIC = "psychiatric"
    ANTIBIOTIC = "antibiotic"
    ANTIVIRAL = "antiviral"
    CANCER = "cancer"
    PAIN = "pain"
    GASTROINTESTINAL = "gastrointestinal"
    RESPIRATORY = "respiratory"
    IMMUNOLOGY = "immunology"
    ENDOCRINE = "endocrine"
    OTHER = "other"


class ApprovalStatus(enum.Enum):
    """Drug approval status"""
    FDA_APPROVED = "fda_approved"
    EMA_APPROVED = "ema_approved"
    INVESTIGATIONAL = "investigational"
    WITHDRAWN = "withdrawn"
    EXPERIMENTAL = "experimental"


class FDADrug(Base):
    """
    FDA-Approved Drugs Table
    Contains 40,000+ approved drugs with clinical data
    """
    __tablename__ = 'fda_drugs'
    
    # Primary identifiers
    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(String(255), nullable=False, index=True)
    generic_name = Column(String(255), index=True)
    brand_names = Column(JSON)  # List of brand names
    
    # Chemical identifiers
    chembl_id = Column(String(50), unique=True, index=True)
    drugbank_id = Column(String(50), index=True)
    pubchem_cid = Column(Integer, index=True)
    cas_number = Column(String(50))
    unii = Column(String(50))
    
    # Structure
    smiles = Column(Text)
    inchi = Column(Text)
    inchi_key = Column(String(50), index=True)
    molecular_formula = Column(String(100))
    
    # Properties
    molecular_weight = Column(Float)
    logp = Column(Float)
    polar_surface_area = Column(Float)
    h_bond_donors = Column(Integer)
    h_bond_acceptors = Column(Integer)
    rotatable_bonds = Column(Integer)
    
    # CNS/BBB properties
    cns_mpo_score = Column(Float)
    bbb_permeability = Column(Float)
    
    # Classification
    therapeutic_category = Column(String(100), index=True)
    drug_class = Column(String(255))
    mechanism_of_action = Column(Text)
    pharmacology = Column(Text)
    
    # Approval status
    approval_status = Column(String(50), default='fda_approved')
    approval_year = Column(Integer)
    first_approval_country = Column(String(100))
    
    # Clinical data counts
    clinical_trial_count = Column(Integer, default=0)
    publication_count = Column(Integer, default=0)
    
    # Safety
    max_phase = Column(Integer, default=4)
    black_box_warning = Column(Boolean, default=False)
    
    # Indications
    approved_indications = Column(JSON)  # List of approved uses
    
    # Timestamps
    created_at = Column(DateTime, default=datetime.utcnow)
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow)
    
    # Relationships
    targets = relationship("DrugTarget", back_populates="drug", cascade="all, delete-orphan")
    
    # Indexes for fast searching
    __table_args__ = (
        Index('idx_fda_drug_category', 'therapeutic_category'),
        Index('idx_fda_drug_search', 'name', 'generic_name'),
        Index('idx_fda_drug_properties', 'molecular_weight', 'logp'),
    )
    
    def to_dict(self):
        """Convert to dictionary for API responses"""
        return {
            'id': self.id,
            'name': self.name,
            'generic_name': self.generic_name,
            'brand_names': self.brand_names or [],
            'chembl_id': self.chembl_id,
            'drugbank_id': self.drugbank_id,
            'smiles': self.smiles,
            'molecular_weight': self.molecular_weight,
            'logp': self.logp,
            'therapeutic_category': self.therapeutic_category,
            'drug_class': self.drug_class,
            'mechanism_of_action': self.mechanism_of_action,
            'approval_status': self.approval_status,
            'clinical_trial_count': self.clinical_trial_count,
            'publication_count': self.publication_count,
            'cns_mpo_score': self.cns_mpo_score,
            'bbb_permeability': self.bbb_permeability
        }


class ChemicalCompound(Base):
    """
    Chemical Compounds Table
    Contains experimental compounds, natural products, and novel molecules
    """
    __tablename__ = 'chemical_compounds'
    
    # Primary identifiers
    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(String(255), index=True)
    
    # Chemical identifiers
    pubchem_cid = Column(Integer, unique=True, index=True)
    chembl_id = Column(String(50), index=True)
    
    # Structure
    smiles = Column(Text, nullable=False)
    inchi = Column(Text)
    inchi_key = Column(String(50), index=True)
    molecular_formula = Column(String(100))
    
    # Properties
    molecular_weight = Column(Float, index=True)
    logp = Column(Float, index=True)
    polar_surface_area = Column(Float)
    h_bond_donors = Column(Integer)
    h_bond_acceptors = Column(Integer)
    rotatable_bonds = Column(Integer)
    aromatic_rings = Column(Integer)
    heavy_atom_count = Column(Integer)
    
    # Drug-likeness
    lipinski_violations = Column(Integer)
    druglikeness_score = Column(Float)
    synthetic_accessibility = Column(Float)
    
    # CNS/BBB properties
    cns_mpo_score = Column(Float)
    bbb_permeability = Column(Float)
    
    # Source
    source = Column(String(100))  # PubChem, ChEMBL, Natural Product, Synthetic
    source_id = Column(String(100))
    
    # Classification
    compound_class = Column(String(255))
    is_natural_product = Column(Boolean, default=False)
    
    # Activity data
    known_targets = Column(JSON)  # List of known protein targets
    bioactivity_data = Column(JSON)  # IC50, Ki, etc.
    
    # Timestamps
    created_at = Column(DateTime, default=datetime.utcnow)
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow)
    
    # Indexes
    __table_args__ = (
        Index('idx_compound_properties', 'molecular_weight', 'logp', 'druglikeness_score'),
        Index('idx_compound_source', 'source'),
    )
    
    def to_dict(self):
        """Convert to dictionary for API responses"""
        return {
            'id': self.id,
            'name': self.name,
            'pubchem_cid': self.pubchem_cid,
            'chembl_id': self.chembl_id,
            'smiles': self.smiles,
            'molecular_weight': self.molecular_weight,
            'logp': self.logp,
            'polar_surface_area': self.polar_surface_area,
            'druglikeness_score': self.druglikeness_score,
            'synthetic_accessibility': self.synthetic_accessibility,
            'lipinski_violations': self.lipinski_violations,
            'compound_class': self.compound_class,
            'is_natural_product': self.is_natural_product,
            'source': self.source
        }


class DrugTarget(Base):
    """
    Drug-Target relationships
    Links drugs to their protein targets
    """
    __tablename__ = 'drug_targets'
    
    id = Column(Integer, primary_key=True, autoincrement=True)
    drug_id = Column(Integer, ForeignKey('fda_drugs.id'), nullable=False, index=True)
    
    # Target info
    target_name = Column(String(255), nullable=False, index=True)
    target_gene = Column(String(50), index=True)
    uniprot_id = Column(String(20), index=True)
    pdb_id = Column(String(10))
    
    # Binding data
    binding_affinity = Column(Float)  # kcal/mol
    ic50 = Column(Float)  # nM
    ki = Column(Float)  # nM
    
    # Mechanism
    action_type = Column(String(100))  # inhibitor, agonist, antagonist, etc.
    
    # Relationship
    drug = relationship("FDADrug", back_populates="targets")
    
    __table_args__ = (
        Index('idx_drug_target_link', 'drug_id', 'target_name'),
    )


class Pathway(Base):
    """
    Biological pathways
    KEGG, Reactome, WikiPathways data
    """
    __tablename__ = 'pathways'
    
    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(String(255), nullable=False, index=True)
    
    # External IDs
    kegg_id = Column(String(50), index=True)
    reactome_id = Column(String(50), index=True)
    wikipathways_id = Column(String(50))
    
    # Description
    description = Column(Text)
    category = Column(String(100))
    
    # Associated genes/proteins
    genes = Column(JSON)  # List of gene symbols
    
    # Disease associations
    diseases = Column(JSON)  # List of associated diseases
    
    created_at = Column(DateTime, default=datetime.utcnow)
