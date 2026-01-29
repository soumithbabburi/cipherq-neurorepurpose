"""
FDA Orange Book Parser
Parses products.txt, patent.txt, exclusivity.txt for drug patent information
"""

import pandas as pd
import logging
from datetime import datetime

logger = logging.getLogger(__name__)

def get_drug_patent_info(drug_name: str) -> dict:
    """
    Get patent and exclusivity information for a drug from Orange Book files
    
    Returns:
        dict with patent_number, expiration_date, exclusivity_info, approval_date
    """
    
    try:
        # Load products file
        df_products = pd.read_csv('products.txt', sep='~', encoding='latin-1')
        
        # Find drug (case insensitive)
        drug_products = df_products[df_products['Ingredient'].str.upper() == drug_name.upper()]
        
        if drug_products.empty:
            # Try partial match
            drug_products = df_products[df_products['Ingredient'].str.contains(drug_name, case=False, na=False)]
        
        if drug_products.empty:
            return {
                'status': 'Not found in Orange Book',
                'patent_number': None,
                'expiration_date': None,
                'approval_date': None
            }
        
        # Get application number
        appl_no = drug_products.iloc[0]['Appl_No']
        approval_date = drug_products.iloc[0]['Approval_Date']
        
        # Load patent file
        try:
            df_patents = pd.read_csv('patent.txt', sep='~', encoding='latin-1')
            
            # Find patents for this application
            drug_patents = df_patents[df_patents['Appl_No'] == appl_no]
            
            if not drug_patents.empty:
                # Get primary patent (earliest expiration usually)
                patent_row = drug_patents.iloc[0]
                patent_number = patent_row['Patent_No']
                patent_expire = patent_row['Patent_Expire_Date_Text']
                
                return {
                    'status': 'Protected',
                    'patent_number': patent_number,
                    'expiration_date': patent_expire,
                    'approval_date': approval_date,
                    'appl_no': appl_no
                }
        except:
            pass
        
        # Load exclusivity file
        try:
            df_exclusivity = pd.read_csv('exclusivity.txt', sep='~', encoding='latin-1')
            
            # Find exclusivity for this application
            drug_excl = df_exclusivity[df_exclusivity['Appl_No'] == appl_no]
            
            if not drug_excl.empty:
                excl_date = drug_excl.iloc[0]['Exclusivity_Date']
                return {
                    'status': 'Market Exclusivity',
                    'patent_number': None,
                    'expiration_date': excl_date,
                    'approval_date': approval_date,
                    'appl_no': appl_no
                }
        except:
            pass
        
        # Has approval but no patent/exclusivity found
        return {
            'status': 'Approved (no active patents)',
            'patent_number': None,
            'expiration_date': None,
            'approval_date': approval_date,
            'appl_no': appl_no
        }
        
    except Exception as e:
        logger.error(f"Orange Book lookup failed for {drug_name}: {e}")
        return {
            'status': 'Data unavailable',
            'patent_number': None,
            'expiration_date': None
        }


__all__ = ['get_drug_patent_info']
