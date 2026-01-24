"""
Setup and environment utilities
"""
from pathlib import Path
import streamlit as st
import os
from config.imports import logger


def setup_local_environment():
    """Setup local environment directories and configurations"""
    local_dirs = ['cache', 'outputs', 'logs', 'data']
    
    for dir_name in local_dirs:
        Path(dir_name).mkdir(exist_ok=True)
    
    if 'workflow_step' not in st.session_state:
        st.session_state.workflow_step = 1
    if 'selected_drugs' not in st.session_state:
        st.session_state.selected_drugs = []
    
    legacy_keys = ['network_data', 'selected_node', 'filtered_drugs', 'target_disease']
    for key in legacy_keys:
        if key in st.session_state:
            del st.session_state[key]


def get_env_var(key: str, default: str = "") -> str:
    """Get environment variable with fallback"""
    return os.getenv(key, default)


def check_api_keys():
    """Check and validate API keys with fallbacks"""
    api_status = {}
    
    nvidia_key = get_env_var('NVIDIA_API_KEY', '')
    api_status['nvidia'] = bool(nvidia_key)
    
    db_url = get_env_var('DATABASE_URL', '')
    api_status['database'] = bool(db_url)
    
    return api_status


def apply_app_styling():
    """Apply professional app styling"""
    from config.imports import STYLING_AVAILABLE, apply_main_theme
    if STYLING_AVAILABLE:
        apply_main_theme()


__all__ = ['setup_local_environment', 'get_env_var', 'check_api_keys', 'apply_app_styling']
