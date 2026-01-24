"""
Configuration Loader - Remove hardcoding from clean_nvidia_app.py
"""
import yaml
import json
import os
from pathlib import Path
from typing import Dict, Any, Optional
import logging

logger = logging.getLogger(__name__)

class ConfigLoader:
    """Load configuration from YAML files and environment variables"""
    
    def __init__(self, config_dir: str = "config", knowledge_dir: str = "knowledge"):
        self.config_dir = Path(config_dir)
        self.knowledge_dir = Path(knowledge_dir)
        self._config_cache = {}
        
    def load_app_config(self) -> Dict[str, Any]:
        """Load main application configuration"""
        return self._load_yaml("app.yaml")
    
    def load_model_config(self) -> Dict[str, Any]:
        """Load ML model registry configuration"""
        return self._load_yaml("models.yaml")
    
    def load_evidence_config(self) -> Dict[str, Any]:
        """Load evidence graph configuration if it exists"""
        try:
            return self._load_yaml("evidence.yaml")
        except FileNotFoundError:
            # Return default evidence configuration
            return {
                "weights": {
                    "drug_target": 0.8,
                    "target_pathway": 0.7, 
                    "pathway_disease": 0.9,
                    "clinical_evidence": 0.6,
                    "publication_evidence": 0.5
                },
                "confidence_thresholds": {
                    "high": 0.8,
                    "medium": 0.6,
                    "low": 0.4
                }
            }
    
    def load_alzheimer_drugs(self) -> Dict[str, Any]:
        """Load Alzheimer's drug exclusion list"""
        file_path = self.knowledge_dir / "alzheimer_drugs.json"
        try:
            with open(file_path, 'r') as f:
                return json.load(f)
        except FileNotFoundError:
            logger.warning(f"Alzheimer drugs file not found: {file_path}")
            return {"existing_ad_drugs": [], "mesh_terms_to_exclude": []}
    
    def get_api_endpoint(self, service: str) -> Optional[str]:
        """Get API endpoint from config or environment"""
        model_config = self.load_model_config()
        
        # Check environment first
        env_key = f"{service.upper()}_API_URL"
        endpoint = os.getenv(env_key)
        if endpoint:
            return endpoint
            
        # Fallback to config
        api_config = model_config.get("api", {})
        return api_config.get(service, {}).get("base_url")
    
    def get_threshold(self, category: str, key: str) -> float:
        """Get threshold value from configuration"""
        app_config = self.load_app_config()
        thresholds = app_config.get("thresholds", {})
        return thresholds.get(category, {}).get(key, 0.5)  # Default fallback
    
    def _load_yaml(self, filename: str) -> Dict[str, Any]:
        """Load and cache YAML configuration file"""
        if filename in self._config_cache:
            return self._config_cache[filename]
            
        file_path = self.config_dir / filename
        try:
            with open(file_path, 'r') as f:
                config = yaml.safe_load(f)
                self._config_cache[filename] = config
                return config
        except FileNotFoundError:
            logger.error(f"Configuration file not found: {file_path}")
            return {}
        except yaml.YAMLError as e:
            logger.error(f"Error parsing YAML file {file_path}: {e}")
            return {}

# Global instance
config = ConfigLoader()