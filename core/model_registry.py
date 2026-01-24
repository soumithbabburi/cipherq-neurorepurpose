"""
Model Registry - Track and display ML models for transparency
"""
import logging
from typing import Dict, Any, List
from datetime import datetime
from core.config_loader import config

logger = logging.getLogger(__name__)

class ModelRegistry:
    """Registry to track all ML models used in the platform"""
    
    def __init__(self):
        self.models = config.load_model_config().get("models", {})
        self.last_checked = datetime.now()
    
    def get_all_models(self) -> List[Dict[str, Any]]:
        """Get list of all models with metadata"""
        model_list = []
        
        for model_key, model_info in self.models.items():
            model_list.append({
                "key": model_key,
                "name": model_info.get("name", "Unknown"),
                "version": model_info.get("version", "Unknown"),
                "type": model_info.get("type", "Unknown"), 
                "description": model_info.get("description", ""),
                "endpoint": model_info.get("endpoint", "N/A"),
                "status": self._check_model_health(model_key)
            })
        
        return model_list
    
    def get_model_info(self, model_key: str) -> Dict[str, Any]:
        """Get detailed info for a specific model"""
        if model_key not in self.models:
            return {}
        
        model_info = self.models[model_key].copy()
        model_info["status"] = self._check_model_health(model_key)
        model_info["last_checked"] = self.last_checked.isoformat()
        return model_info
    
    def _check_model_health(self, model_key: str) -> str:
        """Check if model/service is available"""
        try:
            if model_key == "molecular_docking":
                # Check NVIDIA API key
                import os
                return "Connected" if os.getenv("NVIDIA_API_KEY") else "Demo Mode"
            elif model_key == "conversational_ai":
                # Check Anthropic API key
                import os
                return "Online" if os.getenv("ANTHROPIC_API_KEY") else "Fallback"
            elif model_key == "molecular_properties":
                # Check RDKit availability
                try:
                    from rdkit import Chem
                    return "Available"
                except ImportError:
                    return "Unavailable"
            elif model_key == "graph_analysis":
                # Check if evidence graph builder is available
                return "Active"
            elif model_key == "ml_scoring":
                # Check if ML scoring models are available
                try:
                    from services.ml_scoring_service import ml_scoring_service
                    return "Trained" if ml_scoring_service.models_initialized else "Failed"
                except ImportError:
                    return "Unavailable"
            else:
                return "Unknown"
        except Exception as e:
            logger.warning(f"Health check failed for {model_key}: {e}")
            return "Error"

# Global instance
model_registry = ModelRegistry()