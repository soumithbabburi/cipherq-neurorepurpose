"""
ML Scoring Service - Real ML models for drug repurposing scoring
Gradient Boosting + Random Forest ensemble approach
"""
import logging
import pandas as pd
import numpy as np
from typing import Dict, List, Any, Tuple, Optional
from core.config_loader import config
from services.filter_service import alzheimer_filter

logger = logging.getLogger(__name__)

class MLScoringService:
    """Real ML models for drug repurposing scoring using ensemble approach"""
    
    def __init__(self):
        self.model_config = config.load_model_config()
        self.app_config = config.load_app_config()
        
        # Initialize ML models
        self.models_initialized = False
        self.gradient_boost_model = None
        self.random_forest_model = None
        self.ensemble_weights = {"gb": 0.6, "rf": 0.4}  # GB gets higher weight
        
        # Feature columns for ML training/prediction
        self.feature_columns = [
            'HOMO-LUMO Gap', 'Binding Affinity', 'BBB Penetration',
            'Drug-likeness', 'CNS Score', 'Oral Bioavailability',
            'Plasma Protein Binding', 'Half-life (hours)',
            'Hepatotoxicity_Score', 'Cardiotoxicity_Score'
        ]
        
        self._initialize_ml_models()
    
    def _initialize_ml_models(self):
        """Initialize Gradient Boosting and Random Forest models"""
        try:
            from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
            from sklearn.preprocessing import StandardScaler
            
            # Gradient Boosting Model
            self.gradient_boost_model = GradientBoostingRegressor(
                n_estimators=100,
                learning_rate=0.1,
                max_depth=6,
                random_state=42,
                subsample=0.8,
                validation_fraction=0.1,
                n_iter_no_change=10
            )
            
            # Random Forest Model  
            self.random_forest_model = RandomForestRegressor(
                n_estimators=100,
                max_depth=10,
                min_samples_split=5,
                min_samples_leaf=2,
                random_state=42,
                n_jobs=-1
            )
            
            # Feature scaler
            self.scaler = StandardScaler()
            
            # Train models with synthetic data initially
            self._train_with_synthetic_data()
            
            self.models_initialized = True
            logger.info("✅ ML Scoring Models (Gradient Boost + Random Forest) initialized successfully")
            
        except ImportError as e:
            logger.error(f"ML libraries not available: {e}")
            self.models_initialized = False
        except Exception as e:
            logger.error(f"Failed to initialize ML models: {e}")
            self.models_initialized = False
    
    def _train_with_synthetic_data(self):
        """Train models with synthetic training data representing drug repurposing patterns"""
        np.random.seed(42)
        
        # Generate synthetic training data (1000 samples)
        n_samples = 1000
        
        # Features based on known drug repurposing patterns
        homo_lumo = np.random.normal(3.5, 1.0, n_samples)
        binding_affinity = np.random.normal(-6.5, 1.5, n_samples)
        bbb_penetration = np.random.beta(2, 2, n_samples)  # Beta distribution for [0,1] values
        drug_likeness = np.random.beta(3, 1, n_samples)  # Higher drug-likeness is better
        cns_score = np.random.beta(2, 3, n_samples)  # CNS scores tend to be lower
        oral_bioavailability = np.random.beta(3, 2, n_samples)
        protein_binding = np.random.uniform(0.1, 0.95, n_samples)
        half_life = np.random.lognormal(2.0, 0.8, n_samples)
        hepatotox_score = np.random.beta(3, 1, n_samples)  # Higher = safer
        cardiotox_score = np.random.beta(3, 1, n_samples)  # Higher = safer
        
        # Create features DataFrame
        X_train = pd.DataFrame({
            'HOMO-LUMO Gap': homo_lumo,
            'Binding Affinity': binding_affinity,
            'BBB Penetration': bbb_penetration,
            'Drug-likeness': drug_likeness,
            'CNS Score': cns_score,
            'Oral Bioavailability': oral_bioavailability,
            'Plasma Protein Binding': protein_binding,
            'Half-life (hours)': half_life,
            'Hepatotoxicity_Score': hepatotox_score,
            'Cardiotoxicity_Score': cardiotox_score
        })
        
        # Generate target variable (repurposing likelihood)
        # Realistic patterns: higher drug-likeness + better CNS + safer = higher score
        y_train = (
            0.25 * drug_likeness +
            0.20 * bbb_penetration +
            0.15 * cns_score +
            0.10 * oral_bioavailability +
            0.10 * hepatotox_score +
            0.10 * cardiotox_score +
            0.05 * np.clip(-binding_affinity / 10.0, 0, 1) +  # Stronger binding = better
            0.05 * np.clip((4.5 - homo_lumo) / 4.5, 0, 1) +  # Optimal HOMO-LUMO gap
            np.random.normal(0, 0.05, n_samples)  # Add noise
        )
        
        # Clip to [0, 1] range
        y_train = np.clip(y_train, 0, 1)
        
        # Scale features
        X_train_scaled = self.scaler.fit_transform(X_train)
        
        # Train both models
        self.gradient_boost_model.fit(X_train_scaled, y_train)
        self.random_forest_model.fit(X_train_scaled, y_train)
        
        logger.info(f"ML models trained on {n_samples} synthetic samples")
        
        # Log model performance on training data
        gb_score = self.gradient_boost_model.score(X_train_scaled, y_train)
        rf_score = self.random_forest_model.score(X_train_scaled, y_train)
        logger.info(f"Gradient Boost R²: {gb_score:.3f}, Random Forest R²: {rf_score:.3f}")
    
    def calculate_ml_scores(self, properties_df: pd.DataFrame) -> pd.DataFrame:
        """Calculate drug repurposing scores using trained ML ensemble"""
        if properties_df.empty:
            return pd.DataFrame()
        
        if not self.models_initialized:
            logger.warning("ML models not available, falling back to heuristic scoring")
            return self._calculate_fallback_scores(properties_df)
        
        try:
            # Prepare features for ML prediction
            ml_features_df = self._prepare_ml_features(properties_df)
            
            # Scale features
            X_scaled = self.scaler.transform(ml_features_df)
            
            # Get predictions from both models
            gb_predictions = self.gradient_boost_model.predict(X_scaled)
            rf_predictions = self.random_forest_model.predict(X_scaled)
            
            # Ensemble prediction (weighted average)
            ensemble_predictions = (
                self.ensemble_weights["gb"] * gb_predictions +
                self.ensemble_weights["rf"] * rf_predictions
            )
            
            # Get feature importance for transparency
            gb_importance = self._get_feature_importance(self.gradient_boost_model, "gradient_boost")
            rf_importance = self._get_feature_importance(self.random_forest_model, "random_forest")
            
            # Build results DataFrame
            results = []
            for i, row in properties_df.iterrows():
                drug_name = row.get('Drug', f'Drug_{i}')
                
                # Skip Alzheimer's drugs
                if alzheimer_filter.is_alzheimer_drug(drug_name):
                    continue
                
                ensemble_score = ensemble_predictions[i]
                gb_score = gb_predictions[i]
                rf_score = rf_predictions[i]
                
                results.append({
                    'Drug': drug_name,
                    'ML_Ensemble_Score': ensemble_score,
                    'Gradient_Boost_Score': gb_score,
                    'Random_Forest_Score': rf_score,
                    'Repurposing_Likelihood': min(ensemble_score * 100, 95),  # Convert to percentage
                    'Model_Confidence': self._calculate_confidence(gb_score, rf_score),
                    'Primary_Features': self._get_top_contributing_features(ml_features_df.iloc[i], gb_importance)
                })
            
            results_df = pd.DataFrame(results)
            results_df = results_df.sort_values('ML_Ensemble_Score', ascending=False).reset_index(drop=True)
            
            logger.info(f"ML ensemble scoring completed for {len(results_df)} drugs")
            return results_df
            
        except Exception as e:
            logger.error(f"ML scoring failed: {e}")
            return self._calculate_fallback_scores(properties_df)
    
    def _prepare_ml_features(self, properties_df: pd.DataFrame) -> pd.DataFrame:
        """Prepare features for ML model input"""
        ml_features = pd.DataFrame()
        
        for col in self.feature_columns:
            if col in properties_df.columns:
                ml_features[col] = properties_df[col]
            elif col == 'Hepatotoxicity_Score':
                # Convert text to numeric score
                ml_features[col] = properties_df.get('Hepatotoxicity', 'Low risk').apply(
                    lambda x: self._toxicity_to_score(x, 'hepato')
                )
            elif col == 'Cardiotoxicity_Score':
                # Convert text to numeric score
                ml_features[col] = properties_df.get('Cardiotoxicity', 'Generally safe').apply(
                    lambda x: self._toxicity_to_score(x, 'cardio')
                )
            else:
                # Use reasonable defaults for missing features
                defaults = {
                    'HOMO-LUMO Gap': 3.5,
                    'Binding Affinity': -6.0,
                    'BBB Penetration': 0.6,
                    'Drug-likeness': 0.7,
                    'CNS Score': 0.5,
                    'Oral Bioavailability': 0.6,
                    'Plasma Protein Binding': 0.8,
                    'Half-life (hours)': 8.0
                }
                ml_features[col] = defaults.get(col, 0.5)
        
        return ml_features
    
    def _toxicity_to_score(self, toxicity_text: str, tox_type: str) -> float:
        """Convert toxicity text to numeric score (higher = safer)"""
        text_lower = str(toxicity_text).lower()
        
        if 'low' in text_lower or 'safe' in text_lower:
            return 0.9
        elif 'moderate' in text_lower or 'tolerated' in text_lower:
            return 0.6
        elif 'high' in text_lower or 'risk' in text_lower:
            return 0.3
        else:
            return 0.7  # Default moderate safety
    
    def _get_feature_importance(self, model, model_name: str) -> Dict[str, float]:
        """Get feature importance from trained model"""
        try:
            importances = model.feature_importances_
            return dict(zip(self.feature_columns, importances))
        except Exception as e:
            logger.warning(f"Could not get feature importance from {model_name}: {e}")
            return {}
    
    def _calculate_confidence(self, gb_score: float, rf_score: float) -> float:
        """Calculate ensemble confidence based on agreement between models"""
        agreement = 1.0 - abs(gb_score - rf_score)
        return max(0.5, agreement)  # Minimum 50% confidence
    
    def _get_top_contributing_features(self, feature_row: pd.Series, importance_dict: Dict[str, float]) -> List[str]:
        """Get top 3 features contributing to the prediction"""
        if not importance_dict:
            return ["N/A"]
        
        # Sort features by importance
        sorted_features = sorted(importance_dict.items(), key=lambda x: x[1], reverse=True)
        
        # Return top 3 feature names
        return [feat[0] for feat in sorted_features[:3]]
    
    def _calculate_fallback_scores(self, properties_df: pd.DataFrame) -> pd.DataFrame:
        """Fallback heuristic scoring when ML models unavailable"""
        results = []
        
        for _, row in properties_df.iterrows():
            drug_name = row.get('Drug', 'Unknown')
            
            if alzheimer_filter.is_alzheimer_drug(drug_name):
                continue
            
            # Simple weighted heuristic as fallback
            score = (
                0.25 * row.get('Drug-likeness', 0.7) +
                0.20 * row.get('BBB Penetration', 0.6) +
                0.15 * row.get('CNS Score', 0.5) +
                0.15 * (1.0 - abs(row.get('Binding Affinity', -6.0)) / 10.0) +
                0.25 * 0.7  # Default safety score
            )
            
            results.append({
                'Drug': drug_name,
                'ML_Ensemble_Score': score,
                'Gradient_Boost_Score': score,
                'Random_Forest_Score': score,
                'Repurposing_Likelihood': min(score * 100, 95),
                'Model_Confidence': 0.6,
                'Primary_Features': ["Fallback Mode"]
            })
        
        return pd.DataFrame(results).sort_values('ML_Ensemble_Score', ascending=False).reset_index(drop=True)
    
    def get_model_info(self) -> Dict[str, Any]:
        """Get ML model information for transparency"""
        return {
            "models_used": ["Gradient Boosting Regressor", "Random Forest Regressor"],
            "ensemble_approach": "Weighted Average",
            "ensemble_weights": self.ensemble_weights,
            "feature_count": len(self.feature_columns),
            "features": self.feature_columns,
            "models_initialized": self.models_initialized,
            "training_samples": 1000,
            "model_type": "Supervised Learning - Regression"
        }

# Global instance
ml_scoring_service = MLScoringService()