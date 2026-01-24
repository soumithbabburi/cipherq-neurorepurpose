"""
Configuration Management for CipherQ
Handles environment variables and database connections
"""
import os
import logging
from typing import Dict, Optional

logger = logging.getLogger(__name__)

class Config:
    """Configuration manager for CipherQ application"""
    
    # Database Configuration
    DB_HOST = os.getenv("DB_HOST", "localhost")
    DB_PORT = os.getenv("DB_PORT", "5432")
    DB_NAME = os.getenv("DB_NAME", "cipherq_repurpose")
    DB_USER = os.getenv("DB_USER", "babburisoumith")
    DB_PASSWORD = os.getenv("DB_PASSWORD", "")
    
    # NVIDIA API Configuration
    NVIDIA_API_KEY = os.getenv("NVIDIA_API_KEY", "")
    
    # Application Settings
    DEBUG = os.getenv("DEBUG", "False").lower() == "true"
    LOG_LEVEL = os.getenv("LOG_LEVEL", "INFO")
    
    @classmethod
    def get_db_params(cls) -> Dict[str, str]:
        """Get database connection parameters"""
        return {
            "host": cls.DB_HOST,
            "port": cls.DB_PORT,
            "database": cls.DB_NAME,
            "user": cls.DB_USER,
            "password": cls.DB_PASSWORD
        }
    
    @classmethod
    def validate_config(cls) -> Dict[str, bool]:
        """
        Validate configuration and return status.
        
        Returns:
            Dictionary with validation results
        """
        validation = {
            "database_configured": all([
                cls.DB_HOST,
                cls.DB_NAME,
                cls.DB_USER
            ]),
            "nvidia_api_configured": bool(cls.NVIDIA_API_KEY),
            "all_configured": False
        }
        
        validation["all_configured"] = all([
            validation["database_configured"],
            validation["nvidia_api_configured"]
        ])
        
        return validation
    
    @classmethod
    def print_config_status(cls):
        """Print configuration status to console"""
        validation = cls.validate_config()
        
        print("\n" + "="*50)
        print("CipherQ Configuration Status")
        print("="*50)
        
        # Database
        if validation["database_configured"]:
            print("✅ Database: Configured")
            print(f"   Host: {cls.DB_HOST}")
            print(f"   Database: {cls.DB_NAME}")
            print(f"   User: {cls.DB_USER}")
        else:
            print("❌ Database: Not Configured")
            print("   Set: DB_HOST, DB_NAME, DB_USER, DB_PASSWORD")
        
        # NVIDIA API
        if validation["nvidia_api_configured"]:
            print("✅ NVIDIA API: Configured")
            print(f"   Key: {cls.NVIDIA_API_KEY[:20]}...")
        else:
            print("❌ NVIDIA API: Not Configured")
            print("   Set: NVIDIA_API_KEY")
        
        # Overall status
        print("-"*50)
        if validation["all_configured"]:
            print("✅ All systems configured")
        else:
            print("⚠️  Some systems need configuration")
            print("\nTo configure, set environment variables:")
            print("  export DB_HOST=localhost")
            print("  export DB_NAME=cipherq_repurpose")
            print("  export DB_USER=your_username")
            print("  export DB_PASSWORD=your_password")
            print("  export NVIDIA_API_KEY=your_nvidia_key")
        
        print("="*50 + "\n")
        
        return validation


def load_env_file(env_file: str = ".env"):
    """
    Load environment variables from .env file.
    
    Args:
        env_file: Path to .env file (default: .env)
    """
    if not os.path.exists(env_file):
        logger.warning(f"Environment file {env_file} not found")
        return False
    
    try:
        with open(env_file, 'r') as f:
            for line in f:
                line = line.strip()
                
                # Skip comments and empty lines
                if not line or line.startswith('#'):
                    continue
                
                # Parse KEY=VALUE
                if '=' in line:
                    key, value = line.split('=', 1)
                    key = key.strip()
                    value = value.strip()
                    
                    # Remove quotes if present
                    if value.startswith('"') and value.endswith('"'):
                        value = value[1:-1]
                    elif value.startswith("'") and value.endswith("'"):
                        value = value[1:-1]
                    
                    # Set environment variable if not already set
                    if key not in os.environ:
                        os.environ[key] = value
        
        logger.info(f"Loaded environment variables from {env_file}")
        return True
        
    except Exception as e:
        logger.error(f"Failed to load {env_file}: {e}")
        return False


def create_sample_env_file(filename: str = ".env.sample"):
    """
    Create a sample .env file with all required variables.
    
    Args:
        filename: Name of sample file to create
    """
    sample_content = """# CipherQ Configuration File
# Copy this to .env and fill in your values

# Database Configuration
DB_HOST=localhost
DB_PORT=5432
DB_NAME=cipherq_repurpose
DB_USER=babburisoumith
DB_PASSWORD=

# NVIDIA BioNeMo API
# Get your API key from: https://build.nvidia.com/
NVIDIA_API_KEY=nvapi-your-key-here

# Application Settings
DEBUG=False
LOG_LEVEL=INFO
"""
    
    try:
        with open(filename, 'w') as f:
            f.write(sample_content)
        
        print(f"✅ Created sample environment file: {filename}")
        print(f"   Copy to .env and update with your credentials")
        return True
        
    except Exception as e:
        logger.error(f"Failed to create {filename}: {e}")
        return False


# Auto-load .env if it exists
if os.path.exists(".env"):
    load_env_file()

# Validate and print status on import (optional, can be disabled)
if __name__ == "__main__":
    Config.print_config_status()
    
    # Create sample if requested
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "create-sample":
        create_sample_env_file()
