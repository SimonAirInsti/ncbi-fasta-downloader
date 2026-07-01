#!/usr/bin/env python3
"""
Configuration setup script for NCBI FASTA Downloader.
Prepares the environment and validates configuration.

Configurador para el Descargador NCBI FASTA.
Prepara el entorno y valida la configuración.
"""

import os
import sys
import json
import shutil
from pathlib import Path
from typing import Dict, Any, Optional

class ProjectSetup:
    """Project setup and configuration manager."""
    
    def __init__(self):
        """Initialize the setup manager."""
        self.project_root = Path(__file__).parent
        self.config_file = self.project_root / "config.json"
        self.env_file = self.project_root / ".env"
        
    def create_default_config(self) -> Dict[str, Any]:
        """Create default configuration dictionary."""
        return {
            "email": "",  # User must provide
            "batch_size": 200,
            "max_workers": 4,
            "rate_limit_delay": 0.34,
            "max_retries": 3,
            "retry_delay": 2,
            "save_interval": 300,  # seconds
            "log_level": "INFO",
            "output_dir": "output",
            "files": {
                "input_file": "ncbi_ids_dict.tsv",
                "output_file": "all_sequences.fasta",
                "progress_file": "progress.json",
                "failed_ids_file": "failed_ids.txt",
                "log_file": "ncbi_download.log"
            },
            "ncbi": {
                "database": "protein",
                "rettype": "fasta",
                "retmode": "text",
                "api_base": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
            },
            "uniprot": {
                "max_workers": 35,
                "proteome_csv": "proteomes_Influenza_A_2025_02_05.tsv",
                "fasta_dict": {
                    "Neuraminidase": "uniprotkb_Neuraminidase_2025_02_05.fasta",
                    "Hemagglutinin": "uniprotkb_Hemagglutinin_2025_02_05.fasta",
                    "Matrix protein 2": "uniprotkb_Matrix_protein_2_AND_taxonomy_2025_06_18.fasta",
                    "Matrix protein 1": "uniprotkb_Matrix_protein_1_AND_taxonomy_2025_06_18.fasta",
                    "Nucleoprotein": "uniprotkb_Nucleoprotein_AND_taxonomy_id_2025_06_18.fasta"
                }
            }
        }
    
    def prompt_for_email(self) -> str:
        """Prompt user for email address."""
        print("\n📧 NCBI Email Configuration")
        print("=" * 40)
        print("NCBI requires a valid email address for API access.")
        print("This email is used to identify your requests and contact you if needed.")
        print("Your email will only be stored locally in config.json")
        print()
        
        while True:
            email = input("Enter your email address: ").strip()
            if "@" in email and "." in email:
                return email
            else:
                print("❌ Please enter a valid email address")
    
    def create_config_file(self, force: bool = False) -> bool:
        """Create configuration file with user input."""
        if self.config_file.exists() and not force:
            print(f"✅ Configuration file already exists: {self.config_file}")
            response = input("Overwrite? [y/N]: ").strip().lower()
            if response not in ['y', 'yes']:
                return False
        
        config = self.create_default_config()
        
        # Get email from user
        config["email"] = self.prompt_for_email()
        
        # Create output directory
        output_dir = self.project_root / config["output_dir"]
        output_dir.mkdir(exist_ok=True)
        
        # Update file paths to use output directory
        for file_key in config["files"]:
            if file_key != "input_file":  # Keep input file in project root
                config["files"][file_key] = str(output_dir / config["files"][file_key])
        
        # Save configuration
        with open(self.config_file, 'w') as f:
            json.dump(config, f, indent=2)
        
        print(f"✅ Configuration saved to: {self.config_file}")
        return True
    
    def validate_config(self) -> bool:
        """Validate existing configuration."""
        if not self.config_file.exists():
            print(f"❌ Configuration file not found: {self.config_file}")
            return False
        
        try:
            with open(self.config_file, 'r') as f:
                config = json.load(f)
            
            # Validate required fields
            required_fields = ["email", "batch_size", "max_workers"]
            for field in required_fields:
                if field not in config or not config[field]:
                    print(f"❌ Missing required field: {field}")
                    return False
            
            # Validate email
            if "@" not in config["email"]:
                print("❌ Invalid email address in configuration")
                return False
            
            # Validate numeric fields
            numeric_fields = {
                "batch_size": (1, 1000),
                "max_workers": (1, 10),
                "rate_limit_delay": (0.1, 5.0)
            }
            
            for field, (min_val, max_val) in numeric_fields.items():
                if field in config:
                    value = config[field]
                    if not isinstance(value, (int, float)) or not (min_val <= value <= max_val):
                        print(f"❌ Invalid value for {field}: {value} (should be {min_val}-{max_val})")
                        return False
            
            print("✅ Configuration is valid")
            return True
            
        except json.JSONDecodeError:
            print("❌ Invalid JSON in configuration file")
            return False
        except Exception as e:
            print(f"❌ Error validating configuration: {e}")
            return False
    
    def check_dependencies(self) -> bool:
        """Check if required dependencies are installed."""
        print("\n🔍 Checking Dependencies")
        print("=" * 30)
        
        dependencies = [
            ("biopython", "Bio"),
            ("requests", "requests"),
            ("pandas", "pandas"),
            ("joblib", "joblib"),
            ("tqdm", "tqdm"),
            ("numpy", "numpy")
        ]
        
        missing = []
        for package, import_name in dependencies:
            try:
                __import__(import_name)
                print(f"✅ {package}")
            except ImportError:
                print(f"❌ {package}")
                missing.append(package)
        
        if missing:
            print(f"\n📦 Install missing dependencies:")
            print(f"pip install {' '.join(missing)}")
            return False
        
        return True
    
    def create_sample_input(self) -> None:
        """Create a sample input file for testing."""
        sample_file = self.project_root / "sample_ncbi_ids.tsv"
        
        if sample_file.exists():
            print(f"✅ Sample file already exists: {sample_file}")
            return
        
        # Small sample dataset for testing
        sample_data = {
            "Neuraminidase": ["3024910520", "3024910499", "3024910476"],
            "Hemagglutinin": ["3024910409", "3024910387", "3024910364"]
        }
        
        with open(sample_file, 'w') as f:
            for protein_type, ids in sample_data.items():
                f.write(f"{protein_type}\t{','.join(ids)}\n")
        
        print(f"✅ Sample input file created: {sample_file}")
        print("   Use this file for testing with: python run_download.py --input sample_ncbi_ids.tsv")
    
    def setup_logging_config(self) -> None:
        """Create logging configuration file."""
        logging_config = {
            "version": 1,
            "disable_existing_loggers": False,
            "formatters": {
                "standard": {
                    "format": "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
                },
                "detailed": {
                    "format": "%(asctime)s - %(name)s - %(levelname)s - %(module)s - %(funcName)s - %(message)s"
                }
            },
            "handlers": {
                "console": {
                    "class": "logging.StreamHandler",
                    "level": "INFO",
                    "formatter": "standard",
                    "stream": "ext://sys.stdout"
                },
                "file": {
                    "class": "logging.FileHandler",
                    "level": "DEBUG",
                    "formatter": "detailed",
                    "filename": "output/ncbi_download.log",
                    "mode": "a"
                }
            },
            "loggers": {
                "": {
                    "handlers": ["console", "file"],
                    "level": "DEBUG",
                    "propagate": False
                }
            }
        }
        
        config_file = self.project_root / "logging_config.json"
        with open(config_file, 'w') as f:
            json.dump(logging_config, f, indent=2)
        
        print(f"✅ Logging configuration created: {config_file}")
    
    def run_setup(self) -> bool:
        """Run complete setup process."""
        print("🚀 NCBI FASTA Downloader - Project Setup")
        print("=" * 50)
        print()
        
        # Check dependencies first
        if not self.check_dependencies():
            print("\n❌ Please install missing dependencies before continuing")
            return False
        
        # Create configuration
        print("\n⚙️  Creating Configuration")
        print("=" * 30)
        config_created = self.create_config_file()
        
        if not config_created:
            print("Configuration setup cancelled")
            return False
        
        # Validate configuration
        if not self.validate_config():
            print("❌ Configuration validation failed")
            return False
        
        # Create additional files
        self.create_sample_input()
        self.setup_logging_config()
        
        # Create output directory
        output_dir = self.project_root / "output"
        output_dir.mkdir(exist_ok=True)
        print(f"✅ Output directory ready: {output_dir}")
        
        print(f"\n🎉 Setup Complete!")
        print("=" * 20)
        print("Next steps:")
        print("1. Place your NCBI IDs file as 'ncbi_ids_dict.tsv'")
        print("2. Test with sample data: python run_download.py --input sample_ncbi_ids.tsv")
        print("3. Run full download: python run_download.py")
        print("4. Monitor progress: python check_progress.py")
        
        return True

def main():
    """Main setup function."""
    setup = ProjectSetup()
    
    # Parse command line arguments
    if len(sys.argv) > 1:
        if sys.argv[1] in ["--help", "-h"]:
            print("NCBI FASTA Downloader Setup")
            print("Usage: python setup.py [--force]")
            print("Options:")
            print("  --force    Overwrite existing configuration")
            print("  --validate Validate existing configuration")
            print("  --help     Show this help message")
            return
        elif sys.argv[1] == "--validate":
            setup.validate_config()
            return
        elif sys.argv[1] == "--force":
            setup.run_setup()
            return
    
    # Run interactive setup
    success = setup.run_setup()
    sys.exit(0 if success else 1)

if __name__ == "__main__":
    # Only run the interactive project setup when explicitly requested by the user.
    # Build tools create isolated environments and may execute this file; avoid
    # running interactive checks during those automated builds. To run the
    # interactive setup, call this script with the `--interactive` flag:
    #
    #   python setup.py --interactive
    #
    if "--interactive" in sys.argv:
        main()
