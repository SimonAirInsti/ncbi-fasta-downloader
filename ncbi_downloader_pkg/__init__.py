"""
NCBI FASTA Downloader - A high-performance bioinformatics library.

Download, integrate, and analyze protein sequences from NCBI and UniProt
with batch processing, intelligent rate limiting, and data integration.

Usage:
    from ncbi_downloader_pkg import NCBIDownloaderV2, UniProtDownloader, DataIntegrator

    # Download from NCBI
    downloader = NCBIDownloaderV2("config.json")
    ids = downloader.load_ids_dict("ncbi_ids_dict.tsv")
    downloader.download_all_sequences(ids)

    # Integrate data
    integrator = DataIntegrator("output")
    csv_files = integrator.identify_csv_files()
    fasta_files = integrator.identify_fasta_files()
    ...
"""

__version__ = "2.0.0"
__author__ = "Simon Marino"
__email__ = "marinoperea98@usal.es"
__license__ = "MIT"

# Import main classes from existing modules (sibling files at project root)
import sys
from pathlib import Path

# Ensure the project root is in sys.path so we can import the existing modules
_project_root = str(Path(__file__).resolve().parent.parent)
if _project_root not in sys.path:
    sys.path.insert(0, _project_root)

# Core downloader classes
from ncbi_downloader import (
    NCBIDownloaderV2,
    UniProtDownloader,
    ConfigurationManager,
    DownloadStats,
    main as cli_main,
    create_arg_parser,
)

# Data integration
from data_integrator import DataIntegrator

# Progress checker
from check_progress import check_progress_enhanced

__all__ = [
    # Core classes
    "NCBIDownloaderV2",
    "UniProtDownloader",
    "ConfigurationManager",
    "DownloadStats",
    "DataIntegrator",
    # Utilities
    "check_progress_enhanced",
    # CLI
    "cli_main",
    "create_arg_parser",
]
