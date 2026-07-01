#!/usr/bin/env python3
"""
Example: Use ncbi-fasta-downloader as a library in your own script.
Ejemplo: Usar ncbi-fasta-downloader como libreria en tu propio script.

Shows how to import the package and use its classes programmatically.
"""

# Option A: import the convenience package
from ncbi_downloader_pkg import (
    NCBIDownloaderV2,
    UniProtDownloader,
    DataIntegrator,
    check_progress_enhanced,
)

# Option B: import directly from the modules
# from ncbi_downloader import NCBIDownloaderV2, ConfigurationManager
# from data_integrator import DataIntegrator

print(f"ncbi-fasta-downloader version: {__import__('ncbi_downloader_pkg').__version__}")

# --- Quick integration example ---
integrator = DataIntegrator("output")

# List available files
csv_files = integrator.identify_csv_files(exclude_patterns=["proteome", "unified"])
fasta_files = integrator.identify_fasta_files(exclude_patterns=["unified"])

print(f"CSV  files found: {len(csv_files)}")
print(f"FASTA files found: {len(fasta_files)}")

# Check progress
check_progress_enhanced("output")
