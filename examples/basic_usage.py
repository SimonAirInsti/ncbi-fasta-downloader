#!/usr/bin/env python3
"""
Basic usage example for ncbi-fasta-downloader.
Ejemplo basico de uso de ncbi-fasta-downloader.

Demonstrates:
  - Downloading sequences from NCBI
  - Downloading data from UniProt
  - Checking download progress
"""

from ncbi_downloader import NCBIDownloaderV2, UniProtDownloader
from check_progress import check_progress_enhanced
from pathlib import Path

# --- 1. Configure and download from NCBI ---
# First, create a config.json file (see config_example.json)
# Then:
downloader = NCBIDownloaderV2("config.json")

# Load NCBI IDs from a TSV file
ids_dict = downloader.load_ids_dict("ncbi_ids_dict.tsv")
print(f"Loaded {sum(len(v) for v in ids_dict.values())} sequences across {len(ids_dict)} proteins")

# Download all sequences
downloader.download_all_sequences(ids_dict)


# --- 2. Generate NCBI IDs dynamically ---
downloader.generate_ncbi_ids_file(
    virus_name="Influenza A",
    target_proteins=["Hemagglutinin", "Neuraminidase"],
    max_results_per_protein=500,
)


# --- 3. Download from UniProt ---
downloader.download_uniprot_data(
    virus_name="Influenza A",
    target_proteins=["Hemagglutinin", "Neuraminidase", "Nucleoprotein"],
)


# --- 4. Check progress ---
check_progress_enhanced(output_dir="output", show_details=True)
