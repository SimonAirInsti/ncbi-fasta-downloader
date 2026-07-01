#!/usr/bin/env python3
"""
Complete Influenza A pipeline: download, integrate, and analyze.
Pipeline completo de Influenza A: descarga, integración y análisis.

This example shows how to:
  1. Download sequences from NCBI and UniProt
  2. Integrate data from both sources
  3. Produce unified CSV and FASTA outputs
"""

from ncbi_downloader import NCBIDownloaderV2
from data_integrator import DataIntegrator
from pathlib import Path

# ---------- Configuration ----------
CONFIG_FILE = "config.json"
OUTPUT_DIR = "output"
VIRUS_NAME = "Influenza A"
TARGET_PROTEINS = [
    "Hemagglutinin",
    "Neuraminidase",
    "Matrix protein 1",
    "Matrix protein 2",
    "Nucleoprotein",
]

# ---------- Step 1: Download sequences ----------
print("=" * 60)
print("STEP 1: Downloading sequences")
print("=" * 60)

downloader = NCBIDownloaderV2(CONFIG_FILE)

# Generate NCBI IDs
downloader.generate_ncbi_ids_file(
    virus_name=VIRUS_NAME,
    target_proteins=TARGET_PROTEINS,
    max_results_per_protein=1000,
)

# Load and download
ids_file = Path(OUTPUT_DIR) / "ncbi_ids_dict.tsv"
if ids_file.exists():
    ids_dict = downloader.load_ids_dict(str(ids_file))
    downloader.download_all_sequences(ids_dict)

# UniProt download
downloader.download_uniprot_data(
    virus_name=VIRUS_NAME,
    target_proteins=TARGET_PROTEINS,
)

# ---------- Step 2: Integrate data ----------
print("\n" + "=" * 60)
print("STEP 2: Integrating data")
print("=" * 60)

integrator = DataIntegrator(OUTPUT_DIR)

# Process CSV files (UniProt)
csv_files = integrator.identify_csv_files(
    exclude_patterns=["proteome", "unified"],
    target_proteins=TARGET_PROTEINS,
)
csv_df = integrator.read_and_concatenate_csv_files(csv_files)
print(f"UniProt data: {len(csv_df):,} sequences")

# Process FASTA files (NCBI)
fasta_files = integrator.identify_fasta_files(exclude_patterns=["unified"])
fasta_df = integrator.parse_fasta_files(fasta_files, target_proteins=TARGET_PROTEINS)
print(f"NCBI data: {len(fasta_df):,} sequences")

# Combine
if not csv_df.empty or not fasta_df.empty:
    combined_df = integrator.combine_csv_and_fasta_data(csv_df, fasta_df)

    # Save results
    output_csv = Path(OUTPUT_DIR) / "unified_data.csv"
    combined_df.to_csv(output_csv, index=False)
    print(f"\nUnified CSV saved: {output_csv}  ({len(combined_df):,} sequences)")

    integrator.create_unified_fasta(combined_df, "unified_sequences.fasta")
    print(f"Unified FASTA saved: {Path(OUTPUT_DIR) / 'unified_sequences.fasta'}")

    # Summary
    if "Protein type" in combined_df.columns:
        print("\nSequences per protein:")
        for ptype, count in combined_df["Protein type"].value_counts().items():
            print(f"  {ptype}: {count:,}")

print("\nPipeline completed!")
