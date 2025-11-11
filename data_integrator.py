#!/usr/bin/env python3
"""
Data Integrator for NCBI and UniProt outputs.
Combines FASTA files from NCBI downloads and CSV files from UniProt searches into unified datasets.

This module provides functionality to:
1. Read and identify output files from the output folder
2. Process and concatenate CSV files from UniProt searches
3. Parse FASTA files from NCBI downloads
4. Integrate both data sources into unified CSV and FASTA files
"""

import pandas as pd
import numpy as np
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import logging
import json
from datetime import datetime
from typing import List, Dict, Tuple, Optional
import re

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('data_integration.log'),
        logging.StreamHandler()
    ]
)

class DataIntegrator:
    """
    Integrates NCBI FASTA and UniProt CSV data into unified datasets.
    """
    
    def __init__(self, output_dir: str = "output"):
        """
        Initialize the DataIntegrator.
        
        Args:
            output_dir (str): Path to the output directory containing data files
        """
        self.output_dir = Path(output_dir)
        self.logger = logging.getLogger(__name__)
        
        if not self.output_dir.exists():
            raise FileNotFoundError(f"Output directory not found: {self.output_dir}")
        
        self.logger.info(f"DataIntegrator initialized with output directory: {self.output_dir}")

    def identify_csv_files(self, exclude_patterns: List[str] = None, target_proteins: List[str] = None) -> List[Path]:
        """
        Identify CSV files in the output directory, excluding specific patterns.
        
        Args:
            exclude_patterns (List[str]): Patterns to exclude (default: ["proteome"])
            target_proteins (List[str]): If provided, only include CSV files matching these protein names
            
        Returns:
            List[Path]: List of CSV file paths
        """
        if exclude_patterns is None:
            exclude_patterns = ["proteome"]
        
        csv_files = []
        
        # Normalize target proteins for matching (if provided)
        normalized_target_proteins = None
        if target_proteins:
            # Convert to lowercase and replace spaces with underscores for filename matching
            normalized_target_proteins = [
                protein.lower().replace(' ', '_') 
                for protein in target_proteins
            ]
            self.logger.info(f"Filtering CSV files for target proteins: {target_proteins}")
        
        # Find all CSV files
        for csv_file in self.output_dir.glob("*.csv"):
            # Check if file should be excluded by pattern
            should_exclude = False
            for pattern in exclude_patterns:
                if pattern.lower() in csv_file.name.lower():
                    should_exclude = True
                    self.logger.info(f"Excluding CSV file: {csv_file.name} (matches exclude pattern: {pattern})")
                    break
            
            if should_exclude:
                continue
            
            # If target_proteins is specified, only include matching files
            if normalized_target_proteins:
                filename_lower = csv_file.name.lower()
                
                # Check if any target protein name is in the filename
                matches_target = False
                for target_protein in normalized_target_proteins:
                    if target_protein in filename_lower:
                        matches_target = True
                        self.logger.info(f"Found CSV file: {csv_file.name} (matches target: {target_protein})")
                        break
                
                if matches_target:
                    csv_files.append(csv_file)
                else:
                    self.logger.info(f"Skipping CSV file: {csv_file.name} (doesn't match target proteins)")
            else:
                # No target filter - include all non-excluded files
                csv_files.append(csv_file)
                self.logger.info(f"Found CSV file: {csv_file.name}")
        
        self.logger.info(f"Total CSV files identified: {len(csv_files)}")
        return csv_files
    
    def identify_fasta_files(self, exclude_patterns: List[str] = None) -> List[Path]:
        """
        Identify FASTA files in the output directory.
        
        Args:
            exclude_patterns (List[str]): Patterns to exclude from filename matching
        
        Returns:
            List[Path]: List of FASTA file paths
        """
        if exclude_patterns is None:
            exclude_patterns = ["unified"]
            
        fasta_files = []
        
        # Common FASTA extensions
        fasta_extensions = ["*.fasta", "*.fa", "*.fas", "*.fna", "*.ffn", "*.faa", "*.frn"]
        
        for extension in fasta_extensions:
            for fasta_file in self.output_dir.glob(extension):
                # Check if file should be excluded
                should_exclude = any(pattern.lower() in fasta_file.name.lower() for pattern in exclude_patterns)
                if should_exclude:
                    self.logger.info(f"Excluding FASTA file: {fasta_file.name} (matches pattern: {', '.join([p for p in exclude_patterns if p.lower() in fasta_file.name.lower()])})")
                    continue
                    
                fasta_files.append(fasta_file)
                self.logger.info(f"Found FASTA file: {fasta_file.name}")
        
        # Also check for compressed FASTA files
        for extension in fasta_extensions:
            for fasta_file in self.output_dir.glob(f"{extension}.gz"):
                # Check if file should be excluded
                should_exclude = any(pattern.lower() in fasta_file.name.lower() for pattern in exclude_patterns)
                if should_exclude:
                    self.logger.info(f"Excluding compressed FASTA file: {fasta_file.name} (matches pattern: {', '.join([p for p in exclude_patterns if p.lower() in fasta_file.name.lower()])})")
                    continue
                    
                fasta_files.append(fasta_file)
                self.logger.info(f"Found compressed FASTA file: {fasta_file.name}")
        
        self.logger.info(f"Total FASTA files identified: {len(fasta_files)}")
        return fasta_files
    
    def read_and_concatenate_csv_files(self, csv_files: List[Path]) -> pd.DataFrame:
        """
        Read and concatenate multiple CSV files from UniProt searches.
        
        Args:
            csv_files (List[Path]): List of CSV file paths to concatenate
            
        Returns:
            pd.DataFrame: Concatenated DataFrame with all CSV data
        """
        if not csv_files:
            self.logger.warning("No CSV files provided for concatenation")
            return pd.DataFrame()
        
        dataframes = []
        
        for csv_file in csv_files:
            try:
                self.logger.info(f"Reading CSV file: {csv_file.name}")
                
                # Try different encodings
                encodings = ['utf-8', 'latin-1', 'cp1252']
                df = None
                
                for encoding in encodings:
                    try:
                        df = pd.read_csv(csv_file, encoding=encoding)
                        self.logger.debug(f"Successfully read {csv_file.name} with {encoding} encoding")
                        break
                    except UnicodeDecodeError:
                        continue
                
                if df is None:
                    self.logger.error(f"Could not read {csv_file.name} with any encoding")
                    continue
                
                # Add source file information
                df['Source_File'] = csv_file.name
                df['Data_Source'] = 'UniProt'
                
                # Extract protein type from filename if not present
                if 'Protein type' not in df.columns and 'Protein_Type' not in df.columns:
                    # Try to extract protein type from filename
                    protein_type = self._extract_protein_type_from_filename(csv_file.name)
                    if protein_type:
                        df['Protein_Type'] = protein_type
                
                self.logger.info(f"Read {len(df)} rows from {csv_file.name}")
                dataframes.append(df)
                
            except Exception as e:
                self.logger.error(f"Error reading {csv_file.name}: {e}")
                continue
        
        if not dataframes:
            self.logger.warning("No CSV files could be read successfully")
            return pd.DataFrame()
        
        # Concatenate all DataFrames
        self.logger.info("Concatenating all CSV DataFrames...")
        combined_df = pd.concat(dataframes, ignore_index=True, sort=False)
        
        # Add integration metadata
        combined_df['Integration_Timestamp'] = datetime.now().isoformat()
        combined_df['Total_Source_Files'] = len(csv_files)
        
        self.logger.info(f"Successfully concatenated {len(dataframes)} CSV files into {len(combined_df)} total rows")
        
        return combined_df
    
    def _extract_protein_type_from_filename(self, filename: str) -> Optional[str]:
        """
        Extract protein type from filename.
        
        Args:
            filename (str): Name of the file
            
        Returns:
            Optional[str]: Extracted protein type or None
        """
        # Remove file extension
        name_without_ext = filename.replace('.csv', '').replace('_data', '')
        
        # Common protein names
        protein_patterns = {
            'neuraminidase': 'Neuraminidase',
            'hemagglutinin': 'Hemagglutinin',
            'matrix_protein_1': 'Matrix protein 1',
            'matrix_protein_2': 'Matrix protein 2',
            'nucleoprotein': 'Nucleoprotein'
        }
        
        name_lower = name_without_ext.lower()
        for pattern, protein_type in protein_patterns.items():
            if pattern in name_lower:
                return protein_type
        
        # If no pattern matches, try to clean up the filename
        cleaned_name = name_without_ext.replace('_', ' ').title()
        return cleaned_name if cleaned_name else None
    
    def analyze_csv_structure(self, csv_files: List[Path]) -> Dict:
        """
        Analyze the structure of CSV files to understand their schema.
        
        Args:
            csv_files (List[Path]): List of CSV files to analyze
            
        Returns:
            Dict: Analysis results including column information
        """
        analysis = {
            'files_analyzed': len(csv_files),
            'file_details': [],
            'common_columns': set(),
            'all_columns': set(),
            'column_frequency': {},
            'data_types': {}
        }
        
        for i, csv_file in enumerate(csv_files):
            try:
                df = pd.read_csv(csv_file)
                
                file_info = {
                    'filename': csv_file.name,
                    'rows': len(df),
                    'columns': list(df.columns),
                    'column_count': len(df.columns),
                    'dtypes': df.dtypes.to_dict()
                }
                
                analysis['file_details'].append(file_info)
                
                # Track column frequency
                for col in df.columns:
                    analysis['all_columns'].add(col)
                    if col in analysis['column_frequency']:
                        analysis['column_frequency'][col] += 1
                    else:
                        analysis['column_frequency'][col] = 1
                
                # For the first file, initialize common columns
                if i == 0:
                    analysis['common_columns'] = set(df.columns)
                else:
                    # Find intersection with previous common columns
                    analysis['common_columns'] = analysis['common_columns'].intersection(set(df.columns))
                
                self.logger.info(f"Analyzed {csv_file.name}: {len(df)} rows, {len(df.columns)} columns")
                
            except Exception as e:
                self.logger.error(f"Error analyzing {csv_file.name}: {e}")
        
        # Convert sets to lists for JSON serialization
        analysis['common_columns'] = list(analysis['common_columns'])
        analysis['all_columns'] = list(analysis['all_columns'])
        
        self.logger.info(f"CSV Analysis complete:")
        self.logger.info(f"  - Common columns across all files: {len(analysis['common_columns'])}")
        self.logger.info(f"  - Total unique columns: {len(analysis['all_columns'])}")
        
        return analysis
    
    def save_analysis_report(self, analysis: Dict, output_file: str = "csv_analysis_report.json"):
        """
        Save the CSV analysis report to a JSON file.
        
        Args:
            analysis (Dict): Analysis results from analyze_csv_structure
            output_file (str): Output filename for the report
        """
        output_path = self.output_dir / output_file
        
        try:
            with open(output_path, 'w', encoding='utf-8') as f:
                json.dump(analysis, f, indent=2, default=str)
            
            self.logger.info(f"Analysis report saved to: {output_path}")
            
        except Exception as e:
            self.logger.error(f"Error saving analysis report: {e}")
    
    def parse_fasta_files(self, fasta_files: List[Path]) -> pd.DataFrame:
        """
        Parse FASTA files and convert to DataFrame format with enhanced metadata extraction.
        
        Args:
            fasta_files (List[Path]): List of FASTA file paths
            
        Returns:
            pd.DataFrame: DataFrame with FASTA sequence data and extracted metadata
        """
        if not fasta_files:
            self.logger.warning("No FASTA files provided for parsing")
            return pd.DataFrame()
        
        all_sequences = []
        
        for fasta_file in fasta_files:
            try:
                self.logger.info(f"Parsing FASTA file: {fasta_file.name}")
                
                # Handle compressed files
                if fasta_file.suffix == '.gz':
                    import gzip
                    opener = gzip.open
                    mode = 'rt'
                else:
                    opener = open
                    mode = 'r'
                
                sequences_count = 0
                
                with opener(fasta_file, mode) as handle:
                    for record in SeqIO.parse(handle, "fasta"):
                        # Basic FASTA data
                        sequence_data = {
                            'ID': record.id,
                            'Sequence': str(record.seq),
                            'Description': record.description,
                            'Sequence_Length': len(record.seq),
                            'Source_File': fasta_file.name,
                            'Data_Source': 'NCBI'
                        }
                        
                        # Extract enhanced metadata using the same methods as ncbi_downloader
                        enhanced_metadata = self._extract_enhanced_metadata_from_fasta(record.description, str(record.seq))
                        sequence_data.update(enhanced_metadata)
                        
                        all_sequences.append(sequence_data)
                        sequences_count += 1
                
                self.logger.info(f"Parsed {sequences_count} sequences from {fasta_file.name}")
                
            except Exception as e:
                self.logger.error(f"Error parsing FASTA file {fasta_file.name}: {e}")
                continue
        
        if not all_sequences:
            self.logger.warning("No sequences could be parsed from FASTA files")
            return pd.DataFrame()
        
        # Convert to DataFrame
        fasta_df = pd.DataFrame(all_sequences)
        
        # Add integration metadata
        fasta_df['Integration_Timestamp'] = datetime.now().isoformat()
        fasta_df['Total_Source_Files'] = len(fasta_files)
        
        self.logger.info(f"Successfully parsed {len(fasta_df)} sequences from {len(fasta_files)} FASTA files")
        self.logger.info(f"Enhanced FASTA DataFrame columns: {list(fasta_df.columns)}")
        
        return fasta_df
    
    def _extract_enhanced_metadata_from_fasta(self, description: str, sequence: str) -> Dict:
        """
        Extract enhanced metadata from FASTA description using the same methods as ncbi_downloader.
        This function extracts: Prot ID, Protein type, Organism, and Variants.
        
        Args:
            description (str): FASTA description line
            sequence (str): Sequence string
            
        Returns:
            Dict: Extracted metadata with UniProt-compatible column names
        """
        metadata = {
            'Prot ID': None,
            'Protein type': None,
            'Organism': None,
            'Variants': None
        }
        
        try:
            # Extract Protein ID - similar to ncbi_downloader method
            if '|' in description:
                parts = description.split('|')
                if len(parts) >= 2:
                    # Try to get the protein ID from the standard format
                    metadata['Prot ID'] = parts[1] if len(parts) > 1 else parts[0].replace('>', '')
                else:
                    metadata['Prot ID'] = parts[0].replace('>', '')
            else:
                # If no pipe, take the first part after removing '>'
                first_part = description.split()[0].replace('>', '')
                metadata['Prot ID'] = first_part
            
            # Extract Protein type from metadata in brackets (from ncbi_downloader enhanced headers)
            protein_type_match = re.search(r'\[protein_type=([^\]]+)\]', description)
            if protein_type_match:
                metadata['Protein type'] = protein_type_match.group(1)
            else:
                # Try to infer protein type from description text
                protein_keywords = {
                    'neuraminidase': 'Neuraminidase',
                    'hemagglutinin': 'Hemagglutinin', 
                    'matrix protein 1': 'Matrix protein 1',
                    'matrix protein 2': 'Matrix protein 2',
                    'nucleoprotein': 'Nucleoprotein'
                }
                
                desc_lower = description.lower()
                for keyword, protein_type in protein_keywords.items():
                    if keyword in desc_lower:
                        metadata['Protein type'] = protein_type
                        break
            
            # Extract Organism information - using the same regex pattern as ncbi_downloader
            organism_pattern = re.compile(r'OS=(.*?)(?:\s+OX=|$)')
            organism_match = organism_pattern.search(description)
            if organism_match:
                metadata['Organism'] = organism_match.group(1).strip()
            else:
                # Try alternative patterns for organism extraction
                # Look for patterns like "Influenza A virus (...)"
                virus_pattern = re.compile(r'(Influenza A virus[^,\[\]]*)')
                virus_match = virus_pattern.search(description)
                if virus_match:
                    metadata['Organism'] = virus_match.group(1).strip()
            
            # Extract Variants - using the same approach as ncbi_downloader
            if metadata['Organism']:
                # Extract variants using the same pattern as in ncbi_downloader
                variant_pattern = re.compile(r'Influenza A virus \((.*?)\)')
                variant_match = variant_pattern.search(metadata['Organism'])
                if variant_match:
                    metadata['Variants'] = variant_match.group(1)
                else:
                    # Try to extract variant from other patterns
                    # Look for strain information
                    strain_patterns = [
                        r'strain ([^,\)]+)',
                        r'\(([^)]+strain[^)]*)\)',
                        r'\(([AH][0-9]+[^)]*)\)'  # Pattern for H1N1, H3N2, etc.
                    ]
                    
                    for pattern in strain_patterns:
                        strain_match = re.search(pattern, description, re.IGNORECASE)
                        if strain_match:
                            metadata['Variants'] = strain_match.group(1).strip()
                            break
            
            self.logger.debug(f"Extracted metadata from '{description[:100]}...': {metadata}")
            
        except Exception as e:
            self.logger.warning(f"Error extracting metadata from FASTA description: {e}")
            self.logger.debug(f"Problematic description: {description}")
        
        return metadata
    
    def combine_csv_and_fasta_data(self, csv_df: pd.DataFrame, fasta_df: pd.DataFrame) -> pd.DataFrame:
        """
        Combine CSV and FASTA data into a unified DataFrame.
        
        Args:
            csv_df (pd.DataFrame): DataFrame from CSV files (UniProt)
            fasta_df (pd.DataFrame): DataFrame from FASTA files (NCBI)
            
        Returns:
            pd.DataFrame: Combined DataFrame with unified schema
        """
        self.logger.info("Combining CSV and FASTA data...")
        
        # Standardize column names for merging
        combined_data = []
        
        # Process CSV data (UniProt)
        if not csv_df.empty:
            csv_standardized = csv_df.copy()
            # Ensure consistent column naming for CSV data
            csv_column_mapping = {
                'Prot ID': 'Prot_ID',
                'Protein type': 'Protein_Type'
            }
            csv_standardized = csv_standardized.rename(columns=csv_column_mapping)
            combined_data.append(csv_standardized)
            self.logger.info(f"Added {len(csv_standardized)} rows from CSV data")
        
        # Process FASTA data (NCBI)
        if not fasta_df.empty:
            fasta_standardized = fasta_df.copy()
            # Ensure consistent column naming for FASTA data
            fasta_column_mapping = {
                'Prot ID': 'Prot_ID', 
                'Protein type': 'Protein_Type'
            }
            fasta_standardized = fasta_standardized.rename(columns=fasta_column_mapping)
            combined_data.append(fasta_standardized)
            self.logger.info(f"Added {len(fasta_standardized)} rows from FASTA data")
        
        if not combined_data:
            self.logger.warning("No data to combine")
            return pd.DataFrame()
        
        # Combine all data
        combined_df = pd.concat(combined_data, ignore_index=True, sort=False)
        
        # Calculate sequence length for all entries and ensure it's integer
        self.logger.info("Calculating sequence lengths for all entries...")
        combined_df['Sequence_Length'] = combined_df['Sequence'].str.len()
        # Convert to integer (handles any NaN values by filling with 0)
        combined_df['Sequence_Length'] = combined_df['Sequence_Length'].fillna(0).astype(int)
        self.logger.info(f"Sequence lengths calculated. Range: {combined_df['Sequence_Length'].min()}-{combined_df['Sequence_Length'].max()}")
        
        # Log pre-filtering stats
        pre_filter_count = len(combined_df)
        self.logger.info(f"Pre-filtering total: {pre_filter_count:,} sequences")
        
        # Remove sequences with "fragment" in description
        self.logger.info("Removing sequences containing 'fragment' in description...")
        fragment_mask = combined_df['Description'].str.contains('fragment', case=False, na=False)
        fragments_count = fragment_mask.sum()
        combined_df = combined_df[~fragment_mask]
        
        post_filter_count = len(combined_df)
        self.logger.info(f"Removed {fragments_count:,} fragment sequences")
        self.logger.info(f"Post-fragment-filtering total: {post_filter_count:,} sequences")
        
        # Drop duplicates based on Sequence, keeping first (UniProt entries have priority)
        self.logger.info("Removing duplicate sequences (keeping UniProt entries as priority)...")
        combined_df = combined_df.drop_duplicates(subset="Sequence", keep="first")
        
        post_dedup_count = len(combined_df)
        duplicates_removed = post_filter_count - post_dedup_count
        total_removed = pre_filter_count - post_dedup_count
        self.logger.info(f"Removed {duplicates_removed:,} duplicate sequences")
        self.logger.info(f"Post-deduplication total: {post_dedup_count:,} sequences")
        
        # Add final integration metadata (for internal tracking)
        combined_df['Final_Integration_Timestamp'] = datetime.now().isoformat()
        combined_df['Total_Sequences'] = len(combined_df)
        combined_df['Duplicates_Removed'] = duplicates_removed
        combined_df['Fragments_Removed'] = fragments_count
        combined_df['Total_Removed'] = total_removed
        combined_df['Pre_Filter_Count'] = pre_filter_count
        
        # Calculate some statistics
        csv_count = len(csv_df) if not csv_df.empty else 0
        fasta_count = len(fasta_df) if not fasta_df.empty else 0
        
        self.logger.info(f"Successfully combined, filtered, and deduplicated data:")
        self.logger.info(f"  - CSV sequences: {csv_count:,}")
        self.logger.info(f"  - FASTA sequences: {fasta_count:,}")
        self.logger.info(f"  - Total before filtering: {pre_filter_count:,}")
        self.logger.info(f"  - Fragment sequences removed: {fragments_count:,}")
        self.logger.info(f"  - Duplicates removed: {duplicates_removed:,}")
        self.logger.info(f"  - Final unique sequences: {post_dedup_count:,}")
        
        # Remove metadata columns before returning (keep them only for logging)
        metadata_columns = [
            'Final_Integration_Timestamp', 'Total_Sequences', 'Duplicates_Removed', 
            'Fragments_Removed', 'Total_Removed', 'Pre_Filter_Count',
            'Integration_Timestamp', 'Total_Source_Files'
        ]
        
        # Create clean version without metadata columns
        clean_combined_df = combined_df.drop(columns=[col for col in metadata_columns if col in combined_df.columns])
        
        self.logger.info(f"Removed metadata columns from final dataset")
        self.logger.info(f"Final clean dataset columns: {list(clean_combined_df.columns)}")
        
        return clean_combined_df
    
    def create_unified_fasta(self, combined_df: pd.DataFrame, output_file: str = "unified_sequences.fasta") -> bool:
        """
        Create a unified FASTA file from the combined DataFrame.
        
        Args:
            combined_df (pd.DataFrame): Combined data from CSV and FASTA sources
            output_file (str): Output filename for the unified FASTA
            
        Returns:
            bool: True if successful, False otherwise
        """
        if combined_df.empty:
            self.logger.warning("No data provided for FASTA creation")
            return False
        
        output_path = self.output_dir / output_file
        
        try:
            sequences_written = 0
            
            with open(output_path, 'w') as f:
                for idx, row in combined_df.iterrows():
                    # Create FASTA header
                    header_parts = []
                    
                    # Add ID
                    if pd.notna(row.get('ID')):
                        header_parts.append(str(row['ID']))
                    elif pd.notna(row.get('Prot_ID')):
                        header_parts.append(str(row['Prot_ID']))
                    else:
                        header_parts.append(f"seq_{idx}")
                    
                    # Add protein type
                    if pd.notna(row.get('Protein_Type')):
                        header_parts.append(f"protein_type={row['Protein_Type']}")
                    
                    # Add organism
                    if pd.notna(row.get('Organism')):
                        header_parts.append(f"organism={row['Organism']}")
                    
                    # Add source
                    if pd.notna(row.get('Data_Source')):
                        header_parts.append(f"source={row['Data_Source']}")
                    
                    # Add variants if available
                    if pd.notna(row.get('Variants')):
                        header_parts.append(f"variant={row['Variants']}")
                    
                    # Create header line
                    header = f">{' | '.join(header_parts)}"
                    
                    # Get sequence
                    sequence = row.get('Sequence', '')
                    
                    if sequence and len(sequence.strip()) > 0:
                        f.write(f"{header}\n")
                        # Write sequence in 80-character lines (standard FASTA format)
                        for i in range(0, len(sequence), 80):
                            f.write(f"{sequence[i:i+80]}\n")
                        sequences_written += 1
            
            self.logger.info(f"Successfully created unified FASTA file: {output_path}")
            self.logger.info(f"Sequences written: {sequences_written:,}")
            
            return True
            
        except Exception as e:
            self.logger.error(f"Error creating unified FASTA file: {e}")
            return False
        """
        Test the CSV concatenation functionality.
        
        Args:
            exclude_patterns (List[str]): Patterns to exclude from CSV files
            
        Returns:
            pd.DataFrame: Concatenated DataFrame
        """
        self.logger.info("Starting CSV concatenation test...")
        
        # Step 1: Identify CSV files
        csv_files = self.identify_csv_files(exclude_patterns)
        
        if not csv_files:
            self.logger.warning("No CSV files found for concatenation")
            return pd.DataFrame()
        
        # Step 2: Analyze CSV structure
        analysis = self.analyze_csv_structure(csv_files)
        self.save_analysis_report(analysis)
        
        # Step 3: Concatenate CSV files
        combined_df = self.read_and_concatenate_csv_files(csv_files)
        
        # Step 4: Save concatenated result
        if not combined_df.empty:
            output_file = self.output_dir / "concatenated_uniprot_data.csv"
            combined_df.to_csv(output_file, index=False)
            self.logger.info(f"Concatenated CSV data saved to: {output_file}")
            
            # Print summary
            self.logger.info(f"Concatenation Summary:")
            self.logger.info(f"  - Input files: {len(csv_files)}")
            self.logger.info(f"  - Total rows: {len(combined_df)}")
            self.logger.info(f"  - Total columns: {len(combined_df.columns)}")
            
            # Show column names
            self.logger.info(f"  - Columns: {list(combined_df.columns)}")
            
            # Show data types
            self.logger.info(f"  - Data types:")
            for col, dtype in combined_df.dtypes.items():
                self.logger.info(f"    {col}: {dtype}")
        
        return combined_df


def main():
    """
    Main function to test the complete data integration functionality.
    """
    logger = logging.getLogger(__name__)
    
    try:
        # Initialize the data integrator
        integrator = DataIntegrator("output")
        
        logger.info("="*60)
        logger.info("TESTING COMPLETE DATA INTEGRATION FUNCTIONALITY")
        logger.info("="*60)
        
        # Step 1: Process CSV files
        logger.info("\n[CSV] Step 1: Processing CSV files...")
        csv_files = integrator.identify_csv_files(exclude_patterns=["proteome", "unified"])
        logger.info(f"Found {len(csv_files)} CSV files: {[f.name for f in csv_files]}")
        
        csv_df = integrator.read_and_concatenate_csv_files(csv_files)
        
        if not csv_df.empty:
            logger.info(f"[SUCCESS] CSV processing completed! Shape: {csv_df.shape}")
            logger.info(f"CSV columns: {list(csv_df.columns)}")
        else:
            logger.warning("[WARNING] No CSV data was processed.")
        
        # Step 2: Process FASTA files
        logger.info("\n[FASTA] Step 2: Processing FASTA files...")
        fasta_files = integrator.identify_fasta_files(exclude_patterns=["unified"])
        logger.info(f"Found {len(fasta_files)} FASTA files: {[f.name for f in fasta_files]}")
        
        fasta_df = integrator.parse_fasta_files(fasta_files)
        
        if not fasta_df.empty:
            logger.info(f"[SUCCESS] FASTA processing completed! Shape: {fasta_df.shape}")
            logger.info(f"FASTA columns: {list(fasta_df.columns)}")
            
            # Show sample of extracted metadata
            logger.info("\n[SAMPLE] Sample of extracted FASTA metadata:")
            sample_cols = ['Prot ID', 'Protein type', 'Organism', 'Variants']
            available_cols = [col for col in sample_cols if col in fasta_df.columns]
            if available_cols:
                print(fasta_df[available_cols].head(3).to_string())
            
        else:
            logger.warning("[WARNING] No FASTA data was processed.")
        
        # Step 3: Combine data if both sources have data
        if not csv_df.empty and not fasta_df.empty:
            logger.info("\n[COMBINE] Step 3: Combining CSV and FASTA data...")
            
            combined_df = integrator.combine_csv_and_fasta_data(csv_df, fasta_df)
            
            if not combined_df.empty:
                logger.info(f"[SUCCESS] Data combination completed! Shape: {combined_df.shape}")
                
                # Save combined CSV
                output_csv = integrator.output_dir / "unified_data.csv"
                combined_df.to_csv(output_csv, index=False)
                logger.info(f"[SAVED] Unified CSV saved to: {output_csv}")
                
                # Create unified FASTA
                fasta_success = integrator.create_unified_fasta(combined_df, "unified_sequences.fasta")
                if fasta_success:
                    logger.info(f"[SAVED] Unified FASTA saved to: {integrator.output_dir / 'unified_sequences.fasta'}")
                
                # Show summary statistics
                logger.info(f"\n[SUMMARY] Final Integration Summary:")
                logger.info(f"  - Total sequences: {len(combined_df):,}")
                logger.info(f"  - UniProt source: {len(csv_df):,} sequences")
                logger.info(f"  - NCBI source: {len(fasta_df):,} sequences")
                if 'Protein type' in combined_df.columns:
                    protein_counts = combined_df['Protein type'].value_counts()
                    logger.info(f"  - Protein types: {dict(protein_counts)}")
                
            else:
                logger.warning("[WARNING] Data combination resulted in empty dataset.")
                
        elif not csv_df.empty:
            logger.info("\n[CSV-ONLY] Only CSV data available - saving as unified output...")
            output_csv = integrator.output_dir / "unified_data_csv_only.csv"
            csv_df.to_csv(output_csv, index=False)
            logger.info(f"[SAVED] CSV-only dataset saved to: {output_csv}")
            
        elif not fasta_df.empty:
            logger.info("\n[FASTA-ONLY] Only FASTA data available - saving as unified output...")
            output_csv = integrator.output_dir / "unified_data_fasta_only.csv"
            fasta_df.to_csv(output_csv, index=False)
            logger.info(f"[SAVED] FASTA-only dataset saved to: {output_csv}")
            
            # Also create FASTA output
            fasta_success = integrator.create_unified_fasta(fasta_df, "unified_sequences_fasta_only.fasta")
            if fasta_success:
                logger.info(f"[SAVED] FASTA-only sequences saved to: {integrator.output_dir / 'unified_sequences_fasta_only.fasta'}")
        
        else:
            logger.warning("[WARNING] No data from either CSV or FASTA sources. Check if files exist in the output folder.")
            
        logger.info("\n" + "="*60)
        logger.info("[COMPLETE] DATA INTEGRATION TEST COMPLETED")
        logger.info("="*60)
            
    except FileNotFoundError as e:
        logger.error(f"[ERROR] Error: {e}")
        logger.info("Make sure the 'output' directory exists and contains files from NCBI/UniProt downloads.")
        
    except Exception as e:
        logger.error(f"[ERROR] Unexpected error: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()
