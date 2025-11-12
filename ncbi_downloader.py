#!/usr/bin/env python3
"""
Optimized NCBI FASTA downloader with enhanced configuration and performance.
Enhanced version with additional optimizations and GitHub-ready structure.

Descargador NCBI FASTA optimizado con configuración mejorada y rendimiento.
Versión mejorada con optimizaciones adicionales y estructura lista para GitHub.
"""

from Bio import Entrez, SeqIO
from urllib.error import HTTPError, URLError
import time
import os
import threading
import json
import concurrent.futures
import signal
import sys
import argparse
import logging.config
from typing import List, Dict, Tuple, Optional, Set
from datetime import datetime, timedelta
from pathlib import Path
import hashlib
import gzip
import io
import re
import requests
import urllib.parse
import pandas as pd
import numpy as np
from dataclasses import dataclass
from joblib import Parallel, delayed
from tqdm import tqdm
from dataclasses import dataclass
import requests
import io
import pandas as pd
from joblib import Parallel, delayed
from tqdm import tqdm
import re
import numpy as np
import logging
import urllib.parse

@dataclass
class DownloadStats:
    """Statistics tracking for downloads."""
    total_processed: int = 0
    successful: int = 0
    failed: int = 0
    start_time: Optional[float] = None
    last_save_time: float = 0
    bytes_downloaded: int = 0
    sequences_per_protein: Dict[str, int] = None
    
    def __post_init__(self):
        if self.sequences_per_protein is None:
            self.sequences_per_protein = {}

class UniProtDownloader:
    """UniProt data downloader and processor with dynamic API downloads."""
    
    def __init__(self, output_dir: Path, max_workers: int = 35):
        """Initialize UniProt downloader."""
        self.output_dir = output_dir
        self.max_workers = max_workers
        self.logger = logging.getLogger(f"{__name__}.UniProt")
        
    def download_proteome(self, virus_name: str) -> pd.DataFrame:
        """
        Download proteome data from UniProt API de forma optimizada.
        Solo descarga proteomas específicos y relevantes.
        
        Args:
            virus_name (str): Name of the virus
            
        Returns:
            pd.DataFrame: Proteome data as DataFrame
        """
        try:
            import urllib.parse
            
            # Encode virus name for URL
            encoded_virus = urllib.parse.quote_plus(virus_name)
            
            # URL optimizada - buscar proteomas de referencia y completos
            url = f"https://rest.uniprot.org/proteomes/stream?compressed=true&fields=upid%2Corganism%2Corganism_id%2Cprotein_count%2Cbusco%2Ccpd&format=tsv&query=organism_name%3A%22{encoded_virus}%22+AND+%28proteome_type%3Areference+OR+proteome_type%3Acomplete%29"
            
            self.logger.info(f"Downloading reference proteomes from: {url}")
            
            # Send request
            response = requests.get(url, timeout=60)
            response.raise_for_status()
            
            # Decompress content
            with gzip.GzipFile(fileobj=io.BytesIO(response.content)) as gz:
                content = gz.read().decode('utf-8')
            
            # Convert to DataFrame
            df = pd.read_csv(io.StringIO(content), sep='\t')
            
            if df.empty:
                self.logger.warning(f"No reference proteomes found for '{virus_name}'")
                self.logger.info("Trying broader search...")
                
                # Fallback - búsqueda más específica pero limitada
                url_fallback = f"https://rest.uniprot.org/proteomes/stream?compressed=true&fields=upid%2Corganism%2Corganism_id%2Cprotein_count%2Cbusco%2Ccpd&format=tsv&query=organism_name%3A%22{encoded_virus}%22&size=50"
                
                response = requests.get(url_fallback, timeout=60)
                response.raise_for_status()
                
                with gzip.GzipFile(fileobj=io.BytesIO(response.content)) as gz:
                    content = gz.read().decode('utf-8')
                
                df = pd.read_csv(io.StringIO(content), sep='\t')
            
            if df.empty:
                self.logger.warning(f"No proteomes found for '{virus_name}'")
                return pd.DataFrame()
            
            # Filtrar solo entradas relevantes del virus específico
            if 'Organism' in df.columns:
                virus_filter = df['Organism'].str.contains(virus_name, case=False, na=False)
                df = df[virus_filter]
            
            # Use 'Proteome Id' as index if exists
            if 'Proteome Id' in df.columns:
                df = df.set_index('Proteome Id')
            
            self.logger.info(f"Filtered proteomes downloaded: {len(df)} entries (optimized)")
            return df
            
        except Exception as e:
            self.logger.error(f"Error downloading proteome: {e}")
            return pd.DataFrame()
    
    def download_protein_fasta(self, protein_name: str) -> str:
        """
        Download FASTA file for a specific protein from UniProt.
        
        Args:
            protein_name (str): Name of the protein
            
        Returns:
            str: FASTA content as string
        """
        try:
            import urllib.parse
            
            # Encode protein name for URL
            encoded_protein = urllib.parse.quote_plus(protein_name)
            
            # UniProt API URL for FASTA
            url = f"https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=%28{encoded_protein}%29"
            
            self.logger.info(f"Downloading FASTA for {protein_name}")
            
            # Send request
            response = requests.get(url, timeout=120)
            response.raise_for_status()
            
            # Decompress content
            with gzip.GzipFile(fileobj=io.BytesIO(response.content)) as gz:
                content = gz.read().decode('utf-8')
            
            if not content.strip():
                self.logger.warning(f"Empty FASTA for '{protein_name}'")
                return ""
            
            self.logger.info(f"FASTA downloaded for {protein_name}")
            return content
            
        except Exception as e:
            self.logger.error(f"Error downloading FASTA for {protein_name}: {e}")
            return ""
        
    def fasta_string_to_dataframe(self, fasta_content: str) -> pd.DataFrame:
        """Convert FASTA string content to pandas DataFrame."""
        data = []
        try:
            # Use StringIO to simulate a file
            fasta_io = io.StringIO(fasta_content)
            
            for record in SeqIO.parse(fasta_io, "fasta"):
                data.append((record.id, str(record.seq), str(record.description)))
                
        except Exception as e:
            self.logger.error(f"Error parsing FASTA: {e}")
            
        return pd.DataFrame(data, columns=["ID", "Sequence", "Description"])
        
    def fasta_to_dataframe(self, file_path: str) -> pd.DataFrame:
        """Convert FASTA file to pandas DataFrame."""
        data = []
        try:
            with open(file_path, "r") as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    data.append((record.id, str(record.seq), str(record.description)))
            df = pd.DataFrame(data, columns=["ID", "Sequence", "Description"])
            return df
        except Exception as e:
            self.logger.error(f"Error reading FASTA file {file_path}: {e}")
            return pd.DataFrame(columns=["ID", "Sequence", "Description"])
    
    def extract_data_from_fasta_content(self, fasta_content: str, search_term: str, virus_name: str) -> pd.DataFrame:
        """
        Process FASTA content and filter for specific virus and protein.
        
        Args:
            fasta_content (str): FASTA content as string
            search_term (str): Protein name to search for
            virus_name (str): Virus name for filtering
            
        Returns:
            pd.DataFrame: Filtered and processed protein data
        """
        try:
            # Convert FASTA content to DataFrame
            data = self.fasta_string_to_dataframe(fasta_content)
            
            if data.empty:
                self.logger.warning(f"No data found for {search_term}")
                return pd.DataFrame()

            # Add "Prot ID" column
            data["Prot ID"] = data['ID'].str.split('|').str[1]

            # Add "Protein type" column
            data["Protein type"] = search_term
            
            # Extract organism using regex
            organism_pattern = re.compile(r'OS=(.*?) OX=')
            matched_text = []
            
            for description in data['Description']:
                match = organism_pattern.search(description)
                if match:
                    matched_text.append(match.group(1))
                else:
                    matched_text.append(None)
            
            data['Organism'] = matched_text

            # Filter for specific virus
            virus_filter = data['Organism'].str.contains(virus_name, na=False)
            
            # Filter for protein in description
            protein_filter = data['Description'].str.contains(search_term, case=False, na=False)
            
            # Apply both filters
            filtered_data = data[virus_filter & protein_filter].drop_duplicates(subset=['Sequence'], keep='first')

            # Extract variants
            if not filtered_data.empty:
                variant_pattern = re.compile(rf'{re.escape(virus_name)} \((.*?)\)')
                variants = []
                
                for organism in filtered_data['Organism']:
                    match = variant_pattern.search(organism) if organism else None
                    if match:
                        variants.append(match.group(1))
                    else:
                        variants.append(np.nan)
                
                # Create copy to avoid warnings
                filtered_data = filtered_data.copy()
                filtered_data["Variants"] = variants
            
            self.logger.info(f"Processed {len(filtered_data)} sequences for {search_term}")
            return filtered_data
            
        except Exception as e:
            self.logger.error(f"Error processing FASTA content for {search_term}: {e}")
            return pd.DataFrame()
    
    def extract_data(self, file_path: str, search_term: str) -> pd.DataFrame:
        """
        Filtra el archivo FASTA para encontrar secuencias que contengan un término específico
        en la descripción y pertenezcan a un organismo específico (Influenza A virus).
        """
        try:
            # Utilizamos la función fasta_to_dataframe previamente definida
            data = self.fasta_to_dataframe(file_path)
            
            if data.empty:
                self.logger.warning(f"No data found in {file_path}")
                return pd.DataFrame()

            # Añadimos columna "Prot ID"
            data["Prot ID"] = data['ID'].str.split('|').str[1]

            # Añadimos columna "Protein type"
            data["Protein type"] = search_term
            
            # Definimos el patrón de expresión regular para extraer el organismo
            organism_pattern = re.compile(r'OS=(.*?) OX=')  # Extrae el texto entre "OS=" y "OX="
            
            # Extraemos el nombre del organismo usando el patrón de expresión regular
            matched_text = [organism_pattern.search(description).group(1) if organism_pattern.search(description) else None for description in data['Description']]
            
            # Añadimos la columna 'Organism' al DataFrame
            data['Organism'] = matched_text

            # Filtramos las entradas para que solo contengan "Influenza A virus" en la columna 'Organism'
            influenzas = data['Organism'].str.contains('Influenza A virus', na=False)
            
            # Filtramos las entradas que contienen el término de búsqueda en la columna 'Description'
            search_results = data['Description'].str.contains(search_term, case=False, na=False)
            
            # Aplicamos ambos filtros
            filtered_data = data[influenzas & search_results].drop_duplicates(subset=['Sequence'], keep='first')

            # Creamos columna Variant
            # Regex para capturar variante
            pattern = re.compile(r'Influenza A virus \((.*?)\)')
            # Iniciamos lista
            variants = []
            
            # Extraemos variante de columna organismo
            for organism in filtered_data['Organism']:
                match = pattern.search(organism)
                if match:
                    variants.append(match.group(1))
                else:
                    variants.append(np.nan)
            
            # Crear una copia para evitar warnings
            filtered_data = filtered_data.copy()
            filtered_data["Variants"] = variants
            
            self.logger.info(f"Processed {len(filtered_data)} sequences for {search_term}")
            return filtered_data
            
        except Exception as e:
            self.logger.error(f"Error processing {file_path} for {search_term}: {e}")
            return pd.DataFrame()
    
    def download_and_process_protein(self, protein_name: str, virus_name: str) -> pd.DataFrame:
        """
        Download and process protein data from UniProt API.
        
        Args:
            protein_name (str): Name of the protein
            virus_name (str): Name of the virus for filtering
            
        Returns:
            pd.DataFrame: Processed protein data
        """
        try:
            # Download FASTA content
            fasta_content = self.download_protein_fasta(protein_name)
            
            if not fasta_content:
                self.logger.warning(f"Could not download FASTA for {protein_name}")
                return pd.DataFrame()
            
            # Process FASTA content
            processed_data = self.extract_data_from_fasta_content(fasta_content, protein_name, virus_name)
            
            return processed_data
            
        except Exception as e:
            self.logger.error(f"Error downloading and processing {protein_name}: {e}")
            return pd.DataFrame()
    
    def find_antigenic_entries(self, prot_id: str, keys: List[str]) -> Dict[str, str]:
        """Buscar las proteínas correspondientes a cada clave en UniProt."""
        try:
            # URL de la API UniProt
            url = f"https://rest.uniprot.org/uniprotkb/stream?compressed=true&fields=accession%2Creviewed%2Cid%2Cprotein_name%2Cgene_names%2Corganism_name%2Clength&format=tsv&query=%28%28proteome%3A{prot_id}%29%29"
            
            # Enviar la solicitud GET a la API
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            
            # Descomprimir el contenido gz
            with gzip.GzipFile(fileobj=io.BytesIO(response.content)) as gz:
                content = gz.read().decode('utf-8')
            
            # Leer el contenido en un DataFrame
            df_api = pd.read_csv(io.StringIO(content), sep='\t')
            
            # Crear un diccionario para almacenar los resultados de las claves
            matching_entries = {}
            for key in keys:
                # Filtrar las proteínas que contienen el término de búsqueda (key)
                if 'Protein names' in df_api.columns:
                    prot_indexes = df_api["Protein names"].str.contains(key, case=False, na=False)
                    matching_entries[key] = ", ".join(df_api[prot_indexes]["Entry"].tolist())
                else:
                    matching_entries[key] = "NaN"
            
            return matching_entries
            
        except Exception as e:
            self.logger.error(f"Error fetching UniProt data for {prot_id}: {e}")
            return {key: "NaN" for key in keys}
    
    def add_antigenic_data(self, df: pd.DataFrame, keys: List[str]) -> pd.DataFrame:
        """Añadir datos antigénicos al DataFrame de proteomas de forma optimizada."""
        if df.empty:
            self.logger.warning("DataFrame de proteomas vacío, saltando datos antigénicos")
            return df
            
        try:
            # Limitar a los primeros 10 proteomas más relevantes para acelerar el proceso
            if len(df) > 10:
                self.logger.info(f"Limitando procesamiento a los primeros 10 proteomas más relevantes (de {len(df)})")
                df_to_process = df.head(10)
            else:
                df_to_process = df

            def process_row(index):
                """Procesar cada fila del DataFrame."""
                matching_entries = self.find_antigenic_entries(index, keys)
                return index, matching_entries

            # Procesamiento paralelo optimizado
            n_jobs = min(10, len(df_to_process))  # Máximo 10 workers para evitar sobrecarga
            self.logger.info(f"Processing {len(df_to_process)} proteomes with {n_jobs} workers (optimizado)...")
            
            results = Parallel(n_jobs=n_jobs)(
                delayed(process_row)(index) for index in tqdm(df_to_process.index, desc="Retrieving antigenic proteins")
            )

            # Actualizar el DataFrame con los resultados
            df_copy = df.copy()
            for index, matching_entries in results:
                for key in keys:
                    df_copy.loc[index, f"{key} IDs"] = matching_entries.get(key, "NaN")
            
            self.logger.info("Datos antigénicos añadidos (optimizado)")
            return df_copy
            
        except Exception as e:
            self.logger.error(f"Error adding antigenic data: {e}")
            return df
    
    def process_dynamic_uniprot_data(self, virus_name: str = "Influenza A virus", 
                                    target_proteins: List[str] = None) -> Dict[str, pd.DataFrame]:
        """
        Download and process UniProt data dynamically using APIs.
        
        Args:
            virus_name (str): Name of the virus
            target_proteins (List[str]): List of target proteins
            
        Returns:
            Dict[str, pd.DataFrame]: Dictionary with processed datasets
        """
        try:
            if target_proteins is None:
                target_proteins = [
                    "Neuraminidase", 
                    "Hemagglutinin", 
                    "Matrix protein 2", 
                    "Matrix protein 1", 
                    "Nucleoprotein"
                ]
            
            self.logger.info(f"Starting dynamic UniProt processing for {virus_name}")
            self.logger.info(f"Target proteins: {', '.join(target_proteins)}")
            
            # Download proteome
            self.logger.info("Downloading proteome...")
            proteome_db_incomplete = self.download_proteome(virus_name)
            
            # Add antigenic data
            self.logger.info("Adding antigenic data...")
            proteome_db = self.add_antigenic_data(proteome_db_incomplete, target_proteins)
            
            # Save proteome
            proteome_output = self.output_dir / f"proteome_{virus_name.replace(' ', '_')}.csv"
            proteome_db.to_csv(proteome_output)
            self.logger.info(f"Proteome data saved to: {proteome_output}")
            
            # Process each protein
            results = {"proteome_db": proteome_db}
            
            for protein in target_proteins:
                self.logger.info(f"Processing {protein}...")
                db = self.download_and_process_protein(protein, virus_name)
                results[protein] = db
                
                # Save each DataFrame individually
                output_file = self.output_dir / f"{protein.replace(' ', '_')}_data.csv"
                db.to_csv(output_file, index=False)
                self.logger.info(f"{protein} data saved to: {output_file}")
            
            return results
            
        except Exception as e:
            self.logger.error(f"Error in dynamic UniProt processing: {e}")
            return {}
    
    def process_proteomes_and_fasta(self, proteome_csv: str = None, fasta_dict: Dict[str, str] = None) -> Dict[str, pd.DataFrame]:
        """
        Legacy method for backward compatibility.
        Now uses dynamic downloading if files don't exist.
        """
        try:
            # If files are provided and exist, use legacy method
            if proteome_csv and Path(proteome_csv).exists() and fasta_dict:
                all_exist = all(Path(path).exists() for path in fasta_dict.values())
                if all_exist:
                    self.logger.info("Using legacy file-based processing...")
                    
                    # Read proteome
                    proteome_db_incomplete = pd.read_csv(proteome_csv, sep="\t").set_index("Proteome Id")
                    keys = list(fasta_dict.keys())
                    
                    # Add antigenic data
                    proteome_db = self.add_antigenic_data(proteome_db_incomplete, keys)
                    
                    # Save proteome
                    proteome_output = self.output_dir / "final_proteome_data.csv"
                    proteome_db.to_csv(proteome_output)
                    self.logger.info(f"Proteome data saved to: {proteome_output}")
                    
                    # Process FASTA files
                    results = {"proteome_db": proteome_db}
                    
                    for key, path in fasta_dict.items():
                        self.logger.info(f"Processing FASTA file for {key}: {path}")
                        db = self.extract_data(path, key)
                        results[key] = db
                        
                        # Save each DataFrame individually
                        output_file = self.output_dir / f"{key.replace(' ', '_')}_data.csv"
                        db.to_csv(output_file, index=False)
                        self.logger.info(f"{key} data saved to: {output_file}")
                    
                    return results
            
            # Default to dynamic processing
            self.logger.info("Files not found or not provided, using dynamic API downloading...")
            return self.process_dynamic_uniprot_data()
            
        except Exception as e:
            self.logger.error(f"Error in UniProt processing: {e}")
            return {}

class ConfigurationManager:
    """Manages configuration loading and validation."""
    
    def __init__(self, config_path: str = "config.json"):
        """Initialize configuration manager."""
        self.config_path = Path(config_path)
        self.config = self._load_config()
        
    def _load_config(self) -> Dict:
        """Load configuration from file."""
        if not self.config_path.exists():
            raise FileNotFoundError(f"Configuration file not found: {self.config_path}")
        
        with open(self.config_path, 'r') as f:
            config = json.load(f)
        
        # Validate required fields
        required_fields = ["email", "batch_size", "max_workers"]
        for field in required_fields:
            if field not in config:
                raise ValueError(f"Missing required configuration field: {field}")
        
        return config
    
    def get(self, key: str, default=None):
        """Get configuration value with dot notation support."""
        keys = key.split('.')
        value = self.config
        
        for k in keys:
            if isinstance(value, dict) and k in value:
                value = value[k]
            else:
                return default
        
        return value

class NCBIDownloaderV2:
    """Enhanced NCBI FASTA downloader with additional optimizations."""
    
    def __init__(self, config_path: str = "config.json"):
        """Initialize the enhanced NCBI downloader."""
        # Load configuration
        self.config = ConfigurationManager(config_path)
        
        # Core settings
        self.email = self.config.get('email')
        self.batch_size = self.config.get('batch_size', 200)
        self.max_workers = self.config.get('max_workers', 8)
        self.rate_limit_delay = self.config.get('rate_limit_delay', 0.15)
        self.max_retries = self.config.get('max_retries', 5)
        self.retry_delay = self.config.get('retry_delay', 1)
        self.save_interval = self.config.get('save_interval', 120)
        
        # Adaptive batching settings
        self.adaptive_batching = self.config.get('adaptive_batching', True)
        self.min_batch_size = self.config.get('min_batch_size', 50)
        self.max_batch_size = self.config.get('max_batch_size', 500)
        self.current_batch_size = self.batch_size
        
        # Performance tracking
        self.success_rate = 1.0
        self.avg_response_time = 1.0
        self.consecutive_failures = 0
        
        # File paths
        self.output_dir = Path(self.config.get('output_dir', 'output'))
        self.output_dir.mkdir(exist_ok=True)
        
        self.output_file = self.output_dir / self.config.get('files.output_file', 'all_sequences.fasta')
        self.progress_file = self.output_dir / self.config.get('files.progress_file', 'progress.json')
        self.failed_ids_file = self.output_dir / self.config.get('files.failed_ids_file', 'failed_ids.txt')
        self.log_file = self.output_dir / self.config.get('files.log_file', 'ncbi_download.log')
        
        # Enhanced features
        self.compression_enabled = self.config.get('compression', False)
        self.checksum_validation = self.config.get('checksum_validation', True)
        self.duplicate_detection = self.config.get('duplicate_detection', True)
        
        # Thread-safe locks
        self.file_lock = threading.Lock()
        self.progress_lock = threading.Lock()
        self.rate_limit_lock = threading.Lock()
        self.stats_lock = threading.Lock()
        
        # Statistics and tracking
        self.stats = DownloadStats()
        self.downloaded_hashes: Set[str] = set() if self.duplicate_detection else None
        self.session_start_time = time.time()
        
        # Graceful shutdown
        self.shutdown_requested = False
        signal.signal(signal.SIGINT, self._signal_handler)
        signal.signal(signal.SIGTERM, self._signal_handler)
        
        # Setup logging FIRST (before using self.logger)
        self._setup_logging()
        
        # Initialize logger
        self.logger = logging.getLogger(__name__)
        
        # Log the configured output format
        rettype = self.config.get('ncbi.rettype', 'fasta')
        self.logger.info(f"Configured output format: {rettype.upper()}")
        if rettype in ['gb', 'genbank', 'gp']:
            self.logger.info("GenBank format selected - will download full GenBank records with annotations")
        else:
            self.logger.info("FASTA format selected - will download sequences only")
        
        # Set Entrez configuration
        Entrez.email = self.email
        Entrez.tool = "ncbi-fasta-downloader"
        
        # Initialize UniProt downloader
        self.uniprot_downloader = UniProtDownloader(
            output_dir=self.output_dir,
            max_workers=self.config.get('uniprot.max_workers', 35)
        )
        
        self.logger.info(f"Downloader initialized with batch_size={self.batch_size}, workers={self.max_workers}")
        
    def _setup_logging(self):
        """Setup logging configuration."""
        log_config_file = Path("logging_config.json")
        
        if log_config_file.exists():
            with open(log_config_file, 'r') as f:
                log_config = json.load(f)
                # Update log file path
                if 'handlers' in log_config and 'file' in log_config['handlers']:
                    log_config['handlers']['file']['filename'] = str(self.log_file)
                logging.config.dictConfig(log_config)
        else:
            # Fallback basic configuration
            logging.basicConfig(
                level=getattr(logging, self.config.get('log_level', 'INFO')),
                format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                handlers=[
                    logging.FileHandler(self.log_file),
                    logging.StreamHandler()
                ]
            )
        
    def _signal_handler(self, signum, frame):
        """Handle graceful shutdown on interrupt."""
        self.logger.info("Shutdown signal received. Saving progress...")
        self.shutdown_requested = True
        self._save_progress()
        sys.exit(0)
        
    def _calculate_hash(self, content: str) -> str:
        """Calculate hash for duplicate detection."""
        return hashlib.md5(content.encode('utf-8')).hexdigest()
        
    def _is_duplicate(self, content: str) -> bool:
        """Check if content is a duplicate."""
        if not self.duplicate_detection:
            return False
            
        content_hash = self._calculate_hash(content)
        if content_hash in self.downloaded_hashes:
            return True
        
        self.downloaded_hashes.add(content_hash)
        return False
        
    def fetch_accession_numbers(self, query: str, max_results: int = None) -> List[str]:
        """
        Fetch accession numbers from NCBI based on a search query.
        
        Args:
            query (str): NCBI search query
            max_results (int): Maximum number of results to fetch (None for all)
            
        Returns:
            List[str]: List of accession numbers
        """
        try:
            self.logger.info(f"Searching NCBI with query: {query}")
            
            # Initial search to get total count
            handle = Entrez.esearch(db="protein", term=query)
            record = Entrez.read(handle)
            handle.close()
            
            total_count = int(record["Count"])
            self.logger.info(f"Found {total_count:,} potential sequences")
            
            if total_count == 0:
                self.logger.warning("No sequences found for this query")
                return []
            
            # Limit results if specified
            if max_results and max_results < total_count:
                fetch_count = max_results
                self.logger.info(f"Limiting to {max_results:,} sequences")
            else:
                fetch_count = total_count
            
            ids = []
            batch_size = 10000  # NCBI's maximum for esearch
            
            # Fetch IDs in batches
            with tqdm(total=fetch_count, desc="Fetching IDs", unit="ids") as pbar:
                while len(ids) < fetch_count:
                    try:
                        # Calculate how many to fetch in this batch
                        remaining = fetch_count - len(ids)
                        current_batch_size = min(batch_size, remaining)
                        
                        # Fetch batch
                        handle = Entrez.esearch(
                            db="protein", 
                            term=query, 
                            retmax=current_batch_size, 
                            retstart=len(ids)
                        )
                        record = Entrez.read(handle)
                        handle.close()
                        
                        batch_ids = record["IdList"]
                        if not batch_ids:
                            self.logger.warning("No more IDs returned, stopping")
                            break
                            
                        ids.extend(batch_ids)
                        pbar.update(len(batch_ids))
                        
                        # Rate limiting
                        time.sleep(0.1)
                        
                    except Exception as e:
                        self.logger.error(f"Error fetching batch: {e}")
                        time.sleep(5)  # Wait longer on error
                        continue
            
            # Remove duplicates while preserving order
            unique_ids = list(dict.fromkeys(ids))
            self.logger.info(f"Fetched {len(unique_ids):,} unique IDs")
            
            return unique_ids
            
        except Exception as e:
            self.logger.error(f"Error in fetch_accession_numbers: {e}")
            return []
    
    def generate_ncbi_ids_file(self, virus_name: str, target_proteins: List[str], 
                              output_filename: str = "ncbi_ids_dict.tsv",
                              max_results_per_protein: int = None) -> bool:
        """
        Generate NCBI IDs dictionary file based on virus and target proteins.
        
        Args:
            virus_name (str): Name of the virus (e.g., "Influenza A virus")
            target_proteins (List[str]): List of target protein names
            output_filename (str): Output filename for the TSV file
            max_results_per_protein (int): Maximum results per protein (None for all)
            
        Returns:
            bool: True if successful, False otherwise
        """
        try:
            self.logger.info(f"Generating NCBI IDs file for {virus_name}")
            self.logger.info(f"Target proteins: {', '.join(target_proteins)}")
            
            ncbi_ids_dict = {}
            
            # Process each protein
            for protein in target_proteins:
                self.logger.info(f"Processing {protein}...")
                
                # Build search query similar to your original code
                query = f'("{virus_name}"[Organism] OR ("{virus_name}"[Organism] OR ("{virus_name}"[Organism] OR {virus_name}[All Fields]))) AND "{virus_name}"[porgn] AND "{protein}"[All Fields]'
                
                # Fetch accession numbers
                ids = self.fetch_accession_numbers(query, max_results_per_protein)
                
                if ids:
                    ncbi_ids_dict[protein] = ids
                    self.logger.info(f"✅ {protein}: {len(ids):,} IDs collected")
                else:
                    self.logger.warning(f"⚠️ {protein}: No IDs found")
                    ncbi_ids_dict[protein] = []
                
                # Small delay between proteins to be respectful to NCBI
                time.sleep(1)
            
            # Write to file
            output_path = self.output_dir / output_filename
            with open(output_path, 'w', encoding='utf-8') as f:
                # Write header
                f.write("# NCBI IDs Dictionary\n")
                f.write(f"# Generated: {datetime.now().isoformat()}\n")
                f.write(f"# Virus: {virus_name}\n")
                f.write(f"# Proteins: {', '.join(target_proteins)}\n")
                f.write("# Format: protein_name<TAB>comma_separated_ids\n")
                f.write("#\n")
                
                # Write data
                total_ids = 0
                for protein, ids in ncbi_ids_dict.items():
                    if ids:
                        f.write(f"{protein}\t{','.join(ids)}\n")
                        total_ids += len(ids)
                    else:
                        f.write(f"{protein}\t\n")  # Empty line for proteins with no IDs
                
            self.logger.info(f"✅ NCBI IDs file generated: {output_path}")
            self.logger.info(f"📊 Summary:")
            for protein, ids in ncbi_ids_dict.items():
                self.logger.info(f"  - {protein}: {len(ids):,} IDs")
            self.logger.info(f"📋 Total IDs: {total_ids:,}")
            
            return True
            
        except Exception as e:
            self.logger.error(f"Error generating NCBI IDs file: {e}")
            return False
    
    def load_ids_dict(self, filename: str) -> Dict[str, List[str]]:
        """Load the NCBI IDs dictionary from TSV file with enhanced validation."""
        self.logger.info(f"Loading IDs from {filename}")
        ncbi_ids_dict = {}
        total_ids = 0
        
        try:
            with open(filename, 'r', encoding='utf-8') as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    if not line or line.startswith('#'):  # Skip empty lines and comments
                        continue
                        
                    parts = line.split('\t')
                    if len(parts) >= 2:
                        key = parts[0].strip()
                        ids = [id_str.strip() for id_str in parts[1].split(',') if id_str.strip()]
                        
                        # Validate IDs format (should be numeric for NCBI)
                        valid_ids = []
                        for id_str in ids:
                            if id_str.isdigit():
                                valid_ids.append(id_str)
                            else:
                                self.logger.warning(f"Invalid ID format: {id_str} in line {line_num}")
                        
                        if valid_ids:
                            ncbi_ids_dict[key] = valid_ids
                            total_ids += len(valid_ids)
                            self.logger.info(f"Loaded {key}: {len(valid_ids)} IDs")
                        else:
                            self.logger.warning(f"No valid IDs found for {key} in line {line_num}")
                    else:
                        self.logger.warning(f"Skipping malformed line {line_num}: {line}")
                        
        except FileNotFoundError:
            self.logger.error(f"File {filename} not found")
            raise
        except Exception as e:
            self.logger.error(f"Error loading {filename}: {e}")
            raise
            
        self.logger.info(f"Total valid IDs loaded: {total_ids}")
        return ncbi_ids_dict
        
    def _load_progress(self) -> Dict:
        """Load progress from JSON file with enhanced error handling."""
        if self.progress_file.exists():
            try:
                with open(self.progress_file, 'r') as f:
                    progress = json.load(f)
                    
                # Validate progress format
                if 'processed' in progress:
                    self.logger.info(f"Loaded progress: {progress.get('processed', 0)} sequences processed")
                    
                    # Load downloaded hashes if duplicate detection is enabled
                    if self.duplicate_detection and 'downloaded_hashes' in progress:
                        self.downloaded_hashes = set(progress['downloaded_hashes'])
                        self.logger.info(f"Loaded {len(self.downloaded_hashes)} content hashes")
                    
                    return progress
                else:
                    self.logger.warning("Invalid progress file format")
                    
            except json.JSONDecodeError:
                self.logger.warning("Corrupted progress file, starting fresh")
            except Exception as e:
                self.logger.warning(f"Could not load progress file: {e}")
        
        return {'processed': 0, 'failed_ids': [], 'downloaded_hashes': []}
        
    def _save_progress(self):
        """Save progress to JSON file with enhanced data."""
        with self.progress_lock:
            progress = {
                'processed': self.stats.total_processed,
                'successful': self.stats.successful,
                'failed': self.stats.failed,
                'bytes_downloaded': self.stats.bytes_downloaded,
                'sequences_per_protein': self.stats.sequences_per_protein,
                'session_start': self.session_start_time,
                'timestamp': datetime.now().isoformat(),
                'config_snapshot': {
                    'batch_size': self.batch_size,
                    'max_workers': self.max_workers,
                    'rate_limit_delay': self.rate_limit_delay,
                    'current_batch_size': self.current_batch_size,
                    'success_rate': self.success_rate,
                    'avg_response_time': self.avg_response_time,
                    'consecutive_failures': self.consecutive_failures
                }
            }
            
            # Include hashes if duplicate detection is enabled
            if self.duplicate_detection and self.downloaded_hashes:
                progress['downloaded_hashes'] = list(self.downloaded_hashes)
            
            try:
                # Atomic write using temporary file
                temp_file = self.progress_file.with_suffix('.tmp')
                with open(temp_file, 'w') as f:
                    json.dump(progress, f, indent=2)
                temp_file.replace(self.progress_file)
                
                self.stats.last_save_time = time.time()
                self.logger.debug("Progress saved successfully")
            except Exception as e:
                self.logger.error(f"Could not save progress: {e}")
                
    def _adapt_batch_size(self, success: bool, response_time: float):
        """Dynamically adapt batch size based on performance."""
        if not self.adaptive_batching:
            return
            
        # Update performance metrics
        if success:
            self.consecutive_failures = 0
            self.success_rate = 0.9 * self.success_rate + 0.1 * 1.0
            self.avg_response_time = 0.8 * self.avg_response_time + 0.2 * response_time
            
            # Increase batch size if performing well
            if self.success_rate > 0.95 and self.avg_response_time < 3.0:
                new_size = min(self.max_batch_size, int(self.current_batch_size * 1.1))
                if new_size > self.current_batch_size:
                    self.current_batch_size = new_size
                    self.logger.debug(f"Increased batch size to {self.current_batch_size}")
                    
        else:
            self.consecutive_failures += 1
            self.success_rate = 0.9 * self.success_rate + 0.1 * 0.0
            
            # Decrease batch size on failures
            if self.consecutive_failures >= 2 or self.success_rate < 0.8:
                new_size = max(self.min_batch_size, int(self.current_batch_size * 0.8))
                if new_size < self.current_batch_size:
                    self.current_batch_size = new_size
                    self.logger.info(f"Decreased batch size to {self.current_batch_size} due to failures")
    
    def _get_optimal_rate_limit(self) -> float:
        """Calculate optimal rate limit based on current performance."""
        base_delay = self.rate_limit_delay
        
        if self.consecutive_failures > 0:
            # Increase delay exponentially with consecutive failures
            multiplier = 2 ** min(self.consecutive_failures, 4)
            return base_delay * multiplier
        elif self.success_rate > 0.98 and self.avg_response_time < 2.0:
            # Reduce delay if performing exceptionally well
            return base_delay * 0.7
        else:
            return base_delay
                
    def _save_failed_id(self, accession_id: str, protein_type: str, error: str):
        """Save failed ID to file with enhanced logging."""
        with self.file_lock:
            try:
                with open(self.failed_ids_file, 'a') as f:
                    timestamp = datetime.now().isoformat()
                    f.write(f"{timestamp}\t{protein_type}\t{accession_id}\t{error}\n")
            except Exception as e:
                self.logger.error(f"Could not save failed ID {accession_id}: {e}")
                
    def fetch_sequences_batch_enhanced(self, accession_ids: List[str], protein_type: str) -> List[Tuple[str, str]]:
        """Enhanced batch fetching with adaptive sizing and optimized rate limiting."""
        if not accession_ids:
            return []
            
        results = []
        retry_count = 0
        start_time = time.time()
        
        while retry_count <= self.max_retries:
            # Adaptive rate limiting
            optimal_delay = self._get_optimal_rate_limit()
            with self.rate_limit_lock:
                # Add small jitter to prevent thundering herd
                jitter = optimal_delay * 0.1 * (0.5 - abs(hash(threading.current_thread().name) % 100 - 50) / 100)
                time.sleep(optimal_delay + jitter)
                
            try:
                # Join IDs for batch request
                id_list = ','.join(accession_ids)
                
                # Enhanced Entrez call with timeout
                handle = Entrez.efetch(
                    db=self.config.get('ncbi.database', 'protein'),
                    id=id_list,
                    rettype=self.config.get('ncbi.rettype', 'fasta'),
                    retmode=self.config.get('ncbi.retmode', 'text')
                )
                
                response = handle.read()
                handle.close()
                
                if isinstance(response, bytes):
                    response = response.decode('utf-8')
                    
                # Update bytes downloaded
                with self.stats_lock:
                    self.stats.bytes_downloaded += len(response.encode('utf-8'))
                    
                # Parse and process sequences
                if response.strip():
                    # Choose parser based on rettype configuration
                    rettype = self.config.get('ncbi.rettype', 'fasta')
                    if rettype in ['gb', 'genbank', 'gp']:
                        sequences = self._parse_genbank_batch_enhanced(response, protein_type)
                    else:
                        sequences = self._parse_fasta_batch_enhanced(response, protein_type)
                    
                    # Filter duplicates if enabled
                    if self.duplicate_detection:
                        unique_sequences = []
                        for acc_id, fasta_data in sequences:
                            if not self._is_duplicate(fasta_data):
                                unique_sequences.append((acc_id, fasta_data))
                            else:
                                self.logger.debug(f"Duplicate sequence detected: {acc_id}")
                        sequences = unique_sequences
                    
                    results.extend(sequences)
                    self.logger.debug(f"Fetched {len(sequences)} unique sequences from batch of {len(accession_ids)}")
                    break  # Success, exit retry loop
                else:
                    self.logger.warning(f"Empty response for batch: {accession_ids[:3]}...")
                    break
                    
            except HTTPError as e:
                if "429" in str(e) or "too many requests" in str(e).lower():
                    wait_time = (2 ** retry_count) * 2  # Exponential backoff
                    self.logger.warning(f"Rate limit hit, waiting {wait_time}s (attempt {retry_count + 1})")
                    time.sleep(wait_time)
                elif "400" in str(e):
                    self.logger.error(f"Bad request for batch, skipping: {e}")
                    break
                else:
                    self.logger.error(f"HTTP error for batch (attempt {retry_count + 1}): {e}")
                    if retry_count == self.max_retries:
                        for acc_id in accession_ids:
                            self._save_failed_id(acc_id, protein_type, f"HTTP error: {e}")
                            
            except URLError as e:
                self.logger.error(f"Network error (attempt {retry_count + 1}): {e}")
                if retry_count == self.max_retries:
                    for acc_id in accession_ids:
                        self._save_failed_id(acc_id, protein_type, f"Network error: {e}")
                        
            except Exception as e:
                self.logger.error(f"Unexpected error (attempt {retry_count + 1}): {e}")
                if retry_count == self.max_retries:
                    for acc_id in accession_ids:
                        self._save_failed_id(acc_id, protein_type, f"Unexpected error: {e}")
                        
            retry_count += 1
            if retry_count <= self.max_retries:
                time.sleep(self.retry_delay * retry_count)
        
        # Track performance for adaptive optimization
        response_time = time.time() - start_time
        success = len(results) > 0
        self._adapt_batch_size(success, response_time)
                
        return results
        
    def _parse_fasta_batch_enhanced(self, fasta_text: str, protein_type: str) -> List[Tuple[str, str]]:
        """Enhanced FASTA parsing with validation and metadata."""
        sequences = []
        current_header = None
        current_sequence = []
        
        for line in fasta_text.split('\n'):
            line = line.strip()
            if line.startswith('>'):
                # Save previous sequence if exists
                if current_header and current_sequence:
                    full_sequence = current_header + '\n' + '\n'.join(current_sequence)
                    
                    # Enhanced header with metadata
                    enhanced_header = f"{current_header} [protein_type={protein_type}] [downloaded={datetime.now().isoformat()}]"
                    enhanced_sequence = enhanced_header + '\n' + '\n'.join(current_sequence)
                    
                    # Extract accession ID
                    acc_id = current_header.split()[0][1:]  # Remove '>' and get first part
                    sequences.append((acc_id, enhanced_sequence))
                    
                # Start new sequence
                current_header = line
                current_sequence = []
            elif line and current_header:
                current_sequence.append(line)
                
        # Don't forget the last sequence
        if current_header and current_sequence:
            enhanced_header = f"{current_header} [protein_type={protein_type}] [downloaded={datetime.now().isoformat()}]"
            enhanced_sequence = enhanced_header + '\n' + '\n'.join(current_sequence)
            acc_id = current_header.split()[0][1:]
            sequences.append((acc_id, enhanced_sequence))
            
        return sequences
    
    def _parse_genbank_batch_enhanced(self, genbank_text: str, protein_type: str) -> List[Tuple[str, str]]:
        """
        Enhanced GenBank parsing with validation and metadata extraction.
        
        Args:
            genbank_text (str): GenBank format text
            protein_type (str): Type of protein being processed
            
        Returns:
            List[Tuple[str, str]]: List of (accession_id, genbank_record) tuples
        """
        sequences = []
        
        try:
            # Use BioPython's SeqIO to parse GenBank format
            from io import StringIO
            genbank_io = StringIO(genbank_text)
            
            for record in SeqIO.parse(genbank_io, "genbank"):
                # Get accession ID
                acc_id = record.id
                if not acc_id and record.name:
                    acc_id = record.name
                
                # Add custom annotations for tracking
                record.annotations['protein_type'] = protein_type
                record.annotations['downloaded'] = datetime.now().isoformat()
                
                # Convert record back to GenBank format string
                output_io = StringIO()
                SeqIO.write(record, output_io, "genbank")
                genbank_data = output_io.getvalue()
                
                sequences.append((acc_id, genbank_data))
                
                self.logger.debug(f"Parsed GenBank record: {acc_id} ({len(record.seq)} bp)")
                
        except Exception as e:
            self.logger.error(f"Error parsing GenBank batch: {e}")
            self.logger.debug(f"Problematic GenBank text (first 500 chars): {genbank_text[:500]}")
            
        return sequences
        
    def _append_sequences_to_file(self, sequences: List[Tuple[str, str]], protein_type: str):
        """Enhanced file writing with compression support."""
        if not sequences:
            return
            
        with self.file_lock:
            try:
                file_path = self.output_file
                
                # Add compression if enabled
                if self.compression_enabled:
                    file_path = file_path.with_suffix(file_path.suffix + '.gz')
                    opener = gzip.open
                    mode = 'at'
                else:
                    opener = open
                    mode = 'a'
                
                with opener(file_path, mode, encoding='utf-8') as f:
                    for _, fasta_data in sequences:
                        f.write(fasta_data + '\n')
                
                # Update statistics
                with self.stats_lock:
                    self.stats.successful += len(sequences)
                    if protein_type not in self.stats.sequences_per_protein:
                        self.stats.sequences_per_protein[protein_type] = 0
                    self.stats.sequences_per_protein[protein_type] += len(sequences)
                    
            except Exception as e:
                self.logger.error(f"Error writing to file: {e}")
                for acc_id, _ in sequences:
                    self._save_failed_id(acc_id, protein_type, f"File write error: {e}")
                    
    def process_batch_worker_enhanced(self, batch_data: Tuple[List[str], str]) -> int:
        """Enhanced worker function with better error handling and statistics."""
        accession_batch, protein_type = batch_data
        
        if self.shutdown_requested:
            return 0
            
        batch_start_time = time.time()
        
        try:
            self.logger.info(f"Processing batch of {len(accession_batch)} {protein_type} sequences")
            
            # Fetch sequences with enhanced method
            sequences = self.fetch_sequences_batch_enhanced(accession_batch, protein_type)
            
            # Save successful sequences
            if sequences:
                self._append_sequences_to_file(sequences, protein_type)
                
                batch_duration = time.time() - batch_start_time
                rate = len(sequences) / batch_duration if batch_duration > 0 else 0
                self.logger.info(f"Batch completed: {len(sequences)} sequences in {batch_duration:.1f}s ({rate:.1f} seq/s)")
            
            # Update statistics
            with self.stats_lock:
                self.stats.total_processed += len(accession_batch)
                self.stats.failed += len(accession_batch) - len(sequences)
                
            return len(sequences)
            
        except Exception as e:
            self.logger.error(f"Error processing batch: {e}")
            with self.stats_lock:
                self.stats.total_processed += len(accession_batch)
                self.stats.failed += len(accession_batch)
            return 0
            
    def download_all_sequences(self, ncbi_ids_dict: Dict[str, List[str]]):
        """Enhanced download method with improved progress tracking and optimization."""
        # Load progress
        progress = self._load_progress()
        total_processed = progress.get('processed', 0)
        
        # Create flat list with protein type information
        all_accessions = []
        current_count = 0
        resume_point = 0
        
        for protein_type, accession_list in ncbi_ids_dict.items():
            for accession in accession_list:
                if current_count < total_processed:
                    current_count += 1
                    continue
                if current_count == total_processed:
                    resume_point = len(all_accessions)
                all_accessions.append((accession, protein_type))
                current_count += 1
                
        total_sequences = current_count
        remaining_accessions = all_accessions[resume_point:]
        
        self.logger.info(f"Dataset summary:")
        for protein_type, accession_list in ncbi_ids_dict.items():
            self.logger.info(f"  {protein_type}: {len(accession_list):,} sequences")
        
        self.logger.info(f"Total sequences: {total_sequences:,}")
        self.logger.info(f"Already processed: {total_processed:,}")
        self.logger.info(f"Remaining: {len(remaining_accessions):,}")
        
        if not remaining_accessions:
            self.logger.info("All sequences already processed!")
            return
            
        # Initialize statistics
        self.stats.start_time = time.time()
        self.stats.total_processed = total_processed
        
        # Create optimized batches by protein type
        batches = []
        protein_batches = {}
        
        # Group remaining accessions by protein type
        for accession, protein_type in remaining_accessions:
            if protein_type not in protein_batches:
                protein_batches[protein_type] = []
            protein_batches[protein_type].append(accession)
        
        # Create batches within each protein type with adaptive sizing
        for protein_type, accessions in protein_batches.items():
            current_pos = 0
            while current_pos < len(accessions):
                # Use current adaptive batch size
                batch_size_to_use = self.current_batch_size if self.adaptive_batching else self.batch_size
                batch = accessions[current_pos:current_pos + batch_size_to_use]
                if batch:
                    batches.append((batch, protein_type))
                current_pos += batch_size_to_use
                    
        self.logger.info(f"Created {len(batches)} optimized batches for processing (adaptive batch size: {self.current_batch_size})")
        
        # Process batches with enhanced thread pool
        successful_sequences = 0
        
        try:
            with concurrent.futures.ThreadPoolExecutor(
                max_workers=self.max_workers,
                thread_name_prefix="NCBIWorker"
            ) as executor:
                
                # Submit batches with load balancing
                future_to_batch = {}
                for i, batch_data in enumerate(batches):
                    future = executor.submit(self.process_batch_worker_enhanced, batch_data)
                    future_to_batch[future] = i
                
                # Process completed batches with enhanced monitoring
                completed_batches = 0
                for future in concurrent.futures.as_completed(future_to_batch):
                    if self.shutdown_requested:
                        break
                        
                    try:
                        batch_success_count = future.result(timeout=300)  # 5 minute timeout per batch
                        successful_sequences += batch_success_count
                        completed_batches += 1
                        
                        # Enhanced progress reporting
                        if completed_batches % 5 == 0 or completed_batches == len(batches):
                            elapsed = time.time() - self.stats.start_time
                            processed_pct = (completed_batches / len(batches)) * 100
                            current_rate = self.stats.successful / (elapsed / 60) if elapsed > 0 else 0
                            estimated_remaining = (len(batches) - completed_batches) / (completed_batches / elapsed) if completed_batches > 0 and elapsed > 0 else 0
                            current_optimal_delay = self._get_optimal_rate_limit()
                            
                            self.logger.info(
                                f"Progress: {completed_batches}/{len(batches)} batches ({processed_pct:.1f}%) - "
                                f"Success: {self.stats.successful:,}, Failed: {self.stats.failed:,} - "
                                f"Rate: {current_rate:.1f} seq/min - "
                                f"ETA: {timedelta(seconds=estimated_remaining)} - "
                                f"Size: {self.stats.bytes_downloaded / (1024*1024):.1f} MB - "
                                f"Adaptive: batch={self.current_batch_size}, rate={current_optimal_delay:.3f}s"
                            )
                            
                        # Auto-save progress
                        if time.time() - self.stats.last_save_time > self.save_interval:
                            self._save_progress()
                            
                    except concurrent.futures.TimeoutError:
                        self.logger.error(f"Batch timeout after 5 minutes")
                    except Exception as e:
                        self.logger.error(f"Batch processing error: {e}")
                        
        except KeyboardInterrupt:
            self.logger.info("Download interrupted by user")
        except Exception as e:
            self.logger.error(f"Error in batch processing: {e}")
        finally:
            # Final save and statistics
            self._save_progress()
            
            elapsed = time.time() - self.stats.start_time
            self.logger.info(f"Download session completed!")
            self.logger.info(f"Session statistics:")
            self.logger.info(f"  Duration: {timedelta(seconds=elapsed)}")
            self.logger.info(f"  Successful: {self.stats.successful:,}")
            self.logger.info(f"  Failed: {self.stats.failed:,}")
            self.logger.info(f"  Data downloaded: {self.stats.bytes_downloaded / (1024*1024):.1f} MB")
            self.logger.info(f"  Average rate: {self.stats.successful/(elapsed/60):.1f} sequences/minute")
            
            # Per-protein statistics
            self.logger.info("Per-protein statistics:")
            for protein_type, count in self.stats.sequences_per_protein.items():
                self.logger.info(f"  {protein_type}: {count:,} sequences")
    
    def download_uniprot_data(self, virus_name: str = None, target_proteins: List[str] = None,
                             proteome_csv: str = None, fasta_dict: Dict[str, str] = None):
        """Download and process UniProt data with dynamic parameters."""
        try:
            # Use dynamic approach if virus name and proteins are provided
            if virus_name or target_proteins:
                virus_name = virus_name or "Influenza A virus"
                target_proteins = target_proteins or [
                    "Neuraminidase", 
                    "Hemagglutinin", 
                    "Matrix protein 2", 
                    "Matrix protein 1", 
                    "Nucleoprotein"
                ]
                
                self.logger.info(f"Starting dynamic UniProt download for {virus_name}")
                results = self.uniprot_downloader.process_dynamic_uniprot_data(virus_name, target_proteins)
                
                if results:
                    self.logger.info("Dynamic UniProt processing completed successfully!")
                    self.logger.info(f"Processed datasets: {list(results.keys())}")
                    return results
                else:
                    self.logger.error("Dynamic UniProt processing failed!")
                    return None
            
            # Fallback to legacy method
            if proteome_csv is None:
                proteome_csv = self.config.get('uniprot.proteome_csv', 'proteomes_Influenza_A_2025_02_05.tsv')
            
            if fasta_dict is None:
                fasta_dict = self.config.get('uniprot.fasta_dict', {
                    "Neuraminidase": "uniprotkb_Neuraminidase_2025_02_05.fasta",
                    "Hemagglutinin": "uniprotkb_Hemagglutinin_2025_02_05.fasta",
                    "Matrix protein 2": "uniprotkb_Matrix_protein_2_AND_taxonomy_2025_06_18.fasta",
                    "Matrix protein 1": "uniprotkb_Matrix_protein_1_AND_taxonomy_2025_06_18.fasta",
                    "Nucleoprotein": "uniprotkb_Nucleoprotein_AND_taxonomy_id_2025_06_18.fasta"
                })
            
            self.logger.info("Starting legacy UniProt data download and processing...")
            
            # Process UniProt data (will auto-fallback to dynamic if files don't exist)
            results = self.uniprot_downloader.process_proteomes_and_fasta(proteome_csv, fasta_dict)
            
            if results:
                self.logger.info("UniProt data processing completed successfully!")
                self.logger.info(f"Processed datasets: {list(results.keys())}")
                return results
            else:
                self.logger.error("UniProt data processing failed!")
                return None
                
        except Exception as e:
            self.logger.error(f"Error in UniProt data download: {e}")
            return None
    
    def download_combined_data(self, ncbi_ids_dict: Dict[str, List[str]], 
                              proteome_csv: str = None, 
                              fasta_dict: Dict[str, str] = None):
        """Download both NCBI and UniProt data in parallel."""
        self.logger.info("Starting combined NCBI and UniProt data download...")
        
        # Use ThreadPoolExecutor to run both downloads in parallel
        with concurrent.futures.ThreadPoolExecutor(max_workers=2, thread_name_prefix="CombinedDownload") as executor:
            # Submit both tasks
            ncbi_future = executor.submit(self.download_all_sequences, ncbi_ids_dict)
            uniprot_future = executor.submit(self.download_uniprot_data, proteome_csv, fasta_dict)
            
            # Wait for both to complete
            try:
                # Wait for NCBI download
                ncbi_future.result()
                self.logger.info("NCBI download completed!")
                
                # Wait for UniProt download
                uniprot_results = uniprot_future.result()
                if uniprot_results:
                    self.logger.info("UniProt download completed!")
                else:
                    self.logger.warning("UniProt download failed!")
                
                self.logger.info("Combined download session completed!")
                
            except Exception as e:
                self.logger.error(f"Error in combined download: {e}")

def create_arg_parser() -> argparse.ArgumentParser:
    """Create command line argument parser."""
    parser = argparse.ArgumentParser(
        description="Enhanced NCBI FASTA Downloader",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage with existing NCBI IDs file
  python ncbi_downloader.py                              # Use default config and ncbi_ids_dict.tsv
  python ncbi_downloader.py --input my_ids.tsv           # Use specific input file
  
  # Generate new NCBI IDs file
  python ncbi_downloader.py --generate-ids               # Generate IDs for default proteins
  python ncbi_downloader.py --generate-ids --virus-name "Influenza B virus" # Different virus
  python ncbi_downloader.py --generate-ids --target-proteins Neuraminidase Hemagglutinin # Specific proteins
  python ncbi_downloader.py --generate-ids --max-ids-per-protein 1000 # Limit IDs per protein
  
  # Combined workflows
  python ncbi_downloader.py --generate-ids --ncbi-only   # Generate IDs and download NCBI only
  python ncbi_downloader.py --generate-ids --virus-name "SARS-CoV-2" --target-proteins "spike protein" "nucleocapsid protein"
  
  # Data integration (automatic CSV/FASTA unification)
  python ncbi_downloader.py --integrate                  # Download and integrate automatically
  python ncbi_downloader.py --generate-ids --integrate   # Generate IDs, download, and integrate
  python ncbi_downloader.py --uniprot-only --integrate   # UniProt download with integration
  
  # Configuration and optimization
  python ncbi_downloader.py --config custom.json         # Use custom config
  python ncbi_downloader.py --batch-size 100 --workers 2 # Override settings
  python ncbi_downloader.py --uniprot-only --virus-name "Influenza A virus" # UniProt only
        """
    )
    
    parser.add_argument(
        '--input', '-i',
        type=str,
        default='ncbi_ids_dict.tsv',
        help='Input TSV file with NCBI IDs (default: ncbi_ids_dict.tsv)'
    )
    
    parser.add_argument(
        '--config', '-c',
        type=str,
        default='config.json',
        help='Configuration file path (default: config.json)'
    )
    
    parser.add_argument(
        '--batch-size', '-b',
        type=int,
        help='Override batch size from config'
    )
    
    parser.add_argument(
        '--workers', '-w',
        type=int,
        help='Override number of workers from config'
    )
    
    parser.add_argument(
        '--output-dir', '-o',
        type=str,
        help='Override output directory from config'
    )
    
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Validate configuration and input without downloading'
    )
    
    parser.add_argument(
        '--resume',
        action='store_true',
        help='Force resume from existing progress (default behavior)'
    )
    
    parser.add_argument(
        '--fresh',
        action='store_true',
        help='Start fresh, ignoring existing progress'
    )
    
    parser.add_argument(
        '--uniprot-only',
        action='store_true',
        help='Download only UniProt data (skip NCBI)'
    )
    
    parser.add_argument(
        '--ncbi-only',
        action='store_true',
        help='Download only NCBI data (skip UniProt)'
    )
    
    parser.add_argument(
        '--proteome-csv',
        type=str,
        help='Path to proteome CSV file for UniProt processing'
    )
    
    parser.add_argument(
        '--fasta-dict',
        type=str,
        help='JSON string or file path with FASTA file dictionary for UniProt'
    )
    
    parser.add_argument(
        '--virus-name',
        type=str,
        default='Influenza A virus',
        help='Name of the virus to search for (default: Influenza A virus)'
    )
    
    parser.add_argument(
        '--target-proteins',
        type=str,
        nargs='+',
        help='List of target proteins (space-separated)'
    )
    
    parser.add_argument(
        '--generate-ids',
        action='store_true',
        help='Generate NCBI IDs file instead of using existing one'
    )
    
    parser.add_argument(
        '--max-ids-per-protein',
        type=int,
        help='Maximum number of IDs to fetch per protein (default: all available)'
    )
    
    parser.add_argument(
        '--ids-output',
        type=str,
        default='ncbi_ids_dict.tsv',
        help='Output filename for generated NCBI IDs file (default: ncbi_ids_dict.tsv)'
    )
    
    parser.add_argument(
        '--integrate',
        action='store_true',
        help='Automatically run data integration after successful download to create unified CSV and FASTA files'
    )
    
    return parser

def main():
    """Enhanced main function with command line interface."""
    parser = create_arg_parser()
    args = parser.parse_args()
    
    try:
        # Initialize downloader with configuration
        downloader = NCBIDownloaderV2(args.config)
        
        # Apply command line overrides
        if args.batch_size:
            downloader.batch_size = args.batch_size
        if args.workers:
            downloader.max_workers = args.workers
        if args.output_dir:
            downloader.output_dir = Path(args.output_dir)
            downloader.output_dir.mkdir(exist_ok=True)
        
        # Handle fresh start
        if args.fresh:
            if downloader.progress_file.exists():
                backup_file = downloader.progress_file.with_suffix('.backup')
                downloader.progress_file.rename(backup_file)
                downloader.logger.info(f"Previous progress backed up to {backup_file}")
        
        # Prepare UniProt parameters
        virus_name = args.virus_name
        target_proteins = args.target_proteins
        proteome_csv = args.proteome_csv
        fasta_dict = None
        
        if args.fasta_dict:
            # Try to parse as JSON string or load from file
            try:
                if Path(args.fasta_dict).exists():
                    with open(args.fasta_dict, 'r') as f:
                        fasta_dict = json.load(f)
                else:
                    fasta_dict = json.loads(args.fasta_dict)
            except Exception as e:
                downloader.logger.error(f"Error parsing fasta_dict: {e}")
                sys.exit(1)
        
        # Determine what to download
        download_ncbi = not args.uniprot_only
        download_uniprot = not args.ncbi_only
        
        if args.uniprot_only and args.ncbi_only:
            downloader.logger.error("Cannot specify both --uniprot-only and --ncbi-only")
            sys.exit(1)
        
        # Prepare default target proteins if not specified
        if not target_proteins:
            target_proteins = [
                "Neuraminidase", 
                "Hemagglutinin", 
                "Matrix protein 2", 
                "Matrix protein 1", 
                "Nucleoprotein"
            ]
        
        # Handle NCBI IDs generation or loading
        ncbi_ids_dict = {}
        if download_ncbi:
            if args.generate_ids:
                # Generate new NCBI IDs file
                downloader.logger.info("Generating NCBI IDs file...")
                downloader.logger.info(f"Virus: {virus_name}")
                downloader.logger.info(f"Proteins: {', '.join(target_proteins)}")
                
                success = downloader.generate_ncbi_ids_file(
                    virus_name=virus_name,
                    target_proteins=target_proteins,
                    output_filename=args.ids_output,
                    max_results_per_protein=args.max_ids_per_protein
                )
                
                if not success:
                    downloader.logger.error("Failed to generate NCBI IDs file")
                    sys.exit(1)
                
                # Load the generated file
                ids_file_path = downloader.output_dir / args.ids_output
                if ids_file_path.exists():
                    ncbi_ids_dict = downloader.load_ids_dict(str(ids_file_path))
                else:
                    downloader.logger.error(f"Generated file not found: {ids_file_path}")
                    sys.exit(1)
                    
            else:
                # Load existing NCBI IDs file
                if not Path(args.input).exists():
                    downloader.logger.error(f"Input file not found: {args.input}")
                    downloader.logger.info("Tip: Use --generate-ids to create a new NCBI IDs file")
                    sys.exit(1)
                ncbi_ids_dict = downloader.load_ids_dict(args.input)
        
        if args.dry_run:
            downloader.logger.info("Dry run completed successfully")
            if download_ncbi:
                total_sequences = sum(len(ids) for ids in ncbi_ids_dict.values())
                downloader.logger.info(f"NCBI sequences to process: {total_sequences:,}")
            if download_uniprot:
                downloader.logger.info("UniProt processing would be performed")
            return
        
        # Confirm before starting (unless resuming)
        if download_ncbi and not args.resume and not args.fresh:
            total_sequences = sum(len(ids) for ids in ncbi_ids_dict.values())
            estimated_time_hours = total_sequences / (200 * 60)  # Assuming 200 seq/min
            
            print(f"\nReady to download {total_sequences:,} NCBI sequences")
            if download_uniprot:
                print("UniProt data will also be processed")
            print(f"Estimated time: {estimated_time_hours:.1f} hours")
            print(f"Output directory: {downloader.output_dir}")
            
            response = input("Start download? [y/N]: ").strip().lower()
            if response not in ['y', 'yes']:
                print("Download cancelled.")
                return
        
        # Start downloads based on options
        if download_ncbi and download_uniprot:
            # Combined download with dynamic parameters
            downloader.logger.info("Starting combined NCBI and UniProt download...")
            
            # Use ThreadPoolExecutor for parallel downloads
            with concurrent.futures.ThreadPoolExecutor(max_workers=2, thread_name_prefix="CombinedDownload") as executor:
                # Submit both tasks
                ncbi_future = executor.submit(downloader.download_all_sequences, ncbi_ids_dict)
                uniprot_future = executor.submit(downloader.download_uniprot_data, virus_name, target_proteins, proteome_csv, fasta_dict)
                
                # Wait for both to complete
                try:
                    ncbi_future.result()
                    downloader.logger.info("NCBI download completed!")
                    
                    uniprot_results = uniprot_future.result()
                    if uniprot_results:
                        downloader.logger.info("UniProt download completed!")
                    else:
                        downloader.logger.warning("UniProt download failed!")
                    
                    downloader.logger.info("Combined download session completed!")
                    
                except Exception as e:
                    downloader.logger.error(f"Error in combined download: {e}")
                    
        elif download_ncbi:
            # NCBI only
            downloader.logger.info("Starting NCBI sequence download...")
            downloader.download_all_sequences(ncbi_ids_dict)
        elif download_uniprot:
            # UniProt only with dynamic parameters
            downloader.logger.info(f"Starting UniProt data download for {virus_name}...")
            if target_proteins:
                downloader.logger.info(f"Target proteins: {', '.join(target_proteins)}")
            downloader.download_uniprot_data(virus_name, target_proteins, proteome_csv, fasta_dict)
        
        # Run data integration if requested
        if args.integrate:
            downloader.logger.info("="*60)
            downloader.logger.info("STARTING DATA INTEGRATION")
            downloader.logger.info("="*60)
            
            try:
                # Import and run data integration
                from data_integrator import DataIntegrator
                
                # Initialize integrator with the same output directory
                integrator = DataIntegrator(str(downloader.output_dir))
                
                downloader.logger.info("Processing CSV files...")
                csv_files = integrator.identify_csv_files(exclude_patterns=["proteome", "unified"],
                    target_proteins=args.target_proteins  # Pass the target proteins
                    )
                csv_df = integrator.read_and_concatenate_csv_files(csv_files)
                
                downloader.logger.info("Processing FASTA files...")
                fasta_files = integrator.identify_fasta_files(exclude_patterns=["unified"])
                fasta_df = integrator.parse_fasta_files(fasta_files, target_proteins=args.target_proteins)
                
                if not csv_df.empty or not fasta_df.empty:
                    downloader.logger.info("Combining and processing data...")
                    combined_df = integrator.combine_csv_and_fasta_data(csv_df, fasta_df)
                    
                    if not combined_df.empty:
                        # Save unified outputs
                        output_csv = downloader.output_dir / "unified_data.csv"
                        combined_df.to_csv(output_csv, index=False)
                        downloader.logger.info(f"[INTEGRATION] Unified CSV saved: {output_csv}")
                        
                        # Create unified FASTA
                        fasta_success = integrator.create_unified_fasta(combined_df, "unified_sequences.fasta")
                        if fasta_success:
                            downloader.logger.info(f"[INTEGRATION] Unified FASTA saved: {downloader.output_dir / 'unified_sequences.fasta'}")
                        
                        # Log integration summary
                        downloader.logger.info(f"[INTEGRATION] Integration completed successfully!")
                        downloader.logger.info(f"[INTEGRATION] Final dataset: {len(combined_df):,} unique sequences")
                        downloader.logger.info(f"[INTEGRATION] Output files:")
                        downloader.logger.info(f"[INTEGRATION]   - {output_csv}")
                        downloader.logger.info(f"[INTEGRATION]   - {downloader.output_dir / 'unified_sequences.fasta'}")
                    else:
                        downloader.logger.warning("[INTEGRATION] No data to integrate")
                else:
                    downloader.logger.warning("[INTEGRATION] No source files found for integration")
                    
            except ImportError:
                downloader.logger.error("[INTEGRATION] data_integrator.py not found. Please ensure it's in the same directory.")
            except Exception as e:
                downloader.logger.error(f"[INTEGRATION] Integration failed: {e}")
                import traceback
                traceback.print_exc()
            
            downloader.logger.info("="*60)
            downloader.logger.info("DATA INTEGRATION COMPLETED")
            downloader.logger.info("="*60)
        
    except KeyboardInterrupt:
        print("\nDownload interrupted by user")
    except Exception as e:
        logging.error(f"Fatal error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
