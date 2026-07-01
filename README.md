# NCBI-Uniprot FASTA Downloader

[![Python](https://img.shields.io/badge/Python-3.8%2B-blue)](https://www.python.org)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Code Style](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

**[English](#english) | [Español](#español)**

---

## English

### Overview

A high-performance Python tool for downloading FASTA sequences from NCBI in bulk. Designed to efficiently handle large datasets with intelligent batch processing, parallel execution, and robust error handling.

## Features

- **Dynamic NCBI ID Generation**: Automatically generate NCBI ID files from search queries
- **Data Integration**: Automatically combine CSV and FASTA outputs into unified files
- **Batch Processing**: Download multiple sequences per API call for efficiency
- **Parallel Execution**: Configurable worker threads for optimal performance
- **Smart Rate Limiting**: Respects NCBI's API limits with intelligent timing
- **Robust Error Handling**: Automatic retries with exponential backoff
- **Progress Persistence**: Resume downloads from interruption points
- **Duplicate Detection**: Optional hash-based duplicate filtering
- **Comprehensive Logging**: Detailed logs with performance metrics
- **Graceful Shutdown**: Safe interruption handling with progress preservation
- **Flexible Workflow**: One-step or two-step processes for enhanced usability

## Installation

### Clone and Setup
```bash
git clone https://github.com/SimonAirInsti/ncbi-fasta-downloader.git
cd ncbi-fasta-downloader
pip install -r requirements.txt
python setup.py
```

## Configuration

The setup script will ask for your email address (required by NCBI) and create a `config.json` file:

```json
{
  "email": "your.email@domain.com",
  "batch_size": 200,
  "max_workers": 4,
  "rate_limit_delay": 0.34
}
```

### Prepare Input File (Optional)

**Option 1: Automatic Generation (Recommended)**
Use the new dynamic ID generation to automatically create NCBI ID files:
```bash
# Generate IDs and download in one step
python ncbi_downloader.py --generate-ids --virus "Influenza A virus" --proteins "Neuraminidase,Hemagglutinin,Matrix protein 1" --max-results-per-protein 1000

# Generate IDs file only (two-step process)
python ncbi_downloader.py --generate-ids --virus "Influenza A virus" --proteins "Neuraminidase,Hemagglutinin" --ids-only --output-filename my_ids.tsv
```

**Option 2: Manual TSV File**
Create a TSV file manually with your NCBI IDs:
```
Neuraminidase	3024910520,3024910499,3024910476
Hemagglutinin	3024910409,3024910387,3024910364
Matrix protein 1	164584001,164584002,164584003
```

## Quick Start

### Dynamic ID Generation (New!)
```bash
# Generate NCBI IDs and download sequences in one command
python ncbi_downloader.py --generate-ids --virus "Influenza A virus" --target-proteins "Neuraminidase,Hemagglutinin" --max-results-per-protein 500

# Generate IDs file only for later use
python ncbi_downloader.py --generate-ids --virus "Influenza A virus" --proteins "Neuraminidase,Hemagglutinin" --ids-only

# Use previously generated IDs file
python ncbi_downloader.py --ncbi-ids-file ncbi_ids_dict.tsv
```

### Data Integration (New!)
```bash
# Download and automatically create unified CSV/FASTA files
python ncbi_downloader.py --integrate

# Generate IDs, download, and integrate in one command  
python ncbi_downloader.py --generate-ids --integrate

# Download only UniProt data with integration
python ncbi_downloader.py --uniprot-only --integrate

# Manual integration of existing downloads
python data_integrator.py
```

**Integration Features:**
- Combines CSV and FASTA files from UniProt and NCBI sources
- Removes duplicate sequences (prioritizing UniProt entries)
- Filters out fragment sequences
- Creates unified `unified_data.csv` and `unified_sequences.fasta` files
- Calculates sequence lengths for all entries
- Standardizes column names across data sources

### Traditional Usage
```bash
# Run with default settings
python ncbi_downloader.py

# Use custom input file
python ncbi_downloader.py --input my_ids.tsv

# Override batch size and workers
python ncbi_downloader.py --batch-size 100 --workers 2
```

### Test with Sample Data
```bash
# Test with small dataset
python ncbi_downloader.py --input sample_ncbi_ids.tsv --dry-run
python ncbi_downloader.py --input sample_ncbi_ids.tsv
```

### Monitor Progress
```bash
# Check current status
python check_progress.py

# Follow logs in real-time
tail -f output/ncbi_download.log
```

### 📁 Project Structure

```
ncbi-fasta-downloader/
├── ncbi_downloader.py         # Main downloader script
├── data_integrator.py         # Data integration and unification script (New!)
├── setup.py                   # Interactive setup script
├── config.json                # Configuration (created by setup)
├── requirements.txt           # Python dependencies
├── check_progress.py          # Progress monitoring utility
├── sample_ncbi_ids.tsv        # Sample data for testing
├── test_suite.py              # Comprehensive test suite
├── output/                    # Output directory
│   ├── all_sequences.fasta    # Downloaded sequences
│   ├── unified_data.csv       # Integrated CSV data (with --integrate)
│   ├── unified_sequences.fasta # Integrated FASTA data (with --integrate)
│   ├── progress.json          # Detailed progress tracking
│   ├── failed_ids.txt         # Failed downloads log
│   └── ncbi_download.log      # Execution logs
└── .github/                   # GitHub Actions workflows
```

## 🔄 Data Integration

The `data_integrator.py` script provides powerful data integration capabilities to combine and unify CSV and FASTA outputs from different sources.

### Integration Features

- **🔗 Data Combination**: Merges CSV files (UniProt) and FASTA files (NCBI) into unified outputs
- **🔍 Duplicate Removal**: Eliminates duplicate sequences while prioritizing UniProt entries
- **🧹 Fragment Filtering**: Automatically removes sequences containing "fragment" in their descriptions
- **📏 Sequence Length Calculation**: Computes sequence lengths for all entries (including UniProt data)
- **📊 Column Standardization**: Unifies column names across data sources (`Prot_ID`, `Protein_Type`)
- **🧽 Clean Output**: Removes internal tracking/metadata columns from final dataset
- **📈 Statistics Logging**: Provides detailed logs of processing steps and data transformations

### Integration Usage

#### Automatic Integration (Recommended)
```bash
# Download and integrate automatically
python ncbi_downloader.py --integrate

# Generate IDs, download, and integrate in one command
python ncbi_downloader.py --generate-ids --integrate --virus-name "Influenza A virus" --target-proteins Neuraminidase Hemagglutinin

# UniProt only with integration
python ncbi_downloader.py --uniprot-only --integrate --virus-name "Influenza A virus" --target-proteins Neuraminidase Hemagglutinin
```

#### Manual Integration
```bash
# Integrate existing CSV and FASTA files
python data_integrator.py
```

### Integration Output Files

When integration is enabled, you'll get:

- **`unified_data.csv`**: Combined, deduplicated, and cleaned CSV data with standardized columns
- **`unified_sequences.fasta`**: Combined, deduplicated FASTA sequences
- **Original individual files**: Still created for reference and debugging

### Integration Process

1. **File Identification**: Automatically identifies CSV and FASTA files in the output directory
2. **CSV Processing**: Reads and concatenates CSV files from UniProt downloads
3. **FASTA Processing**: Parses FASTA files and extracts metadata (ID, protein type, organism, variants)
4. **Data Combination**: Merges CSV and FASTA data with consistent schema
5. **Sequence Length Calculation**: Computes lengths for all sequences (especially UniProt entries)
6. **Fragment Filtering**: Removes sequences with "fragment" in description
7. **Deduplication**: Removes duplicate sequences by sequence content (UniProt priority)
8. **Column Standardization**: Ensures consistent column names (`Prot_ID`, `Protein_Type`)
9. **Cleanup**: Removes internal metadata columns from final output
10. **File Generation**: Creates unified CSV and FASTA files

### Integration Statistics

The integration process provides detailed logging:

```
[SUCCESS] CSV processing completed! Shape: (290444, 11)
[SUCCESS] FASTA processing completed! Shape: (769743, 12)
[SUCCESS] Data combination completed! Shape: (215005, 10)

[SUMMARY] Final Integration Summary:
  - Total sequences: 215,005
  - UniProt source: 290,444 sequences  
  - NCBI source: 769,743 sequences
  - Duplicates removed: 760,111
  - Fragments removed: 85,071
```

### 🔧 Advanced Configuration

#### Performance Tuning
```bash
# High-performance setup (stable connection required)
python ncbi_downloader.py --batch-size 500 --workers 6

# Conservative setup (unstable connection)
python ncbi_downloader.py --batch-size 50 --workers 2
```

#### Command Line Options
```bash
python ncbi_downloader.py --help
```

**Core Options:**
- `--generate-ids`: Generate NCBI IDs file dynamically
- `--virus`: Target virus name (e.g., "Influenza A virus")
- `--proteins`: Comma-separated protein names
- `--max-results-per-protein`: Limit results per protein
- `--ids-only`: Generate IDs file without downloading
- `--output-filename`: Custom output filename for IDs file

**Traditional Options:**
- `--input`: Input TSV file path
- `--config`: Configuration file path  
- `--batch-size`: Sequences per API call
- `--workers`: Number of parallel workers
- `--dry-run`: Test mode without downloading

### 🔍 Monitoring and Troubleshooting

#### Check Progress
```bash
python check_progress.py
```

#### View Logs
```bash
# Recent activity
tail -50 output/ncbi_download.log

# Follow in real-time
tail -f output/ncbi_download.log

# Filter errors
grep ERROR output/ncbi_download.log
```

#### Resume Downloads
The downloader automatically resumes from the last checkpoint if interrupted.

## 🎯 Complete Workflow Examples

### Example 1: Quick Start for Influenza A Research
```bash
# One command to generate IDs, download, and create unified dataset
python ncbi_downloader.py --generate-ids --integrate \
  --virus-name "Influenza A virus" \
  --target-proteins Neuraminidase Hemagglutinin "Matrix protein 1" \
  --max-ids-per-protein 1000

# Result: unified_data.csv and unified_sequences.fasta with clean, deduplicated data
```

### Example 2: UniProt Only Workflow
```bash
# Download UniProt data with automatic integration
python ncbi_downloader.py --uniprot-only --integrate \
  --virus-name "Influenza A virus" \
  --target-proteins Neuraminidase Hemagglutinin "Matrix protein 2" "Matrix protein 1" Nucleoprotein

# Result: Clean CSV and FASTA files from UniProt only
```

### Example 3: Large-Scale Data Collection
```bash
# Step 1: Generate comprehensive ID list
python ncbi_downloader.py --generate-ids \
  --virus-name "Influenza A virus" \
  --target-proteins Neuraminidase Hemagglutinin "Matrix protein 2" "Matrix protein 1" Nucleoprotein \
  --max-ids-per-protein 10000 \
  --ids-output large_influenza_ids.tsv

# Step 2: Download with optimized settings
python ncbi_downloader.py --input large_influenza_ids.tsv \
  --batch-size 300 --workers 8 --integrate

# Result: Large unified dataset with optimized performance
```

### Example 4: Two-Step Process for Different Viruses
```bash
# Step 1: Download SARS-CoV-2 data
python ncbi_downloader.py --generate-ids --integrate \
  --virus-name "SARS-CoV-2" \
  --target-proteins "spike protein" "nucleocapsid protein" \
  --max-ids-per-protein 2000

# Step 2: Add Influenza data to same output directory
python ncbi_downloader.py --generate-ids \
  --virus-name "Influenza A virus" \
  --target-proteins Neuraminidase Hemagglutinin \
  --max-ids-per-protein 1000

# Step 3: Integrate all data manually
python data_integrator.py

# Result: Multi-virus unified dataset
```

### Example 5: Testing and Development
```bash
# Test configuration without downloading
python ncbi_downloader.py --generate-ids --dry-run \
  --virus-name "Influenza A virus" \
  --target-proteins Neuraminidase \
  --max-ids-per-protein 100

# Small test download with integration
python ncbi_downloader.py --generate-ids --integrate \
  --virus-name "Influenza A virus" \
  --target-proteins Neuraminidase \
  --max-ids-per-protein 50 \
  --batch-size 25

# Result: Small test dataset to verify workflow
```

### Example 6: Advanced Custom Configuration
```bash
# Create custom config file first
cat > custom_config.json << EOF
{
    "email": "researcher@university.edu",
    "batch_size": 250,
    "max_workers": 6,
    "rate_limit_delay": 0.1,
    "output_dir": "influenza_study_2024"
}
EOF

# Use custom configuration with integration
python ncbi_downloader.py --config custom_config.json --integrate \
  --generate-ids --virus-name "Influenza A virus" \
  --target-proteins Neuraminidase Hemagglutinin "Matrix protein 1" \
  --max-ids-per-protein 2000

# Result: Customized workflow with specific settings
```

### Integration Workflow Comparison

| Approach | Command | Result | Use Case |
|----------|---------|---------|----------|
| **Automatic** | `--integrate` | Unified files immediately | Most users, streamlined workflow |
| **Manual** | `python data_integrator.py` | Flexible integration timing | Advanced users, custom processing |
| **No integration** | Default behavior | Individual source files | Legacy workflows, debugging |

### Output File Structure

After running with `--integrate`:
```
output/
├── unified_data.csv           # 🎯 Main result: clean, unified CSV
├── unified_sequences.fasta    # 🎯 Main result: clean, unified FASTA
├── Neuraminidase_data.csv     # Individual UniProt CSV
├── Hemagglutinin_data.csv     # Individual UniProt CSV  
├── all_sequences.fasta        # Individual NCBI FASTA
├── proteome_*.csv             # Proteome reference data
└── ncbi_download.log          # Detailed processing logs
```

**Key Benefits:**
- ✅ **unified_data.csv**: Ready for analysis (no duplicates, no fragments, standardized columns)
- ✅ **unified_sequences.fasta**: Ready for bioinformatics tools (clean, deduplicated sequences)
- ✅ **Individual files**: Available for debugging and reference
- ✅ **Detailed logs**: Complete audit trail of all processing steps

## 🔧 Troubleshooting

### Common Issues and Solutions

#### Empty or Small Results
**Problem**: Few or no sequences downloaded
```bash
# Check if virus name is correct
python ncbi_downloader.py --generate-ids --dry-run \
  --virus-name "Influenza A virus" --target-proteins Neuraminidase

# Try broader search terms
python ncbi_downloader.py --generate-ids \
  --virus-name "Influenza" --target-proteins "neuraminidase"
```

#### Network/Rate Limiting Issues
**Problem**: Frequent timeouts or API errors
```bash
# Reduce batch size and increase delays
python ncbi_downloader.py --batch-size 50 --workers 2 \
  --config config_conservative.json
```

Create `config_conservative.json`:
```json
{
    "rate_limit_delay": 1.0,
    "timeout": 30,
    "max_retries": 5
}
```

#### Unicode/Encoding Errors
**Problem**: Character encoding issues in logs
- **Solution**: Already fixed in current version (ASCII-only log messages)
- Check that your terminal supports UTF-8 output

#### Integration Issues
**Problem**: Integration produces unexpected results
```bash
# Run integration manually to see detailed logs
python data_integrator.py

# Check integration statistics in log file
grep "Integration Statistics" output/ncbi_download.log
```

#### Memory Issues with Large Datasets
**Problem**: Out of memory with large downloads
```bash
# Process in smaller batches
python ncbi_downloader.py --batch-size 100 --max-ids-per-protein 5000

# Or split by protein
for protein in "Neuraminidase" "Hemagglutinin"; do
    python ncbi_downloader.py --generate-ids \
        --virus-name "Influenza A virus" \
        --target-proteins "$protein" \
        --max-ids-per-protein 10000
done
python data_integrator.py  # Integrate all results
```

### Log Analysis

Check `output/ncbi_download.log` for:
- **ID Generation**: Number of IDs found per protein
- **Download Progress**: Success/failure rates
- **Integration Statistics**: Duplicates removed, fragments filtered
- **Errors**: Detailed error messages with timestamps

Example log analysis:
```bash
# Check download statistics
grep "Successfully downloaded" output/ncbi_download.log | wc -l

# Check integration results
grep "Integration Statistics" output/ncbi_download.log

# Find errors
grep "ERROR" output/ncbi_download.log
```

### 🤝 Contributing

1. Fork the repository
2. Create a feature branch: `git checkout -b feature-name`
3. Make changes with tests
4. Submit a pull request

### 📝 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

### 🙏 Acknowledgments

- NCBI for providing the E-utilities API
- BioPython community for the excellent Bio library
- Contributors and users providing feedback

---

## Español

### Descripción General

Herramienta Python de alto rendimiento para descargar secuencias FASTA de NCBI en masa. Diseñada para manejar eficientemente grandes conjuntos de datos con procesamiento inteligente por lotes, ejecución paralela y manejo robusto de errores.

### 🚀 Características Principales

- **Generación Dinámica de IDs de NCBI**: Genera automáticamente archivos de IDs de NCBI desde consultas de búsqueda
- **Procesamiento por Lotes**: Descarga múltiples secuencias por llamada API para mayor eficiencia
- **Ejecución Paralela**: Hilos de trabajo configurables para rendimiento óptimo
- **Limitación Inteligente de Velocidad**: Respeta los límites de la API de NCBI con temporización inteligente
- **Manejo Robusto de Errores**: Reintentos automáticos con retroceso exponencial
- **Persistencia de Progreso**: Reanuda descargas desde puntos de interrupción
- **Detección de Duplicados**: Filtrado opcional de duplicados basado en hash
- **Registro Completo**: Logs detallados con métricas de rendimiento
- **Apagado Elegante**: Manejo seguro de interrupciones con preservación del progreso
- **Flujo de Trabajo Flexible**: Procesos de uno o dos pasos para mayor usabilidad

### 🛠️ Instalación y Configuración

#### 1. Clonar y Configurar
```bash
git clone https://github.com/SimonAirInsti/ncbi-fasta-downloader.git
cd ncbi-fasta-downloader
pip install -r requirements.txt
python setup.py
```

#### 2. Configuración
El script de configuración solicitará tu dirección de email (requerida por NCBI) y creará un archivo `config.json`:

```json
{
  "email": "tu.email@dominio.com",
  "batch_size": 200,
  "max_workers": 4,
  "rate_limit_delay": 0.34
}
```

#### 3. Preparar Archivo de Entrada (Opcional)

**Opción 1: Generación Automática (Recomendado)**
Usa la nueva generación dinámica de IDs para crear automáticamente archivos de IDs de NCBI:
```bash
# Generar IDs y descargar en un solo paso
python ncbi_downloader.py --generate-ids --virus "Influenza A virus" --proteins "Neuraminidase,Hemagglutinin,Matrix protein 1" --max-results-per-protein 1000

# Generar solo archivo de IDs (proceso de dos pasos)
python ncbi_downloader.py --generate-ids --virus "Influenza A virus" --proteins "Neuraminidase,Hemagglutinin" --ids-only --output-filename mis_ids.tsv
```

**Opción 2: Archivo TSV Manual**
Crea un archivo TSV manualmente con tus IDs de NCBI:
```
Neuraminidase	3024910520,3024910499,3024910476
Hemagglutinin	3024910409,3024910387,3024910364
Matrix protein 1	164584001,164584002,164584003
```

### 🎯 Inicio Rápido

#### Generación Dinámica de IDs (¡Nuevo!)
```bash
# Generar IDs de NCBI y descargar secuencias en un comando
python ncbi_downloader.py --generate-ids --virus "Influenza A virus" --proteins "Neuraminidase,Hemagglutinin" --max-results-per-protein 500

# Generar solo archivo de IDs para uso posterior
python ncbi_downloader.py --generate-ids --virus "Influenza A virus" --proteins "Neuraminidase,Hemagglutinin" --ids-only

# Usar archivo de IDs generado previamente
python ncbi_downloader.py --ncbi-ids-file ncbi_ids_dict.tsv
```

#### Uso Tradicional
```bash
# Ejecutar con configuración predeterminada
python ncbi_downloader.py

# Usar archivo de entrada personalizado
python ncbi_downloader.py --input mis_ids.tsv

# Sobrescribir tamaño de lote y trabajadores
python ncbi_downloader.py --batch-size 100 --workers 2
```

#### Probar con Datos de Muestra
```bash
# Probar con dataset pequeño
python ncbi_downloader.py --input sample_ncbi_ids.tsv --dry-run
python ncbi_downloader.py --input sample_ncbi_ids.tsv
```

#### Monitorear Progreso
```bash
# Verificar estado actual
python check_progress.py

# Seguir logs en tiempo real
tail -f output/ncbi_download.log
```

### 🔧 Configuración Avanzada

#### Ajuste de Rendimiento
```bash
# Configuración de alto rendimiento (conexión estable requerida)
python ncbi_downloader.py --batch-size 500 --workers 6

# Configuración conservadora (conexión inestable)
python ncbi_downloader.py --batch-size 50 --workers 2
```

### 🔍 Monitoreo y Solución de Problemas

#### Verificar Progreso
```bash
python check_progress.py
```

#### Ver Logs
```bash
# Actividad reciente
tail -50 output/ncbi_download.log

# Seguir en tiempo real
tail -f output/ncbi_download.log

# Filtrar errores
grep ERROR output/ncbi_download.log
```

### 🤝 Contribuir

1. Hacer fork del repositorio
2. Crear una rama de características: `git checkout -b nombre-caracteristica`
3. Realizar cambios con pruebas
4. Enviar un pull request

### 📝 Licencia

Este proyecto está licenciado bajo la Licencia MIT - ver el archivo [LICENSE](LICENSE) para detalles.

---

**⭐ Si encuentras útil este proyecto, ¡considera darle una estrella en GitHub!**
