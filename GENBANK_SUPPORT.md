# GenBank Format Support

## Overview

The NCBI downloader now supports downloading sequences in **GenBank format** in addition to FASTA format. GenBank format provides richer annotation data including features, references, and detailed metadata.

## Features Added

### 1. GenBank Download Support
- **New parser**: `_parse_genbank_batch_enhanced()` in `ncbi_downloader.py`
- Automatically detects format based on `ncbi.rettype` configuration
- Preserves GenBank annotations and features
- Adds custom metadata (protein_type, download timestamp)

### 2. GenBank Integration Support
- Data integrator now recognizes GenBank files (`.gb`, `.gbk`, `.genbank`)
- Extracts metadata from GenBank annotations and features
- Supports both FASTA and GenBank in the same integration workflow

### 3. Enhanced Metadata Extraction
- **From GenBank records**:
  - Protein ID from record.id
  - Organism from annotations
  - Protein type from CDS/Protein features or custom annotations
  - Variants/strain from source features

## Configuration

### To Use GenBank Format:

Edit your config file (e.g., `config_gb.json`):

```json
{
  "email": "your.email@example.com",
  "files": {
    "output_file": "all_sequences.gb",
    "...": "..."
  },
  "ncbi": {
    "database": "protein",
    "rettype": "gb",
    "retmode": "text",
    "api_base": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
  }
}
```

**Supported `rettype` values**:
- `"fasta"` - FASTA format (default)
- `"gb"` - GenBank format
- `"genbank"` - GenBank format (alternative)
- `"gp"` - GenPept format (protein GenBank)

### To Use FASTA Format (default):

```json
{
  "files": {
    "output_file": "all_sequences.fasta"
  },
  "ncbi": {
    "rettype": "fasta",
    "retmode": "text"
  }
}
```

## Usage Examples

### Download in GenBank Format

```bash
# Using GenBank config
python ncbi_downloader.py --config config_gb.json --input ncbi_ids_dict.tsv

# With integration (will handle GenBank files)
python ncbi_downloader.py --config config_gb.json \
  --generate-ids --integrate \
  --virus-name "Influenza A virus" \
  --target-proteins Neuraminidase Hemagglutinin
```

### Download in FASTA Format

```bash
# Using default config (FASTA)
python ncbi_downloader.py --input ncbi_ids_dict.tsv

# Or specify FASTA config
python ncbi_downloader.py --config config.json --input ncbi_ids_dict.tsv
```

### Manual Integration with Mixed Formats

```bash
# Integration works with both FASTA and GenBank files
python data_integrator.py
```

The integrator automatically detects file format based on extension and processes accordingly.

## File Extensions

| Format | Extensions |
|--------|-----------|
| FASTA | `.fasta`, `.fa`, `.fas`, `.fna`, `.ffn`, `.faa`, `.frn` |
| GenBank | `.gb`, `.gbk`, `.genbank` |
| Compressed | Add `.gz` to any extension (e.g., `.fasta.gz`, `.gb.gz`) |

## Differences Between Formats

### FASTA Format
- **Pros**:
  - Smaller file size
  - Faster to download and parse
  - Standard format for sequence analysis tools
  - Sufficient for most sequence-based analyses

- **Cons**:
  - Limited metadata (only ID, description, sequence)
  - No annotation information

### GenBank Format
- **Pros**:
  - Rich annotation data (features, CDS, genes)
  - Complete metadata (organism, references, taxonomy)
  - Strain and isolate information
  - Publication references
  - Feature qualifiers (protein_id, product, etc.)

- **Cons**:
  - Larger file size (3-5x larger than FASTA)
  - Slower to download and parse
  - More complex structure

## When to Use Each Format

### Use FASTA When:
- You only need sequence data
- File size is a concern
- Speed is important
- Running sequence alignment or similarity searches
- Building phylogenetic trees (sequence-only)

### Use GenBank When:
- You need detailed annotations
- Extracting gene features or CDS regions
- Analyzing protein domains or motifs
- Need complete taxonomic information
- Require strain/isolate metadata
- Want publication references

## Implementation Details

### Parser Selection Logic

```python
# In ncbi_downloader.py
rettype = self.config.get('ncbi.rettype', 'fasta')
if rettype in ['gb', 'genbank', 'gp']:
    sequences = self._parse_genbank_batch_enhanced(response, protein_type)
else:
    sequences = self._parse_fasta_batch_enhanced(response, protein_type)
```

### Format Detection in Integrator

```python
# In data_integrator.py
file_format = "fasta"
if fasta_file.suffix.lower() in ['.gb', '.gbk', '.genbank']:
    file_format = "genbank"

# Parse with appropriate format
for record in SeqIO.parse(handle, file_format):
    # Process record...
```

### Metadata Extraction Hierarchy

**For GenBank files**:
1. Custom annotations (added during download)
2. Record annotations (organism, source)
3. Feature qualifiers (product, protein_id, strain)
4. Description parsing (fallback)

**For FASTA files**:
1. Enhanced headers (if present)
2. Description pattern matching
3. Keyword-based inference

## Logging

When using GenBank format, you'll see:

```
INFO - Configured output format: GB
INFO - GenBank format selected - will download full GenBank records with annotations
INFO - Parsing GenBank file: all_sequences.gb
INFO - Parsed GenBank record: XP_123456 (1500 bp)
```

When using FASTA format:

```
INFO - Configured output format: FASTA
INFO - FASTA format selected - will download sequences only
INFO - Parsing FASTA file: all_sequences.fasta
```

## Troubleshooting

### Issue: Empty GenBank files
**Cause**: Parser is trying to read GenBank data as FASTA
**Solution**: Ensure `ncbi.rettype` in config matches file extension

### Issue: Integration fails with GenBank files
**Cause**: File extension not recognized
**Solution**: Use `.gb`, `.gbk`, or `.genbank` extension

### Issue: Missing protein types in GenBank data
**Cause**: GenBank records don't have protein_type in features
**Solution**: Parser will try to infer from description as fallback

## Testing

Test GenBank download:
```bash
# Create test config
cat > test_gb.json << EOF
{
  "email": "test@example.com",
  "batch_size": 50,
  "files": {
    "output_file": "test_sequences.gb"
  },
  "ncbi": {
    "rettype": "gb"
  }
}
EOF

# Download small test set
python ncbi_downloader.py --config test_gb.json \
  --generate-ids --max-ids-per-protein 10 \
  --virus-name "Influenza A virus" \
  --target-proteins Neuraminidase
```

## Compatibility

- ✅ Works with `--integrate` flag
- ✅ Supports compressed files (`.gb.gz`)
- ✅ Compatible with target protein filtering
- ✅ Works with manual integration via `data_integrator.py`
- ✅ Unified output includes both FASTA and GenBank sources

## Performance Considerations

| Aspect | FASTA | GenBank | Impact |
|--------|-------|---------|--------|
| File size | 1x | 3-5x | Disk space |
| Download time | 1x | 2-3x | Network bandwidth |
| Parse time | 1x | 1.5-2x | Processing time |
| Memory usage | 1x | 1.5-2x | RAM requirement |

For large-scale downloads (>10,000 sequences), consider:
- Using smaller batch sizes for GenBank
- Enabling compression
- Processing in chunks
- Using FASTA unless annotations are specifically needed
