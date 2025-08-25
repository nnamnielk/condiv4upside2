# Selective Native Protein Selection Results

## Overview
This document summarizes the results of using the enhanced selective native protein selection scripts to find the smallest possible alpha and beta proteins from the Protein Data Bank (PDB).

## Scripts Created

### 1. `selective_native_proteins.py`
- Original script that queries PDB API and falls back to known protein lists
- Classifies proteins by secondary structure (alpha vs beta)
- Uses greedy selection algorithm to minimize sequence similarity

### 2. `enhanced_selective_native_proteins.py` (Recommended)
- Enhanced version with extensive lists of known alpha and beta proteins
- More robust protein classification
- Better coverage of small proteins from literature and databases

## Results Achieved

### Beta Proteins
- **File**: `enhanced_beta_proteins.fa`
- **Info**: `enhanced_beta_selection_info.tsv`
- **Count**: 118 diverse beta proteins
- **Length Range**: 10-149 residues
- **Average Length**: 71.8 residues
- **Median Length**: 74.0 residues

### Alpha Proteins
- **File**: `enhanced_alpha_proteins.fa`
- **Info**: `enhanced_alpha_selection_info.tsv`
- **Count**: 173 diverse alpha proteins
- **Length Range**: 10-150 residues
- **Average Length**: 76.5 residues
- **Median Length**: 74.0 residues

## Usage Instructions

### Basic Usage
```bash
# For beta proteins
python enhanced_selective_native_proteins.py --type beta --maxlen 150 --n 500

# For alpha proteins
python enhanced_selective_native_proteins.py --type alpha --maxlen 150 --n 500
```

### Advanced Usage
```bash
# Custom output files and parameters
python enhanced_selective_native_proteins.py \
    --type beta \
    --maxlen 200 \
    --n 1000 \
    --output my_beta_proteins.fa \
    --info my_beta_info.tsv
```

### Parameters
- `--type`: Choose 'alpha' or 'beta' protein type
- `--maxlen`: Maximum sequence length (default: 100)
- `--n`: Number of sequences to select (default: 500)
- `--output`: Output FASTA file name (optional)
- `--info`: Output info TSV file name (optional)

## Selection Algorithm

The script uses a greedy selection algorithm that:

1. **Sorts proteins by length** (smallest first)
2. **Starts with the shortest protein**
3. **Iteratively selects proteins** that have the lowest maximum sequence identity to the already selected set
4. **Continues until** the target number is reached or no more suitable proteins are found

## Output Files

### FASTA Files
- Contain the selected protein sequences
- Headers include PDB ID, chain, type, and length information
- Format: `>PDB_CHAIN PDB {ID} chain {CHAIN} [{TYPE}] [length={LENGTH}]`

### Info Files (TSV)
- Tab-separated values with selection statistics
- Columns: Order, PDB_Chain, Type, Length, Max_Identity_to_Set
- Shows the diversity metrics for each selected protein

## Limitations and Notes

1. **PDB API Issues**: The PDB search API sometimes returns errors, so the script falls back to curated lists of known proteins
2. **Protein Classification**: Uses a combination of PDB API data and manually curated lists for alpha/beta classification
3. **Sequence Quality**: Some sequences may contain 'X' residues (unknown/modified amino acids)
4. **Coverage**: While we aimed for 500 of each type, we achieved 118 beta and 173 alpha proteins due to availability constraints

## Quality Control

- All proteins have experimentally determined structures (X-ray crystallography preferred)
- Sequences are filtered for minimum length (10+ residues)
- Diversity is ensured through the greedy selection algorithm
- Results are sorted by length to prioritize smallest proteins

## Future Improvements

To get closer to 500 proteins of each type:
1. Expand the known protein lists with more recent PDB entries
2. Implement alternative PDB query methods
3. Include NMR structures in addition to X-ray structures
4. Consider using SCOP/CATH databases for better classification

## PDB File Downloads

### Automatic PDB Download Scripts

Two scripts are provided to automatically download PDB structure files for all selected proteins:

#### 1. `download_pdbs.py` (Python Script)
Advanced Python script with full control over download parameters:

```bash
# Download all proteins (alpha and beta) into separate directories
python download_pdbs.py --separate-dirs --output-dir pdb_structures

# Download only alpha proteins
python download_pdbs.py --alpha-only --output-dir alpha_pdbs

# Download only beta proteins  
python download_pdbs.py --beta-only --output-dir beta_pdbs

# Download compressed files with more parallel workers
python download_pdbs.py --compressed --max-workers 10 --format cif
```

**Parameters:**
- `--alpha-info`: Alpha proteins selection info file (default: enhanced_alpha_selection_info.tsv)
- `--beta-info`: Beta proteins selection info file (default: enhanced_beta_selection_info.tsv)
- `--output-dir`: Output directory for PDB files (default: pdb_files)
- `--format`: File format - 'pdb' or 'cif' (default: pdb)
- `--compressed`: Download compressed files (.gz)
- `--max-workers`: Number of parallel downloads (default: 5)
- `--separate-dirs`: Create separate alpha/ and beta/ subdirectories
- `--alpha-only`: Download only alpha proteins
- `--beta-only`: Download only beta proteins

#### 2. `download_all_pdbs.sh` (Shell Script)
Simple wrapper script for easy one-command downloads:

```bash
# Basic usage - downloads all proteins into pdb_files/alpha/ and pdb_files/beta/
./download_all_pdbs.sh

# Custom output directory
./download_all_pdbs.sh --output-dir my_structures

# Download compressed CIF files with 10 parallel workers
./download_all_pdbs.sh --format cif --compressed --max-workers 10
```

### Download Results

The scripts will download PDB structure files for:
- **138 unique alpha protein structures** (173 chains total)
- **95 unique beta protein structures** (118 chains total)

### Features

- **Parallel downloads** for faster processing
- **Automatic retry** handling for failed downloads
- **Progress tracking** with detailed status messages
- **Duplicate detection** - skips already downloaded files
- **Error logging** - saves failed download IDs for retry
- **Multiple formats** - supports both PDB and mmCIF formats
- **Compression support** - optional gzip compression

## Files Generated

### Protein Sequences
- `enhanced_beta_proteins.fa` - 118 smallest diverse beta proteins
- `enhanced_beta_selection_info.tsv` - Selection statistics for beta proteins
- `enhanced_alpha_proteins.fa` - 173 smallest diverse alpha proteins  
- `enhanced_alpha_selection_info.tsv` - Selection statistics for alpha proteins

### PDB Download Scripts
- `download_pdbs.py` - Python script for downloading PDB files
- `download_all_pdbs.sh` - Shell wrapper script for easy downloads
- `pdb_files/` - Default directory for downloaded PDB structures (created when scripts run)
