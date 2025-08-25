#!/bin/bash

# Download all PDB files for selected alpha and beta proteins
# This script downloads PDB files into organized directories

echo "=========================================="
echo "PDB Download Script for Selected Proteins"
echo "=========================================="

# Set default parameters
OUTPUT_DIR="pdb_files"
MAX_WORKERS=5
FORMAT="pdb"

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --max-workers)
            MAX_WORKERS="$2"
            shift 2
            ;;
        --format)
            FORMAT="$2"
            shift 2
            ;;
        --compressed)
            COMPRESSED="--compressed"
            shift
            ;;
        --help)
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --output-dir DIR     Output directory (default: pdb_files)"
            echo "  --max-workers N      Number of parallel downloads (default: 5)"
            echo "  --format FORMAT      File format: pdb or cif (default: pdb)"
            echo "  --compressed         Download compressed files (.gz)"
            echo "  --help               Show this help message"
            echo ""
            echo "This script will download PDB files for all selected alpha and beta proteins"
            echo "into separate subdirectories: alpha/ and beta/"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

echo "Configuration:"
echo "  Output directory: $OUTPUT_DIR"
echo "  Max workers: $MAX_WORKERS"
echo "  Format: $FORMAT"
echo "  Compressed: ${COMPRESSED:-no}"
echo ""

# Check if Python script exists
if [ ! -f "download_pdbs.py" ]; then
    echo "Error: download_pdbs.py not found in current directory"
    exit 1
fi

# Check if selection info files exist
if [ ! -f "enhanced_alpha_selection_info.tsv" ]; then
    echo "Warning: enhanced_alpha_selection_info.tsv not found"
fi

if [ ! -f "enhanced_beta_selection_info.tsv" ]; then
    echo "Warning: enhanced_beta_selection_info.tsv not found"
fi

# Run the download script
echo "Starting PDB downloads..."
python download_pdbs.py \
    --alpha-info enhanced_alpha_selection_info.tsv \
    --beta-info enhanced_beta_selection_info.tsv \
    --output-dir "$OUTPUT_DIR" \
    --format "$FORMAT" \
    --max-workers "$MAX_WORKERS" \
    --separate-dirs \
    $COMPRESSED

# Check results
if [ $? -eq 0 ]; then
    echo ""
    echo "=========================================="
    echo "Download completed successfully!"
    echo "=========================================="
    
    if [ -d "$OUTPUT_DIR/alpha" ]; then
        ALPHA_COUNT=$(ls -1 "$OUTPUT_DIR/alpha" | wc -l)
        echo "Alpha proteins downloaded: $ALPHA_COUNT files"
    fi
    
    if [ -d "$OUTPUT_DIR/beta" ]; then
        BETA_COUNT=$(ls -1 "$OUTPUT_DIR/beta" | wc -l)
        echo "Beta proteins downloaded: $BETA_COUNT files"
    fi
    
    echo "Files saved to: $OUTPUT_DIR"
    
    # Show directory structure
    echo ""
    echo "Directory structure:"
    tree "$OUTPUT_DIR" 2>/dev/null || find "$OUTPUT_DIR" -type d | head -10
    
else
    echo ""
    echo "=========================================="
    echo "Download failed with errors"
    echo "=========================================="
    exit 1
fi
