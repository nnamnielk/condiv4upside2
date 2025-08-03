#!/bin/bash

# Wrapper script to download SCOP data and select diverse protein domains

set -e  # Exit on any error

# Default parameters
CLASSES="A B A/B A+B"
MAXLEN=70
N=30
OUTPUT="scope/selected_domains.fa"
INFO="scope/selection_info.tsv"

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --classes)
            CLASSES="$2"
            shift 2
            ;;
        --maxlen)
            MAXLEN="$2"
            shift 2
            ;;
        --n)
            N="$2"
            shift 2
            ;;
        --output)
            OUTPUT="$2"
            shift 2
            ;;
        --info)
            INFO="$2"
            shift 2
            ;;
        --help|-h)
            echo "Usage: $0 [options]"
            echo "Options:"
            echo "  --classes    SCOP classes to include (default: 'A B A/B A+B')"
            echo "  --maxlen     Maximum sequence length (default: 70)"
            echo "  --n          Number of sequences to select (default: 30)"
            echo "  --output     Output FASTA file (default: scope/selected_domains.fa)"
            echo "  --info       Output info file (default: scope/selection_info.tsv)"
            echo "  --help, -h   Show this help message"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

echo "=== SCOP Domain Selection Pipeline ==="
echo "Classes: $CLASSES"
echo "Max length: $MAXLEN"
echo "Number to select: $N"
echo "Output FASTA: $OUTPUT"
echo "Output info: $INFO"
echo

# Step 1: Download data if not present
if [ ! -f "scope/astral-40.fa" ] || [ ! -f "scope/dir.hie.scope.2.08-stable.txt" ]; then
    echo "=== Step 1: Downloading SCOP data ==="
    chmod +x scope/download_fasta.sh
    ./scope/download_fasta.sh
    echo
else
    echo "=== Step 1: SCOP data already present ==="
    echo
fi

# Step 2: Run domain selection
echo "=== Step 2: Selecting diverse domains ==="
python3 scope/select_domains.py \
    --fasta scope/astral-40.fa \
    --hie scope/dir.hie.scope.2.08-stable.txt \
    --classes $CLASSES \
    --maxlen $MAXLEN \
    --n $N \
    --output $OUTPUT \
    --info $INFO

echo
echo "=== Pipeline completed successfully ==="
echo "Selected domains: $OUTPUT"
echo "Selection info: $INFO"
