#!/bin/bash

# Condiv training launch script
# Generated automatically

echo "Starting condiv training setup..."

# Set environment variables
export UPSIDE_HOME="/project2/trsosnic/upside2"
export PYTHONPATH="$UPSIDE_HOME/py:$PYTHONPATH"

# Check if UPSIDE_HOME exists
if [ ! -d "$UPSIDE_HOME" ]; then
    echo "Error: UPSIDE_HOME directory not found: $UPSIDE_HOME"
    echo "Please set the correct path to your Upside installation"
    exit 1
fi

# Parameters
INIT_DIR="init_param"
PROTEIN_DIR="condiv_training_proteins"
PROTEIN_LIST="condiv_training_proteins/protein_list.txt"
OUTPUT_DIR="condiv_training_results"

echo "Configuration:"
echo "  Init parameters: $INIT_DIR"
echo "  Protein directory: $PROTEIN_DIR"
echo "  Protein list: $PROTEIN_LIST"
echo "  Output directory: $OUTPUT_DIR"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Initialize condiv training
echo "Initializing condiv training..."
python condiv2.py initialize "$INIT_DIR" "$PROTEIN_DIR" "$PROTEIN_LIST" "$OUTPUT_DIR"

if [ $? -eq 0 ]; then
    echo "Initialization successful!"
    echo ""
    echo "To start training, run:"
    echo "  python condiv2.py restart $OUTPUT_DIR/initial_checkpoint.pkl 100"
    echo ""
    echo "Training files are in: $OUTPUT_DIR"
else
    echo "Initialization failed!"
    exit 1
fi
