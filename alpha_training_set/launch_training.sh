#!/bin/bash

# Condiv training script for 0 smallest alpha proteins
# Generated automatically

echo "Starting condiv training for 0 alpha proteins"
echo "Training directory: alpha_training_set"

# List of proteins
PROTEINS=()

echo "Proteins to train:"
for protein in "${PROTEINS[@]}"; do
    echo "  $protein"
done

# Check if .up files exist
echo ""
echo "Checking .up files..."
missing_up=0
for protein in "${PROTEINS[@]}"; do
    up_file="alpha_training_set/$protein.up"
    if [ -f "$up_file" ]; then
        echo "  ✓ $protein.up"
    else
        echo "  ✗ $protein.up (missing)"
        missing_up=$((missing_up + 1))
    fi
done

if [ $missing_up -gt 0 ]; then
    echo ""
    echo "Error: $missing_up .up files are missing!"
    echo "Please run the preparation script first."
    exit 1
fi

echo ""
echo "All .up files found. Ready to start training!"
echo ""

# Launch condiv training
echo "Launching condiv training..."
python condiv2.py \
    --input-dir alpha_training_set \
    --proteins  \
    --output-dir results_alpha_training \
    --epochs 1000 \
    --batch-size 32

echo "Training completed!"
