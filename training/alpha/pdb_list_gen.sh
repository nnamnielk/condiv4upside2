#!/bin/bash

if [[ "$1" == "-h" ]] || [[ "$1" == "--help" ]]; then
    echo "Usage: $0 [DIRECTORY]"
    echo "Extracts unique PDB names (before '.') from the specified directory"
    echo "DIRECTORY: Directory to search (default: PDB_train)"
    echo "Pipe the output to a file if needed: $0 [DIRECTORY] > output.txt"
    exit 0
fi

# Use provided directory or default to PDB_train
DIR="${1:-PDB_train}"

ls "$DIR" | cut -d'.' -f1 | sort -u
