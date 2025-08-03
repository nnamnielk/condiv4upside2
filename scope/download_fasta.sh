#!/bin/bash

# Download ASTRAL SCOPe 2.08 non-redundant FASTA (40% identity cutoff)
echo "Downloading ASTRAL SCOPe 2.08 FASTA..."
curl -k -o scope/astral-40.fa "https://scop.berkeley.edu/downloads/scopeseq-2.08/astral-scopedom-seqres-gd-sel-gs-bib-40-2.08.fa"

# Download SCOP hierarchy file for class mapping
echo "Downloading SCOP hierarchy file..."
curl -k -o scope/dir.hie.scope.2.08-stable.txt "https://scop.berkeley.edu/downloads/parse/dir.hie.scope.2.08-stable.txt"

echo "Downloads completed."
