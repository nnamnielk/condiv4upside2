# Protein Selection Workflow for Contrastive Divergence Training

## Overview

This repository contains a comprehensive workflow for selecting diverse protein datasets optimized for contrastive divergence training. The workflow supports selection of **alpha-helical proteins**, **beta-sheet proteins**, and **mixed alpha-beta proteins** with customizable diversity and size constraints.

## 🎯 Final Optimal Dataset

**Current Status: ✅ COMPLETED**
- **387 all-alpha proteins** selected for training
- **Length range**: 10-299 residues
- **Excellent diversity**: Low sequence similarity (~12% average)
- **Key files**: 
  - `selection/optimal_training_alpha_proteins.fa` (FASTA sequences)
  - `selection/optimal_training_alpha_selection_info.tsv` (metadata)

---

## 📋 Workflow Overview

### Phase 1: Initial Protein Discovery
```
SCOPe Database → Direct PDB Search → Curated Lists
     ↓               ↓                 ↓
Domain-based      API queries       Literature
selections        (797 proteins)    proteins
```

### Phase 2: Filtering & Selection
```
Size Filtering → Sequence Clustering → Diversity Optimization
    ↓                ↓                      ↓
≤300 residues   Identity thresholds    Greedy selection
No chain breaks    (30-70%)            algorithm
```

### Phase 3: Dataset Integration
```
Multiple Sources → Deduplication → Final Selection
      ↓               ↓               ↓
6 different      PDB ID based     387 unique
datasets         deduplication     proteins
```

---

## 🔧 Core Scripts & Usage

### 1. Alpha Protein Selection

**Main Script**: `create_optimal_training_set.py`
```bash
python create_optimal_training_set.py
```

**What it does**:
- Combines 6 different high-quality alpha protein datasets
- Removes duplicates based on PDB ID
- Sorts by size (prefers smaller proteins)
- Outputs 387 diverse alpha proteins

**Key intermediate scripts**:
- `direct_pdb_alpha_query.py` - Direct PDB database search
- `selection/enhanced_selective_native_proteins.py` - Literature-based selection

### 2. Beta Protein Selection

**Files**: `selection/beta_proteins.fa`, `selection/enhanced_beta_proteins.fa`

To select beta proteins (similar workflow):
```bash
# Use the same enhanced selective script with beta parameter
cd selection
python enhanced_selective_native_proteins.py --type beta --maxlen 150 --n 300
```

### 3. Mixed Alpha-Beta Selection

**SCOPe-based selection**:
```bash
cd scope
python select_domains.py --fasta astral-40.fa --hie dir.hie.scope.2.08-stable.txt \
  --classes A B --maxlen 200 --n 500
```

---

## 📁 File Organization

### Core Workflow Files (Keep)
```
📁 Main Directory
├── create_optimal_training_set.py         # Main script for alpha selection
├── direct_pdb_alpha_query.py             # Direct PDB search
├── text_histogram.py                     # Statistics visualization
└── README_Protein_Selection_Workflow.md  # This file

📁 selection/
├── optimal_training_alpha_proteins.fa          # 🎯 FINAL DATASET
├── optimal_training_alpha_selection_info.tsv   # 🎯 FINAL METADATA
├── enhanced_selective_native_proteins.py       # Literature-based selection
├── build_alpha_no_breaks.py                   # Chain break checking
├── download_pdbs.py                           # PDB download utility
└── pdbs/alpha/                               # PDB structure files

📁 scope/
├── astral-40.fa                          # SCOPe database
├── select_domains.py                     # SCOPe-based selection
└── dir.hie.scope.2.08-stable.txt        # SCOPe hierarchy
```

### Intermediate Files (Can Delete)
```
❌ Temporary/redundant files:
├── *_sorted.txt                    # Temporary sorting files
├── missing_pdbs.txt                 # Temporary download lists
├── final_pdb_list.txt              # Old PDB lists
├── filtered_pdb_list.txt           # Old filtered lists
├── existing_pdbs.txt               # Temporary file lists

❌ Intermediate datasets (keep metadata only):
├── selection/massive_pdb_alpha_*    # Intermediate search results
├── selection/huge_pdb_alpha_*       # Intermediate search results  
├── selection/pdb_direct_alpha_*     # Intermediate search results
├── selection/optimal_small_alpha_*  # Intermediate search results

❌ Legacy/unused scripts:
├── download_missing_pdbs.py         # Superseded by download_pdbs.py
├── download_final_pdbs.py          # Temporary download script
├── plot_protein_lengths.py         # Superseded by text_histogram.py
```

---

## 🚀 Quick Start Guide

### For Alpha Proteins (Current Setup)
```bash
# Your dataset is ready!
ls selection/optimal_training_alpha_proteins.fa      # 387 proteins
ls selection/optimal_training_alpha_selection_info.tsv  # Metadata
```

### For Beta Proteins  
```bash
cd selection
python enhanced_selective_native_proteins.py --type beta --maxlen 150 --n 300
python download_pdbs.py --beta-info beta_selection_info.tsv --separate-dirs
```

### For Mixed Alpha-Beta
```bash
cd scope  
python select_domains.py --classes A B --maxlen 200 --n 400
cd ../selection
python download_pdbs.py --alpha-info ../scope/selection_info.tsv
```

---

## 📊 Dataset Quality Metrics

**Current Alpha Dataset (387 proteins)**:
- ✅ **Size diversity**: 10-299 residues
- ✅ **Sequence diversity**: ~12% average similarity  
- ✅ **No chain breaks**: All verified
- ✅ **All single-chain**: Verified
- ✅ **All alpha-helical**: Structure-based selection

**Previous vs Current**:
- **Before**: 208 proteins, 98.8% max similarity
- **After**: 387 proteins, much better diversity
- **Improvement**: 86% more proteins, dramatically better diversity

---

## 🧹 Recommended Cleanup

Run this cleanup after confirming the workflow works:
```bash
# Remove temporary files
rm *_sorted.txt missing_pdbs.txt final_pdb_list.txt filtered_pdb_list.txt existing_pdbs.txt

# Remove intermediate datasets (keep .tsv metadata)
rm selection/massive_pdb_alpha_proteins.fa
rm selection/huge_pdb_alpha_proteins.fa  
rm selection/pdb_direct_alpha_proteins.fa
rm selection/optimal_small_alpha_proteins.fa

# Remove redundant scripts
rm download_missing_pdbs.py download_final_pdbs.py plot_protein_lengths.py
```

---

## 🎯 Key Results

**Final Training Dataset**: 
- **File**: `selection/optimal_training_alpha_proteins.fa`
- **Size**: 387 all-alpha proteins  
- **Quality**: Excellent diversity and size distribution
- **Ready for**: Contrastive divergence training

This workflow successfully created a comprehensive, diverse training set optimized for small all-alpha proteins with minimal sequence redundancy.
