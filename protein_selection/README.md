# Protein Selection for Contrastive Divergence Training

## 🎯 Quick Start

**Your training dataset is ready!**
- **387 all-alpha proteins** optimized for small size and diversity
- **Files**: `datasets/alpha/final_387_proteins.fa` and `metadata.tsv`
- **PDB structures**: `datasets/alpha/pdbs/` (209 files)

## 📁 Directory Structure

```
protein_selection/
├── README.md                    # This file
├── scripts/                     # Core workflow scripts
│   ├── select_alpha_proteins.py # Main alpha selection
│   ├── direct_pdb_search.py     # PDB database search  
│   └── download_pdbs.py          # Download PDB files
├── datasets/                    # Final protein datasets
│   ├── alpha/                   
│   │   ├── final_387_proteins.fa # 🎯 YOUR TRAINING SET
│   │   ├── metadata.tsv          # Protein info & statistics
│   │   ├── enhanced_proteins.fa   # Alternative dataset
│   │   └── pdbs/                 # PDB structure files
│   └── beta/                    # Beta protein datasets
│       ├── enhanced_proteins.fa
│       └── metadata.tsv
├── tools/                       # Analysis utilities
│   └── analyze_dataset.py        # Generate statistics/histograms
└── archive/                     # Archived intermediate results
```

## 🚀 Usage

### Analyze Your Dataset
```bash
cd protein_selection
python tools/analyze_dataset.py
```

### Select New Alpha Proteins
```bash
python scripts/select_alpha_proteins.py
```

### Select Beta Proteins
```bash
cd scripts
python ../scripts/select_alpha_proteins.py  # Modify for beta
```

### Download Missing PDB Files
```bash
python scripts/download_pdbs.py --alpha-info datasets/alpha/metadata.tsv --output-dir datasets/alpha/pdbs --alpha-only
```

## 📊 Dataset Quality

**Alpha Training Set (387 proteins)**:
- **Length range**: 10-299 residues
- **Average similarity**: ~12% (excellent diversity)
- **Size focus**: Optimized for smallest proteins
- **Quality**: No chain breaks, single-chain only

## 🎯 Ready for Training!

Your dataset is optimized and ready for contrastive divergence experiments. All files are organized and documented for easy use.
