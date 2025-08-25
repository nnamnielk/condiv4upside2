# Diverse Native Protein Selection

## Overview
This pipeline selects diverse native protein structures with experimentally determined coordinates, suitable for training with contrastive divergence.

## Results
- **18 diverse protein chains** selected from known small proteins (<70 residues)
- **Greedy algorithm** minimizes pairwise sequence similarity
- All proteins have **native state PDB structures** available

## Key Features
- Native state structures only (no domains/subsets)
- Length filtered: 6-69 residues
- Diverse sequence space coverage
- Maximum pairwise identity ranges from 0.00 to 0.75

## Files
- `selected_proteins.fa`: Selected protein sequences in FASTA format
- `selection_info.tsv`: Selection order and diversity metrics
- `select_native_proteins.py`: Main selection script

## Usage
```bash
python select_native_proteins.py --n 30 --maxlen 70
```

## Selected Proteins Summary
| Order | PDB Chain | Length | Max Identity | Description |
|-------|-----------|--------|--------------|-------------|
| 1     | 1BDD      | 60     | 0.0000      | Hirudin variant-2 |
| 2     | 1FE1      | 40     | 0.0000      | Iron-sulfur cluster |
| 3     | 1AGL      | 6      | 0.0312      | Hexapeptide |
| 4     | 1AVP      | 11     | 0.1333      | Vasopressin |
| 5     | 1CRN      | 46     | 0.1648      | Crambin |
| ...   | ...       | ...    | ...         | ... |

Note: Some sequences contain 'X' residues (unknown/modified) and should be filtered for training purposes.
