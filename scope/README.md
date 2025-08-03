# SCOP Domain Selection Pipeline

This pipeline filters SCOP protein domains by structural class and sequence length, then uses a greedy algorithm to select a diverse set with minimal sequence similarity.

## Files

- `download_fasta.sh`: Downloads ASTRAL SCOP 2.08 FASTA and hierarchy files
- `select_domains.py`: Main Python script for filtering and selection
- `run_selection.sh`: Wrapper script that runs the complete pipeline

## Usage

### Quick start (default parameters)
```bash
./scope/run_selection.sh
```

### Custom parameters
```bash
./scope/run_selection.sh --classes A B --maxlen 50 --n 20
```

### Parameters
- `--classes`: SCOP classes to include (A=all-alpha, B=all-beta, A/B=alpha/beta, A+B=alpha+beta)
- `--maxlen`: Maximum sequence length (default: 70)
- `--n`: Number of sequences to select (default: 30)
- `--output`: Output FASTA file (default: scope/selected_domains.fa)
- `--info`: Output info file (default: scope/selection_info.tsv)

## Algorithm

1. **Download**: Fetches ASTRAL SCOP 2.08 non-redundant sequences (<40% identity)
2. **Filter**: Selects domains by structural class and sequence length
3. **Greedy Selection**: Iteratively adds sequences with minimal maximum identity to the current set

## Output

- `selected_domains.fa`: FASTA file with selected sequences
- `selection_info.tsv`: Tab-separated file with selection order, domain IDs, classes, lengths, and identity scores

## Example Results

From 15,177 total sequences, 731 match filtering criteria (<70 residues, all classes):
- Class a (all-alpha): 388 domains
- Class b (all-beta): 192 domains  
- Class c (alpha/beta): 7 domains
- Class d (alpha+beta): 144 domains

The greedy algorithm successfully minimizes sequence similarity while maintaining structural diversity across SCOP classes.
