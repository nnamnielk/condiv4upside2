#!/usr/bin/env python3
"""
Build a large alpha-protein training set (170–500) with no chain breaks.

- Source PDBs: selection/pdbs/alpha/*.pdb (already present locally)
- Selection info: selection/enhanced_alpha_selection_info.tsv (curated, low similarity)
- Outputs (per protein) in condiv_training_proteins/:
    <PDB>.initial.npy  (CA coordinates, N x 3, float32, Angstrom)
    <PDB>.fasta
    <PDB>.chi          (placeholder/dummy)
- Also writes condiv_training_proteins/protein_list.txt with required "prot" header.

Env overrides:
  TARGET_COUNT            default 500    (max proteins to prepare)
  MAX_PROTEIN_LENGTH      default 300    (skip longer ones)
  MIN_PROTEIN_LENGTH      default 40     (skip very small ones)
  CHAIN_BREAK_THRESHOLD   default 6.0    (Angstrom; max allowed gap between consecutive CA)
"""

import os
import sys
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser

SRC_PDB_DIR = "selection/pdbs/alpha"
SELECTION_TSV = "selection/enhanced_alpha_selection_info.tsv"
OUT_DIR = "condiv_training_proteins"

def has_no_chain_break(ca_coords: np.ndarray, threshold: float = 6.0) -> bool:
    """
    Return True if no chain break is detected based on consecutive CA-CA distances.
    ca_coords: (N, 3) array in Angstrom
    threshold: maximum allowed gap distance for a continuous chain
    """
    if ca_coords.ndim != 2 or ca_coords.shape[1] != 3 or ca_coords.shape[0] < 2:
        return False
    diffs = np.diff(ca_coords, axis=0)
    dists = np.sqrt((diffs ** 2).sum(axis=1))
    return float(dists.max(initial=0.0)) <= threshold

def extract_ca_coords(pdb_file: str) -> np.ndarray:
    """
    Extract CA coordinates from a PDB file as an (N,3) float32 array in Angstrom.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)
    ca = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_id()[0] == ' ':  # standard residue
                    if residue.has_id('CA'):
                        ca.append(residue['CA'].get_coord())
    if not ca:
        raise ValueError(f"No CA atoms found in {pdb_file}")
    return np.asarray(ca, dtype=np.float32)

def create_fasta_from_pdb(pdb_file: str) -> str:
    """
    Create FASTA sequence string from PDB residues (standard residues only).
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)
    aa_codes = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
        'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G',
        'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
        'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
        'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }
    seq = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_id()[0] == ' ':
                    res = residue.get_resname()
                    seq.append(aa_codes.get(res, 'X'))
    if not seq:
        raise ValueError(f"No valid residues for FASTA in {pdb_file}")
    return ''.join(seq)

def write_fasta(path: str, pdb_id: str, seq: str) -> None:
    with open(path, 'w') as f:
        f.write(f">{pdb_id}\n")
        f.write(seq + "\n")

def write_dummy_chi(path: str, pdb_id: str) -> None:
    with open(path, 'w') as f:
        f.write(f"# Chi angles for {pdb_id}\n")
        f.write("# Placeholder file\n")

def main() -> int:
    target_count = int(os.environ.get("TARGET_COUNT", "500"))
    max_len = int(os.environ.get("MAX_PROTEIN_LENGTH", "300"))
    min_len = int(os.environ.get("MIN_PROTEIN_LENGTH", "40"))
    chain_break_thr = float(os.environ.get("CHAIN_BREAK_THRESHOLD", "6.0"))

    if not os.path.isdir(SRC_PDB_DIR):
        print(f"Error: source PDB directory not found: {SRC_PDB_DIR}", file=sys.stderr)
        return 1
    if not os.path.isfile(SELECTION_TSV):
        print(f"Error: selection TSV not found: {SELECTION_TSV}", file=sys.stderr)
        return 1

    os.makedirs(OUT_DIR, exist_ok=True)

    # Load selection table and build list of unique PDB IDs sorted by length
    df = pd.read_csv(SELECTION_TSV, sep="\t")
    if "PDB_Chain" not in df.columns or "Length" not in df.columns:
        print("Error: selection TSV missing required columns PDB_Chain and/or Length", file=sys.stderr)
        return 1
    df = df.sort_values("Length")
    # Filter by length window
    df = df[(df["Length"] >= min_len) & (df["Length"] <= max_len)]

    unique_ids = []
    seen = set()
    for pdb_chain in df["PDB_Chain"]:
        pdb_id = str(pdb_chain).split("_")[0]
        if pdb_id not in seen:
            unique_ids.append(pdb_id)
            seen.add(pdb_id)

    print(f"Candidates after length filter [{min_len}, {max_len}]: {len(unique_ids)}")

    prepared = []
    for pdb_id in unique_ids:
        if len(prepared) >= target_count:
            break
        pdb_path = os.path.join(SRC_PDB_DIR, f"{pdb_id}.pdb")
        if not os.path.isfile(pdb_path):
            # silently skip missing PDB files
            continue
        try:
            ca = extract_ca_coords(pdb_path)
            if not has_no_chain_break(ca, threshold=chain_break_thr):
                # enforce no chain breaks
                continue

            base = os.path.join(OUT_DIR, pdb_id)
            # save .initial.npy
            np.save(base + ".initial.npy", ca.astype(np.float32))
            # save .fasta
            seq = create_fasta_from_pdb(pdb_path)
            write_fasta(base + ".fasta", pdb_id, seq)
            # save .chi placeholder
            write_dummy_chi(base + ".chi", pdb_id)

            prepared.append(pdb_id)
            if len(prepared) % 25 == 0:
                print(f"Prepared {len(prepared)} proteins...")

        except Exception as e:
            # skip problematic entries
            continue

    if not prepared:
        print("No proteins prepared. Consider relaxing thresholds or checking source data.", file=sys.stderr)
        return 2

    # Write protein_list.txt with required header
    list_path = os.path.join(OUT_DIR, "protein_list.txt")
    with open(list_path, "w") as f:
        f.write("prot\n")
        for pdb_id in prepared:
            f.write(f"{pdb_id}\n")

    print("==============================================")
    print("SETUP COMPLETE")
    print(f"Total prepared: {len(prepared)} (target={target_count})")
    print(f"Output dir: {OUT_DIR}")
    print(f"Protein list: {list_path}")
    print("Controls:")
    print(f"  CHAIN_BREAK_THRESHOLD={chain_break_thr} Å  MIN_LEN={min_len}  MAX_LEN={max_len}")
    print("Re-run initialization:")
    print("  sbatch init.sbatch")
    print("Then run training:")
    print("  sbatch srun.sbatch")
    return 0

if __name__ == "__main__":
    sys.exit(main())
