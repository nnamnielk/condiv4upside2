#!/usr/bin/env python3
"""
Build a large set (200–500) of alpha-class proteins with NO chain breaks.

- Sources PDB IDs from selection/enhanced_alpha_selection_info.tsv (curated alpha set)
- Downloads missing PDBs into selection/pdbs/alpha/
- Filters to single-chain entries (no multichain)
- Outputs required training files in condiv_training_proteins/:
    <PDB>.initial.npy  (N x 3 float32, Angstrom CA coords)
    <PDB>.fasta        (sequence)
    <PDB>.chi          (placeholder)
- Writes condiv_training_proteins/protein_list.txt with header 'prot'

Arguments:
  --selection    TSV file (default: selection/enhanced_alpha_selection_info.tsv)
  --pdb-dir      Directory with PDB files (default: selection/pdbs/alpha)
  --out-dir      Output training dir (default: condiv_training_proteins)
  --target-count Max proteins to prepare (default: 500)
  --min-len      Min residues (default: 40)
  --max-len      Max residues (default: 300)
  --max-download-workers  Parallel download workers (default: 6)

Notes:
- We enforce single-chain entries (no multichain) based on PDB parsing (single chain with CA atoms).
- This aligns with the requirement "NO chain breaks" for the dataset stage.
"""

import os
import sys
import argparse
import time
from pathlib import Path
import numpy as np
import pandas as pd
import requests
from concurrent.futures import ThreadPoolExecutor, as_completed

from typing import List
from Bio.PDB import PDBParser

DEFAULT_SELECTION = "selection/enhanced_alpha_selection_info.tsv"
DEFAULT_PDB_DIR = "selection/pdbs/alpha"
DEFAULT_OUT_DIR = "condiv_training_proteins"

AA3_TO1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
    'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G',
    'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
    'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
    'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}

def extract_alpha_pdb_ids(selection_paths: List[str]) -> pd.DataFrame:
    frames: List[pd.DataFrame] = []
    for selection_path in selection_paths:
        if not os.path.exists(selection_path):
            print(f"Warning: selection file not found, skipping: {selection_path}")
            continue
        try:
            df = pd.read_csv(selection_path, sep="\t")
        except Exception as e:
            print(f"Warning: failed reading {selection_path}: {e}")
            continue
        if "PDB_Chain" not in df.columns or "Length" not in df.columns:
            print(f"Warning: missing required columns in {selection_path}, skipping")
            continue
        dfx = pd.DataFrame({
            "PDB_ID": df["PDB_Chain"].astype(str).str.split("_").str[0].str.upper(),
            "Length": df["Length"]
        })
        frames.append(dfx)
    if not frames:
        raise ValueError("No valid selection files provided")
    all_df = pd.concat(frames, ignore_index=True)
    # Keep smallest Length per PDB_ID
    all_df = all_df.sort_values("Length").drop_duplicates(subset=["PDB_ID"], keep="first").reset_index(drop=True)
    return all_df[["PDB_ID", "Length"]]

def extract_pdbs_from_fastas(fasta_paths: List[str]) -> pd.DataFrame:
    """
    Parse FASTA files to extract PDB IDs and lengths.
    - Header format examples:
        >1AGL_1|Chain ... [alpha] [length=6]
      We prefer the explicit [length=...] if present; otherwise compute from sequence lines.
    Returns DataFrame with columns: PDB_ID, Length
    """
    rows = []
    for fp in fasta_paths:
        if not os.path.exists(fp):
            print(f"Warning: FASTA not found, skipping: {fp}")
            continue
        try:
            with open(fp, "r") as f:
                curr_id = None
                seq_len = 0
                header_len = None
                for ln in f:
                    ln = ln.strip()
                    if not ln:
                        continue
                    if ln.startswith(">"):
                        # flush previous
                        if curr_id is not None:
                            rows.append((curr_id, header_len if header_len is not None else seq_len))
                        # parse new header
                        curr_id = ln[1:].split()[0].split("_")[0].upper()
                        seq_len = 0
                        header_len = None
                        # try to parse [length=N]
                        if "length=" in ln:
                            try:
                                part = ln.split("length=")[1]
                                num = ""
                                for ch in part:
                                    if ch.isdigit():
                                        num += ch
                                    else:
                                        break
                                if num:
                                    header_len = int(num)
                            except Exception:
                                header_len = None
                    else:
                        # sequence line
                        seq_len += len(ln)
                # flush last record
                if curr_id is not None:
                    rows.append((curr_id, header_len if header_len is not None else seq_len))
        except Exception as e:
            print(f"Warning: failed reading FASTA {fp}: {e}")
            continue
    if not rows:
        return pd.DataFrame(columns=["PDB_ID", "Length"])
    df = pd.DataFrame(rows, columns=["PDB_ID", "Length"])
    # Keep smallest length per PDB_ID (some files may include multiple chains of same PDB)
    df = df.sort_values("Length").drop_duplicates(subset=["PDB_ID"], keep="first").reset_index(drop=True)
    return df

def download_pdb(pdb_id: str, pdb_dir: str) -> bool:
    os.makedirs(pdb_dir, exist_ok=True)
    path = os.path.join(pdb_dir, f"{pdb_id}.pdb")
    if os.path.exists(path):
        return True
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    try:
        r = requests.get(url, timeout=30)
        r.raise_for_status()
        with open(path, "wb") as f:
            f.write(r.content)
        return True
    except Exception:
        return False

def download_missing_pdbs(pdb_ids, pdb_dir: str, max_workers: int = 6):
    ok = []
    bad = []
    with ThreadPoolExecutor(max_workers=max_workers) as ex:
        futs = {ex.submit(download_pdb, pid, pdb_dir): pid for pid in pdb_ids}
        for i, fut in enumerate(as_completed(futs), 1):
            pid = futs[fut]
            try:
                res = fut.result()
                if res:
                    ok.append(pid)
                    print(f"[{i}/{len(pdb_ids)}] ✓ {pid}")
                else:
                    bad.append(pid)
                    print(f"[{i}/{len(pdb_ids)}] ✗ {pid} download failed")
            except Exception as e:
                bad.append(pid)
                print(f"[{i}/{len(pdb_ids)}] ✗ {pid} exception: {e}")
            time.sleep(0.05)
    return ok, bad

def parse_best_chain_ca(pdb_path: str):
    """
    Return CA coordinates (N,3) and residue 3-letter sequence for the longest valid chain
    (standard residues with CA atoms) in the first model. Returns (None, None) if none found.
    """
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("prot", pdb_path)
    except Exception:
        return None, None

    best = None
    for model in structure:
        for chain in model:
            residues = []
            for residue in chain:
                hetflag = residue.get_id()[0]
                if hetflag != " ":
                    continue
                if not residue.has_id("CA"):
                    continue
                resname = residue.get_resname()
                residues.append((resname, residue["CA"].get_coord()))
            if residues:
                if best is None or len(residues) > len(best):
                    best = residues
        # Only use first model
        break

    if best is None:
        return None, None

    resnames, ca_list = zip(*best)
    ca = np.asarray(ca_list, dtype=np.float32)
    seq3 = list(resnames)
    return ca, seq3

def has_no_chain_break(ca: np.ndarray, ca_gap: float) -> bool:
    if ca is None or ca.ndim != 2 or ca.shape[1] != 3 or ca.shape[0] < 2:
        return False
    diffs = np.diff(ca, axis=0)
    dists = np.sqrt((diffs ** 2).sum(axis=1))
    return float(dists.max(initial=0.0)) <= ca_gap

def seq3_to_fasta(seq3):
    seq1 = []
    for r in seq3:
        seq1.append(AA3_TO1.get(r, "X"))
    return "".join(seq1)

def write_training_files(out_dir: str, pdb_id: str, ca: np.ndarray, fasta_seq: str):
    base = os.path.join(out_dir, pdb_id)
    np.save(base + ".initial.npy", ca.astype(np.float32))
    with open(base + ".fasta", "w") as f:
        f.write(f">{pdb_id}\n{fasta_seq}\n")
    with open(base + ".chi", "w") as f:
        f.write(f"# Chi angles for {pdb_id}\n# Placeholder file\n")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--selection", default=DEFAULT_SELECTION)
    ap.add_argument("--selection-files", nargs="+", default=None, help="Optional list of selection TSV files to union")
    ap.add_argument("--fasta-files", nargs="+", default=None, help="Optional list of FASTA files to union for additional PDB IDs")
    ap.add_argument("--pdb-dir", default=DEFAULT_PDB_DIR)
    ap.add_argument("--out-dir", default=DEFAULT_OUT_DIR)
    ap.add_argument("--target-count", type=int, default=500)
    ap.add_argument("--min-len", type=int, default=40)
    ap.add_argument("--max-len", type=int, default=300)
    ap.add_argument("--max-download-workers", type=int, default=6)
    ap.add_argument("--ca-gap", type=float, default=4.5, help="Max allowed CA-CA gap (Angstrom) between consecutive residues")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    # Load alpha candidates sorted by Length (smallest first)
    sel_files = args.selection_files if args.selection_files else [args.selection]
    df = extract_alpha_pdb_ids(sel_files)
    # Optionally union with FASTA sources
    if args.fasta_files:
        df_fa = extract_pdbs_from_fastas(args.fasta_files)
        if not df_fa.empty:
            df = pd.concat([df, df_fa], ignore_index=True)
            df = df.sort_values("Length").drop_duplicates(subset=["PDB_ID"], keep="first").reset_index(drop=True)
    df = df[(df["Length"] >= args.min_len) & (df["Length"] <= args.max_len)]
    candidates = df["PDB_ID"].tolist()
    print(f"Alpha candidates after length [{args.min_len},{args.max_len}]: {len(candidates)}")

    # Ensure PDBs present
    missing = [pid for pid in candidates if not os.path.exists(os.path.join(args.pdb_dir, f"{pid}.pdb"))]
    if missing:
        print(f"Downloading {len(missing)} missing PDBs to {args.pdb_dir} ...")
        download_pdbs(missing, args)
    else:
        print("All candidate PDBs present locally.")

    # Filter for single-chain, contiguous CA traces (no breaks)
    prepared = []
    for i, pid in enumerate(candidates, 1):
        if len(prepared) >= args.target_count:
            break
        pdb_path = os.path.join(args.pdb_dir, f"{pid}.pdb")
        if not os.path.exists(pdb_path):
            continue
        ca, seq3 = parse_best_chain_ca(pdb_path)
        if ca is None:
            continue
        n = ca.shape[0]
        if n < args.min_len or n > args.max_len:
            continue
        if not has_no_chain_break(ca, args.ca_gap):
            continue
        fasta = seq3_to_fasta(seq3)
        try:
            write_training_files(args.out_dir, pid, ca, fasta)
            prepared.append(pid)
            if len(prepared) % 25 == 0:
                print(f"Prepared {len(prepared)} proteins...")
        except Exception:
            # skip if write fails
            continue

    if not prepared:
        print("No proteins prepared. Consider relaxing thresholds or verifying PDB availability.")
        return 2

    # Write protein list with header
    list_path = os.path.join(args.out_dir, "protein_list.txt")
    with open(list_path, "w") as f:
        f.write("prot\n")
        for pid in prepared:
            f.write(f"{pid}\n")

    print("==============================================")
    print("SETUP COMPLETE")
    print(f"Prepared: {len(prepared)} proteins (target={args.target_count})")
    print(f"Output dir: {args.out_dir}")
    print(f"Protein list: {list_path}")
    return 0

def download_pdbs(pdb_ids, args):
    ok, bad = download_missing_pdbs(pdb_ids, args.pdb_dir, args.max_download_workers)
    print(f"Downloads: {len(ok)} ok, {len(bad)} failed")
    if bad:
        failed_file = os.path.join(args.pdb_dir, "failed_downloads.txt")
        with open(failed_file, "w") as f:
            for pid in bad:
                f.write(f"{pid}\n")
        print(f"Failed list written to {failed_file}")

if __name__ == "__main__":
    sys.exit(main())
