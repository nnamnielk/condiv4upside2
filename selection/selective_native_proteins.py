#!/usr/bin/env python3

import argparse
import sys
import json
import requests
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import pairwise2
import numpy as np
import os
import time

def get_protein_secondary_structure_info(pdb_id):
    """
    Get secondary structure information from RCSB PDB API.
    Returns classification as 'alpha', 'beta', 'mixed', or 'other'.
    """
    try:
        # Query RCSB PDB API for secondary structure information
        api_url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
        response = requests.get(api_url, timeout=10)
        
        if response.status_code == 200:
            data = response.json()
            
            # Check for SCOP classification if available
            if 'rcsb_entry_info' in data:
                entry_info = data['rcsb_entry_info']
                
                # Look for structure classification
                if 'structure_determination_methodology' in entry_info:
                    method = entry_info['structure_determination_methodology']
                    if method and 'X-RAY DIFFRACTION' not in method:
                        return 'other'  # Skip non-X-ray structures
            
            # Try to get polymer entity information for better classification
            polymer_url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/1"
            polymer_response = requests.get(polymer_url, timeout=10)
            
            if polymer_response.status_code == 200:
                polymer_data = polymer_response.json()
                
                # Check for secondary structure content
                if 'rcsb_polymer_entity_feature_summary' in polymer_data:
                    features = polymer_data['rcsb_polymer_entity_feature_summary']
                    
                    alpha_content = 0
                    beta_content = 0
                    
                    for feature in features:
                        if feature.get('type') == 'HELIX':
                            alpha_content += feature.get('count', 0)
                        elif feature.get('type') == 'SHEET':
                            beta_content += feature.get('count', 0)
                    
                    # Classify based on secondary structure content
                    total_ss = alpha_content + beta_content
                    if total_ss > 0:
                        alpha_ratio = alpha_content / total_ss
                        beta_ratio = beta_content / total_ss
                        
                        if alpha_ratio > 0.7:
                            return 'alpha'
                        elif beta_ratio > 0.7:
                            return 'beta'
                        else:
                            return 'mixed'
            
            # Fallback: use known classifications for common small proteins
            known_alpha = {
                '1BDD', '1CRN', '1PGB', '1VII', '1UBQ', '1SHG', '1PIN', '1WRP', 
                '1ENH', '1CTF', '1BBA', '1MJC', '1ERC', '1APS', '1HTM', '1VCC',
                '1TOP', '1AVP', '1PEF', '1DSG', '1OSO', '1AGL', '1CQB', '2TRX',
                '1ROP', '1ZIF', '1FSD', '1SN3', '1FIS', '1RGS', '1P68', '1POU',
                '1GBS', '1EYE', '1OVA'
            }
            
            known_beta = {
                '1FE1', '1TEN', '1SHF', '1BF4', '1QRE', '1WIT', '1GB1', '1FNF',
                '1CLB', '1TIT', '1CSP', '1BPI', '1PTQ', '1SRL', '1TUL', '1BEE',
                '1CEY', '1DFN', '1FBV', '1GDO', '1HFH', '1IGD', '1JMX', '1KNG',
                '1L2Y', '1MFI', '1NLS', '1ORC', '1PDA', '1QDD', '1RIS', '1SJT',
                '1TKV', '1UCS', '1VFB', '1WGE', '1XNB', '1YCS', '1ZAA'
            }
            
            if pdb_id.upper() in known_alpha:
                return 'alpha'
            elif pdb_id.upper() in known_beta:
                return 'beta'
                
        return 'other'
        
    except Exception as e:
        print(f"Error getting secondary structure for {pdb_id}: {e}")
        return 'other'

def query_proteins_by_type(protein_type='alpha', max_length=100, max_entries=1000):
    """
    Query RCSB PDB for protein structures of specific type (alpha or beta).
    """
    chains = []
    
    print(f"Querying PDB for {protein_type} protein structures (<{max_length} residues)")
    
    # Use a broader search to get high-quality structures
    search_url = "https://search.rcsb.org/rcsbsearch/v2/query"
    search_query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "exptl.method",
                        "value": "X-RAY DIFFRACTION"
                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entry_info.resolution_combined",
                        "operator": "less_or_equal",
                        "value": 3.0
                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "entity_poly.rcsb_entity_polymer_type",
                        "value": "Protein"
                    }
                }
            ]
        },
        "return_type": "entry",
        "request_options": {
            "results_verbosity": "minimal",
            "return_all_hits": True
        }
    }
    
    try:
        response = requests.post(search_url, json=search_query)
        response.raise_for_status()
        result = response.json()
        pdb_ids = result.get("result_set", [])
        print(f"Found {len(pdb_ids)} high-resolution protein structures")
        
        # Process entries
        processed = 0
        type_found = 0
        
        for pdb_id in pdb_ids:
            if len(chains) >= max_entries:
                break
                
            try:
                # Get secondary structure classification
                ss_type = get_protein_secondary_structure_info(pdb_id)
                
                if ss_type != protein_type:
                    processed += 1
                    continue
                
                # Get FASTA sequence
                fasta_url = f"https://www.rcsb.org/fasta/entry/{pdb_id}"
                fasta_response = requests.get(fasta_url, timeout=15)
                
                if fasta_response.status_code == 200:
                    fasta_content = fasta_response.text.strip()
                    if fasta_content and not fasta_content.startswith('<!DOCTYPE'):
                        from io import StringIO
                        fasta_io = StringIO(fasta_content)
                        
                        try:
                            for record in SeqIO.parse(fasta_io, 'fasta'):
                                seq_len = len(record.seq)
                                if 20 <= seq_len <= max_length:  # At least 20 residues
                                    chain_id = record.id.split('_')[-1] if '_' in record.id else 'A'
                                    record.id = f"{pdb_id}_{chain_id}"
                                    record.description = f"PDB {pdb_id} chain {chain_id} [{protein_type}] [length={seq_len}]"
                                    chains.append(record)
                                    type_found += 1
                                    
                                    if len(chains) >= max_entries:
                                        break
                        except Exception as parse_error:
                            print(f"Parse error for {pdb_id}: {parse_error}")
                            continue
                
                processed += 1
                if processed % 50 == 0:
                    print(f"Processed {processed} entries, found {type_found} {protein_type} proteins, collected {len(chains)} chains")
                
                # Add small delay to avoid overwhelming the API
                time.sleep(0.1)
                    
            except Exception as e:
                print(f"Error processing {pdb_id}: {e}")
                processed += 1
                continue
                
    except Exception as e:
        print(f"Error querying PDB: {e}")
        
        # Fallback to known proteins of the specified type
        if protein_type == 'alpha':
            known_proteins = [
                '1BDD', '1CRN', '1PGB', '1VII', '1UBQ', '1SHG', '1PIN', '1WRP', 
                '1ENH', '1CTF', '1BBA', '1MJC', '1ERC', '1APS', '1HTM', '1VCC',
                '1TOP', '1AVP', '1PEF', '1DSG', '1OSO', '1AGL', '1CQB', '2TRX',
                '1ROP', '1ZIF', '1FSD', '1SN3', '1FIS', '1RGS', '1P68', '1POU'
            ]
        else:  # beta
            known_proteins = [
                '1FE1', '1TEN', '1SHF', '1BF4', '1QRE', '1WIT', '1GB1', '1FNF',
                '1CLB', '1TIT', '1CSP', '1BPI', '1PTQ', '1SRL', '1TUL', '1BEE',
                '1CEY', '1DFN', '1FBV', '1GDO', '1HFH', '1IGD', '1JMX', '1KNG'
            ]
        
        print(f"Falling back to known {protein_type} proteins: {len(known_proteins)} proteins")
        
        for pdb_id in known_proteins:
            if len(chains) >= max_entries:
                break
                
            try:
                fasta_url = f"https://www.rcsb.org/fasta/entry/{pdb_id}"
                fasta_response = requests.get(fasta_url, timeout=10)
                
                if fasta_response.status_code == 200:
                    fasta_content = fasta_response.text.strip()
                    if fasta_content:
                        from io import StringIO
                        fasta_io = StringIO(fasta_content)
                        
                        for record in SeqIO.parse(fasta_io, 'fasta'):
                            seq_len = len(record.seq)
                            if seq_len <= max_length:
                                chain_id = record.id.split('_')[-1] if '_' in record.id else 'A'
                                record.id = f"{pdb_id}_{chain_id}"
                                record.description = f"PDB {pdb_id} chain {chain_id} [{protein_type}] [length={seq_len}]"
                                chains.append(record)
            except:
                continue
    
    # Sort by length to get smallest proteins first
    chains.sort(key=lambda x: len(x.seq))
    
    print(f"Found {len(chains)} {protein_type} protein chains with length <= {max_length}")
    return chains

def calculate_sequence_identity(seq1, seq2):
    """Calculate sequence identity between two sequences."""
    if len(seq1) == 0 or len(seq2) == 0:
        return 0.0
    
    alignments = pairwise2.align.globalxx(seq1, seq2, one_alignment_only=True)
    if not alignments:
        return 0.0
    
    alignment = alignments[0]
    seq1_aligned = alignment.seqA
    seq2_aligned = alignment.seqB
    
    matches = sum(1 for a, b in zip(seq1_aligned, seq2_aligned) if a == b and a != '-' and b != '-')
    aligned_length = sum(1 for a, b in zip(seq1_aligned, seq2_aligned) if a != '-' or b != '-')
    
    return matches / aligned_length if aligned_length > 0 else 0.0

def greedy_selection(sequences, n):
    """
    Greedy algorithm to select n sequences with minimal pairwise similarity.
    """
    n_select = min(n, len(sequences))
    selected_indices = []
    selection_info = []
    
    if len(sequences) == 0:
        return [], []
    
    # Start with shortest sequence
    selected_indices.append(0)
    selection_info.append((0, 0.0))
    print(f"Selected sequence 1: {sequences[0].id} (length: {len(sequences[0].seq)})")
    
    if len(sequences) == 1 or n_select == 1:
        return selected_indices, selection_info
    
    # Greedy selection for remaining sequences
    for step in range(1, n_select):
        best_idx = -1
        min_max_identity = 1.0
        
        for i, seq in enumerate(sequences):
            if i in selected_indices:
                continue
            
            max_identity = 0.0
            for selected_idx in selected_indices:
                identity = calculate_sequence_identity(str(seq.seq), str(sequences[selected_idx].seq))
                max_identity = max(max_identity, identity)
            
            if max_identity < min_max_identity:
                min_max_identity = max_identity
                best_idx = i
        
        if best_idx != -1:
            selected_indices.append(best_idx)
            selection_info.append((best_idx, min_max_identity))
            print(f"Selected sequence {step+1}: {sequences[best_idx].id} (length: {len(sequences[best_idx].seq)}, max identity: {min_max_identity:.3f})")
        else:
            break
    
    return selected_indices, selection_info

def main():
    parser = argparse.ArgumentParser(description='Select smallest diverse alpha or beta protein structures')
    parser.add_argument('--type', choices=['alpha', 'beta'], required=True,
                       help='Protein type to select (alpha or beta)')
    parser.add_argument('--maxlen', type=int, default=100, help='Maximum sequence length')
    parser.add_argument('--n', type=int, default=500, help='Number of sequences to select')
    parser.add_argument('--output', help='Output FASTA file (default: {type}_proteins.fa)')
    parser.add_argument('--info', help='Output info file (default: {type}_selection_info.tsv)')
    
    args = parser.parse_args()
    
    # Set default output filenames
    if not args.output:
        args.output = f"{args.type}_proteins.fa"
    if not args.info:
        args.info = f"{args.type}_selection_info.tsv"
    
    print(f"Searching for {args.type} proteins")
    print(f"Maximum length: {args.maxlen} residues")
    print(f"Target selection: {args.n} sequences")
    
    # Query PDB for protein chains of specified type
    sequences = query_proteins_by_type(args.type, args.maxlen, max_entries=2000)
    
    if len(sequences) == 0:
        print(f"No {args.type} sequences found matching criteria!")
        sys.exit(1)
    
    print(f"Found {len(sequences)} {args.type} protein chains matching criteria")
    
    # Apply greedy selection
    print(f"\nSelecting {min(args.n, len(sequences))} diverse sequences using greedy algorithm...")
    selected_indices, selection_info = greedy_selection(sequences, args.n)
    
    # Write selected sequences to FASTA
    selected_sequences = [sequences[i] for i in selected_indices]
    with open(args.output, 'w') as f:
        SeqIO.write(selected_sequences, f, 'fasta')
    
    # Write selection info
    with open(args.info, 'w') as f:
        f.write("Order\tPDB_Chain\tType\tLength\tMax_Identity_to_Set\n")
        for order, (idx, max_identity) in enumerate(selection_info):
            seq = sequences[idx]
            f.write(f"{order+1}\t{seq.id}\t{args.type}\t{len(seq.seq)}\t{max_identity:.4f}\n")
    
    print(f"\nSelected {len(selected_sequences)} {args.type} sequences saved to: {args.output}")
    print(f"Selection info saved to: {args.info}")
    
    # Print summary statistics
    lengths = [len(seq.seq) for seq in selected_sequences]
    print(f"\nLength statistics:")
    print(f"  Minimum: {min(lengths)} residues")
    print(f"  Maximum: {max(lengths)} residues")
    print(f"  Average: {np.mean(lengths):.1f} residues")
    print(f"  Median: {np.median(lengths):.1f} residues")

if __name__ == '__main__':
    main()
