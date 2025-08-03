#!/usr/bin/env python3

import argparse
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import pairwise2
import numpy as np

def parse_classification_from_fasta(fasta_file):
    """Parse SCOP classifications directly from FASTA headers."""
    domain_to_class = {}
    
    for record in SeqIO.parse(fasta_file, 'fasta'):
        header = record.description
        domain_id = record.id
        
        # Extract classification from header (e.g., "d1dlwa_ a.1.1.1 (A:) ...")
        parts = header.split()
        if len(parts) >= 2:
            classification = parts[1]  # e.g., "a.1.1.1"
            if '.' in classification:
                class_code = classification.split('.')[0]  # Extract 'a' from 'a.1.1.1'
                domain_to_class[domain_id] = class_code
    
    return domain_to_class

def calculate_sequence_identity(seq1, seq2):
    """Calculate sequence identity between two sequences."""
    if len(seq1) == 0 or len(seq2) == 0:
        return 0.0
    
    # Use simple character matching for speed
    min_len = min(len(seq1), len(seq2))
    max_len = max(len(seq1), len(seq2))
    
    # Align sequences using global alignment
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
    Returns indices of selected sequences and their selection order info.
    """
    if len(sequences) <= n:
        return list(range(len(sequences))), []
    
    selected_indices = []
    selection_info = []
    
    # Start with first sequence
    selected_indices.append(0)
    selection_info.append((0, 0.0))  # (index, max_identity_to_selected_set)
    
    print(f"Selected sequence 1: {sequences[0].id}")
    
    # Greedy selection for remaining n-1 sequences
    for step in range(1, n):
        best_idx = -1
        min_max_identity = 1.0
        
        # For each remaining sequence, find max identity to selected set
        for i, seq in enumerate(sequences):
            if i in selected_indices:
                continue
            
            max_identity = 0.0
            for selected_idx in selected_indices:
                identity = calculate_sequence_identity(str(seq.seq), str(sequences[selected_idx].seq))
                max_identity = max(max_identity, identity)
            
            # Choose sequence with minimum max_identity
            if max_identity < min_max_identity:
                min_max_identity = max_identity
                best_idx = i
        
        if best_idx != -1:
            selected_indices.append(best_idx)
            selection_info.append((best_idx, min_max_identity))
            print(f"Selected sequence {step+1}: {sequences[best_idx].id} (max identity to set: {min_max_identity:.3f})")
    
    return selected_indices, selection_info

def main():
    parser = argparse.ArgumentParser(description='Select diverse protein domains from SCOP')
    parser.add_argument('--fasta', required=True, help='Path to ASTRAL FASTA file')
    parser.add_argument('--hie', required=True, help='Path to SCOP hierarchy file')
    parser.add_argument('--classes', nargs='+', choices=['A', 'B', 'A/B', 'A+B'], 
                       default=['A', 'B', 'A/B', 'A+B'], help='SCOP classes to include')
    parser.add_argument('--maxlen', type=int, default=70, help='Maximum sequence length')
    parser.add_argument('--n', type=int, default=30, help='Number of sequences to select')
    parser.add_argument('--output', default='selected_domains.fa', help='Output FASTA file')
    parser.add_argument('--info', default='selection_info.tsv', help='Output info file')
    
    args = parser.parse_args()
    
    # Map user classes to SCOP class codes
    class_mapping = {'A': 'a', 'B': 'b', 'A/B': 'c', 'A+B': 'd'}
    target_classes = [class_mapping[cls] for cls in args.classes]
    
    print(f"Parsing classifications from FASTA file: {args.fasta}")
    domain_to_class = parse_classification_from_fasta(args.fasta)
    print(f"Loaded {len(domain_to_class)} domain classifications")
    
    print(f"Reading FASTA file for filtering: {args.fasta}")
    sequences = []
    total_sequences = 0
    
    for record in SeqIO.parse(args.fasta, 'fasta'):
        total_sequences += 1
        domain_id = record.id
        
        # Check if domain has classification
        if domain_id not in domain_to_class:
            continue
        
        domain_class = domain_to_class[domain_id]
        
        # Filter by class and length
        if domain_class in target_classes and len(record.seq) < args.maxlen:
            sequences.append(record)
    
    print(f"Found {len(sequences)} domains matching criteria from {total_sequences} total sequences")
    print(f"Class distribution:")
    class_counts = {}
    for seq in sequences:
        cls = domain_to_class[seq.id]
        class_counts[cls] = class_counts.get(cls, 0) + 1
    for cls, count in sorted(class_counts.items()):
        print(f"  {cls}: {count}")
    
    if len(sequences) == 0:
        print("No sequences match the filtering criteria!")
        sys.exit(1)
    
    print(f"\nSelecting {min(args.n, len(sequences))} diverse sequences using greedy algorithm...")
    selected_indices, selection_info = greedy_selection(sequences, args.n)
    
    # Write selected sequences to FASTA
    selected_sequences = [sequences[i] for i in selected_indices]
    with open(args.output, 'w') as f:
        SeqIO.write(selected_sequences, f, 'fasta')
    
    # Write selection info
    with open(args.info, 'w') as f:
        f.write("Order\tDomain_ID\tClass\tLength\tMax_Identity_to_Set\n")
        for order, (idx, max_identity) in enumerate(selection_info):
            seq = sequences[idx]
            domain_class = domain_to_class[seq.id]
            f.write(f"{order+1}\t{seq.id}\t{domain_class}\t{len(seq.seq)}\t{max_identity:.4f}\n")
    
    print(f"\nSelected {len(selected_sequences)} sequences saved to: {args.output}")
    print(f"Selection info saved to: {args.info}")

if __name__ == '__main__':
    main()
