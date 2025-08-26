#!/usr/bin/env python3

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def read_selection_info(filename):
    """Read selection info file and return as list of records."""
    records = []
    with open(filename, 'r') as f:
        lines = f.readlines()
        header = lines[0].strip().split('\t') if lines else []
        
        for line in lines[1:]:  # Skip header
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                order = parts[0]
                pdb_chain = parts[1] 
                prot_type = parts[2]
                length = int(parts[3])
                
                # Check if this file has similarity/identity column
                if len(parts) > 4 and 'Source' not in header:
                    try:
                        similarity = float(parts[4])
                    except ValueError:
                        similarity = 0.0  # Default if parsing fails
                else:
                    similarity = 0.0  # No similarity column (e.g., PDB direct)
                
                records.append({
                    'order': order,
                    'pdb_chain': pdb_chain, 
                    'type': prot_type,
                    'length': length,
                    'similarity': similarity,
                    'source': filename
                })
    return records

def combine_datasets():
    """Combine all available alpha protein datasets into optimal training set."""
    all_records = []
    
    # Priority 1: Direct PDB results (highest quality)
    try:
        pdb_records = read_selection_info('selection/pdb_direct_alpha_selection_info.tsv')
        print(f"Loaded {len(pdb_records)} direct PDB proteins")
        all_records.extend(pdb_records)
    except FileNotFoundError:
        print("Direct PDB results not found")
    
    # Priority 2: Optimal small results  
    try:
        optimal_records = read_selection_info('selection/optimal_small_alpha_selection_info.tsv')
        print(f"Loaded {len(optimal_records)} optimal small proteins")
        all_records.extend(optimal_records)
    except FileNotFoundError:
        print("Optimal small results not found")
    
    # Priority 3: Filtered results (add ALL filtered proteins)
    try:
        filtered_records = read_selection_info('selection/filtered_alpha_selection_info.tsv')
        # Take ALL filtered proteins for maximum coverage
        print(f"Loaded {len(filtered_records)} filtered proteins (all available)")
        all_records.extend(filtered_records)
    except FileNotFoundError:
        print("Filtered results not found")
    
    # Priority 4: Enhanced results (add medium proteins for more diversity)
    try:
        enhanced_records = read_selection_info('selection/enhanced_alpha_selection_info.tsv')
        # Take ALL enhanced proteins for maximum coverage
        print(f"Loaded {len(enhanced_records)} enhanced proteins (all available)")
        all_records.extend(enhanced_records)
    except FileNotFoundError:
        print("Enhanced results not found")
    
    # Priority 5: Massive PDB search results (latest and largest)
    try:
        massive_records = read_selection_info('selection/massive_pdb_alpha_selection_info.tsv')
        print(f"Loaded {len(massive_records)} massive PDB search proteins")
        all_records.extend(massive_records)
    except FileNotFoundError:
        print("Massive PDB search results not found")
    
    # Priority 6: HUGE PDB search results (ultimate and largest!)
    try:
        huge_records = read_selection_info('selection/huge_pdb_alpha_selection_info.tsv')
        print(f"Loaded {len(huge_records)} HUGE PDB search proteins")
        all_records.extend(huge_records)
    except FileNotFoundError:
        print("Huge PDB search results not found")
    
    # Remove duplicates based on PDB chain ID
    seen_chains = set()
    unique_records = []
    
    for record in all_records:
        pdb_id = record['pdb_chain'].split('_')[0].split('|')[0]
        if pdb_id not in seen_chains:
            seen_chains.add(pdb_id)
            unique_records.append(record)
    
    print(f"After deduplication: {len(unique_records)} unique proteins")
    
    # Sort by length (prefer smaller proteins) then by similarity
    unique_records.sort(key=lambda x: (x['length'], x['similarity']))
    
    # Keep all unique proteins for maximum training set size
    print(f"Keeping all {len(unique_records)} unique proteins for maximum training set")
    
    return unique_records

def create_combined_fasta(records, output_file):
    """Create combined FASTA file from records."""
    sequences = []
    
    # Try to get sequences from existing FASTA files
    fasta_sources = [
        'selection/pdb_direct_alpha_proteins.fa',
        'selection/optimal_small_alpha_proteins.fa', 
        'selection/enhanced_alpha_proteins.fa',
        'selection/massive_pdb_alpha_proteins.fa',
        'selection/huge_pdb_alpha_proteins.fa'
    ]
    
    # Load all available sequences
    all_seqs = {}
    for fasta_file in fasta_sources:
        try:
            for record in SeqIO.parse(fasta_file, 'fasta'):
                pdb_id = record.id.split('_')[0] 
                all_seqs[pdb_id] = record
        except FileNotFoundError:
            continue
    
    # Create sequences for selected records
    missing_seqs = []
    for i, record in enumerate(records, 1):
        pdb_id = record['pdb_chain'].split('_')[0].split('|')[0]
        
        if pdb_id in all_seqs:
            seq_record = all_seqs[pdb_id]
            # Update ID and description for consistency
            seq_record.id = f"{pdb_id}_{i}"
            seq_record.description = f"PDB {pdb_id} [alpha] [length={record['length']}] [similarity={record['similarity']:.3f}]"
            sequences.append(seq_record)
        else:
            missing_seqs.append(pdb_id)
    
    if missing_seqs:
        print(f"Warning: Missing sequences for {len(missing_seqs)} proteins: {missing_seqs[:10]}...")
    
    # Write sequences to FASTA file
    if sequences:
        with open(output_file, 'w') as f:
            SeqIO.write(sequences, f, 'fasta')
        print(f"Wrote {len(sequences)} sequences to {output_file}")
    
    return sequences

def main():
    print("Creating optimal alpha protein training set...")
    
    # Combine datasets
    records = combine_datasets()
    
    if len(records) == 0:
        print("No records found!")
        return
    
    # Create combined FASTA
    sequences = create_combined_fasta(records, 'selection/optimal_training_alpha_proteins.fa')
    
    # Write selection info
    with open('selection/optimal_training_alpha_selection_info.tsv', 'w') as f:
        f.write("Order\tPDB_Chain\tType\tLength\tSimilarity\tSource\n")
        for i, record in enumerate(records, 1):
            f.write(f"{i}\t{record['pdb_chain']}\t{record['type']}\t{record['length']}\t{record['similarity']:.4f}\t{record['source']}\n")
    
    # Print statistics
    lengths = [r['length'] for r in records]
    similarities = [r['similarity'] for r in records]
    
    print(f"\nOptimal Training Set Statistics:")
    print(f"Total proteins: {len(records)}")
    print(f"Length range: {min(lengths)}-{max(lengths)} residues")
    print(f"Average length: {sum(lengths)/len(lengths):.1f} residues")
    print(f"Similarity range: {min(similarities):.3f}-{max(similarities):.3f}")
    print(f"Average similarity: {sum(similarities)/len(similarities):.3f}")
    
    # Size distribution
    size_ranges = [(0, 50), (50, 100), (100, 150), (150, 200), (200, 300)]
    for start, end in size_ranges:
        count = len([l for l in lengths if start <= l < end])
        if count > 0:
            print(f"{start}-{end-1} residues: {count} proteins ({count/len(lengths)*100:.1f}%)")

if __name__ == '__main__':
    main()
