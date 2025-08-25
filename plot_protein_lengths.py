#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np

def read_protein_lengths(filename):
    """Read protein lengths from the filtered selection info file."""
    lengths = []
    with open(filename, 'r') as f:
        lines = f.readlines()
        # Skip header
        for line in lines[1:]:
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                length = int(parts[3])  # Length is in column 4
                lengths.append(length)
    return lengths

def create_histogram(lengths, output_file='protein_lengths_histogram.png'):
    """Create and save a histogram of protein lengths."""
    plt.figure(figsize=(10, 6))
    
    # Create histogram with appropriate bins
    bins = np.arange(0, max(lengths) + 20, 20)  # 20-residue bins
    n, bins, patches = plt.hist(lengths, bins=bins, edgecolor='black', alpha=0.7, color='skyblue')
    
    plt.xlabel('Protein Length (residues)', fontsize=12)
    plt.ylabel('Number of Proteins', fontsize=12)
    plt.title(f'Distribution of Protein Lengths in Alpha Training Set\n(Total: {len(lengths)} proteins)', fontsize=14)
    plt.grid(True, alpha=0.3)
    
    # Add statistics text
    mean_length = np.mean(lengths)
    median_length = np.median(lengths)
    min_length = min(lengths)
    max_length = max(lengths)
    
    stats_text = f'Min: {min_length}\nMax: {max_length}\nMean: {mean_length:.1f}\nMedian: {median_length:.1f}'
    plt.text(0.7, 0.9, stats_text, transform=plt.gca().transAxes, 
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
             verticalalignment='top', fontsize=10)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Histogram saved as {output_file}")
    
    # Print summary statistics
    print(f"\nProtein Length Statistics:")
    print(f"Total proteins: {len(lengths)}")
    print(f"Min length: {min_length} residues")
    print(f"Max length: {max_length} residues")
    print(f"Mean length: {mean_length:.1f} residues")
    print(f"Median length: {median_length:.1f} residues")
    print(f"Standard deviation: {np.std(lengths):.1f} residues")
    
    # Show distribution by ranges
    print(f"\nLength Distribution:")
    ranges = [(0, 50), (50, 100), (100, 150), (150, 200), (200, 250)]
    for start, end in ranges:
        count = len([l for l in lengths if start <= l < end])
        print(f"{start}-{end-1} residues: {count} proteins ({count/len(lengths)*100:.1f}%)")

def main():
    lengths = read_protein_lengths('selection/filtered_alpha_selection_info.tsv')
    if lengths:
        create_histogram(lengths)
    else:
        print("No protein length data found!")

if __name__ == '__main__':
    main()
