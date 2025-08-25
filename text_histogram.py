#!/usr/bin/env python3

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

def create_text_histogram(lengths):
    """Create a text-based histogram of protein lengths."""
    print(f"\nProtein Length Distribution (Total: {len(lengths)} proteins)")
    print("=" * 60)
    
    # Calculate statistics
    min_length = min(lengths)
    max_length = max(lengths)
    mean_length = sum(lengths) / len(lengths)
    sorted_lengths = sorted(lengths)
    median_length = sorted_lengths[len(lengths) // 2]
    std_dev = (sum((x - mean_length) ** 2 for x in lengths) / len(lengths)) ** 0.5
    
    print(f"Min length: {min_length} residues")
    print(f"Max length: {max_length} residues")
    print(f"Mean length: {mean_length:.1f} residues")
    print(f"Median length: {median_length:.1f} residues")
    print(f"Standard deviation: {std_dev:.1f} residues")
    print()
    
    # Create bins and count frequencies
    bin_size = 20
    min_bin = (min_length // bin_size) * bin_size
    max_bin = ((max_length // bin_size) + 1) * bin_size
    
    bins = {}
    for i in range(min_bin, max_bin, bin_size):
        bins[i] = 0
    
    for length in lengths:
        bin_start = (length // bin_size) * bin_size
        bins[bin_start] += 1
    
    # Find max count for scaling
    max_count = max(bins.values())
    scale = 50 / max_count if max_count > 0 else 1
    
    print("Length Range (residues)  Count    Percentage  Histogram")
    print("-" * 60)
    
    for bin_start in sorted(bins.keys()):
        bin_end = bin_start + bin_size - 1
        count = bins[bin_start]
        percentage = (count / len(lengths)) * 100
        bar_length = int(count * scale)
        bar = "█" * bar_length
        
        print(f"{bin_start:3d} - {bin_end:3d}            {count:3d}     {percentage:5.1f}%     {bar}")
    
    print()
    print("Detailed breakdown by size categories:")
    print("-" * 40)
    
    categories = [
        ("Very small", 0, 50),
        ("Small", 50, 100), 
        ("Medium", 100, 150),
        ("Large", 150, 200),
        ("Very large", 200, 250)
    ]
    
    for name, start, end in categories:
        count = len([l for l in lengths if start <= l < end])
        percentage = (count / len(lengths)) * 100
        print(f"{name:12} ({start:3d}-{end-1:3d}): {count:3d} proteins ({percentage:5.1f}%)")

def main():
    lengths = read_protein_lengths('selection/filtered_alpha_selection_info.tsv')
    if lengths:
        create_text_histogram(lengths)
    else:
        print("No protein length data found!")

if __name__ == '__main__':
    main()
