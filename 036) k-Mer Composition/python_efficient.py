def kmer_composition_efficient(dna, k=4):
    """Efficient k-mer composition using integer encoding."""
    # Map bases to numbers: A=0, C=1, G=2, T=3
    base_to_num = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    
    # Pre-calculate powers of 4 for position weights
    powers = [4 ** (k - 1 - i) for i in range(k)]
    
    # Initialize array for all possible k-mers (4^k elements)
    total_kmers = 4 ** k
    counts = [0] * total_kmers
    
    # If sequence is shorter than k, return zeros
    if len(dna) < k:
        return counts
    
    # Calculate first k-mer index
    current_index = 0
    for i in range(k):
        current_index = current_index * 4 + base_to_num[dna[i]]
    counts[current_index] = 1
    
    # Slide window through the sequence
    for i in range(k, len(dna)):
        # Remove leftmost base contribution
        current_index %= (total_kmers // 4)
        # Shift left and add new base
        current_index = current_index * 4 + base_to_num[dna[i]]
        counts[current_index] += 1
    
    return counts

def main():
    # Read FASTA file
    with open("Dataset.txt", "r") as f:
        lines = f.readlines()
    
    # Extract DNA sequence
    dna = ""
    for line in lines:
        if not line.startswith(">"):
            dna += line.strip().upper()
    
    # Calculate 4-mer composition efficiently
    composition = kmer_composition_efficient(dna, 4)
    
    # Output result
    print(" ".join(map(str, composition)))

if __name__ == "__main__":
    main()