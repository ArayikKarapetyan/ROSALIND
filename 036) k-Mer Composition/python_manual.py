def generate_all_kmers(k, alphabet='ACGT'):
    """Generate all k-mers in lexicographic order."""
    from itertools import product
    return [''.join(p) for p in product(alphabet, repeat=k)]

def kmer_composition(dna, k=4):
    """Calculate k-mer composition for DNA string."""
    # Generate all possible k-mers in lexicographic order
    all_kmers = generate_all_kmers(k)
    
    # Create dictionary for counting
    kmer_counts = {kmer: 0 for kmer in all_kmers}
    
    # Count occurrences of each k-mer
    for i in range(len(dna) - k + 1):
        kmer = dna[i:i+k]
        if kmer in kmer_counts:
            kmer_counts[kmer] += 1
    
    # Return counts in lexicographic order
    return [kmer_counts[kmer] for kmer in all_kmers]

def main():
    # Read FASTA file
    with open("Dataset.txt", "r") as f:
        lines = f.readlines()
    
    # Extract DNA sequence
    dna = ""
    for line in lines:
        if not line.startswith(">"):
            dna += line.strip()
    
    # Calculate 4-mer composition
    composition = kmer_composition(dna, 4)
    
    # Output result
    print(" ".join(map(str, composition)))

if __name__ == "__main__":
    main()