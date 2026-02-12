from Bio import SeqIO
from itertools import product

def kmer_composition_biopython(seq, k=4):
    """Calculate k-mer composition using BioPython."""
    # Generate all k-mers in lexicographic order
    alphabet = ['A', 'C', 'G', 'T']
    all_kmers = [''.join(p) for p in product(alphabet, repeat=k)]
    
    # Initialize count dictionary
    kmer_counts = {kmer: 0 for kmer in all_kmers}
    
    # Count k-mers
    seq_str = str(seq)
    for i in range(len(seq_str) - k + 1):
        kmer = seq_str[i:i+k]
        if kmer in kmer_counts:
            kmer_counts[kmer] += 1
    
    # Return counts in lexicographic order
    return [kmer_counts[kmer] for kmer in all_kmers]

def main():
    # Read sequence from FASTA file
    record = SeqIO.read("Dataset.txt", "fasta")
    
    # Calculate 4-mer composition
    composition = kmer_composition_biopython(record.seq, 4)
    
    # Output result
    print(" ".join(map(str, composition)))

if __name__ == "__main__":
    main()