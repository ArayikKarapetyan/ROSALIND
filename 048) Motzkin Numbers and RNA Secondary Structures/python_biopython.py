from Bio import SeqIO
from functools import lru_cache

def count_rna_matchings_biopython(seq):
    """
    Count RNA noncrossing matchings using BioPython sequence.
    """
    rna = str(seq)
    n = len(rna)
    MOD = 1000000
    
    # Define valid base pairs
    valid_pairs = {('A', 'U'), ('U', 'A'), ('C', 'G'), ('G', 'C')}
    
    @lru_cache(maxsize=None)
    def dp(i, j):
        if i >= j:  # empty or single base
            return 1
        
        total = dp(i + 1, j)  # i unpaired
        
        for k in range(i + 1, j + 1):
            if (rna[i], rna[k]) in valid_pairs:
                total = (total + dp(i + 1, k - 1) * dp(k + 1, j)) % MOD
        
        return total
    
    return dp(0, n - 1) % MOD

def main():
    # Read sequence using BioPython
    record = SeqIO.read("Dataset.txt", "fasta")
    
    # Count matchings
    result = count_rna_matchings_biopython(record.seq)
    
    # Output result
    print(result)

if __name__ == "__main__":
    main()