from functools import lru_cache

def count_rna_matchings_efficient(rna):
    """
    Efficient computation using memoization (top-down DP).
    Counts all noncrossing matchings (not necessarily perfect).
    """
    n = len(rna)
    MOD = 1000000
    
    # Precompute pairing matrix for faster lookup
    can_pair = {
        ('A', 'U'): True, ('U', 'A'): True,
        ('C', 'G'): True, ('G', 'C'): True,
    }
    
    @lru_cache(maxsize=None)
    def count(i, j):
        """Count matchings for substring rna[i:j+1]"""
        if i > j:  # empty substring
            return 1
        
        # Start with case where i is unpaired
        total = count(i + 1, j)
        
        # Try pairing i with every possible k
        for k in range(i + 1, j + 1):
            if (rna[i], rna[k]) in can_pair:
                # Split: (i+1..k-1) and (k+1..j)
                left = count(i + 1, k - 1)
                right = count(k + 1, j)
                total = (total + left * right) % MOD
        
        return total % MOD
    
    return count(0, n - 1)

def main():
    # Read FASTA file
    with open("Dataset.txt", "r") as f:
        lines = f.readlines()
    
    rna = ""
    for line in lines:
        if not line.startswith(">"):
            rna += line.strip().upper()
    
    result = count_rna_matchings_efficient(rna)
    print(result)

if __name__ == "__main__":
    main()