def count_noncrossing_matchings(rna):
    """
    Count noncrossing matchings (not necessarily perfect) in RNA string.
    Only canonical base pairs: A-U and C-G.
    Uses dynamic programming similar to Motzkin numbers with constraints.
    """
    n = len(rna)
    MOD = 1000000
    
    # DP table: dp[i][j] = number of matchings for substring rna[i:j+1]
    dp = [[0] * n for _ in range(n)]
    
    # Initialize: empty substring has 1 matching (empty matching)
    for i in range(n + 1):
        if i < n:
            dp[i][i - 1] = 1  # empty substring
    
    # Fill DP table for increasing lengths
    for length in range(n):
        for i in range(n - length):
            j = i + length
            
            # Case 1: base i is unpaired
            total = dp[i + 1][j] if i + 1 <= j else 1
            
            # Case 2: base i is paired with some k
            for k in range(i + 1, j + 1):
                # Check if bases can pair
                if (rna[i] == 'A' and rna[k] == 'U') or \
                   (rna[i] == 'U' and rna[k] == 'A') or \
                   (rna[i] == 'C' and rna[k] == 'G') or \
                   (rna[i] == 'G' and rna[k] == 'C'):
                    
                    # Split into two independent regions
                    left = dp[i + 1][k - 1] if i + 1 <= k - 1 else 1
                    right = dp[k + 1][j] if k + 1 <= j else 1
                    total = (total + left * right) % MOD
            
            dp[i][j] = total % MOD
    
    return dp[0][n - 1] % MOD

def main():
    # Read FASTA file
    with open("Dataset.txt", "r") as f:
        lines = f.readlines()
    
    # Extract RNA sequence
    rna = ""
    for line in lines:
        if not line.startswith(">"):
            rna += line.strip()
    
    # Count noncrossing matchings
    result = count_noncrossing_matchings(rna)
    print(result)

if __name__ == "__main__":
    main()