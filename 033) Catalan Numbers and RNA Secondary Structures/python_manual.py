def count_noncrossing_matchings(rna):
    """
    Count noncrossing perfect matchings in RNA string using dynamic programming.
    Only canonical base pairs: A-U and C-G.
    """
    n = len(rna)
    if n % 2 != 0:
        return 0
    
    # DP table: dp[i][j] = number of matchings for substring rna[i:j+1]
    dp = [[0] * n for _ in range(n)]
    
    # Base case: empty substring has 1 matching
    for i in range(n + 1):
        if i < n:
            dp[i][i - 1] = 1  # empty substring
    
    # Fill DP table for increasing lengths
    for length in range(2, n + 1, 2):  # only even lengths
        for i in range(n - length + 1):
            j = i + length - 1
            
            # Initialize with matching where i is not paired
            total = 0
            
            # Try pairing i with k (i+1, i+3, i+5, ...)
            for k in range(i + 1, j + 1, 2):
                # Check if bases can pair
                if (rna[i] == 'A' and rna[k] == 'U') or \
                   (rna[i] == 'U' and rna[k] == 'A') or \
                   (rna[i] == 'C' and rna[k] == 'G') or \
                   (rna[i] == 'G' and rna[k] == 'C'):
                    
                    # Split into two independent regions
                    left = dp[i + 1][k - 1] if i + 1 <= k - 1 else 1
                    right = dp[k + 1][j] if k + 1 <= j else 1
                    total = (total + left * right) % 1000000
            
            dp[i][j] = total % 1000000
    
    return dp[0][n - 1] % 1000000


def main():
    # Read input
    with open("Dataset.txt", "r") as f:
        lines = f.readlines()
    
    # Skip fasta header if present
    rna = ""
    for line in lines:
        if not line.startswith(">"):
            rna += line.strip()
    
    result = count_noncrossing_matchings(rna)
    print(result)


if __name__ == "__main__":
    main()