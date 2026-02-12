def longest_common_subsequence(s, t):
    """
    Find longest common subsequence of two strings using dynamic programming.
    Returns one LCS (if multiple exist).
    """
    m, n = len(s), len(t)
    
    # DP table: dp[i][j] = length of LCS for s[:i] and t[:j]
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    
    # Fill DP table
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if s[i-1] == t[j-1]:
                dp[i][j] = dp[i-1][j-1] + 1
            else:
                dp[i][j] = max(dp[i-1][j], dp[i][j-1])
    
    # Backtrack to find one LCS
    lcs_chars = []
    i, j = m, n
    while i > 0 and j > 0:
        if s[i-1] == t[j-1]:
            lcs_chars.append(s[i-1])
            i -= 1
            j -= 1
        elif dp[i-1][j] > dp[i][j-1]:
            i -= 1
        else:
            j -= 1
    
    # Reverse to get correct order
    return ''.join(reversed(lcs_chars))

def main():
    # Read FASTA file
    sequences = []
    with open("Dataset.txt", "r") as f:
        current_seq = ""
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_seq:
                    sequences.append(current_seq)
                    current_seq = ""
            else:
                current_seq += line
        if current_seq:
            sequences.append(current_seq)
    
    # Get the two sequences
    s, t = sequences[0], sequences[1]
    
    # Find LCS
    lcs = longest_common_subsequence(s, t)
    
    # Output result
    print(lcs)

if __name__ == "__main__":
    main()