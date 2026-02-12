def shortest_common_supersequence_efficient(s, t):
    """
    Efficient SCS with LCS-based reconstruction.
    Length of SCS = len(s) + len(t) - len(LCS)
    """
    m, n = len(s), len(t)
    
    # First compute LCS length table
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if s[i-1] == t[j-1]:
                dp[i][j] = dp[i-1][j-1] + 1
            else:
                dp[i][j] = max(dp[i-1][j], dp[i][j-1])
    
    # Reconstruct SCS using LCS table
    scs_chars = []
    i, j = m, n
    
    while i > 0 and j > 0:
        if s[i-1] == t[j-1]:
            # Character is part of LCS, add once
            scs_chars.append(s[i-1])
            i -= 1
            j -= 1
        elif dp[i-1][j] > dp[i][j-1]:
            # Character from s is not in LCS
            scs_chars.append(s[i-1])
            i -= 1
        else:
            # Character from t is not in LCS
            scs_chars.append(t[j-1])
            j -= 1
    
    # Add remaining characters
    while i > 0:
        scs_chars.append(s[i-1])
        i -= 1
    while j > 0:
        scs_chars.append(t[j-1])
        j -= 1
    
    # Reverse to get correct order
    return ''.join(reversed(scs_chars))

def main():
    # Read input
    with open("Dataset.txt", "r") as f:
        lines = f.readlines()
    
    s = lines[0].strip()
    t = lines[1].strip()
    
    scs = shortest_common_supersequence_efficient(s, t)
    print(scs)

if __name__ == "__main__":
    main()