def shortest_common_supersequence(s, t):
    """
    Find shortest common supersequence of two strings using dynamic programming.
    """
    m, n = len(s), len(t)
    
    # DP table: dp[i][j] = length of SCS for s[:i] and t[:j]
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    
    # Initialize first row and column
    for i in range(m + 1):
        dp[i][0] = i  # Need to add all characters of s
    for j in range(n + 1):
        dp[0][j] = j  # Need to add all characters of t
    
    # Fill DP table
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if s[i-1] == t[j-1]:
                dp[i][j] = dp[i-1][j-1] + 1
            else:
                dp[i][j] = min(dp[i-1][j], dp[i][j-1]) + 1
    
    # Reconstruct SCS from DP table
    scs_chars = []
    i, j = m, n
    
    while i > 0 and j > 0:
        if s[i-1] == t[j-1]:
            scs_chars.append(s[i-1])
            i -= 1
            j -= 1
        elif dp[i-1][j] < dp[i][j-1]:
            scs_chars.append(s[i-1])
            i -= 1
        else:
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
    # Read input from Dataset.txt
    with open("Dataset.txt", "r") as f:
        lines = f.readlines()
    
    # Get the two DNA strings
    s = lines[0].strip()
    t = lines[1].strip()
    
    # Find shortest common supersequence
    scs = shortest_common_supersequence(s, t)
    
    # Output result
    print(scs)

if __name__ == "__main__":
    main()