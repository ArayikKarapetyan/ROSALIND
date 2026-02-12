from Bio import Seq

def shortest_common_supersequence_biopython(seq1, seq2):
    """
    SCS using BioPython sequences.
    """
    s = str(seq1)
    t = str(seq2)
    
    m, n = len(s), len(t)
    
    # Compute LCS length table
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if s[i-1] == t[j-1]:
                dp[i][j] = dp[i-1][j-1] + 1
            else:
                dp[i][j] = max(dp[i-1][j], dp[i][j-1])
    
    # Reconstruct SCS
    result = []
    i, j = m, n
    
    while i > 0 and j > 0:
        if s[i-1] == t[j-1]:
            result.append(s[i-1])
            i -= 1
            j -= 1
        elif dp[i-1][j] >= dp[i][j-1]:
            result.append(s[i-1])
            i -= 1
        else:
            result.append(t[j-1])
            j -= 1
    
    # Add remaining
    while i > 0:
        result.append(s[i-1])
        i -= 1
    while j > 0:
        result.append(t[j-1])
        j -= 1
    
    return ''.join(reversed(result))

def main():
    # Read input (simple format for this problem)
    with open("Dataset.txt", "r") as f:
        lines = f.readlines()
    
    s = lines[0].strip()
    t = lines[1].strip()
    
    scs = shortest_common_supersequence_biopython(s, t)
    print(scs)

if __name__ == "__main__":
    main()