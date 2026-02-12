def edit_distance_efficient(s, t):
    """
    Edit distance with full matrix for clarity.
    Still O(mn) time but O(mn) space.
    """
    m, n = len(s), len(t)
    
    # Create DP table
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    
    # Initialize first row and column
    for i in range(m + 1):
        dp[i][0] = i  # i deletions
    for j in range(n + 1):
        dp[0][j] = j  # j insertions
    
    # Fill DP table
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if s[i-1] == t[j-1]:
                sub_cost = dp[i-1][j-1]
            else:
                sub_cost = dp[i-1][j-1] + 1
            
            del_cost = dp[i-1][j] + 1
            ins_cost = dp[i][j-1] + 1
            
            dp[i][j] = min(sub_cost, del_cost, ins_cost)
    
    return dp[m][n]

def main():
    # Read input
    with open("Dataset.txt", "r") as f:
        lines = f.readlines()
    
    sequences = []
    current_seq = ""
    for line in lines:
        line = line.strip()
        if line.startswith(">"):
            if current_seq:
                sequences.append(current_seq)
                current_seq = ""
        else:
            current_seq += line
    if current_seq:
        sequences.append(current_seq)
    
    s, t = sequences[0], sequences[1]
    distance = edit_distance_efficient(s, t)
    print(distance)

if __name__ == "__main__":
    main()