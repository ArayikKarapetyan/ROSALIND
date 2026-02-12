# Despite the name, this might not be so efficient version of solution, you can also get RecursionError, but still we can learn some concepts from this code I think 

from functools import lru_cache

def lcs_efficient(s, t):
    """
    Efficient LCS using memoization and recursion.
    Returns one LCS.
    """
    
    @lru_cache(maxsize=None)
    def dp(i, j):
        """Length of LCS for s[:i] and t[:j]"""
        if i == 0 or j == 0:
            return 0
        if s[i-1] == t[j-1]:
            return dp(i-1, j-1) + 1
        return max(dp(i-1, j), dp(i, j-1))
    
    # Backtrack using DP function
    def backtrack(i, j):
        if i == 0 or j == 0:
            return ""
        if s[i-1] == t[j-1]:
            return backtrack(i-1, j-1) + s[i-1]
        if dp(i-1, j) >= dp(i, j-1):
            return backtrack(i-1, j)
        return backtrack(i, j-1)
    
    return backtrack(len(s), len(t))

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
    
    s, t = sequences[0], sequences[1]
    result = lcs_efficient(s, t)
    print(result)

if __name__ == "__main__":
    main()