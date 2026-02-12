# May cause a RecursionError due to excessive recursion depth

from Bio import SeqIO
from functools import lru_cache

def lcs_biopython(seq1, seq2):
    """
    LCS using BioPython sequences and memoization.
    """
    s = str(seq1)
    t = str(seq2)
    
    @lru_cache(maxsize=None)
    def lcs_length(i, j):
        if i == 0 or j == 0:
            return 0
        if s[i-1] == t[j-1]:
            return lcs_length(i-1, j-1) + 1
        return max(lcs_length(i-1, j), lcs_length(i, j-1))
    
    # Backtrack to reconstruct LCS
    result = []
    i, j = len(s), len(t)
    
    while i > 0 and j > 0:
        if s[i-1] == t[j-1]:
            result.append(s[i-1])
            i -= 1
            j -= 1
        elif lcs_length(i-1, j) >= lcs_length(i, j-1):
            i -= 1
        else:
            j -= 1
    
    return ''.join(reversed(result))

def main():
    # Read sequences from FASTA file using BioPython
    records = list(SeqIO.parse("Dataset.txt", "fasta"))
    seq1 = records[0].seq
    seq2 = records[1].seq
    
    # Find LCS
    lcs = lcs_biopython(seq1, seq2)
    
    # Output result
    print(lcs)

if __name__ == "__main__":
    main()