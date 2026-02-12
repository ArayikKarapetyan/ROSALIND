from math import comb
from collections import Counter

def count_maximum_matchings(rna):
    """
    Count maximum matchings in RNA considering Watson-Crick base pairs.
    For maximum matching, we pair as many bases as possible.
    """
    # Count each base
    counts = Counter(rna)
    A = counts['A']
    U = counts['U']
    C = counts['C']
    G = counts['G']
    
    # Maximum matching size is limited by min of complementary pairs
    # We can pair min(A, U) A-U pairs and min(C, G) C-G pairs
    
    # Count ways to choose which specific bases to pair
    # For A-U pairs: choose min(A, U) bases from max(A, U) to pair
    au_pairs = min(A, U)
    au_ways = 1
    
    if A >= U:
        # Choose which U's to pair with A's
        # Number of ways: U! (ways to pair specific U's with chosen A's)
        # Actually: number of injections from smaller set to larger set
        au_ways = comb(A, U) * comb(U, U)  # Choose which A's to use, then match
        # More precisely: P(A, U) = A! / (A-U)!
        from math import factorial
        au_ways = factorial(A) // factorial(A - U)
    else:
        from math import factorial
        au_ways = factorial(U) // factorial(U - A)
    
    # For C-G pairs: similar logic
    cg_pairs = min(C, G)
    cg_ways = 1
    
    if C >= G:
        from math import factorial
        cg_ways = factorial(C) // factorial(C - G)
    else:
        from math import factorial
        cg_ways = factorial(G) // factorial(G - C)
    
    # Total ways: multiply independent pairings
    return au_ways * cg_ways

def main():
    # Read FASTA file
    with open("Dataset.txt", "r") as f:
        lines = f.readlines()
    
    # Extract RNA sequence
    rna = ""
    for line in lines:
        if not line.startswith(">"):
            rna += line.strip()
    
    # Count maximum matchings
    result = count_maximum_matchings(rna)
    print(result)

if __name__ == "__main__":
    main()