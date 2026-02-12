from math import factorial
from collections import Counter

def count_maximum_matchings_efficient(rna):
    """
    Efficient counting of maximum RNA matchings using combinatorics.
    Maximum matching size = min(A,U) + min(C,G)
    Counting: permutations of choosing specific bases to pair.
    """
    counts = Counter(rna)
    A = counts.get('A', 0)
    U = counts.get('U', 0)
    C = counts.get('C', 0)
    G = counts.get('G', 0)
    
    # Calculate number of ways for A-U pairs
    if A >= U:
        au_ways = factorial(A) // factorial(A - U)
    else:
        au_ways = factorial(U) // factorial(U - A)
    
    # Calculate number of ways for C-G pairs
    if C >= G:
        cg_ways = factorial(C) // factorial(C - G)
    else:
        cg_ways = factorial(G) // factorial(G - C)
    
    return au_ways * cg_ways

def main():
    # Read input
    with open("Dataset.txt", "r") as f:
        lines = f.readlines()
    
    rna = ""
    for line in lines:
        if not line.startswith(">"):
            rna += line.strip().upper()
    
    result = count_maximum_matchings_efficient(rna)
    print(result)

if __name__ == "__main__":
    main()