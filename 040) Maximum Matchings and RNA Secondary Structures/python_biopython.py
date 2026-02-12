from Bio import SeqIO
from math import factorial
from collections import Counter

def count_maximum_matchings_biopython(seq):
    """
    Count maximum matchings using BioPython sequence object.
    """
    rna = str(seq)
    counts = Counter(rna)
    
    A = counts.get('A', 0)
    U = counts.get('U', 0)
    C = counts.get('C', 0)
    G = counts.get('G', 0)
    
    # A-U pairings
    if A >= U:
        au_ways = factorial(A) // factorial(A - U)
    else:
        au_ways = factorial(U) // factorial(U - A)
    
    # C-G pairings
    if C >= G:
        cg_ways = factorial(C) // factorial(C - G)
    else:
        cg_ways = factorial(G) // factorial(G - C)
    
    return au_ways * cg_ways

def main():
    # Read sequence from FASTA file
    record = SeqIO.read("Dataset.txt", "fasta")
    
    # Count maximum matchings
    result = count_maximum_matchings_biopython(record.seq)
    
    # Output result
    print(result)

if __name__ == "__main__":
    main()