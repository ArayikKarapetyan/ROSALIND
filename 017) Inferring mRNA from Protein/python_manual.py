def count_rna_strings_manual(protein):
    """Count possible RNA strings that could translate to given protein."""
    # RNA codon table with number of codons for each amino acid
    codon_counts = {
        'A': 4,  # Ala: GCU, GCC, GCA, GCG
        'C': 2,  # Cys: UGU, UGC
        'D': 2,  # Asp: GAU, GAC
        'E': 2,  # Glu: GAA, GAG
        'F': 2,  # Phe: UUU, UUC
        'G': 4,  # Gly: GGU, GGC, GGA, GGG
        'H': 2,  # His: CAU, CAC
        'I': 3,  # Ile: AUU, AUC, AUA
        'K': 2,  # Lys: AAA, AAG
        'L': 6,  # Leu: UUA, UUG, CUU, CUC, CUA, CUG
        'M': 1,  # Met: AUG (also start)
        'N': 2,  # Asn: AAU, AAC
        'P': 4,  # Pro: CCU, CCC, CCA, CCG
        'Q': 2,  # Gln: CAA, CAG
        'R': 6,  # Arg: CGU, CGC, CGA, CGG, AGA, AGG
        'S': 6,  # Ser: UCU, UCC, UCA, UCG, AGU, AGC
        'T': 4,  # Thr: ACU, ACC, ACA, ACG
        'V': 4,  # Val: GUU, GUC, GUA, GUG
        'W': 1,  # Trp: UGG
        'Y': 2,  # Tyr: UAU, UAC
        # Stop codons: 3 possibilities (UAA, UAG, UGA)
    }
    
    mod = 1000000
    
    # Start with 1 (empty product)
    total = 1
    
    # Multiply by codon counts for each amino acid
    for aa in protein:
        if aa not in codon_counts:
            return 0  # Invalid amino acid
        total = (total * codon_counts[aa]) % mod
    
    # Multiply by 3 for stop codon
    total = (total * 3) % mod
    
    return total

with open("Dataset.txt", "r") as file:
    protein = file.read().strip()

print(f"Python Manual: {count_rna_strings_manual(protein)}")
