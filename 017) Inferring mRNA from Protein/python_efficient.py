def count_rna_strings_list(protein):
    """Count possible RNA strings using list indexing."""
    # Create a mapping from amino acid letter to codon count
    # Using ASCII codes for efficiency
    counts = [0] * 256  # ASCII table size
    
    # Fill with codon counts
    counts[ord('A')] = 4  # Alanine
    counts[ord('C')] = 2  # Cysteine
    counts[ord('D')] = 2  # Aspartic acid
    counts[ord('E')] = 2  # Glutamic acid
    counts[ord('F')] = 2  # Phenylalanine
    counts[ord('G')] = 4  # Glycine
    counts[ord('H')] = 2  # Histidine
    counts[ord('I')] = 3  # Isoleucine
    counts[ord('K')] = 2  # Lysine
    counts[ord('L')] = 6  # Leucine
    counts[ord('M')] = 1  # Methionine
    counts[ord('N')] = 2  # Asparagine
    counts[ord('P')] = 4  # Proline
    counts[ord('Q')] = 2  # Glutamine
    counts[ord('R')] = 6  # Arginine
    counts[ord('S')] = 6  # Serine
    counts[ord('T')] = 4  # Threonine
    counts[ord('V')] = 4  # Valine
    counts[ord('W')] = 1  # Tryptophan
    counts[ord('Y')] = 2  # Tyrosine
    
    mod = 1000000
    total = 1
    
    for aa in protein:
        count = counts[ord(aa)]
        if count == 0:
            return 0  # Invalid amino acid
        total = (total * count) % mod
    
    # Multiply by stop codon possibilities
    total = (total * 3) % mod
    
    return total

with open("Dataset.txt", "r") as file:
    protein = file.read().strip()

print(f"Python List: {count_rna_strings_list(protein)}")
