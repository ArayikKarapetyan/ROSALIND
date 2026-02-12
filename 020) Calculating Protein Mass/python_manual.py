def calculate_protein_mass_manual(protein):
    """Calculate protein mass using monoisotopic mass table."""
    # Monoisotopic mass table for amino acids
    mass_table = {
        'A': 71.03711,   # Alanine
        'C': 103.00919,  # Cysteine
        'D': 115.02694,  # Aspartic acid
        'E': 129.04259,  # Glutamic acid
        'F': 147.06841,  # Phenylalanine
        'G': 57.02146,   # Glycine
        'H': 137.05891,  # Histidine
        'I': 113.08406,  # Isoleucine
        'K': 128.09496,  # Lysine
        'L': 113.08406,  # Leucine
        'M': 131.04049,  # Methionine
        'N': 114.04293,  # Asparagine
        'P': 97.05276,   # Proline
        'Q': 128.05858,  # Glutamine
        'R': 156.10111,  # Arginine
        'S': 87.03203,   # Serine
        'T': 101.04768,  # Threonine
        'V': 99.06841,   # Valine
        'W': 186.07931,  # Tryptophan
        'Y': 163.06333   # Tyrosine
    }
    
    total_mass = 0.0
    
    for aa in protein:
        if aa in mass_table:
            total_mass += mass_table[aa]
        else:
            # Handle invalid amino acids if needed
            raise ValueError(f"Invalid amino acid: {aa}")
    
    return total_mass

with open("Dataset.txt", "r") as file:
    protein = file.read().strip()

print(f"Python Manual: {calculate_protein_mass_manual(protein):.3f}")
