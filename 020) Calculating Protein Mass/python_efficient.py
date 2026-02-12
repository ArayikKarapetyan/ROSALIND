# Alternative using list for faster lookup
def calculate_protein_mass_array(protein):
    """Calculate mass using array lookup for efficiency."""
    # Create array for ASCII codes (A=65, Z=90)
    mass_array = [0.0] * 128  # ASCII extended
    
    # Fill with masses
    mass_array[ord('A')] = 71.03711
    mass_array[ord('C')] = 103.00919
    mass_array[ord('D')] = 115.02694
    mass_array[ord('E')] = 129.04259
    mass_array[ord('F')] = 147.06841
    mass_array[ord('G')] = 57.02146
    mass_array[ord('H')] = 137.05891
    mass_array[ord('I')] = 113.08406
    mass_array[ord('K')] = 128.09496
    mass_array[ord('L')] = 113.08406
    mass_array[ord('M')] = 131.04049
    mass_array[ord('N')] = 114.04293
    mass_array[ord('P')] = 97.05276
    mass_array[ord('Q')] = 128.05858
    mass_array[ord('R')] = 156.10111
    mass_array[ord('S')] = 87.03203
    mass_array[ord('T')] = 101.04768
    mass_array[ord('V')] = 99.06841
    mass_array[ord('W')] = 186.07931
    mass_array[ord('Y')] = 163.06333
    
    total_mass = 0.0
    
    for aa in protein:
        mass = mass_array[ord(aa)]
        if mass == 0.0:
            raise ValueError(f"Invalid amino acid: {aa}")
        total_mass += mass
    
    return total_mass

with open("Dataset.txt", "r") as file:
    protein = file.read().strip()

print(f"Python Array: {calculate_protein_mass_array(protein):.3f}")
