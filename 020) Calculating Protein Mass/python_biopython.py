from Bio.SeqUtils import molecular_weight

def calculate_protein_mass_biopython(protein):
    """Calculate protein mass using BioPython."""
    # BioPython's molecular_weight calculates average mass by default
    # We need monoisotopic mass

    mass = molecular_weight(protein, seq_type='protein', monoisotopic=True)
    return mass

with open("Dataset.txt", "r") as file:
    protein = file.read().strip()


print(f"Python BioPython: {calculate_protein_mass_biopython(protein):.3f}")