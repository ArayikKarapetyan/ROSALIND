from Bio.Seq import Seq

# Method with using special library: Biopython
def transcribe_dna_to_rna_biopython(dna_string):
    """Alternative BioPython approach."""
    return str(Seq(dna_string).transcribe())


path = "Dataset.txt"    # !!! Make sure the path is correct !!! 

with open(path, 'r') as file:
    data = file.read().strip()


# Method 2: BioPython (requires: pip install biopython)
print("BioPython:", transcribe_dna_to_rna_biopython(data))