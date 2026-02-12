from Bio.Seq import Seq

def reverse_complement_biopython(dna_string):
    """Compute reverse complement using BioPython."""
    dna_seq = Seq(dna_string)
    return str(dna_seq.reverse_complement())


path = "Dataset.txt"    # !!! Make sure the path is correct !!! 

with open(path, 'r') as file:
    data = file.read().strip()

print("Python BioPython:", reverse_complement_biopython(data))