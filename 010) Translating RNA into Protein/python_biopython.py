from Bio.Seq import Seq

def rna_to_protein_biopython(rna_string):
    """Translate RNA to protein using BioPython."""
    rna_seq = Seq(rna_string)
    protein_seq = rna_seq.translate(to_stop=True)  # Stop at stop codon
    return str(protein_seq)


with open("Dataset.txt", "r") as file:
    data = file.read().strip()

print("Python BioPython:", rna_to_protein_biopython(data))