from Bio import SeqIO
from Bio.Seq import Seq
from io import StringIO

def rna_splicing_biopython(fasta_string):
    """RNA splicing using BioPython."""
    # Parse FASTA
    fasta_io = StringIO(fasta_string)
    records = list(SeqIO.parse(fasta_io, "fasta"))
    
    # First record is DNA, rest are introns
    dna_seq = str(records[0].seq)
    introns = [str(record.seq) for record in records[1:]]
    
    # Remove introns
    exons_only = dna_seq
    for intron in introns:
        exons_only = exons_only.replace(intron, "")
    
    # Transcribe and translate
    exons_seq = Seq(exons_only)
    rna_seq = exons_seq.transcribe()
    protein_seq = rna_seq.translate(to_stop=True)
    
    return str(protein_seq)


with open("Dataset.txt", "r") as file:
    fasta_data = file.read()

print("\nPython BioPython Solution:")
print(rna_splicing_biopython(fasta_data))
