from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

# Alternative with biopython

def find_highest_gc_biopython(fasta_string):
    """Find sequence with highest GC content using BioPython."""
    from io import StringIO
    
    highest_id = ""
    highest_gc = 0.0
    
    # Parse FASTA from string
    fasta_io = StringIO(fasta_string)
    for record in SeqIO.parse(fasta_io, "fasta"):
        gc_content = gc_fraction(record.seq) * 100
        if gc_content > highest_gc:
            highest_gc = gc_content
            highest_id = record.id
    
    return highest_id, highest_gc



path = "Dataset.txt"    # !!! Make sure the path is correct !!! 

with open(path, 'r') as file:
    data = file.read().strip()


id2, gc2 = find_highest_gc_biopython(data)
print(f"Python BioPython: {id2}\n{gc2:.6f}")