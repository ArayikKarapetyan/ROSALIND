from Bio import SeqIO
from io import StringIO

def find_subsequence_biopython(fasta_string):
    """Find subsequence indices using BioPython."""
    # Parse FASTA using BioPython
    fasta_io = StringIO(fasta_string)
    records = list(SeqIO.parse(fasta_io, "fasta"))
    
    s = str(records[0].seq)
    t = str(records[1].seq)
    
    # Use greedy algorithm (BioPython doesn't have built-in subsequence finder)
    indices = []
    i = j = 0
    
    while i < len(s) and j < len(t):
        if s[i] == t[j]:
            indices.append(i + 1)
            j += 1
        i += 1
    
    if j == len(t):
        return indices
    else:
        raise ValueError("t is not a subsequence of s")
    

with open("Dataset.txt", "r") as file:
    fasta_data = file.read()
indices_biopython = find_subsequence_biopython(fasta_data)
print("BioPython:", " ".join(map(str, indices_biopython)))