from Bio import SeqIO
from io import StringIO

def calculate_tt_ratio_biopython(fasta_string):
    """Calculate transition/transversion ratio using BioPython."""
    fasta_io = StringIO(fasta_string)
    records = list(SeqIO.parse(fasta_io, "fasta"))
    
    s1 = str(records[0].seq)
    s2 = str(records[1].seq)
    
    transitions = 0
    transversions = 0
    
    # Define purines and pyrimidines
    purines = {'A', 'G'}
    pyrimidines = {'C', 'T'}
    
    for b1, b2 in zip(s1, s2):
        if b1 == b2:
            continue
        
        # Check if both are purines or both are pyrimidines (transition)
        if (b1 in purines and b2 in purines) or (b1 in pyrimidines and b2 in pyrimidines):
            transitions += 1
        else:
            transversions += 1
    
    if transversions == 0:
        return float('inf')
    
    return transitions / transversions

with open("Dataset.txt", "r") as file:
    fasta_data = file.read()

ratio_biopython = calculate_tt_ratio_biopython(fasta_data)
print(f"BioPython: {ratio_biopython:.11f}")