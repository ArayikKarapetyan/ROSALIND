from Bio.Seq import Seq
from Bio import pairwise2

def hamming_distance_biopython(s, t):
    """Calculate Hamming distance using alignment tools."""
    
    # For Hamming distance, sequences must be same length
    if len(s) != len(t):
        raise ValueError("Strings must be of equal length")
    
    # Create a simple alignment and count mismatches
    alignments = pairwise2.align.globalxx(s, t)
    best_alignment = alignments[0]
    
    # Count mismatches in aligned sequences
    aligned_s = best_alignment.seqA
    aligned_t = best_alignment.seqB
    
    mismatches = sum(1 for a, b in zip(aligned_s, aligned_t) if a != b)
    return mismatches



# print(f"Python BioPython: {hamming_distance_biopython(s, t)}")        

# !!! At the time I executed Biopython version I got this warning, so in your case it might  work correctly. !!! 
# ------------------------------------------------------------------------------------------------------------------------------------------
# BiopythonDeprecationWarning: Bio.pairwise2 has been deprecated, and we intend to remove it in a future release of Biopython. As an alternative, please consider using Bio.Align.PairwiseAligner as a replacement, and contact the Biopython developers if you still need the Bio.pairwise2 module.
# ------------------------------------------------------------------------------------------------------------------------------------------