from Bio import SeqIO
import numpy as np

def p_distance_biopython(seq1, seq2):
    """
    Calculate p-distance for BioPython Seq objects.
    """
    s1 = str(seq1)
    s2 = str(seq2)
    
    if len(s1) != len(s2):
        raise ValueError("Sequences must have equal length")
    
    differences = sum(c1 != c2 for c1, c2 in zip(s1, s2))
    return differences / len(s1)

def compute_distance_matrix_biopython(records):
    """
    Compute distance matrix from BioPython SeqRecords.
    """
    n = len(records)
    matrix = np.zeros((n, n))
    
    for i in range(n):
        for j in range(i, n):
            if i == j:
                matrix[i, j] = 0.0
            else:
                dist = p_distance_biopython(records[i].seq, records[j].seq)
                matrix[i, j] = dist
                matrix[j, i] = dist
    
    return matrix

def main():
    # Read sequences using BioPython
    records = list(SeqIO.parse("Dataset.txt", "fasta"))
    
    # Compute distance matrix
    matrix = compute_distance_matrix_biopython(records)
    
    # Output matrix
    for i in range(len(records)):
        row = " ".join(f"{matrix[i, j]:.5f}" for j in range(len(records)))
        print(row)

if __name__ == "__main__":
    main()