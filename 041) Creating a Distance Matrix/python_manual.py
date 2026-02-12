def p_distance(s1, s2):
    """
    Calculate p-distance between two strings of equal length.
    p-distance = proportion of differing characters.
    """
    if len(s1) != len(s2):
        raise ValueError("Strings must have equal length")
    
    differences = sum(c1 != c2 for c1, c2 in zip(s1, s2))
    return differences / len(s1)

def compute_distance_matrix(sequences):
    """
    Compute p-distance matrix for a list of sequences.
    """
    n = len(sequences)
    matrix = [[0.0] * n for _ in range(n)]
    
    for i in range(n):
        for j in range(i, n):
            if i == j:
                matrix[i][j] = 0.0
            else:
                dist = p_distance(sequences[i], sequences[j])
                matrix[i][j] = dist
                matrix[j][i] = dist
    
    return matrix

def main():
    # Read FASTA file
    sequences = []
    with open("Dataset.txt", "r") as f:
        current_seq = ""
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_seq:
                    sequences.append(current_seq)
                    current_seq = ""
            else:
                current_seq += line
        if current_seq:
            sequences.append(current_seq)
    
    # Compute distance matrix
    matrix = compute_distance_matrix(sequences)
    
    # Output matrix with 5 decimal places
    for row in matrix:
        print(" ".join(f"{value:.5f}" for value in row))

if __name__ == "__main__":
    main()