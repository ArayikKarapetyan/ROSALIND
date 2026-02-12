import numpy as np

def p_distance_fast(s1, s2):
    """
    Fast p-distance calculation using numpy.
    """
    # Convert strings to numpy arrays of characters
    arr1 = np.array(list(s1))
    arr2 = np.array(list(s2))
    
    # Count mismatches
    mismatches = np.sum(arr1 != arr2)
    return mismatches / len(s1)

def compute_distance_matrix_efficient(sequences):
    """
    Compute distance matrix efficiently.
    """
    n = len(sequences)
    length = len(sequences[0])
    
    # Convert sequences to 2D numpy array
    seq_array = np.array([list(seq) for seq in sequences])
    
    # Initialize distance matrix
    dist_matrix = np.zeros((n, n))
    
    # Compute distances
    for i in range(n):
        for j in range(i + 1, n):
            dist = np.sum(seq_array[i] != seq_array[j]) / length
            dist_matrix[i, j] = dist
            dist_matrix[j, i] = dist
    
    return dist_matrix

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
    
    # Verify all sequences have same length
    length = len(sequences[0])
    for i, seq in enumerate(sequences):
        if len(seq) != length:
            print(f"Error: Sequence {i} has different length")
            return
    
    # Compute distance matrix
    matrix = compute_distance_matrix_efficient(sequences)
    
    # Output with 5 decimal places
    for row in matrix:
        print(" ".join(f"{val:.5f}" for val in row))

if __name__ == "__main__":
    main()