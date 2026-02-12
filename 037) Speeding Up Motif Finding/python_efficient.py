def failure_array_efficient(s):
    """
    Efficient computation of failure array (prefix function).
    O(n) time complexity.
    """
    n = len(s)
    P = [0] * n
    
    for i in range(1, n):
        j = P[i-1]
        # Use while loop to find matching prefix
        while j > 0 and s[i] != s[j]:
            j = P[j-1]
        # If characters match, extend the prefix
        if s[i] == s[j]:
            P[i] = j + 1
        else:
            P[i] = 0  # This is redundant since already 0, but explicit
    
    return P

def main():
    # Read FASTA file
    with open("Dataset.txt", "r") as f:
        lines = f.readlines()
    
    # Extract DNA sequence
    dna = ""
    for line in lines:
        if not line.startswith(">"):
            dna += line.strip()
    
    # Compute failure array
    result = failure_array_efficient(dna)
    
    # Output
    with open("output.txt", "w") as out_file:
        out_file.write(" ".join(map(str, failure_array_efficient)))


if __name__ == "__main__":
    main()