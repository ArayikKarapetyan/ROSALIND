def compute_failure_array(s):
    """
    Compute failure array (prefix function) for string s.
    P[i] = length of longest proper prefix that is also a suffix for s[0:i+1]
    """
    n = len(s)
    P = [0] * n
    P[0] = 0  # By convention
    
    for i in range(1, n):
        j = P[i-1]
        while j > 0 and s[i] != s[j]:
            j = P[j-1]
        if s[i] == s[j]:
            j += 1
        P[i] = j
    
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
    failure_array = compute_failure_array(dna)
    
    # Output result
    with open("output.txt", "w") as out_file:
        out_file.write(" ".join(map(str, failure_array)))




if __name__ == "__main__":
    main()