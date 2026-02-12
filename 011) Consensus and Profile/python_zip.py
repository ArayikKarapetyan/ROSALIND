def parse_fasta_manual(fasta_string):
    """Parse FASTA format manually."""
    sequences = {}
    current_id = ""
    
    for line in fasta_string.strip().split('\n'):
        if line.startswith('>'):
            current_id = line[1:].strip()
            sequences[current_id] = ""
        else:
            sequences[current_id] += line.strip()
    
    return sequences
# Alternative using numpy-like approach
def compute_profile_and_consensus_zip(sequences):
    """Compute using zip for transposition."""
    seq_list = list(sequences.values())
    
    if not seq_list:
        return "", {}
    
    n = len(seq_list[0])
    
    # Transpose: get list of columns
    columns = list(zip(*seq_list))
    
    # Initialize profile
    profile = {
        'A': [0] * n,
        'C': [0] * n,
        'G': [0] * n,
        'T': [0] * n
    }
    
    # Fill profile
    for j, column in enumerate(columns):
        for nucleotide in column:
            profile[nucleotide][j] += 1
    
    # Compute consensus
    consensus_chars = []
    for j in range(n):
        max_nuc = max(['A', 'C', 'G', 'T'], key=lambda x: profile[x][j])
        consensus_chars.append(max_nuc)
    
    consensus = ''.join(consensus_chars)
    
    return consensus, profile


with open("Dataset.txt", "r") as file:
    fasta_data = file.read()

sequences = parse_fasta_manual(fasta_data)
consensus, profile = compute_profile_and_consensus_zip(sequences)

print("Consensus:", consensus)
print("\nProfile Matrix:")
for nuc in ['A', 'C', 'G', 'T']:
    print(f"{nuc}: {' '.join(map(str, profile[nuc]))}")