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

def compute_profile_and_consensus_manual(sequences):
    """Compute profile matrix and consensus string manually."""
    # Get sequences as list
    seq_list = list(sequences.values())
    
    if not seq_list:
        return "", {}
    
    n = len(seq_list[0])  # Sequence length
    m = len(seq_list)     # Number of sequences
    
    # Initialize profile matrix (A, C, G, T)
    profile = {
        'A': [0] * n,
        'C': [0] * n,
        'G': [0] * n,
        'T': [0] * n
    }
    
    # Fill profile matrix
    for seq in seq_list:
        for j in range(n):
            nucleotide = seq[j]
            profile[nucleotide][j] += 1
    
    # Compute consensus string
    consensus_chars = []
    for j in range(n):
        max_count = -1
        max_nuc = 'A'
        
        for nuc in ['A', 'C', 'G', 'T']:
            if profile[nuc][j] > max_count:
                max_count = profile[nuc][j]
                max_nuc = nuc
        
        consensus_chars.append(max_nuc)
    
    consensus = ''.join(consensus_chars)
    
    return consensus, profile



with open("Dataset.txt", "r") as file:
    fasta_data = file.read()

sequences = parse_fasta_manual(fasta_data)
consensus, profile = compute_profile_and_consensus_manual(sequences)

print("Consensus:", consensus)
print("\nProfile Matrix:")
for nuc in ['A', 'C', 'G', 'T']:
    print(f"{nuc}: {' '.join(map(str, profile[nuc]))}")