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
    
    return list(sequences.values())

def longest_common_substring_manual(sequences):
    """Find longest common substring using naive approach."""
    if not sequences:
        return ""
    
    # Start with first sequence
    first_seq = sequences[0]
    n = len(first_seq)
    
    longest = ""
    
    # Try all possible substrings of first sequence
    for i in range(n):
        for j in range(i + 1, n + 1):
            substr = first_seq[i:j]
            
            # Check if substring is in all sequences
            if all(substr in seq for seq in sequences[1:]):
                if len(substr) > len(longest):
                    longest = substr
    
    return longest

with open("Dataset.txt", "r") as file:
    fasta_data = file.read()
sequences = parse_fasta_manual(fasta_data)

print(f"Python Manual: {longest_common_substring_manual(sequences)}")
