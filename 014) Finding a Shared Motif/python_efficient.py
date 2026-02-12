
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

# More efficient using suffix-based approach
def longest_common_substring_suffix(sequences):
    """Find LCS using suffix comparison."""
    if not sequences:
        return ""
    
    first_seq = sequences[0]
    n = len(first_seq)
    
    longest = ""
    
    # Generate all suffixes of first sequence
    for i in range(n):
        # For each starting position, try increasing lengths
        for length in range(1, n - i + 1):
            if length <= len(longest):
                continue  # No need to check shorter substrings
            
            substr = first_seq[i:i + length]
            
            if all(substr in seq for seq in sequences[1:]):
                longest = substr
            else:
                break  # If this length fails, longer will also fail
    
    return longest


with open("Dataset.txt", "r") as file:
    fasta_data = file.read()
sequences = parse_fasta_manual(fasta_data)


print(f"Python Suffix: {longest_common_substring_suffix(sequences)}")
