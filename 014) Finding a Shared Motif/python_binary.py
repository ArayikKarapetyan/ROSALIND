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

# Using binary search for length
def longest_common_substring_binary(sequences):
    """Find LCS using binary search on length."""
    if not sequences:
        return ""
    
    first_seq = sequences[0]
    n = len(first_seq)
    
    # Binary search on possible lengths
    low, high = 0, n
    result = ""
    
    while low <= high:
        mid = (low + high) // 2
        
        # Check if common substring of length mid exists
        found = False
        current = ""
        
        for i in range(n - mid + 1):
            substr = first_seq[i:i + mid]
            if all(substr in seq for seq in sequences[1:]):
                found = True
                current = substr
                break
        
        if found:
            result = current
            low = mid + 1  # Try longer
        else:
            high = mid - 1  # Try shorter
    
    return result


with open("Dataset.txt", "r") as file:
    fasta_data = file.read()
sequences = parse_fasta_manual(fasta_data)

print(f"Python Binary: {longest_common_substring_binary(sequences)}")
