def parse_fasta(fasta_string):
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

def overlap_length(a, b, min_overlap=1):
    """Calculate maximum overlap where end of a matches start of b."""
    max_len = min(len(a), len(b))
    
    # Try larger overlaps first (problem says > half length)
    for overlap in range(max_len, min_overlap - 1, -1):
        if a[-overlap:] == b[:overlap]:
            return overlap
    return 0

def find_max_overlap_pair(strings):
    """Find pair with maximum overlap."""
    max_overlap = 0
    best_pair = (-1, -1)  # (i, j) where end of i overlaps start of j
    
    for i in range(len(strings)):
        for j in range(len(strings)):
            if i == j:
                continue
            
            overlap = overlap_length(strings[i], strings[j])
            if overlap > max_overlap:
                max_overlap = overlap
                best_pair = (i, j)
    
    return best_pair, max_overlap

def merge_strings(a, b, overlap):
    """Merge two strings given overlap length."""
    return a + b[overlap:]

def shortest_superstring_manual(strings):
    """Find shortest superstring using greedy overlap merging."""
    # Make a copy to work with
    current_strings = strings[:]
    
    # Continue until only one string remains
    while len(current_strings) > 1:
        # Find best pair to merge
        (i, j), overlap = find_max_overlap_pair(current_strings)
        
        if overlap == 0:
            # No overlap found, just concatenate
            merged = current_strings[i] + current_strings[j]
        else:
            # Merge with overlap
            merged = merge_strings(current_strings[i], current_strings[j], overlap)
        
        # Remove the two strings and add the merged one
        # Remove the larger index first to avoid index shifting issues
        idx1, idx2 = max(i, j), min(i, j)
        del current_strings[idx1]
        del current_strings[idx2]
        current_strings.append(merged)
    
    return current_strings[0] if current_strings else ""



with open("Dataset.txt", "r") as file:
    strings = parse_fasta(file.read())

result = shortest_superstring_manual(strings)
print(result)

