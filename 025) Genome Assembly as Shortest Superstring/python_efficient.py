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

def merge_strings(a, b, overlap):
    """Merge two strings given overlap length."""
    return a + b[overlap:]

# More efficient version with precomputed overlaps
def shortest_superstring_efficient(strings):
    """Efficient shortest superstring using precomputed overlaps."""
    n = len(strings)
    
    # Precompute all overlaps
    overlaps = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j:
                overlaps[i][j] = overlap_length(strings[i], strings[j])
    
    # Greedy merging
    current_strings = strings[:]
    current_indices = list(range(n))
    
    while len(current_strings) > 1:
        # Find maximum overlap among remaining strings
        max_overlap = -1
        best_i = -1
        best_j = -1
        
        for idx_i, i in enumerate(current_indices):
            for idx_j, j in enumerate(current_indices):
                if i == j:
                    continue
                if overlaps[i][j] > max_overlap:
                    max_overlap = overlaps[i][j]
                    best_i = idx_i
                    best_j = idx_j
        
        if max_overlap <= 0:
            # No overlaps, just concatenate all
            result = ''.join(current_strings)
            return result
        
        # Merge the best pair
        i_val = current_indices[best_i]
        j_val = current_indices[best_j]
        
        merged = merge_strings(current_strings[best_i], 
                               current_strings[best_j], 
                               overlaps[i_val][j_val])
        
        # Update lists
        # Remove larger index first
        if best_i > best_j:
            del current_strings[best_i]
            del current_indices[best_i]
            del current_strings[best_j]
            del current_indices[best_j]
        else:
            del current_strings[best_j]
            del current_indices[best_j]
            del current_strings[best_i]
            del current_indices[best_i]
        
        # Add merged string with new index
        current_strings.append(merged)
        current_indices.append(len(strings))  # New unique index
        strings.append(merged)  # Add to global list for overlap updates
        
        # Update overlap matrix for new string
        new_idx = len(strings) - 1
        # Add row and column for new string
        for row in overlaps:
            row.append(0)
        overlaps.append([0] * (len(strings)))
        
        # Calculate overlaps with new string
        for k in range(len(strings) - 1):
            overlaps[new_idx][k] = overlap_length(strings[new_idx], strings[k])
            overlaps[k][new_idx] = overlap_length(strings[k], strings[new_idx])
    
    return current_strings[0] if current_strings else ""




with open("Dataset.txt", "r") as file:
    strings = parse_fasta(file.read())


result2 = shortest_superstring_efficient(strings)
print(result2)
