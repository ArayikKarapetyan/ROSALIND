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

# More efficient version using dictionary of suffixes
def build_overlap_graph_efficient(sequences, k=3):
    """Build overlap graph more efficiently."""
    edges = []
    
    # Create dictionaries for quick lookup
    prefix_map = {}  # prefix -> list of sequence IDs
    suffix_map = {}  # suffix -> list of sequence IDs
    
    for seq_id, sequence in sequences.items():
        if len(sequence) >= k:
            prefix = sequence[:k]
            suffix = sequence[-k:]
            
            if prefix not in prefix_map:
                prefix_map[prefix] = []
            prefix_map[prefix].append(seq_id)
            
            if suffix not in suffix_map:
                suffix_map[suffix] = []
            suffix_map[suffix].append(seq_id)
    
    # Find edges where suffix matches prefix
    for suffix, s_ids in suffix_map.items():
        if suffix in prefix_map:
            t_ids = prefix_map[suffix]
            for s_id in s_ids:
                for t_id in t_ids:
                    if s_id != t_id:  # Avoid self-loops
                        edges.append((s_id, t_id))
    
    return edges


with open("Dataset.txt", "r") as file:
    data = file.read()

edges2 = build_overlap_graph_efficient(parse_fasta_manual(data))
print("\nPython Efficient:")
for edge in edges2:
    print(f"{edge[0]} {edge[1]}")

