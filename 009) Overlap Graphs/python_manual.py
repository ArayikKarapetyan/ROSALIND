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

def build_overlap_graph_manual(sequences, k=3):
    """Build overlap graph manually."""
    edges = []
    ids = list(sequences.keys())
    
    for i in range(len(ids)):
        for j in range(len(ids)):
            if i == j:
                continue  # Skip self-loops
            
            s_id = ids[i]
            t_id = ids[j]
            s_seq = sequences[s_id]
            t_seq = sequences[t_id]
            
            # Check if suffix of s matches prefix of t
            if s_seq[-k:] == t_seq[:k]:
                edges.append((s_id, t_id))
    
    return edges


with open("Dataset.txt", "r") as file:
    data = file.read()

edges1 = build_overlap_graph_manual(parse_fasta_manual(data))
print("Python Manual:")
for edge in edges1:
    print(f"{edge[0]} {edge[1]}")

