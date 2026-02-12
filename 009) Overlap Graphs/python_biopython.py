from Bio import SeqIO
from io import StringIO



def build_overlap_graph_biopython(fasta_string, k=3):
    """Build overlap graph using BioPython."""
    edges = []
    
    # Parse FASTA from string
    fasta_io = StringIO(fasta_string)
    records = list(SeqIO.parse(fasta_io, "fasta"))
    
    # Convert to dictionary
    sequences = {record.id: str(record.seq) for record in records}
    
    # Build graph
    ids = list(sequences.keys())
    for i in range(len(ids)):
        for j in range(len(ids)):
            if i == j:
                continue
            
            s_id = ids[i]
            t_id = ids[j]
            s_seq = sequences[s_id]
            t_seq = sequences[t_id]
            
            if s_seq[-k:] == t_seq[:k]:
                edges.append((s_id, t_id))
    
    return edges

with open("Dataset.txt", "r") as file:
    data = file.read()

edges3 = build_overlap_graph_biopython(data)
print("\nPython BioPython:")
for edge in edges3:
    print(f"{edge[0]} {edge[1]}")
