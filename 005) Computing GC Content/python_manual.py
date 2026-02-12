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

def calculate_gc_content(dna_string):
    """Calculate GC content percentage."""
    if not dna_string:
        return 0.0
    
    gc_count = dna_string.count('G') + dna_string.count('C')
    return (gc_count / len(dna_string)) * 100

def find_highest_gc_manual(fasta_string):
    """Find sequence with highest GC content manually."""
    sequences = parse_fasta_manual(fasta_string)
    
    highest_id = ""
    highest_gc = 0.0
    
    for seq_id, dna in sequences.items():
        gc_content = calculate_gc_content(dna)
        if gc_content > highest_gc:
            highest_gc = gc_content
            highest_id = seq_id
    
    return highest_id, highest_gc



path = "Dataset.txt"    # !!! Make sure the path is correct !!! 

with open(path, 'r') as file:
    data = file.read().strip()

id1, gc1 = find_highest_gc_manual(data)
print(f"Python Manual: {id1}\n{gc1:.6f}")
