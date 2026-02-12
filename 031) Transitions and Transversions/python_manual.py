def calculate_transition_transversion_manual(s1, s2):
    """
    Calculate transition/transversion ratio manually.
    
    Transitions: A↔G or C↔T (purine↔purine or pyrimidine↔pyrimidine)
    Transversions: All other substitutions
    """
    transitions = 0
    transversions = 0
    
    # Define transition pairs
    transition_pairs = {('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')}
    
    for base1, base2 in zip(s1, s2):
        if base1 == base2:
            continue  # No mutation
        
        if (base1, base2) in transition_pairs:
            transitions += 1
        else:
            transversions += 1
    
    # Avoid division by zero
    if transversions == 0:
        return float('inf')  # Only transitions occurred
    
    return transitions / transversions

def parse_fasta_manual(fasta_string):
    sequences = {}
    current_id = ""
    for line in fasta_string.strip().split('\n'):
        if line.startswith('>'):
            current_id = line[1:].strip()
            sequences[current_id] = ""
        else:
            sequences[current_id] += line.strip()
    return sequences

with open("Dataset.txt", "r") as file:
    fasta_data = file.read()

sequences = parse_fasta_manual(fasta_data)
s1 = list(sequences.values())[0]
s2 = list(sequences.values())[1]

ratio_manual = calculate_transition_transversion_manual(s1, s2)
print(f"Manual: {ratio_manual:.11f}")

