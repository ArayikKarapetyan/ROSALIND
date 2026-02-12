def parse_fasta(fasta_string):
    """Parse FASTA format manually."""
    lines = fasta_string.strip().split('\n')
    sequences = {}
    current_id = ""
    
    for line in lines:
        if line.startswith('>'):
            current_id = line[1:].strip()
            sequences[current_id] = ""
        else:
            sequences[current_id] += line.strip()
    
    return list(sequences.values())[0]  # Return first sequence

def count_perfect_matchings_manual(rna):
    """Count perfect matchings for RNA using combinatorial approach."""
    from math import factorial
    
    # Count occurrences of each base
    a_count = rna.count('A')
    u_count = rna.count('U')
    c_count = rna.count('C')
    g_count = rna.count('G')
    
    # Verify the condition
    if a_count != u_count or c_count != g_count:
        return 0
    
    # Perfect matchings are independent for AU and CG pairs
    # For n AU pairs: number of perfect matchings = n!
    # For m CG pairs: number of perfect matchings = m!
    # Total = n! × m!
    
    au_matchings = factorial(a_count)
    cg_matchings = factorial(c_count)
    
    return au_matchings * cg_matchings

with open("Dataset.txt", "r") as file:
    rna = parse_fasta(file.read())
print(f"Perfect matchings: {count_perfect_matchings_manual(rna)}")

