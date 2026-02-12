from math import factorial
from collections import Counter

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

def count_perfect_matchings_counter(rna):
    """Count perfect matchings using Counter for base counting."""
    counts = Counter(rna)
    
    if counts['A'] != counts['U'] or counts['C'] != counts['G']:
        return 0
    
    # Calculate using factorial
    return factorial(counts['A']) * factorial(counts['C'])

# Using product of ranges for factorial
def count_perfect_matchings_product(rna):
    """Calculate factorial using product."""
    from functools import reduce
    from operator import mul
    
    a_count = rna.count('A')
    c_count = rna.count('C')
    
    if rna.count('U') != a_count or rna.count('G') != c_count:
        return 0
    
    # au_matchings = a_count!
    au_matchings = reduce(mul, range(1, a_count + 1), 1) if a_count > 0 else 1
    
    # cg_matchings = c_count!
    cg_matchings = reduce(mul, range(1, c_count + 1), 1) if c_count > 0 else 1
    
    return au_matchings * cg_matchings

# For very large numbers (though n ≤ 80 so factorial(40) fits in Python int)
def count_perfect_matchings_bigint(rna):
    """Handle potentially large numbers (though not needed for n≤80)."""
    a_count = rna.count('A')
    c_count = rna.count('C')
    
    if rna.count('U') != a_count or rna.count('G') != c_count:
        return 0
    
    # Python automatically handles big integers
    return factorial(a_count) * factorial(c_count)

with open("Dataset.txt", "r") as file:
    rna = parse_fasta(file.read())

print("\nPython Counter Solution:")
print(f"Perfect matchings: {count_perfect_matchings_counter(rna)}")

print("\nPython Product Solution:")
print(f"Perfect matchings: {count_perfect_matchings_product(rna)}")