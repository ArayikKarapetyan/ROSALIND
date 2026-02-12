def count_nucleotides_dna(dna_string):
    """
    Count occurrences of A, C, G, T in a DNA string
    
    Args:
        dna_string (str): DNA sequence
    
    Returns:
        tuple: Counts of A, C, G, T in that order
    """
    # Simple counting using dictionary
    counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    
    for nucleotide in dna_string:
        if nucleotide in counts:
            counts[nucleotide] += 1
    
    return counts['A'], counts['C'], counts['G'], counts['T']

# Alternative concise version
def count_nucleotides_dna_concise(dna_string):
    """Concise version using count() method"""
    return (
        dna_string.count('A'),
        dna_string.count('C'),
        dna_string.count('G'),
        dna_string.count('T')
    )


# Usage

path = "Dataset.txt"    # !!! Make sure the path is correct !!! 

with open(path, 'r') as file:
        data = file.read().strip()
a_count, c_count, g_count, t_count = count_nucleotides_dna(data)
print(f"{a_count} {c_count} {g_count} {t_count}")


# result = count_nucleotides_dna_concise(dna)
# print(f"{result[0]} {result[1]} {result[2]} {result[3]}")