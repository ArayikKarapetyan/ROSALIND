def reverse_complement_manual(dna_string):
    """Compute reverse complement manually."""
    complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    # Reverse the string and complement each character
    reverse_comp = ''.join(complement_map[base] for base in reversed(dna_string))
    return reverse_comp

# Alternative using translation table
def reverse_complement_translate(dna_string):
    """Compute reverse complement using translation table."""
    complement_table = str.maketrans('ATCG', 'TAGC')
    # Reverse, then translate
    return dna_string[::-1].translate(complement_table)




path = "Dataset.txt"    # !!! Make sure the path is correct !!! 

with open(path, 'r') as file:
    data = file.read().strip()

print("Python Manual:", reverse_complement_manual(data))