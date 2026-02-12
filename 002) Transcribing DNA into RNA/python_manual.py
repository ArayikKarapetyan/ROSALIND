def transcribe_dna_to_rna(dna_string):
    """Transcribe DNA to RNA by replacing T with U."""
    # Simple string replacement
    return dna_string.replace('T', 'U')

path = "Dataset.txt"    # !!! Make sure the path is correct !!! 

with open(path, 'r') as file:
    data = file.read().strip()
# Method 1: Python from scratch
print("Python (from scratch):", transcribe_dna_to_rna(data))
