import numpy as np

def calculate_log_probabilities_numpy(dna_string, gc_contents):
    """Calculate using NumPy for vectorized operations."""
    # Convert to NumPy array
    gc_array = np.array(gc_contents)
    
    # Count nucleotides
    dna_np = np.array(list(dna_string))
    at_count = np.sum((dna_np == 'A') | (dna_np == 'T'))
    cg_count = len(dna_string) - at_count
    
    # Calculate log probabilities
    log_probs = at_count * np.log10((1 - gc_array) / 2) + cg_count * np.log10(gc_array / 2)
    
    return log_probs.tolist()

with open("Dataset.txt", "r") as file:
    lines = file.read().strip().split('\n')
    dna_string = lines[0].strip()
    gc_contents = list(map(float, lines[1].split()))
    
print("Python NumPy:", *calculate_log_probabilities_numpy(dna_string, gc_contents))