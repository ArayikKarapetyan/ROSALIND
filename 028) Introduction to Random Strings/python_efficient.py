import math

# More efficient version using counts
def calculate_log_probabilities_efficient(dna_string, gc_contents):
    """Efficient calculation using nucleotide counts."""
    results = []
    
    # Pre-count nucleotides
    at_count = 0
    cg_count = 0
    for nucleotide in dna_string:
        if nucleotide in ['A', 'T']:
            at_count += 1
        else:  # 'C' or 'G'
            cg_count += 1
    
    for gc_content in gc_contents:
        prob_at = (1 - gc_content) / 2
        prob_cg = gc_content / 2
        
        # Log probability = at_count * log(prob_at) + cg_count * log(prob_cg)
        log_prob = at_count * math.log10(prob_at) + cg_count * math.log10(prob_cg)
        results.append(log_prob)
    
    return results

with open("Dataset.txt", "r") as file:
    lines = file.read().strip().split('\n')
    dna_string = lines[0].strip()
    gc_contents = list(map(float, lines[1].split()))
    
print("Python Efficient:", *calculate_log_probabilities_efficient(dna_string, gc_contents))
