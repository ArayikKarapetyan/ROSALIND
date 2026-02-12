import math

def calculate_log_probabilities_manual(dna_string, gc_contents):
    """Calculate log10 probabilities for random DNA string matching."""
    results = []
    
    # Count nucleotides in the DNA string
    a_t_count = dna_string.count('A') + dna_string.count('T')
    c_g_count = dna_string.count('C') + dna_string.count('G')
    total_length = len(dna_string)
    
    for gc_content in gc_contents:
        # Calculate probabilities
        # For GC-content = x:
        # P(G) = P(C) = x/2
        # P(A) = P(T) = (1-x)/2
        prob_at = (1 - gc_content) / 2
        prob_cg = gc_content / 2
        
        # Probability of exact match = ∏ P(nucleotide)
        # Use logarithms to avoid underflow
        log_prob = 0.0
        
        for nucleotide in dna_string:
            if nucleotide in ['A', 'T']:
                log_prob += math.log10(prob_at)
            else:  # 'C' or 'G'
                log_prob += math.log10(prob_cg)
        
        results.append(log_prob)
    
    return results

with open("Dataset.txt", "r") as file:
    lines = file.read().strip().split('\n')
    dna_string = lines[0].strip()
    gc_contents = list(map(float, lines[1].split()))
    
print("Python Manual:", *calculate_log_probabilities_manual(dna_string, gc_contents))
