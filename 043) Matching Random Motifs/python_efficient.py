import math

def probability_at_least_one_match_efficient(N, gc_content, motif):
    """
    Efficient calculation using logarithms to avoid underflow for large N.
    """
    # Calculate log probability of a single match
    log_prob_single = 0.0
    for base in motif:
        if base in 'GC':
            log_prob_single += math.log(gc_content / 2)
        else:  # A or T
            log_prob_single += math.log((1 - gc_content) / 2)
    
    # Convert back to probability
    prob_single = math.exp(log_prob_single)
    
    # For large N, use logarithms to compute (1-p)^N
    if prob_single == 1.0:
        return 1.0
    elif prob_single == 0.0:
        return 0.0
    else:
        # log(1-p) for p small: use math.log1p for better precision
        log_prob_not_single = math.log1p(-prob_single)
        
        # (1-p)^N = exp(N * log(1-p))
        log_prob_none_match = N * log_prob_not_single
        prob_none_match = math.exp(log_prob_none_match)
        
        return 1 - prob_none_match

def main():
    # Read input
    with open("Dataset.txt", "r") as f:
        lines = f.readlines()
    
    # Parse input
    N_str, gc_str = lines[0].strip().split()
    N = int(N_str)
    gc_content = float(gc_str)
    motif = lines[1].strip()
    
    # Calculate probability
    prob = probability_at_least_one_match_efficient(N, gc_content, motif)
    
    # Output
    print(f"{prob:.3f}")

if __name__ == "__main__":
    main()