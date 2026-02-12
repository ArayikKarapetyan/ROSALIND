def probability_at_least_one_match(N, gc_content, motif):
    """
    Calculate probability that at least one of N random strings equals the motif.
    """
    # Calculate probability of a single random string matching the motif
    prob_single = 1.0
    for base in motif:
        if base in 'GC':
            prob_single *= gc_content / 2
        else:  # A or T
            prob_single *= (1 - gc_content) / 2
    
    # Probability that a single random string does NOT match
    prob_not_single = 1 - prob_single
    
    # Probability that N independent strings all do NOT match
    prob_none_match = prob_not_single ** N
    
    # Probability that at least one matches
    prob_at_least_one = 1 - prob_none_match
    
    return prob_at_least_one

def main():
    # Read input from Dataset.txt
    with open("Dataset.txt", "r") as f:
        lines = f.readlines()
    
    # Parse N and GC content from first line
    N_str, gc_str = lines[0].strip().split()
    N = int(N_str)
    gc_content = float(gc_str)
    
    # Parse motif from second line
    motif = lines[1].strip()
    
    # Calculate probability
    prob = probability_at_least_one_match(N, gc_content, motif)
    
    # Output with reasonable precision
    print(f"{prob:.3f}")

if __name__ == "__main__":
    main()