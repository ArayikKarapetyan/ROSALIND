import math

def probability_at_least_one_match_biopython(N, gc_content, motif):
    """
    BioPython-style implementation for probability calculation.
    """
    # Count GC and AT bases in motif
    gc_count = motif.count('G') + motif.count('C')
    at_count = len(motif) - gc_count
    
    # Probability of generating the exact motif in one trial
    # P(motif) = (gc_content/2)^gc_count * ((1-gc_content)/2)^at_count
    prob_single = ((gc_content / 2) ** gc_count) * (((1 - gc_content) / 2) ** at_count)
    
    # Probability that none of N trials succeed
    if prob_single == 0:
        return 0.0
    elif prob_single == 1:
        return 1.0
    else:
        prob_none = (1 - prob_single) ** N
        return 1 - prob_none

def main():
    # Read input
    with open("Dataset.txt", "r") as f:
        lines = f.readlines()
    
    # Parse
    N_str, gc_str = lines[0].strip().split()
    N = int(N_str)
    gc_content = float(gc_str)
    motif = lines[1].strip()
    
    # Calculate
    prob = probability_at_least_one_match_biopython(N, gc_content, motif)
    
    # Output
    print(f"{prob:.3f}")

if __name__ == "__main__":
    main()