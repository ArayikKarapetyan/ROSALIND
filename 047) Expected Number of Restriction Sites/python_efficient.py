def expected_occurrences_efficient(n, s, gc_contents):
    """
    Efficient calculation with precomputation.
    """
    results = []
    length = len(s)
    
    # Pre-count GC and AT bases in s
    gc_count = s.count('G') + s.count('C')
    at_count = length - gc_count
    
    # For each GC-content
    for gc_content in gc_contents:
        # Probability of generating s at specific position
        # (gc_content/2)^gc_count * ((1-gc_content)/2)^at_count
        prob_s = ((gc_content / 2) ** gc_count) * (((1 - gc_content) / 2) ** at_count)
        
        # Expected occurrences
        expected = prob_s * (n - length + 1)
        results.append(expected)
    
    return results

def main():
    # Read input
    with open("Dataset.txt", "r") as f:
        lines = [line.strip() for line in f if line.strip()]
    
    n = int(lines[0])
    s = lines[1]
    gc_contents = list(map(float, lines[2].split()))
    
    results = expected_occurrences_efficient(n, s, gc_contents)
    print(" ".join(f"{x:.3f}" for x in results))

if __name__ == "__main__":
    main()