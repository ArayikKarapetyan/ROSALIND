def expected_occurrences(n, s, gc_contents):
    """
    Calculate expected number of times substring s appears in random DNA of length n
    for each GC-content value.
    """
    results = []
    length = len(s)
    
    # Calculate probability of generating s for each GC-content
    for gc_content in gc_contents:
        # Probability of generating s at a specific position
        prob_s = 1.0
        for base in s:
            if base in 'GC':
                prob_s *= gc_content / 2
            else:  # A or T
                prob_s *= (1 - gc_content) / 2
        
        # Expected number of occurrences
        # For each possible starting position (n - length + 1 positions)
        expected = prob_s * (n - length + 1)
        results.append(expected)
    
    return results

def main():
    # Read input from Dataset.txt
    with open("Dataset.txt", "r") as f:
        lines = [line.strip() for line in f if line.strip()]
    
    # Parse input
    n = int(lines[0])           # Length of random DNA string
    s = lines[1]               # Substring to search for
    gc_contents = list(map(float, lines[2].split()))  # GC-content values
    
    # Calculate expected occurrences
    results = expected_occurrences(n, s, gc_contents)
    
    # Output with 3 decimal places
    print(" ".join(f"{x:.3f}" for x in results))

if __name__ == "__main__":
    main()