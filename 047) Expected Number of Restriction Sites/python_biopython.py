def expected_occurrences_biopython(n, s, gc_contents):
    """
    BioPython-style implementation for expected occurrences.
    """
    results = []
    length = len(s)
    
    # Count base types
    from collections import Counter
    counts = Counter(s)
    gc_count = counts.get('G', 0) + counts.get('C', 0)
    at_count = length - gc_count
    
    for gc_content in gc_contents:
        # Probability formula
        prob_s = ((gc_content / 2) ** gc_count) * (((1 - gc_content) / 2) ** at_count)
        
        # Expected value: linearity of expectation
        expected = prob_s * (n - length + 1)
        results.append(expected)
    
    return results

def main():
    # Read input (simpler format for this problem)
    with open("Dataset.txt", "r") as f:
        lines = [line.strip() for line in f if line.strip()]
    
    n = int(lines[0])
    s = lines[1]
    gc_contents = list(map(float, lines[2].split()))
    
    results = expected_occurrences_biopython(n, s, gc_contents)
    print(" ".join(f"{x:.3f}" for x in results))

if __name__ == "__main__":
    main()