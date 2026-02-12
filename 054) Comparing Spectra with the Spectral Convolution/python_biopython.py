from collections import Counter

def spectral_convolution_biopython(S1, S2):
    """
    BioPython-style spectral convolution calculation.
    """
    # Round differences to handle floating-point precision
    PRECISION = 5
    
    conv = Counter()
    
    for s1 in S1:
        for s2 in S2:
            diff = round(s1 - s2, PRECISION)
            conv[diff] += 1
    
    return conv

def analyze_convolution(convolution):
    """
    Analyze convolution to find maximum multiplicity.
    """
    if not convolution:
        return 0, 0.0
    
    # Find maximum multiplicity
    max_multiplicity = max(convolution.values())
    
    # Collect all x with this multiplicity
    max_x_values = [abs(x) for x, mult in convolution.items() 
                   if mult == max_multiplicity]
    
    # Return any (first) x value with max multiplicity
    return max_multiplicity, max_x_values[0]

def main():
    # Read input
    with open("Dataset.txt", "r") as f:
        lines = f.readlines()
    
    # Parse multisets
    S1 = [float(x) for x in lines[0].strip().split()]
    S2 = [float(x) for x in lines[1].strip().split()]
    
    # Compute convolution
    conv = spectral_convolution_biopython(S1, S2)
    
    # Analyze results
    max_mult, x_val = analyze_convolution(conv)
    
    # Output
    print(max_mult)
    print(f"{x_val:.5f}")

if __name__ == "__main__":
    main()