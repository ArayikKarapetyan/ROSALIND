def spectral_convolution(S1, S2):
    """
    Compute spectral convolution S1 ⊖ S2.
    Returns dictionary of differences and their multiplicities.
    """
    convolution = {}
    
    for s1 in S1:
        for s2 in S2:
            diff = round(s1 - s2, 5)  # Round to 5 decimal places
            if diff not in convolution:
                convolution[diff] = 0
            convolution[diff] += 1
    
    return convolution

def find_max_multiplicity(convolution):
    """
    Find the largest multiplicity and corresponding x value.
    Returns (max_multiplicity, x_value).
    """
    max_multiplicity = 0
    best_x = 0
    
    for x, multiplicity in convolution.items():
        if multiplicity > max_multiplicity:
            max_multiplicity = multiplicity
            best_x = abs(x)  # Return absolute value
        elif multiplicity == max_multiplicity:
            # If tie, choose larger absolute value? Or any...
            # Problem says "you may return any such value"
            if abs(x) > abs(best_x):
                best_x = abs(x)
    
    return max_multiplicity, best_x

def main():
    # Read input from Dataset.txt
    with open("Dataset.txt", "r") as f:
        lines = f.readlines()
    
    # Parse the two multisets
    S1 = list(map(float, lines[0].strip().split()))
    S2 = list(map(float, lines[1].strip().split()))
    
    # Compute spectral convolution
    conv = spectral_convolution(S1, S2)
    
    # Find maximum multiplicity and corresponding x
    max_mult, x_val = find_max_multiplicity(conv)
    
    # Output results
    print(max_mult)
    print(f"{x_val:.5f}")

if __name__ == "__main__":
    main()