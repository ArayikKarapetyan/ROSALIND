from collections import Counter
import itertools

def spectral_convolution_efficient(S1, S2):
    """
    Efficient spectral convolution using Counter.
    """
    # Generate all differences
    differences = (round(s1 - s2, 5) for s1 in S1 for s2 in S2)
    
    # Count multiplicities
    return Counter(differences)

def find_max_multiplicity_efficient(convolution):
    """
    Find maximum multiplicity efficiently.
    """
    if not convolution:
        return 0, 0.0
    
    # Find max multiplicity
    max_mult = max(convolution.values())
    
    # Find all x with this multiplicity
    candidates = [abs(x) for x, mult in convolution.items() if mult == max_mult]
    
    # Return any candidate (first one)
    return max_mult, candidates[0]

def main():
    # Read input
    with open("Dataset.txt", "r") as f:
        lines = f.readlines()
    
    S1 = list(map(float, lines[0].strip().split()))
    S2 = list(map(float, lines[1].strip().split()))
    
    # Compute convolution
    conv = spectral_convolution_efficient(S1, S2)
    
    # Find results
    max_mult, x_val = find_max_multiplicity_efficient(conv)
    
    # Output
    print(max_mult)
    print(f"{x_val:.5f}")

if __name__ == "__main__":
    main()