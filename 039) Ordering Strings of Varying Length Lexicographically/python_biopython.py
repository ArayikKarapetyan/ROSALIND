def generate_strings_biopython_style(alphabet, n):
    """
    Generate strings with BioPython-like approach.
    Not actually using BioPython, but structured similarly.
    """
    from itertools import product
    
    results = []
    
    # Generate strings of each length from 1 to n
    for length in range(1, n + 1):
        # Generate all combinations of given length
        for combo in product(alphabet, repeat=length):
            results.append(''.join(combo))
    
    # Now we need to sort according to the special lexicographic rule
    # where shorter strings come before longer strings with same prefix
    
    # Custom comparator function
    def compare_strings(a, b):
        """Compare two strings according to problem's rules."""
        min_len = min(len(a), len(b))
        
        # Compare character by character
        for i in range(min_len):
            if alphabet.index(a[i]) < alphabet.index(b[i]):
                return -1
            elif alphabet.index(a[i]) > alphabet.index(b[i]):
                return 1
        
        # If all compared characters are equal
        if len(a) < len(b):
            return -1
        elif len(a) > len(b):
            return 1
        else:
            return 0
    
    # Sort using custom comparator
    results.sort(key=lambda x: [alphabet.index(c) for c in x] + [len(x)])
    
    return results

def main():
    # Read input
    with open("Dataset.txt", "r") as f:
        lines = f.readlines()
    
    alphabet = lines[0].strip().split()
    n = int(lines[1].strip())
    
    # Generate strings
    strings = generate_strings_biopython_style(alphabet, n)
    
    # Output result
    with open("output.txt", "w") as out_file:
        out_file.write(" ".join(map(str, strings)))

if __name__ == "__main__":
    main()