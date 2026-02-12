def generate_strings_varying_length(alphabet, n):
    """
    Generate all strings of length at most n from given alphabet,
    ordered lexicographically according to custom alphabet order.
    """
    results = []
    
    def backtrack(current):
        if current:
            results.append(current)
        if len(current) < n:
            for char in alphabet:
                backtrack(current + char)
    
    backtrack("")
    return results

def main():
    # Read input from Dataset.txt
    with open("Dataset.txt", "r") as f:
        lines = f.readlines()
    
    # Parse alphabet (split by whitespace)
    alphabet = lines[0].strip().split()
    
    # Parse n
    n = int(lines[1].strip())
    
    # Generate strings
    strings = generate_strings_varying_length(alphabet, n)
    
    # Output result
    with open("output.txt", "w") as out_file:
        out_file.write(" ".join(map(str, strings)))

if __name__ == "__main__":
    main()