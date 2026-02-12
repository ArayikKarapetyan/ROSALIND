def generate_strings_iterative(alphabet, n):
    """
    Generate strings iteratively using stack instead of recursion.
    More memory efficient for large alphabets.
    """
    results = []
    stack = [""]
    
    while stack:
        current = stack.pop()
        if current:
            results.append(current)
        
        # Add next characters in reverse order to maintain lex order
        if len(current) < n:
            # Add in reverse to maintain correct order when popping from stack
            for char in reversed(alphabet):
                stack.append(current + char)
    
    return results

def main():
    # Read input
    with open("Dataset.txt", "r") as f:
        lines = f.readlines()
    
    alphabet = lines[0].strip().split()
    n = int(lines[1].strip())
    
    # Generate strings
    strings = generate_strings_iterative(alphabet, n)
    
    # Output result
    with open("output.txt", "w") as out_file:
        out_file.write(" ".join(map(str, strings)))

if __name__ == "__main__":
    main()