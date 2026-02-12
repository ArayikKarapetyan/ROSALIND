def parse_set(set_str):
    """Parse set string like '{1, 2, 3, 4, 5}' into Python set."""
    # Remove braces and split by comma
    if set_str.startswith('{') and set_str.endswith('}'):
        set_str = set_str[1:-1]
    
    if not set_str.strip():
        return set()
    
    # Split and convert to integers
    elements = [int(x.strip()) for x in set_str.split(',')]
    return set(elements)

def format_set(s):
    """Format set as string with sorted elements."""
    if not s:
        return "{}"
    return "{" + ", ".join(str(x) for x in sorted(s)) + "}"

def set_operations(n, A, B):
    """Perform all set operations."""
    # Universal set U = {1, 2, ..., n}
    U = set(range(1, n + 1))
    
    # Perform operations
    union = A | B
    intersection = A & B
    A_minus_B = A - B
    B_minus_A = B - A
    A_complement = U - A
    B_complement = U - B
    
    return union, intersection, A_minus_B, B_minus_A, A_complement, B_complement

def main():
    # Read input from Dataset.txt
    with open("Dataset.txt", "r") as f:
        lines = [line.strip() for line in f if line.strip()]
    
    # Parse input
    n = int(lines[0])
    A = parse_set(lines[1])
    B = parse_set(lines[2])
    
    # Perform set operations
    results = set_operations(n, A, B)
    
    # Output results
    with open("output.txt", "w") as f:
        for result in results:
            f.write(format_set(result) + "\n")

if __name__ == "__main__":
    main()