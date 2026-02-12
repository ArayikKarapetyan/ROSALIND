def parse_set_biopython(set_str):
    """Parse set string, BioPython style."""
    if set_str.startswith('{') and set_str.endswith('}'):
        set_str = set_str[1:-1]
    
    elements = set()
    if set_str.strip():
        for item in set_str.split(','):
            elements.add(int(item.strip()))
    
    return elements

def set_operations_biopython(n, A, B):
    """Perform set operations in BioPython style."""
    U = set(range(1, n + 1))
    
    results = [
        A | B,          # union
        A & B,          # intersection
        A - B,          # A minus B
        B - A,          # B minus A
        U - A,          # complement of A
        U - B           # complement of B
    ]
    
    return results

def format_set_biopython(s):
    """Format set with sorted elements."""
    return "{" + ", ".join(str(x) for x in sorted(s)) + "}"

def main():
    # Read input
    with open("Dataset.txt", "r") as f:
        lines = [line.strip() for line in f if line.strip()]
    
    n = int(lines[0])
    A = parse_set_biopython(lines[1])
    B = parse_set_biopython(lines[2])
    
    results = set_operations_biopython(n, A, B)
    
    with open("output.txt", "w") as f:
        for mask in results:
            f.write(format_set_biopython(mask, n) + "\n")

if __name__ == "__main__":
    main()