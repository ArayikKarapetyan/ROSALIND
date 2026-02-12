def parse_set_bits(set_str, n):
    """Parse set into bitmask representation."""
    if set_str.startswith('{') and set_str.endswith('}'):
        set_str = set_str[1:-1]
    
    if not set_str.strip():
        return 0
    
    mask = 0
    for x in set_str.split(','):
        num = int(x.strip())
        if 1 <= num <= n:
            mask |= 1 << (num - 1)  # Set bit at position (num-1)
    
    return mask

def format_set_from_mask(mask, n):
    """Convert bitmask to formatted set string."""
    elements = []
    for i in range(n):
        if mask & (1 << i):
            elements.append(str(i + 1))
    
    if not elements:
        return "{}"
    return "{" + ", ".join(elements) + "}"

def set_operations_bits(n, A_mask, B_mask):
    """Perform set operations using bitwise operations."""
    # Universal set mask: all bits from 0 to n-1 set
    U_mask = (1 << n) - 1
    
    # Bitwise operations
    union = A_mask | B_mask
    intersection = A_mask & B_mask
    A_minus_B = A_mask & ~B_mask
    B_minus_A = B_mask & ~A_mask
    A_complement = U_mask & ~A_mask  # Mask with U to ensure only n bits
    B_complement = U_mask & ~B_mask
    
    return union, intersection, A_minus_B, B_minus_A, A_complement, B_complement

def main():
    # Read input
    with open("Dataset.txt", "r") as f:
        lines = [line.strip() for line in f if line.strip()]
    
    n = int(lines[0])
    
    # Parse sets as bitmasks
    A_mask = parse_set_bits(lines[1], n)
    B_mask = parse_set_bits(lines[2], n)
    
    # Perform operations
    results = set_operations_bits(n, A_mask, B_mask)
    
    # Output results
    with open("output.txt", "w") as f:
        for mask in results:
            f.write(format_set_from_mask(mask, n) + "\n")

if __name__ == "__main__":
    main()