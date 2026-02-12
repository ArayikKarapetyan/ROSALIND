def generate_lexicographic_strings_manual(alphabet, n):
    """Generate all strings of length n from alphabet in lex order."""
    from itertools import product
    
    # Generate all combinations
    strings = [''.join(combo) for combo in product(alphabet, repeat=n)]
    
    # Sort lexicographically (itertools.product already generates in order
    # for sorted input, but we'll sort to be sure)
    strings.sort()
    
    return strings


with open("Dataset.txt", "r") as file:
    lines = file.read().strip().split('\n')
    alphabet = lines[0].split()
    n = int(lines[1])


strings = generate_lexicographic_strings_manual(alphabet, n)
for s in strings:
    print(s)
