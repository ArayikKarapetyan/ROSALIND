# Iterative implementation
def generate_lexicographic_iterative(alphabet, n):
    """Generate strings iteratively."""
    results = []
    
    # Start with list of empty strings
    current = [""]
    
    for position in range(n):
        next_level = []
        for string in current:
            for char in alphabet:
                next_level.append(string + char)
        current = next_level
    
    return current


with open("Dataset.txt", "r") as file:
    lines = file.read().strip().split('\n')
    alphabet = lines[0].split()
    n = int(lines[1])


strings_rec = generate_lexicographic_iterative(alphabet, n)
for s in strings_rec:
    print(s)
