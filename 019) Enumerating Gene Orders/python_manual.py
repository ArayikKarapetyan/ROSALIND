from itertools import permutations

def generate_permutations_manual(n):
    """Generate all permutations of 1..n manually."""
    
    # Generate all permutations
    perms = list(permutations(range(1, n + 1)))
    
    # Format output
    result = []
    result.append(str(len(perms)))
    
    for perm in perms:
        result.append(' '.join(map(str, perm)))
    
    return '\n'.join(result)

with open("Dataset.txt", "r") as file:
    n = int(file.read().strip())

print("Python Manual (itertools):")
print(generate_permutations_manual(n))

