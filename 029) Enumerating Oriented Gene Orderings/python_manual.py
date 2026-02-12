import itertools

def generate_signed_permutations_manual(n):
    """Generate all signed permutations manually."""
    if n <= 0:
        return 0, []
    
    numbers = list(range(1, n + 1))
    total_count = 0
    permutations = []
    
    # Generate all permutations of numbers
    for perm in itertools.permutations(numbers):
        # Generate all sign combinations (2^n possibilities)
        for signs in itertools.product([-1, 1], repeat=n):
            # Apply signs to permutation
            signed_perm = tuple(perm[i] * signs[i] for i in range(n))
            permutations.append(signed_perm)
            total_count += 1
    
    return total_count, permutations

with open("Dataset.txt", "r") as file:
    n = int(file.read().strip())

count, permutations = generate_signed_permutations_manual(n)
print(count)
for perm in permutations:
    print(" ".join(str(x) for x in perm))



