import numpy as np
import itertools
import math

def generate_signed_permutations_numpy(n):
    """Generate signed permutations using NumPy."""
    numbers = np.arange(1, n + 1)
    permutations = []
    
    # Generate all permutations
    for perm in itertools.permutations(numbers):
        perm_array = np.array(perm)
        
        # Generate all sign combinations
        for signs in itertools.product([-1, 1], repeat=n):
            signed = perm_array * signs
            permutations.append(tuple(signed))
    
    total = math.factorial(n) * (2 ** n)
    return total, permutations


with open("Dataset.txt", "r") as file:
    n = int(file.read().strip())

count, permutations = generate_signed_permutations_numpy(n)
print(count)
for perm in permutations:
    print(" ".join(str(x) for x in perm))


