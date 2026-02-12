import itertools
import math

def generate_signed_permutations(n):
    numbers = list(range(1, n + 1))
    all_permutations = []
    
    for perm in itertools.permutations(numbers):
        for signs in itertools.product([-1, 1], repeat=n):
            signed = tuple(p * s for p, s in zip(perm, signs))
            all_permutations.append(signed)
    
    total = math.factorial(n) * (2 ** n)
    return total, all_permutations


with open("Dataset.txt", "r") as file:
    n = int(file.read().strip())

count, permutations = generate_signed_permutations(n)

print(count)
for perm in permutations:
    print(" ".join(map(str, perm)))
