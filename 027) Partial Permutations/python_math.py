import math

def partial_permutations_math(n, k, mod=1000000):
    """Calculate using math library (handles small numbers)."""
    if k > n:
        return 0
    
    # For small n, we can use math.factorial
    # But for n up to 100, we need modular arithmetic
    result = 1
    for i in range(n, n - k, -1):
        result = (result * i) % mod
    
    return result

with open("Dataset.txt", "r") as file:
    n, k = map(int, file.read().strip().split())
print(f"Python Math: {partial_permutations_math(n, k)}")