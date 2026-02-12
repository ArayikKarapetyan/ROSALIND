# Iterative implementation using heap's algorithm
def generate_permutations_iterative(n):
    """Generate permutations using Heap's algorithm."""
    # Initialize array with numbers 1..n
    arr = list(range(1, n + 1))
    results = [arr[:]]  # Copy of initial array
    
    # Heap's algorithm for generating permutations
    c = [0] * n
    i = 0
    
    while i < n:
        if c[i] < i:
            if i % 2 == 0:
                arr[0], arr[i] = arr[i], arr[0]
            else:
                arr[c[i]], arr[i] = arr[i], arr[c[i]]
            results.append(arr[:])
            c[i] += 1
            i = 0
        else:
            c[i] = 0
            i += 1
    
    output = [str(len(results))]
    for perm in results:
        output.append(' '.join(map(str, perm)))
    
    return '\n'.join(output)

with open("Dataset.txt", "r") as file:
    n = int(file.read().strip())



print("\nPython Iterative (Heap's):")
print(generate_permutations_iterative(n))