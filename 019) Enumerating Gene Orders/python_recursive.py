# Using recursion with generator
def generate_permutations_recursive(n):
    """Generate permutations using recursion with yield."""
    def permute(arr, start=0):
        if start == len(arr):
            yield arr[:]
        else:
            for i in range(start, len(arr)):
                arr[start], arr[i] = arr[i], arr[start]
                yield from permute(arr, start + 1)
                arr[start], arr[i] = arr[i], arr[start]
    
    arr = list(range(1, n + 1))
    perms = list(permute(arr))
    
    output = [str(len(perms))]
    for perm in perms:
        output.append(' '.join(map(str, perm)))
    
    return '\n'.join(output)

with open("Dataset.txt", "r") as file:
    n = int(file.read().strip())



print("\nPython Recursive:")
print(generate_permutations_recursive(n))
