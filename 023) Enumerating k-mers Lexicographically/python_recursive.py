# Manual implementation using recursion
def generate_lexicographic_recursive(alphabet, n):
    """Generate strings using recursion."""
    results = []
    
    def backtrack(current):
        if len(current) == n:
            results.append(current)
            return
        
        for char in alphabet:
            backtrack(current + char)
    
    backtrack("")
    return results

# Using recursion with yield
def generate_lexicographic_recursive_yield(alphabet, n, prefix=""):
    """Recursive generator version."""
    if n == 0:
        yield prefix
    else:
        for char in alphabet:
            yield from generate_lexicographic_recursive_yield(alphabet, n-1, prefix + char)



with open("Dataset.txt", "r") as file:
    lines = file.read().strip().split('\n')
    alphabet = lines[0].split()
    n = int(lines[1])

strings_rec = generate_lexicographic_recursive(alphabet, n)
for s in strings_rec:
    print(s)
