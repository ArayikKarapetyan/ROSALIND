def partial_permutations_manual(n, k, mod=1000000):
    """Calculate partial permutations P(n,k) manually."""
    if k > n:
        return 0
    
    result = 1
    # P(n,k) = n! / (n-k)! = n × (n-1) × ... × (n-k+1)
    for i in range(n, n - k, -1):
        result = (result * i) % mod
    
    return result

with open("Dataset.txt", "r") as file:
    n, k = map(int, file.read().strip().split())
print(f"Python Manual: {partial_permutations_manual(n, k)}")
