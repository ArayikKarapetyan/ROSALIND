# More efficient version
def partial_permutations_efficient(n, k, mod=1000000):
    """Efficient calculation of partial permutations."""
    if k > n:
        return 0
    
    # P(n,k) = n × (n-1) × ... × (n-k+1)
    # Calculate this directly without computing full factorials
    result = 1
    for i in range(n - k + 1, n + 1):
        result = (result * i) % mod
    
    return result

with open("Dataset.txt", "r") as file:
    n, k = map(int, file.read().strip().split())
print(f"Python Efficient: {partial_permutations_efficient(n, k)}")
