def rabbit_pairs_manual(n, k):
    """Calculate rabbit pairs using dynamic programming."""
    # Base cases
    if n == 1 or n == 2:
        return 1
    
    # DP approach: Fn = Fn-1 + k*Fn-2
    # F1 = 1, F2 = 1
    prev1 = 1  # Fn-1 (initially F2)
    prev2 = 1  # Fn-2 (initially F1)
    
    for i in range(3, n + 1):
        current = prev1 + k * prev2
        prev2, prev1 = prev1, current
    
    return prev1




with open("Dataset.txt", "r") as file:
    n, k = map(int, file.read().strip().split())
    
print(f"Python Manual: {rabbit_pairs_manual(n, k)}")
