import numpy as np

# Alternative using list for DP
def rabbit_pairs_dp(n, k):
    """Calculate rabbit pairs using DP list."""
    if n <= 2:
        return 1
    
    dp = [0] * (n + 1)
    dp[1] = 1
    dp[2] = 1
    
    for i in range(3, n + 1):
        dp[i] = dp[i-1] + k * dp[i-2]
    
    return dp[n]


def rabbit_pairs_numpy(n, k):
    """Calculate rabbit pairs using numpy for large n."""
    if n <= 2:
        return 1
    
    # Using matrix exponentiation for O(log n) solution
    # This is more efficient but not necessary for n ≤ 40
    # Simple DP is sufficient
    
    # DP with numpy array
    dp = np.zeros(n + 1, dtype=np.int64)
    dp[1] = 1
    dp[2] = 1
    
    for i in range(3, n + 1):
        dp[i] = dp[i-1] + k * dp[i-2]
    
    return int(dp[n])


with open("Dataset.txt", "r") as file:
    n, k = map(int, file.read().strip().split())


print(f"Python DP: {rabbit_pairs_dp(n, k)}")