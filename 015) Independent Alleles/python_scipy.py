from scipy.stats import binom

def probability_at_least_n_scipy(k, N):
    """Calculate using scipy's binomial distribution."""
    from scipy.stats import binom
    n = 2 ** k
    p = 0.25
    # 1 - CDF(N-1) = P(X >= N)
    probability = 1 - binom.cdf(N - 1, n, p)
    return probability

with open("Dataset.txt", "r") as file:
    k, N = map(int, file.read().strip().split())

print(f"Python NumPy: {probability_at_least_n_scipy(k, N):.6f}")