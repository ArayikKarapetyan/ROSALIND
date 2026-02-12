from math import comb

def probability_at_least_n_manual(k, N):
    """Calculate probability of at least N AaBb in generation k."""
    # Total organisms in generation k
    total_organisms = 2 ** k
    
    # Probability that a single offspring is AaBb
    # When AaBb mates with AaBb for two independent traits:
    # For each trait: Aa x Aa → 0.5 Aa, aa x bb → 0.5 Bb
    # Since independent: 0.5 * 0.5 = 0.25 probability of AaBb
    p_single = 0.25
    
    # Use binomial distribution: X ~ Binomial(n=total_organisms, p=0.25)
    # P(X >= N) = 1 - P(X < N) = 1 - Σ_{i=0}^{N-1} C(n,i) * p^i * (1-p)^{n-i}
    
    probability_less_than_n = 0.0
    n = total_organisms
    
    for i in range(N):
        # Binomial probability for exactly i successes
        prob_i = comb(n, i) * (p_single ** i) * ((1 - p_single) ** (n - i))
        probability_less_than_n += prob_i
    
    probability_at_least_n = 1 - probability_less_than_n
    return probability_at_least_n


with open("Dataset.txt", "r") as file:
    k, N = map(int, file.read().strip().split())
    
print(f"Python Manual: {probability_at_least_n_manual(k, N):.6f}")
