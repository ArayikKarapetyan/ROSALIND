# Log-space calculation to handle very large numbers
import math
from math import comb

def probability_at_least_n_logspace(k, N):
    """Calculate in log space to avoid numerical issues."""
    n = 2 ** k
    p = 0.25
    q = 1 - p
    
    # Calculate log probabilities
    log_q_n = n * math.log(q)
    
    # Initialize log probability of 0 successes
    log_prob_sum = float('-inf')  # log(0)
    
    # Calculate sum of probabilities for i < N
    for i in range(N):
        if i == 0:
            log_prob_i = log_q_n
        else:
            # log(C(n,i)) + i*log(p) + (n-i)*log(q)
            # Use log of binomial coefficient: log(C(n,i)) = log(n!) - log(i!) - log((n-i)!)
            # Or use recurrence: log(C(n,i+1)) = log(C(n,i)) + log(n-i) - log(i+1)
            if i == 1:
                log_comb = math.log(n)
            else:
                # This is simplified; in practice would need proper log comb calculation
                log_comb = math.log(comb(n, i))
            
            log_prob_i = log_comb + i * math.log(p) + (n - i) * math.log(q)
        
        # Sum probabilities in log space: log(a+b) = log(a) + log(1 + exp(log(b)-log(a)))
        if log_prob_sum == float('-inf'):
            log_prob_sum = log_prob_i
        else:
            # log_sum = log(a+b) = log(a) + log(1 + exp(log(b)-log(a)))
            diff = log_prob_i - log_prob_sum
            if diff > 30:  # If b >> a
                log_prob_sum = log_prob_i
            elif diff < -30:  # If a >> b
                pass  # Keep log_prob_sum
            else:
                log_prob_sum = log_prob_sum + math.log(1 + math.exp(diff))
    
    # Probability at least N = 1 - sum(P(i) for i < N)
    if log_prob_sum > 0:  # Shouldn't happen with probabilities
        prob_at_least = 0.0
    else:
        prob_less_than_n = math.exp(log_prob_sum)
        prob_at_least = 1 - prob_less_than_n
    
    return max(0.0, min(1.0, prob_at_least))

with open("Dataset.txt", "r") as file:
    k, N = map(int, file.read().strip().split())

print(f"Python Logspace: {probability_at_least_n_logspace(k, N):.6f}")
