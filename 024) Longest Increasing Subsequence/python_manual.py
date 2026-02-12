def parse_input_manual(input_string):
    """Parse input: first line n, second line permutation."""
    lines = input_string.strip().split('\n')
    n = int(lines[0].strip())
    permutation = list(map(int, lines[1].strip().split()))
    return n, permutation

def longest_increasing_subsequence_manual(seq):
    """Find longest increasing subsequence using DP."""
    n = len(seq)
    if n == 0:
        return []
    
    # DP arrays
    dp = [1] * n  # Length of LIS ending at i
    prev = [-1] * n  # Previous index in LIS
    
    # Compute DP
    for i in range(n):
        for j in range(i):
            if seq[j] < seq[i] and dp[j] + 1 > dp[i]:
                dp[i] = dp[j] + 1
                prev[i] = j
    
    # Find max length and its index
    max_len = max(dp)
    max_idx = dp.index(max_len)
    
    # Reconstruct LIS
    lis = []
    curr = max_idx
    while curr != -1:
        lis.append(seq[curr])
        curr = prev[curr]
    
    return lis[::-1]  # Reverse to get increasing order

def longest_decreasing_subsequence_manual(seq):
    """Find longest decreasing subsequence using DP."""
    n = len(seq)
    if n == 0:
        return []
    
    # DP arrays
    dp = [1] * n  # Length of LDS ending at i
    prev = [-1] * n  # Previous index in LDS
    
    # Compute DP
    for i in range(n):
        for j in range(i):
            if seq[j] > seq[i] and dp[j] + 1 > dp[i]:
                dp[i] = dp[j] + 1
                prev[i] = j
    
    # Find max length and its index
    max_len = max(dp)
    max_idx = dp.index(max_len)
    
    # Reconstruct LDS
    lds = []
    curr = max_idx
    while curr != -1:
        lds.append(seq[curr])
        curr = prev[curr]
    
    return lds[::-1]  # Reverse to get decreasing order

def process_permutation_manual(input_string):
    """Main processing function."""
    n, permutation = parse_input_manual(input_string)
    
    lis = longest_increasing_subsequence_manual(permutation)
    lds = longest_decreasing_subsequence_manual(permutation)
    
    return lis, lds

with open("Dataset.txt", "r") as file:
    input_data = file.read()

lis, lds = process_permutation_manual(input_data)

print(' '.join(map(str, lis)))
print(' '.join(map(str, lds)))

