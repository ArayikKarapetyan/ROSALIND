from math import comb

def sum_combinations_biopython(n, m):
    """
    Uses Python's math.comb for combination calculation.
    Note: This may be inefficient for large n, but works for n<=2000.
    """
    MOD = 1000000
    total = 0
    
    for k in range(m, n + 1):
        total = (total + comb(n, k)) % MOD
    
    return total

def main():
    # Read input
    with open("Dataset.txt", "r") as f:
        n_str, m_str = f.readline().strip().split()
        n = int(n_str)
        m = int(m_str)
    
    result = sum_combinations_biopython(n, m)
    print(result)

if __name__ == "__main__":
    main()