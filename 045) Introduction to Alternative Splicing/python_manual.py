def sum_combinations(n, m):
    """
    Calculate sum of C(n,k) for k from m to n, modulo 1,000,000.
    Uses Pascal's triangle with modular arithmetic.
    """
    MOD = 1000000
    
    # Initialize Pascal's triangle row
    row = [0] * (n + 1)
    row[0] = 1  # C(0,0) = 1
    
    # Build Pascal's triangle row by row
    for i in range(1, n + 1):
        # Compute new row from previous row
        new_row = [0] * (n + 1)
        new_row[0] = 1
        for j in range(1, i + 1):
            new_row[j] = (row[j-1] + row[j]) % MOD
        row = new_row
    
    # Sum values from m to n
    total = 0
    for k in range(m, n + 1):
        total = (total + row[k]) % MOD
    
    return total

def main():
    # Read input from Dataset.txt
    with open("Dataset.txt", "r") as f:
        n_str, m_str = f.readline().strip().split()
        n = int(n_str)
        m = int(m_str)
    
    # Calculate sum
    result = sum_combinations(n, m)
    print(result)

if __name__ == "__main__":
    main()