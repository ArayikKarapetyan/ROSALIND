def sum_combinations_efficient(n, m):
    """
    Efficient calculation using Pascal's triangle with O(n^2) time and O(n) space.
    Also uses symmetry: C(n,k) = C(n,n-k)
    """
    MOD = 1000000
    
    # If m <= n/2, we can use the complementary range
    if m <= n // 2:
        # Compute sum from 0 to m-1, then subtract from 2^n
        # But wait, we need sum from m to n
        
        # Actually, let's compute directly using Pascal's triangle
        
        # Initialize current row of Pascal's triangle
        row = [1] + [0] * n
        
        # We'll accumulate sum as we build the triangle
        # But we need the final row for all k
        
        # Build Pascal's triangle up to row n
        for i in range(1, n + 1):
            # Update row in-place from right to left
            for j in range(i, 0, -1):
                row[j] = (row[j-1] + row[j]) % MOD
            row[0] = 1
        
        # Sum from m to n
        total = 0
        for k in range(m, n + 1):
            total = (total + row[k]) % MOD
        
        return total
    else:
        # If m > n/2, compute from m to n directly
        # This is smaller range
        row = [1] + [0] * n
        
        for i in range(1, n + 1):
            for j in range(i, 0, -1):
                row[j] = (row[j-1] + row[j]) % MOD
            row[0] = 1
        
        total = 0
        for k in range(m, n + 1):
            total = (total + row[k]) % MOD
        
        return total

def main():
    # Read input
    with open("Dataset.txt", "r") as f:
        n_str, m_str = f.readline().strip().split()
        n = int(n_str)
        m = int(m_str)
    
    # Calculate
    result = sum_combinations_efficient(n, m)
    print(result)

if __name__ == "__main__":
    main()