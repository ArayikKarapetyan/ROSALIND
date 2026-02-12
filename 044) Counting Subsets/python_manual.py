def count_subsets(n):
    """
    Count all subsets of a set with n elements.
    Total subsets = 2^n.
    """
    MOD = 1000000
    
    # 2^n modulo MOD using fast modular exponentiation
    result = 1
    base = 2
    exp = n
    
    while exp > 0:
        if exp % 2 == 1:
            result = (result * base) % MOD
        base = (base * base) % MOD
        exp //= 2
    
    return result

def main():
    # Read input from Dataset.txt
    with open("Dataset.txt", "r") as f:
        n = int(f.readline().strip())
    
    # Count subsets
    result = count_subsets(n)
    print(result)

if __name__ == "__main__":
    main()