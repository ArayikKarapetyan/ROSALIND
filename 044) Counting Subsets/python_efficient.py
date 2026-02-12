def count_subsets_efficient(n):
    """
    Efficient calculation of 2^n mod 1,000,000.
    Uses Python's pow with modular exponentiation.
    """
    MOD = 1000000
    return pow(2, n, MOD)

def main():
    # Read input
    with open("Dataset.txt", "r") as f:
        n = int(f.readline().strip())
    
    # Calculate result
    result = count_subsets_efficient(n)
    print(result)

if __name__ == "__main__":
    main()