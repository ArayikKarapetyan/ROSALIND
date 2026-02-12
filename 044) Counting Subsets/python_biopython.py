def count_subsets_biopython_style(n):
    """
    BioPython-style implementation (though not really BioPython specific).
    Demonstrates modular arithmetic approach.
    """
    MOD = 1000000
    
    # Each element can be either in or out of subset: 2 choices per element
    # Total: 2^n mod MOD
    
    # Using iterative multiplication
    result = 1
    for _ in range(n):
        result = (result * 2) % MOD
    
    return result

def main():
    # Read input
    with open("Dataset.txt", "r") as f:
        n = int(f.readline().strip())
    
    # Count subsets
    result = count_subsets_biopython_style(n)
    print(result)

if __name__ == "__main__":
    main()