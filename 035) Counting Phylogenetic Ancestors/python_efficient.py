def count_internal_nodes_formula(n):
    """
    Efficient calculation using direct formula derivation.
    
    Derivation:
    In an unrooted binary tree:
    - Leaves: degree 1 (n leaves)
    - Internal nodes: degree 3 (i internal nodes)
    - Handshake lemma: sum(degrees) = 2 * edges
    - Tree property: edges = nodes - 1 = (n + i) - 1
    
    Equation: n*1 + i*3 = 2 * (n + i - 1)
    Solving: n + 3i = 2n + 2i - 2
           3i - 2i = 2n - n - 2
               i = n - 2
    """
    return n - 2

def main():
    # Read input
    with open("Dataset.txt", "r") as f:
        n = int(f.readline().strip())
    
    # Calculate and output
    print(count_internal_nodes_formula(n))

if __name__ == "__main__":
    main()