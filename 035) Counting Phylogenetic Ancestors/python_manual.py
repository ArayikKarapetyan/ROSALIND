def count_internal_nodes(n_leaves):
    """
    Calculate number of internal nodes in an unrooted binary tree.
    
    For an unrooted binary tree:
    - All leaves have degree 1
    - All internal nodes have degree 3
    - Using handshake lemma: sum(degrees) = 2 * edges
    - Number of edges = (total_degree_sum) / 2
    """
    # For an unrooted binary tree:
    # Let L = number of leaves (degree 1)
    # Let I = number of internal nodes (degree 3)
    # Total nodes N = L + I
    # Total edges E = (L*1 + I*3) / 2
    # Also, in any tree: E = N - 1
    # So: (L + 3I)/2 = (L + I) - 1
    # Solving: L + 3I = 2L + 2I - 2
    # => 3I = L + 2I - 2
    # => I = L - 2
    
    return n_leaves - 2

def main():
    # Read input
    with open("Dataset.txt", "r") as f:
        n = int(f.readline().strip())
    
    # Calculate internal nodes
    result = count_internal_nodes(n)
    print(result)

if __name__ == "__main__":
    main()