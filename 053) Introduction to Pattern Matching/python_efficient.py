def build_trie_efficient(strings):
    """
    Build trie efficiently without explicit node objects.
    Returns list of edges as (parent_id, child_id, char).
    """
    # Node storage: id -> dict of children {char: child_id}
    nodes = {1: {}}  # root has id 1
    edges = []
    next_id = 2
    
    for s in strings:
        current_id = 1  # start at root
        
        for char in s:
            # Check if current node has child with this char
            if char in nodes[current_id]:
                # Follow existing edge
                current_id = nodes[current_id][char]
            else:
                # Create new node and edge
                nodes[current_id][char] = next_id
                nodes[next_id] = {}  # new empty node
                edges.append((current_id, next_id, char))
                current_id = next_id
                next_id += 1
    
    return edges

def main():
    # Read input
    with open("Dataset.txt", "r") as f:
        strings = [line.strip() for line in f if line.strip()]
    
    # Build trie
    edges = build_trie_efficient(strings)
    
    # Output
    with open("output.txt", "w") as f:
        for parent, child, char in edges:
            f.write(f"{parent} {child} {char}" + "\n")
        


if __name__ == "__main__":
    main()