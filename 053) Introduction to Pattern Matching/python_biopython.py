def build_trie_biopython(strings):
    """
    BioPython-style trie construction.
    Could be useful for biological sequence patterns.
    """
    class TrieNode:
        __slots__ = ('id', 'children')
        
        def __init__(self, id):
            self.id = id
            self.children = {}
    
    root = TrieNode(1)
    next_id = 2
    edges = []
    
    for s in strings:
        current = root
        
        for char in s:
            if char not in current.children:
                # Create new node
                new_node = TrieNode(next_id)
                current.children[char] = new_node
                edges.append((current.id, new_node.id, char))
                next_id += 1
                current = new_node
            else:
                # Follow existing path
                current = current.children[char]
    
    return edges

def main():
    # Read input
    with open("Dataset.txt", "r") as f:
        strings = [line.strip() for line in f if line.strip()]
    
    # Build trie
    edges = build_trie_biopython(strings)
    
    # Output
    with open("output.txt", "w") as f:
        for edge in edges:
            f.write(f"{edge[0]} {edge[1]} {edge[2]}" + "\n")

if __name__ == "__main__":
    main()