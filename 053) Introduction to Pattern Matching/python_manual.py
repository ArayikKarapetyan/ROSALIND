class TrieNode:
    def __init__(self, id):
        self.id = id
        self.children = {}  # char -> TrieNode
    
    def add_child(self, char, child_node):
        self.children[char] = child_node

def build_trie(strings):
    """Build a trie from a list of strings."""
    root = TrieNode(1)
    next_id = 2
    edges = []
    
    for s in strings:
        current = root
        for char in s:
            if char not in current.children:
                # Create new node
                new_node = TrieNode(next_id)
                current.add_child(char, new_node)
                edges.append((current.id, new_node.id, char))
                next_id += 1
                current = new_node
            else:
                # Follow existing edge
                current = current.children[char]
    
    return edges

def main():
    # Read input from Dataset.txt
    with open("Dataset.txt", "r") as f:
        strings = [line.strip() for line in f if line.strip()]
    
    # Build trie and get edges
    edges = build_trie(strings)
    
    # Output edges
    with open("output.txt", "w") as f:
        for edge in edges:
            f.write(f"{edge[0]} {edge[1]} {edge[2]}" + "\n")
            
if __name__ == "__main__":
    main()