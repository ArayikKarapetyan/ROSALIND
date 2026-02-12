class TreeNode:
    def __init__(self, name):
        self.name = name
        self.children = []
        self.parent = None
    
    def add_child(self, child):
        child.parent = self
        self.children.append(child)

def parse_newick(newick):
    """Parse Newick format string into a tree."""
    newick = newick.strip()
    if newick.endswith(';'):
        newick = newick[:-1]
    
    stack = []
    current = TreeNode("")  # root
    i = 0
    name = ""
    
    while i < len(newick):
        c = newick[i]
        
        if c == '(':
            # Start new internal node
            new_node = TreeNode("")
            current.add_child(new_node)
            stack.append(current)
            current = new_node
            i += 1
            
        elif c == ')':
            # End of internal node
            if name:
                current.name = name
                name = ""
            current = stack.pop()
            i += 1
            
            # Read node name if present
            while i < len(newick) and newick[i] not in ',()':
                name += newick[i]
                i += 1
            if name:
                current.name = name
                name = ""
                
        elif c == ',':
            # End of current child
            if name:
                current.name = name
                name = ""
            current = stack[-1]  # Go back to parent
            # Start new sibling
            new_node = TreeNode("")
            current.add_child(new_node)
            current = new_node
            i += 1
            
        else:
            # Accumulate name
            name += c
            i += 1
    
    # Handle root name
    if name:
        current.name = name
    
    return current

def find_node(node, target_name):
    """Find node with given name in tree."""
    if node.name == target_name:
        return node
    
    for child in node.children:
        result = find_node(child, target_name)
        if result:
            return result
    
    return None

def find_lca(node1, node2):
    """Find lowest common ancestor of two nodes."""
    # Get paths to root
    path1 = set()
    current = node1
    while current:
        path1.add(current)
        current = current.parent
    
    # Find first common node
    current = node2
    while current:
        if current in path1:
            return current
        current = current.parent
    
    return None

def distance_between_nodes(node1, node2):
    """Calculate distance between two nodes in tree."""
    lca = find_lca(node1, node2)
    
    # Distance = depth(node1) + depth(node2) - 2*depth(lca)
    def depth(node):
        d = 0
        while node and node.parent:
            d += 1
            node = node.parent
        return d
    
    d1 = depth(node1)
    d2 = depth(node2)
    dlca = depth(lca) if lca else 0
    
    return d1 + d2 - 2 * dlca

def main():
    # Read input from Dataset.txt
    with open("Dataset.txt", "r") as f:
        lines = [line.strip() for line in f if line.strip()]
    
    results = []
    i = 0
    
    while i < len(lines):
        # Parse tree
        tree_str = lines[i]
        i += 1
        
        # Parse node pair
        if i >= len(lines):
            break
        node1_name, node2_name = lines[i].split()
        i += 1
        
        # Build tree
        root = parse_newick(tree_str)
        
        # Find nodes
        node1 = find_node(root, node1_name)
        node2 = find_node(root, node2_name)
        
        if node1 and node2:
            dist = distance_between_nodes(node1, node2)
            results.append(str(dist))
        else:
            results.append("0")  # One or both nodes not found
    
    # Output results
    print(" ".join(results))

if __name__ == "__main__":
    main()