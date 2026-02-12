class TreeNode:
    __slots__ = ('name', 'parent', 'children', 'depth')
    
    def __init__(self, name):
        self.name = name
        self.parent = None
        self.children = []
        self.depth = 0

def parse_newick_efficient(newick):
    """Parse Newick format efficiently."""
    newick = newick.strip().rstrip(';')
    
    root = TreeNode("")
    current = root
    stack = []
    name_buffer = []
    node_map = {}
    
    i = 0
    while i < len(newick):
        c = newick[i]
        
        if c == '(':
            # New internal node
            new_node = TreeNode("")
            current.children.append(new_node)
            new_node.parent = current
            stack.append(current)
            current = new_node
            i += 1
            
        elif c == ')':
            # End internal node
            if name_buffer:
                current.name = ''.join(name_buffer)
                if current.name:
                    node_map[current.name] = current
                name_buffer = []
            
            current = stack.pop()
            i += 1
            
            # Read node name after ')'
            while i < len(newick) and newick[i] not in ',()':
                name_buffer.append(newick[i])
                i += 1
            if name_buffer:
                current.name = ''.join(name_buffer)
                if current.name:
                    node_map[current.name] = current
                name_buffer = []
                
        elif c == ',':
            # New sibling
            if name_buffer:
                current.name = ''.join(name_buffer)
                if current.name:
                    node_map[current.name] = current
                name_buffer = []
            
            current = stack[-1]
            new_node = TreeNode("")
            current.children.append(new_node)
            new_node.parent = current
            current = new_node
            i += 1
            
        else:
            name_buffer.append(c)
            i += 1
    
    # Handle any remaining name
    if name_buffer:
        current.name = ''.join(name_buffer)
        if current.name:
            node_map[current.name] = current
    
    # Calculate depths
    def set_depths(node, depth):
        node.depth = depth
        for child in node.children:
            set_depths(child, depth + 1)
    
    set_depths(root, 0)
    
    return root, node_map

def distance_between_nodes_efficient(node1, node2):
    """Efficient distance calculation using depths and LCA."""
    # Move up to same depth
    a, b = node1, node2
    while a.depth > b.depth:
        a = a.parent
    while b.depth > a.depth:
        b = b.parent
    
    # Find LCA
    while a != b:
        a = a.parent
        b = b.parent
    
    lca = a
    
    # Distance = d1 + d2 - 2*dlca
    return node1.depth + node2.depth - 2 * lca.depth

def main():
    # Read input
    with open("Dataset.txt", "r") as f:
        lines = [line.strip() for line in f if line.strip()]
    
    results = []
    i = 0
    
    while i < len(lines):
        tree_str = lines[i]
        i += 1
        
        if i >= len(lines):
            break
            
        node1_name, node2_name = lines[i].split()
        i += 1
        
        # Parse tree with node mapping
        root, node_map = parse_newick_efficient(tree_str)
        
        node1 = node_map.get(node1_name)
        node2 = node_map.get(node2_name)
        
        if node1 and node2:
            dist = distance_between_nodes_efficient(node1, node2)
            results.append(str(dist))
        else:
            results.append("0")
    
    print(" ".join(results))

if __name__ == "__main__":
    main()