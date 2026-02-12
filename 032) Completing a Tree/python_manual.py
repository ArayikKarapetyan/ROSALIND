def count_edges_to_complete_tree_manual(n, edges_list):
    """Calculate minimum edges to complete a tree manually."""
    # Build adjacency list
    adj_list = {i: [] for i in range(1, n + 1)}
    
    for u, v in edges_list:
        adj_list[u].append(v)
        adj_list[v].append(u)
    
    # Count connected components using DFS
    visited = set()
    components = 0
    
    for node in range(1, n + 1):
        if node not in visited:
            components += 1
            # DFS to mark all nodes in this component
            stack = [node]
            visited.add(node)
            
            while stack:
                current = stack.pop()
                for neighbor in adj_list[current]:
                    if neighbor not in visited:
                        visited.add(neighbor)
                        stack.append(neighbor)
    
    # Formula: For a tree with n nodes and k components,
    # need to add (k - 1) edges to connect all components
    return components - 1

def read_input(filename):
    with open(filename, 'r') as file:
        lines = file.read().strip().split('\n')
    
    n = int(lines[0].strip())
    edges_list = []
    
    for line in lines[1:]:
        if line.strip():
            u, v = map(int, line.strip().split())
            edges_list.append((u, v))
    
    return n, edges_list

n, edges_list = read_input("Dataset.txt")
print(f"Python Manual: {count_edges_to_complete_tree_manual(n, edges_list)}")
