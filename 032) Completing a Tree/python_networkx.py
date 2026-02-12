import networkx as nx

def count_edges_to_complete_tree_networkx(n, edges_list):
    """Calculate minimum edges to complete a tree using NetworkX."""
    # Create graph
    G = nx.Graph()
    G.add_nodes_from(range(1, n + 1))
    G.add_edges_from(edges_list)
    
    # Count connected components
    components = nx.number_connected_components(G)
    
    # Minimum edges to add = components - 1
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

print(f"Python NetworkX: {count_edges_to_complete_tree_networkx(n, edges_list)}")