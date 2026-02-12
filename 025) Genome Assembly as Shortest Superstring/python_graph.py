def parse_fasta(fasta_string):
    """Parse FASTA format manually."""
    sequences = {}
    current_id = ""
    
    for line in fasta_string.strip().split('\n'):
        if line.startswith('>'):
            current_id = line[1:].strip()
            sequences[current_id] = ""
        else:
            sequences[current_id] += line.strip()
    
    return list(sequences.values())

def overlap_length(a, b, min_overlap=1):
    """Calculate maximum overlap where end of a matches start of b."""
    max_len = min(len(a), len(b))
    
    # Try larger overlaps first (problem says > half length)
    for overlap in range(max_len, min_overlap - 1, -1):
        if a[-overlap:] == b[:overlap]:
            return overlap
    return 0

def build_overlap_graph(strings, min_overlap=1):
    """Build directed graph where edge weight = overlap length."""
    n = len(strings)
    graph = [[0] * n for _ in range(n)]
    
    for i in range(n):
        for j in range(n):
            if i != j:
                overlap = overlap_length(strings[i], strings[j], min_overlap)
                graph[i][j] = overlap
    
    return graph

def greedy_superstring_from_graph(strings, graph):
    """Greedy assembly from overlap graph."""
    n = len(strings)
    used = [False] * n
    result = ""
    
    # Start with arbitrary string
    current = 0
    used[current] = True
    result = strings[current]
    
    while True:
        # Find best unused string that overlaps with current
        best_next = -1
        best_overlap = 0
        
        for i in range(n):
            if not used[i] and graph[current][i] > best_overlap:
                best_overlap = graph[current][i]
                best_next = i
        
        if best_next == -1:
            # No more overlaps, add remaining strings
            for i in range(n):
                if not used[i]:
                    result += strings[i]
                    used[i] = True
            break
        
        # Merge with best overlap
        result += strings[best_next][best_overlap:]
        used[best_next] = True
        current = best_next
    
    return result

def shortest_superstring_graph(strings):
    """Shortest superstring using graph approach."""
    # Build overlap graph with minimum overlap > half length
    # As per problem statement: overlap > half length
    min_len = min(len(s) for s in strings)
    min_overlap = min_len // 2 + 1
    
    graph = build_overlap_graph(strings, min_overlap)
    
    # Try starting from each string
    best_result = None
    best_length = float('inf')
    
    for start in range(len(strings)):
        result = greedy_superstring_from_graph(strings, graph)
        if len(result) < best_length:
            best_length = len(result)
            best_result = result
    
    return best_result


with open("Dataset.txt", "r") as file:
    strings = parse_fasta(file.read())

result3 = shortest_superstring_graph(strings)
print(result3)