def parse_fasta_manual(fasta_string):
    sequences = {}
    current_id = ""
    for line in fasta_string.strip().split('\n'):
        if line.startswith('>'):
            current_id = line[1:].strip()
            sequences[current_id] = ""
        else:
            sequences[current_id] += line.strip()
    return sequences


def find_subsequence_indices_manual(s, t):
    """Find indices where t appears as subsequence in s."""
    indices = []
    i = 0  # pointer in s
    j = 0  # pointer in t
    
    while i < len(s) and j < len(t):
        if s[i] == t[j]:
            indices.append(i + 1)  # 1-based indexing
            j += 1
        i += 1
    
    # If we found all characters in t
    if j == len(t):
        return indices
    else:
        return None

# Alternative greedy approach
def find_subsequence_greedy(s, t):
    """Find subsequence indices using greedy algorithm."""
    indices = []
    s_pos = 0
    t_pos = 0
    
    while t_pos < len(t) and s_pos < len(s):
        if s[s_pos] == t[t_pos]:
            indices.append(s_pos + 1)
            t_pos += 1
        s_pos += 1
    
    return indices if t_pos == len(t) else None


with open("Dataset.txt", "r") as file:
    fasta_data = file.read()
sequences = parse_fasta_manual(fasta_data)
s = list(sequences.values())[0]
t = list(sequences.values())[1]
indices_manual = find_subsequence_indices_manual(s, t)
print("Manual:", " ".join(map(str, indices_manual)))