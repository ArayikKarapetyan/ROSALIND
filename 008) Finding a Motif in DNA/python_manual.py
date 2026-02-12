def find_substring_positions_manual(s, t):
    """Find all positions of substring t in s manually."""
    positions = []
    t_len = len(t)
    
    # Slide window through s
    for i in range(len(s) - t_len + 1):
        if s[i:i + t_len] == t:
            # Convert to 1-based indexing as per problem
            positions.append(i + 1)
    
    return positions

with open("Dataset.txt", "r") as file:
    lines = file.read().strip().split('\n')
    s = lines[0].strip()
    t = lines[1].strip()

print("Python Manual:", *find_substring_positions_manual(s, t))
