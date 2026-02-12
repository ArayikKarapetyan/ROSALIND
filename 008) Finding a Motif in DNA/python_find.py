# Python has find() but it doesn't find all occurrences easily
def find_substring_positions_find(s, t):
    """Find positions using find() method."""
    positions = []
    start = 0
    
    while True:
        pos = s.find(t, start)
        if pos == -1:
            break
        positions.append(pos + 1)
        start = pos + 1  # Move past current match (allows overlapping)
    
    return positions

with open("Dataset.txt", "r") as file:
    lines = file.read().strip().split('\n')
    s = lines[0].strip()
    t = lines[1].strip()

print("Python Find:", *find_substring_positions_find(s, t))
