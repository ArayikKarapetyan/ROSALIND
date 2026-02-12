import re

def find_substring_positions_regex(s, t):
    """Find positions using regular expressions."""
    positions = []
    
    # Use lookahead to find overlapping matches
    pattern = f'(?={re.escape(t)})'
    matches = re.finditer(pattern, s)
    
    for match in matches:
        positions.append(match.start() + 1)
    
    return positions

with open("Dataset.txt", "r") as file:
    lines = file.read().strip().split('\n')
    s = lines[0].strip()
    t = lines[1].strip()

print("Python Regex:", *find_substring_positions_regex(s, t))