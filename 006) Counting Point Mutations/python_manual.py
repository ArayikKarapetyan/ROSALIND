def hamming_distance_manual(s, t):
    """Calculate Hamming distance manually."""
    if len(s) != len(t):
        raise ValueError("Strings must be of equal length")
    
    distance = 0
    for i in range(len(s)):
        if s[i] != t[i]:
            distance += 1
    
    return distance


path = "Dataset.txt"    # !!! Make sure the path is correct !!! 

with open(path, 'r') as file:
    data = file.read().strip()

s, t = data.split("\n")
print(f"Python Manual: {hamming_distance_manual(s, t)}")
