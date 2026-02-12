# Alternative using zip
def hamming_distance_zip(s, t):
    """Calculate Hamming distance using zip."""
    if len(s) != len(t):
        raise ValueError("Strings must be of equal length")
    
    return sum(1 for a, b in zip(s, t) if a != b)


path = "Dataset.txt"    # !!! Make sure the path is correct !!! 

with open(path, 'r') as file:
    data = file.read().strip()

s, t = data.split("\n")
print(f"Python Zip: {hamming_distance_zip(s, t)}")