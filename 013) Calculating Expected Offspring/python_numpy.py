# Using numpy for array operations
import numpy as np

def expected_dominant_offspring_numpy(counts):
    """Calculate using NumPy."""
    counts_array = np.array(counts, dtype=np.float64)
    probs = np.array([1.0, 1.0, 1.0, 0.75, 0.5, 0.0])
    
    expected = np.sum(counts_array * probs * 2)
    return expected


with open("Dataset.txt", "r") as file:
    counts = list(map(int, file.read().strip().split()))
    
print(f"Python NumPy: {expected_dominant_offspring_numpy(counts)}")