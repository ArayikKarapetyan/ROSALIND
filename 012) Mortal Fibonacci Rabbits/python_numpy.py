import numpy as np

def mortal_rabbits_numpy(n, m):
    """Calculate mortal rabbits using NumPy."""
    if n <= 0:
        return 0
    if n == 1:
        return 1
    
    # Initialize age array
    ages = np.zeros(m, dtype=np.int64)
    ages[0] = 1
    
    for month in range(2, n + 1):
        # Newborns = sum of rabbits age 1 to m-1
        newborns = np.sum(ages[1:])
        
        # Shift ages: age 1 gets previous age 0, etc.
        # Age m-1 rabbits die
        ages[1:] = ages[:-1]
        ages[0] = newborns
    
    return np.sum(ages)

# Matrix-based approach (more mathematical)
def mortal_rabbits_matrix(n, m):
    """Calculate mortal rabbits using matrix exponentiation."""
    if n <= 0:
        return 0
    if n == 1:
        return 1
    
    # State vector: ages 0 to m-1
    state = np.zeros(m, dtype=np.int64)
    state[0] = 1
    
    # Transition matrix
    transition = np.zeros((m, m), dtype=np.int64)
    
    # Newborns column: all 1s except first (newborns don't reproduce)
    transition[0, 1:] = 1
    
    # Aging: shift down by one
    for i in range(1, m):
        transition[i, i-1] = 1
    
    # Apply transition n-1 times
    for _ in range(n-1):
        state = transition @ state
    
    return np.sum(state)


with open("Dataset.txt", "r") as file:
    n, m = map(int, file.read().strip().split())

print(f"Python NumPy: {mortal_rabbits_numpy(n, m)}")

