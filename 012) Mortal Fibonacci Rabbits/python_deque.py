def mortal_rabbits_queue(n, m):
    """Calculate mortal rabbits using deque."""
    from collections import deque
    
    if n <= 0:
        return 0
    if n == 1:
        return 1
    
    # Use deque for easy rotation
    # Each position represents rabbits of that age
    ages = deque([0] * m)
    ages[0] = 1
    
    for month in range(2, n + 1):
        # Newborns = sum of reproducing ages (1 to m-1)
        newborns = sum(list(ages)[1:])
        
        # Age rabbits: shift right
        ages.pop()  # Oldest die
        ages.appendleft(newborns)  # Newborns become age 0
    
    return sum(ages)

with open("Dataset.txt", "r") as file:
    n, m = map(int, file.read().strip().split())
    
print(f"Python Manual: {mortal_rabbits_queue(n, m)}")

