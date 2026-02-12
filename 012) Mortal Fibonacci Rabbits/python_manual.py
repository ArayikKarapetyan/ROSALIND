def mortal_rabbits_manual(n, m):
    """Calculate mortal rabbits using dynamic programming."""
    if n <= 0:
        return 0
    if n == 1:
        return 1
    
    # DP array: rabbit pairs by age (0-indexed, age 0 = newborns)
    # ages range from 0 to m-1 (rabbits die after m months)
    ages = [0] * m
    ages[0] = 1  # Start with 1 newborn pair
    
    for month in range(2, n + 1):
        # Newborns = sum of reproducing rabbits (age 1 to m-1)
        newborns = sum(ages[1:])
        
        # Shift ages: rabbits age by 1 month
        # Rabbits at age m-1 die (don't carry over)
        for age in range(m-1, 0, -1):
            ages[age] = ages[age-1]
        
        # Set newborns
        ages[0] = newborns
    
    return sum(ages)

with open("Dataset.txt", "r") as file:
    n, m = map(int, file.read().strip().split())
    
print(f"Python Manual: {mortal_rabbits_manual(n, m)}")

