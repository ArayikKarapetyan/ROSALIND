# For this problem, we can use itertools for combinations
import itertools

def dominant_probability_itertools(k, m, n):
    """Calculate probability using itertools."""
    # Create population list
    population = ['AA'] * k + ['Aa'] * m + ['aa'] * n
    total = len(population)
    
    # Generate all possible pairs
    pairs = list(itertools.combinations(population, 2))
    
    dominant_count = 0
    total_pairs = len(pairs)
    
    # Punnett square probabilities
    # AA x AA: 100% dominant
    # AA x Aa: 100% dominant
    # AA x aa: 100% dominant
    # Aa x Aa: 75% dominant
    # Aa x aa: 50% dominant
    # aa x aa: 0% dominant
    
    for parent1, parent2 in pairs:
        if parent1 == 'AA' or parent2 == 'AA':
            # At least one homozygous dominant
            dominant_count += 1
        elif parent1 == 'Aa' and parent2 == 'Aa':
            # Both heterozygous
            dominant_count += 0.75
        elif (parent1 == 'Aa' and parent2 == 'aa') or (parent1 == 'aa' and parent2 == 'Aa'):
            # One heterozygous, one homozygous recessive
            dominant_count += 0.5
        else:
            # Both homozygous recessive
            dominant_count += 0
    
    return dominant_count / total_pairs


with open("./Dataset.txt", "r") as file:
    k, m, n = map(int, file.read().strip().split())


print(f"Python Exact: {dominant_probability_itertools(k, m, n):.5f}")