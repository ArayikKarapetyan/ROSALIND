def expected_dominant_offspring_manual(counts):
    """Calculate expected dominant offspring manually."""
    # Unpack counts for each mating pair type
    AA_AA, AA_Aa, AA_aa, Aa_Aa, Aa_aa, aa_aa = counts
    
    # Probabilities of dominant offspring for each mating pair type
    # Based on Mendelian inheritance
    p_AA_AA = 1.0      # AA x AA: always dominant
    p_AA_Aa = 1.0      # AA x Aa: always dominant
    p_AA_aa = 1.0      # AA x aa: always dominant
    p_Aa_Aa = 0.75     # Aa x Aa: 75% dominant
    p_Aa_aa = 0.5      # Aa x aa: 50% dominant
    p_aa_aa = 0.0      # aa x aa: never dominant
    
    # Each couple has exactly 2 offspring
    offspring_per_couple = 2
    
    # Calculate expected dominant offspring for each pair type
    expected = (
        AA_AA * p_AA_AA * offspring_per_couple +
        AA_Aa * p_AA_Aa * offspring_per_couple +
        AA_aa * p_AA_aa * offspring_per_couple +
        Aa_Aa * p_Aa_Aa * offspring_per_couple +
        Aa_aa * p_Aa_aa * offspring_per_couple +
        aa_aa * p_aa_aa * offspring_per_couple
    )
    
    return expected

# More compact version
def expected_dominant_offspring_compact(counts):
    """Compact calculation."""
    AA_AA, AA_Aa, AA_aa, Aa_Aa, Aa_aa, aa_aa = counts
    
    # Dominant probabilities for each pair type
    probs = [1.0, 1.0, 1.0, 0.75, 0.5, 0.0]
    
    # Each couple has 2 offspring
    expected = sum(count * prob * 2 for count, prob in zip(counts, probs))
    
    return expected

with open("Dataset.txt", "r") as file:
    counts = list(map(int, file.read().strip().split()))
    
print(f"Python Manual: {expected_dominant_offspring_manual(counts)}")
print(f"Python Compact: {expected_dominant_offspring_compact(counts)}")
