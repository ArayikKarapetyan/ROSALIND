# Using fractions for exact calculation
from fractions import Fraction

def expected_dominant_offspring_exact(counts):
    """Calculate exact expected value using fractions."""
    AA_AA, AA_Aa, AA_aa, Aa_Aa, Aa_aa, aa_aa = counts
    
    # Exact fractions for probabilities
    p_AA_AA = Fraction(1, 1)      # 1
    p_AA_Aa = Fraction(1, 1)      # 1
    p_AA_aa = Fraction(1, 1)      # 1
    p_Aa_Aa = Fraction(3, 4)      # 0.75
    p_Aa_aa = Fraction(1, 2)      # 0.5
    p_aa_aa = Fraction(0, 1)      # 0
    
    offspring_per_couple = 2
    
    # Calculate using fractions
    expected = Fraction(0, 1)
    expected += Fraction(AA_AA, 1) * p_AA_AA * offspring_per_couple
    expected += Fraction(AA_Aa, 1) * p_AA_Aa * offspring_per_couple
    expected += Fraction(AA_aa, 1) * p_AA_aa * offspring_per_couple
    expected += Fraction(Aa_Aa, 1) * p_Aa_Aa * offspring_per_couple
    expected += Fraction(Aa_aa, 1) * p_Aa_aa * offspring_per_couple
    expected += Fraction(aa_aa, 1) * p_aa_aa * offspring_per_couple
    
    return float(expected)

with open("Dataset.txt", "r") as file:
    counts = list(map(int, file.read().strip().split()))
    
print(f"Python Exact: {expected_dominant_offspring_exact(counts)}")
