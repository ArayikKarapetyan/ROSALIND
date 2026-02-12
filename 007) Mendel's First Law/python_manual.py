def dominant_probability_manual(k, m, n):
    """Calculate probability of dominant phenotype manually."""
    total = k + m + n
    
    # Calculate probabilities for all possible pairings
    # Probability both homozygous dominant
    p_kk = (k / total) * ((k - 1) / (total - 1))
    
    # Probability one homozygous dominant, one heterozygous
    p_km = (k / total) * (m / (total - 1))
    p_mk = (m / total) * (k / (total - 1))
    p_km_total = p_km + p_mk  # 2 * (k * m) / (total * (total - 1))
    
    # Probability one homozygous dominant, one homozygous recessive
    p_kn = (k / total) * (n / (total - 1))
    p_nk = (n / total) * (k / (total - 1))
    p_kn_total = p_kn + p_nk  # 2 * (k * n) / (total * (total - 1))
    
    # Probability both heterozygous
    p_mm = (m / total) * ((m - 1) / (total - 1))
    
    # Probability one heterozygous, one homozygous recessive
    p_mn = (m / total) * (n / (total - 1))
    p_nm = (n / total) * (m / (total - 1))
    p_mn_total = p_mn + p_nm  # 2 * (m * n) / (total * (total - 1))
    
    # Probability both homozygous recessive
    p_nn = (n / total) * ((n - 1) / (total - 1))
    
    # Calculate probability of dominant phenotype for each pairing
    # KK: always dominant = 1
    # KM or MK: always dominant = 1
    # KN or NK: always dominant = 1
    # MM: 3/4 probability of dominant
    # MN or NM: 1/2 probability of dominant
    # NN: 0 probability of dominant
    
    dominant_prob = (
        p_kk * 1 +                    # KK
        p_km_total * 1 +              # KM or MK
        p_kn_total * 1 +              # KN or NK
        p_mm * 0.75 +                 # MM
        p_mn_total * 0.5 +            # MN or NM
        p_nn * 0                      # NN
    )
    
    return dominant_prob

# More compact version
def dominant_probability_compact(k, m, n):
    """Compact calculation of dominant probability."""
    total = k + m + n
    
    # Calculate all probabilities in one go
    # Total pairs: total * (total - 1)
    total_pairs = total * (total - 1)
    
    # Calculate dominant offspring probability
    dominant = (
        k * (k - 1) * 1 +               # KK
        k * m * 1 * 2 +                 # KM or MK
        k * n * 1 * 2 +                 # KN or NK
        m * (m - 1) * 0.75 +            # MM
        m * n * 0.5 * 2 +               # MN or NM
        n * (n - 1) * 0                 # NN
    )
    
    return dominant / total_pairs

with open("./Dataset.txt", "r") as file:
    k, m, n = map(int, file.read().strip().split())

print(f"Python Manual: {dominant_probability_manual(k, m, n):.5f}")
print(f"Python Compact: {dominant_probability_compact(k, m, n):.5f}")
