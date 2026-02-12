MASS_TABLE = {
    'A': 71.03711, 'C': 103.00919, 'D': 115.02694, 'E': 129.04259,
    'F': 147.06841, 'G': 57.02146, 'H': 137.05891, 'I': 113.08406,
    'K': 128.09496, 'L': 113.08406, 'M': 131.04049, 'N': 114.04293,
    'P': 97.05276, 'Q': 128.05858, 'R': 156.10111, 'S': 87.03203,
    'T': 101.04768, 'V': 99.06841, 'W': 186.07931, 'Y': 163.06333
}

EPS = 0.01

def find_protein_from_prefix_spectrum(L):
    protein = []

    for i in range(1, len(L)):
        diff = L[i] - L[i-1]

        for aa, mass in MASS_TABLE.items():
            if abs(diff - mass) < EPS:
                protein.append(aa)
                break

    return "".join(protein)

def main():
    with open("Dataset.txt") as f:
        L = [float(line.strip()) for line in f if line.strip()]

    protein = find_protein_from_prefix_spectrum(L)
    print(protein)

if __name__ == "__main__":
    main()
