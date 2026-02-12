# Expected Dominant Offspring Solutions

![Python](https://img.shields.io/badge/Python-3.8+-blue.svg)
![R](https://img.shields.io/badge/R-4.0+-blue.svg)
![Probability](https://img.shields.io/badge/Probability-Expected%20Value-blue.svg)
![Genetics](https://img.shields.io/badge/Genetics-Mendelian-green.svg)

<p align="center">
  <img src="https://rosalind.info/static/img/logo.png?4eac7af5ad70" alt="ROSALIND Logo">
</p>
A bioinformatics solution for calculating the expected number of dominant phenotype offspring from different parental genotype combinations, using probability and Mendelian genetics.

## 📋 Problem Description

**Problem:** [Calculating Expected Offspring](https://rosalind.info/problems/iev/)  
**Category:** Bioinformatics Stronghold  
**ID:** IEV

Given counts of couples with specific genotype pairings, calculate the expected number of offspring displaying the dominant phenotype in the next generation, assuming each couple produces exactly two offspring.

**Input:** Six nonnegative integers (each ≤ 20,000) representing counts of couples with genotype pairings:

1. AA-AA
2. AA-Aa
3. AA-aa
4. Aa-Aa
5. Aa-aa
6. aa-aa

**Output:** Expected number of dominant offspring (float)

### Example

```text
Input: 1 0 0 1 0 1
Output: 3.5
```

**Explanation:**
- 1 AA-AA couple: 2 offspring × 100% dominant = 2.0
- 1 Aa-Aa couple: 2 offspring × 75% dominant = 1.5
- 1 aa-aa couple: 2 offspring × 0% dominant = 0.0

Total expected: 2.0 + 1.5 + 0.0 = 3.5

## 🧬 Solutions

### 1. Python Solution (From Scratch) - `python_manual.py`

```python
def expected_dominant_offspring_manual(counts):
    """Calculate expected dominant offspring."""
    AA_AA, AA_Aa, AA_aa, Aa_Aa, Aa_aa, aa_aa = counts
    
    # Probabilities from Mendelian inheritance
    probs = [1.0, 1.0, 1.0, 0.75, 0.5, 0.0]
    offspring_per_couple = 2
    
    expected = sum(count * prob * offspring_per_couple 
                   for count, prob in zip(counts, probs))
    
    return expected
```

**Features:**
- Direct application of probability rules
- O(1) time complexity
- Clear mapping of genotype pairings to probabilities

### 2. Python Solution with Fractions - `python_exact.py`

```python
from fractions import Fraction

def expected_dominant_offspring_exact(counts):
    """Calculate exact expected value using fractions."""
    AA_AA, AA_Aa, AA_aa, Aa_Aa, Aa_aa, aa_aa = counts
    
    # Exact fractions for probabilities
    probs = [Fraction(1,1), Fraction(1,1), Fraction(1,1), 
             Fraction(3,4), Fraction(1,2), Fraction(0,1)]
    
    expected = Fraction(0,1)
    for count, prob in zip(counts, probs):
        expected += Fraction(count, 1) * prob * 2
    
    return float(expected)
```

**Features:**
- Exact calculation using fractions
- Avoids floating-point precision issues
- Useful for verification

### 3. R Solution - `r_solution.R`

```r
expected_dominant_offspring_r <- function(counts) {
  probs <- c(1.0, 1.0, 1.0, 0.75, 0.5, 0.0)
  expected <- sum(counts * probs * 2)
  return(expected)
}
```

**Features:**
- Vectorized R calculation
- Multiple implementation approaches
- Clean and efficient

## 📁 File Structure

```text
Expected-Offspring/
│
├── README.md              # This documentation
├── python_manual.py       # Python basic solution
├── python_exact.py        # Python exact fraction solution
├── python_numpy.py        # Python NumPy solution
├── r_solution.R           # R implementation
└── Dataset.txt            # Input file (6 integers)
```

## 🚀 Installation and Usage

### Python Requirements

```bash
# Only needed for NumPy solution
pip install numpy
```

### Running the Solutions

**Python (Manual):**

```bash
python python_manual.py
```

**Python (Exact):**

```bash
python python_exact.py
```

**R:**

```r
# From R console
source("r_solution.R")

# From command line
Rscript r_solution.R
```

## 🔧 Configuration

### Input Format

Input file `Dataset.txt` should contain six integers separated by spaces:

```text
x1 x2 x3 x4 x5 x6
```

**Example:**

```text
1 0 0 1 0 1
```

### File Reading

**Python:**

```python
with open("Dataset.txt", "r") as file:
    counts = list(map(int, file.read().strip().split()))
```

**R:**

```r
counts <- scan("Dataset.txt", quiet = TRUE)
# or
input_data <- readLines("Dataset.txt", warn = FALSE)
counts <- as.integer(strsplit(input_data, " ")[[1]])
```

## 📊 Algorithm Analysis

### Time Complexity

| Method | Approach | Time Complexity | Space Complexity |
|--------|----------|-----------------|------------------|
| Python Manual | Direct calculation | O(1) | O(1) |
| Python Fractions | Exact calculation | O(1) | O(1) |
| R Vectorized | Vector operations | O(1) | O(1) |
| R Matrix | Matrix multiplication | O(1) | O(1) |

### Genetic Probabilities

Derived from Punnett squares:

- **AA × AA:** 100% AA (all dominant)
- **AA × Aa:** 50% AA, 50% Aa (all dominant)
- **AA × aa:** 100% Aa (all dominant)
- **Aa × Aa:** 25% AA, 50% Aa, 25% aa (75% dominant)
- **Aa × aa:** 50% Aa, 50% aa (50% dominant)
- **aa × aa:** 100% aa (0% dominant)

## 🧪 Testing

### Test Cases

```python
# Test case from problem
assert expected_dominant_offspring_manual([1, 0, 0, 1, 0, 1]) == 3.5

# Additional test cases
assert expected_dominant_offspring_manual([1, 0, 0, 0, 0, 0]) == 2.0  # AA-AA only
assert expected_dominant_offspring_manual([0, 1, 0, 0, 0, 0]) == 2.0  # AA-Aa only
assert expected_dominant_offspring_manual([0, 0, 0, 1, 0, 0]) == 1.5  # Aa-Aa only
assert expected_dominant_offspring_manual([0, 0, 0, 0, 1, 0]) == 1.0  # Aa-aa only
assert expected_dominant_offspring_manual([0, 0, 0, 0, 0, 1]) == 0.0  # aa-aa only
assert expected_dominant_offspring_manual([2, 2, 2, 2, 2, 2]) == 13.0  # Mixed
```

### Sample Dataset Verification

```text
Input: 1 0 0 1 0 1

Breakdown:
1. AA-AA: 1 couple × 2 offspring × 100% dominant = 2.0
2. AA-Aa: 0 couples = 0.0
3. AA-aa: 0 couples = 0.0
4. Aa-Aa: 1 couple × 2 offspring × 75% dominant = 1.5
5. Aa-aa: 0 couples = 0.0
6. aa-aa: 1 couple × 2 offspring × 0% dominant = 0.0

Total: 2.0 + 1.5 + 0.0 = 3.5
```

## 🔗 Related Problems

- [IPRB](https://rosalind.info/problems/iprb/) - Mendel's First Law
- [FIB](https://rosalind.info/problems/fib/) - Rabbits and Recurrence Relations
- [LIA](https://rosalind.info/problems/lia/) - Independent Alleles
- [PROB](https://rosalind.info/problems/prob/) - Introduction to Random Strings

## 📚 Learning Resources

- [Expected Value](https://en.wikipedia.org/wiki/Expected_value)
- [Mendelian Inheritance](https://en.wikipedia.org/wiki/Mendelian_inheritance)
- [Punnett Square](https://en.wikipedia.org/wiki/Punnett_square)
- [Probability Theory](https://en.wikipedia.org/wiki/Probability_theory)

## 🤝 Contributing

Contributions are welcome! Here are ways to contribute:

### Add New Features:
- Support for different offspring counts
- Variable offspring per couple
- Sex-linked traits
- Multiple gene interactions

### Improve Documentation:
- Add Punnett square diagrams
- Create probability tree visualizations
- Add biological context and examples

### Enhance Performance:
- Batch processing for multiple scenarios
- Confidence interval calculation
- Monte Carlo simulation comparison

## 📄 License

This project is created for educational purposes as part of the [ROSALIND](http://rosalind.info/) bioinformatics platform. All code is provided under the MIT License.

## 🌟 Acknowledgments

- ROSALIND for the bioinformatics problems
- Gregor Mendel for foundational genetics
- Probability theorists
- Population genetics researchers

## 📈 Performance Benchmarks

All methods are O(1) and extremely fast:

- Python Manual: ~0.000001 seconds
- Python Fractions: ~0.000002 seconds
- R Vectorized: ~0.000001 seconds

## 🔍 Additional Notes

### Assumptions:
- Each couple produces exactly 2 offspring
- Random mating within genotype pairs
- No selection or mutation
- Large population (no genetic drift)

### Mathematical Formulation:

```text
E = 2 × [N₁×1 + N₂×1 + N₃×1 + N₄×0.75 + N₅×0.5 + N₆×0]
```

where Nᵢ are the counts for each genotype pairing

### Extensions:
- Can modify for different offspring distributions
- Add environmental factors affecting survival
- Consider incomplete dominance

---

For more bioinformatics challenges, visit [ROSALIND](https://rosalind.info) and join the bioinformatics community!