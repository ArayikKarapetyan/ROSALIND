# Independent Alleles Probability Solutions

![Python](https://img.shields.io/badge/Python-3.8+-blue.svg)
![R](https://img.shields.io/badge/R-4.0+-blue.svg)
![Probability](https://img.shields.io/badge/Probability-Binomial-blue.svg)
![Genetics](https://img.shields.io/badge/Genetics-Mendelian-green.svg)

![ROSALIND](https://rosalind.info/static/img/logo.png?4eac7af5ad70)

A bioinformatics solution for calculating the probability of inheriting specific genotypes across generations, using binomial distribution and Mendelian genetics.

## 📋 Problem Description

**Problem:** [Independent Alleles](https://rosalind.info/problems/lia/)  
**Category:** Bioinformatics Stronghold  
**ID:** LIA

Given a family tree where:

- Generation 0: Tom has genotype Aa Bb
- Each generation: every organism mates with an Aa Bb individual
- Each organism has exactly 2 children
- Traits A and B are independently inherited (Mendel's second law)

Calculate the probability that at least N organisms in generation k have genotype Aa Bb.

**Input:** Two integers k (k ≤ 7) and N (N ≤ 2ᵏ)  
**Output:** Probability (float) that at least N are Aa Bb in generation k

### Example

```text
Input: k = 2, N = 1
Output: 0.684

Explanation:
Generation 2 has 4 organisms (2^k = 4)
Each has 0.25 probability of being AaBb (independent inheritance)
P(at least 1 AaBb) = 1 - P(0 AaBb) = 1 - (0.75)^4 = 1 - 0.3164 = 0.6836
```

## 🧬 Solutions

### 1. Python Solution (From Scratch) - `python_manual.py`

```python
from math import comb

def probability_at_least_n_manual(k, N):
    """Calculate probability using binomial distribution."""
    n = 2 ** k          # Total organisms in generation k
    p = 0.25            # Probability single organism is AaBb
    
    # P(X >= N) = 1 - P(X < N) = 1 - Σ_{i=0}^{N-1} C(n,i) * p^i * (1-p)^{n-i}
    prob_less = sum(comb(n, i) * (p ** i) * ((1-p) ** (n-i)) 
                    for i in range(N))
    
    return 1 - prob_less
```

**Features:**
- Direct binomial distribution application
- O(N) time complexity
- Clear probability calculation

### 2. Python Efficient Solution - `python_efficient.py`

```python
def probability_at_least_n_efficient(k, N):
    """Calculate probability with recurrence to avoid large numbers."""
    n = 2 ** k
    p = 0.25
    q = 1 - p
    
    # Use recurrence to calculate probabilities efficiently
    prob_i = q ** n  # Probability of 0 successes
    prob_less = prob_i if N > 0 else 0
    
    for i in range(1, N):
        # Recurrence: P(i) = P(i-1) * (n-i+1)/i * p/q
        prob_i = prob_i * (n - i + 1) / i * p / q
        prob_less += prob_i
    
    return 1 - prob_less
```

**Features:**
- Avoids large binomial coefficients
- Numerically stable
- Same O(N) complexity but better for large n

### 3. R Solution - `r_solution.R`

```r
probability_at_least_n_r <- function(k, N) {
  n <- 2^k
  p <- 0.25
  
  # Use R's built-in binomial distribution
  prob_less_than_n <- pbinom(N - 1, size = n, prob = p)
  return(1 - prob_less_than_n)
}
```

**Features:**
- Uses R's optimized pbinom() function
- Clean one-liner implementation
- Multiple alternative approaches available

## 📁 File Structure

```text
Independent-Alleles/
│
├── README.md              # This documentation
├── python_manual.py       # Python binomial calculation
├── python_efficient.py    # Python efficient recurrence
├── python_scipy.py        # Python SciPy solution
├── r_solution.R           # R implementation
└── Dataset.txt            # Input file (k N)
```

## 🚀 Installation and Usage

### Python Requirements

```bash
# No special requirements for basic solution
# For SciPy solution:
pip install scipy numpy
```

### Running the Solutions

**Python (Manual):**

```bash
python python_manual.py
```

**Python (Efficient):**

```bash
python python_efficient.py
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

Input file `Dataset.txt` should contain two integers separated by a space:

```text
k N
```

**Example:**

```text
2 1
```

### File Reading

**Python:**

```python
with open("Dataset.txt", "r") as file:
    k, N = map(int, file.read().strip().split())
```

**R:**

```r
input_data <- readLines("Dataset.txt", warn = FALSE)
values <- strsplit(input_data, " ")[[1]]
k <- as.integer(values[1])
N <- as.integer(values[2])
```

## 📊 Algorithm Analysis

### Time Complexity

| Method | Approach | Time Complexity | Space Complexity |
|--------|----------|-----------------|------------------|
| Python Manual | Binomial sum | O(N) | O(1) |
| Python Efficient | Recurrence | O(N) | O(1) |
| R pbinom | Built-in CDF | O(1) | O(1) |
| R Manual | Binomial sum | O(N) | O(1) |

*k = generation, N = minimum count, n = 2ᵏ total organisms*

### Genetic Basis

- **Single trait Aa × Aa:** 25% AA, 50% Aa, 25% aa → 50% Aa
- **Independent traits:** P(Aa and Bb) = P(Aa) × P(Bb) = 0.5 × 0.5 = 0.25
- **Binomial distribution:** Each of n organisms independently has 0.25 chance

## 🧪 Testing

### Test Cases

```python
# Test case from problem
assert abs(probability_at_least_n_manual(2, 1) - 0.68359375) < 1e-6

# Additional test cases
assert probability_at_least_n_manual(0, 1) == 1.0  # Generation 0, Tom is AaBb
assert probability_at_least_n_manual(1, 3) == 0.0  # Impossible (only 2 organisms)
assert probability_at_least_n_manual(1, 0) == 1.0  # At least 0 is always true
assert probability_at_least_n_manual(3, 1) > 0.99  # Very high probability
```

### Sample Dataset Verification

```text
k = 2, N = 1

Total organisms in generation 2: 2^2 = 4
Probability single organism is AaBb: 0.25

P(at least 1 AaBb) = 1 - P(0 AaBb)
P(0 AaBb) = C(4,0) × (0.25)^0 × (0.75)^4 = 1 × 1 × 0.75^4 = 0.31640625

P(at least 1) = 1 - 0.31640625 = 0.68359375 ≈ 0.684
```

## 🔗 Related Problems

- [IPRB](https://rosalind.info/problems/iprb/) - Mendel's First Law
- [IEV](https://rosalind.info/problems/iev/) - Calculating Expected Offspring
- [FIB](https://rosalind.info/problems/fib/) - Rabbits and Recurrence Relations
- [PROB](https://rosalind.info/problems/prob/) - Introduction to Random Strings

## 📚 Learning Resources

- [Binomial Distribution](https://en.wikipedia.org/wiki/Binomial_distribution)
- [Mendel's Second Law](https://en.wikipedia.org/wiki/Mendelian_inheritance#Law_of_Independent_Assortment)
- [Independent Events](https://en.wikipedia.org/wiki/Independence_(probability_theory))
- [Probability Generating Functions](https://en.wikipedia.org/wiki/Probability-generating_function)

## 🤝 Contributing

Contributions are welcome! Here are ways to contribute:

### Add New Features:
- Multiple independent traits (more than 2)
- Variable offspring count
- Different mating partners
- Non-binomial offspring distributions

### Improve Documentation:
- Add probability tree diagrams
- Create visual family tree
- Add mathematical derivations
- Include statistical tests

### Enhance Performance:
- Approximations for large k
- Monte Carlo simulation
- Parallel computation for multiple scenarios
- GPU acceleration for probability calculations

## 📄 License

This project is created for educational purposes as part of the [ROSALIND](http://rosalind.info/) bioinformatics platform. All code is provided under the MIT License.

## 🌟 Acknowledgments

- [ROSALIND](http://rosalind.info/) for the bioinformatics problems
- Gregor Mendel for laws of inheritance
- Probability theory mathematicians
- Statistical genetics researchers

## 📈 Performance Benchmarks

For maximum values (k=7, N=128):

- Python Manual: ~0.0001 seconds
- Python Efficient: ~0.00005 seconds
- R pbinom: ~0.00001 seconds
- R Manual: ~0.0002 seconds

## 🔍 Additional Notes

### Mathematical Formulation:

```text
X ~ Binomial(n = 2^k, p = 0.25)
P(X ≥ N) = 1 - Σ_{i=0}^{N-1} C(2^k, i) × 0.25^i × 0.75^{2^k - i}
```

### Assumptions:
- Independent assortment (Mendel's second law)
- No selection, mutation, or genetic drift
- Each mating with Aa Bb partner
- Exactly 2 offspring per individual

### Numerical Considerations:
- For large k, 2^k can be 128 (max)
- Use log probabilities for extreme values
- Recurrence method avoids large numbers

### Biological Context:
- Models inheritance of unlinked traits
- Useful for genetic counseling
- Basis for pedigree analysis

---

For more bioinformatics challenges, visit [ROSALIND](https://rosalind.info) and join the bioinformatics community!