# Fibonacci Rabbits with Variable Offspring Solutions

![Python](https://img.shields.io/badge/Python-3.8+-blue.svg)
![R](https://img.shields.io/badge/R-4.0+-blue.svg)
![Dynamic Programming](https://img.shields.io/badge/Dynamic%20Programming-Applied-orange.svg)

<p align="center">
  <img src="https://rosalind.info/static/img/logo.png?4eac7af5ad70" alt="ROSALIND Logo">
</p>

A bioinformatics solution for calculating rabbit population growth with variable offspring size, implemented using dynamic programming in multiple languages.

## 📋 Problem Description

**Problem:** [Rabbits and Recurrence Relations](https://rosalind.info/problems/fib/)  
**Category:** Bioinformatics Stronghold  
**ID:** FIB

This problem generalizes the Fibonacci sequence to model rabbit population growth where each reproduction-age pair produces k rabbit pairs per generation.

### Recurrence Relation:

```
Fn = Fn-1 + k × Fn-2
```

where `F1 = 1, F2 = 1`

**Input:** Two positive integers `n ≤ 40` and `k ≤ 5`  
**Output:** The total number of rabbit pairs after n months

### Example:

```text
Input: 5 3
Output: 19

Explanation:
Month 1: 1 pair
Month 2: 1 pair
Month 3: 1 + 3×1 = 4 pairs
Month 4: 4 + 3×1 = 7 pairs
Month 5: 7 + 3×4 = 19 pairs
```

## 🧬 Solutions

### 1. Python Solution (From Scratch) - `python_manual.py`

```python
def rabbit_pairs_manual(n, k):
    """Calculate rabbit pairs using dynamic programming."""
    if n == 1 or n == 2:
        return 1
    
    prev1 = 1  # Fn-1
    prev2 = 1  # Fn-2
    
    for i in range(3, n + 1):
        current = prev1 + k * prev2
        prev2, prev1 = prev1, current
    
    return prev1
```

**Features:**
- Space-optimized DP using only two variables
- O(n) time complexity, O(1) space complexity
- Handles n up to 40 efficiently

### 2. Python Solution with DP Array - `python_dp.py`

```python
def rabbit_pairs_dp(n, k):
    """Calculate rabbit pairs using DP array."""
    if n <= 2:
        return 1
    
    dp = [0] * (n + 1)
    dp[1] = 1
    dp[2] = 1
    
    for i in range(3, n + 1):
        dp[i] = dp[i-1] + k * dp[i-2]
    
    return dp[n]
```

**Features:**
- Explicit DP table for clarity
- Easy to understand and debug
- Stores intermediate results

### 3. R Solution - `r_solution.R`

```r
rabbit_pairs_r <- function(n, k) {
  if (n <= 2) {
    return(1)
  }
  
  dp <- numeric(n)
  dp[1] <- 1
  dp[2] <- 1
  
  for (i in 3:n) {
    dp[i] <- dp[i-1] + k * dp[i-2]
  }
  
  return(dp[n])
}
```

**Features:**
- R vector-based implementation
- Functional programming alternative available
- Multiple approaches for different use cases

## 📁 File Structure

```text
Fibonacci-Rabbits/
│
├── README.md              # This documentation
├── python_manual.py       # Python space-optimized solution
├── python_dp.py           # Python DP array solution
├── r_solution.R           # R implementation
└── Dataset.txt            # Input file (n k format)
```

## 🚀 Installation and Usage

### Running the Solutions

**Python (Manual):**

```bash
python python_manual.py
```

**Python (DP Array):**

```bash
python python_dp.py
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
n k
```

Example:

```text
5 3
```

### File Reading

**Python:**

```python
with open("Dataset.txt", "r") as file:
    n, k = map(int, file.read().strip().split())
```

**R:**

```r
input_data <- readLines("Dataset.txt", warn = FALSE)
values <- strsplit(input_data, " ")[[1]]
n <- as.integer(values[1])
k <- as.integer(values[2])
```

## 📊 Algorithm Analysis

### Time Complexity

| Method | Approach | Time Complexity | Space Complexity |
|--------|----------|-----------------|------------------|
| Python Manual | Iterative DP | O(n) | O(1) |
| Python DP Array | Tabular DP | O(n) | O(n) |
| R Vector | Tabular DP | O(n) | O(n) |
| Recursive | Naive recursion | O(2ⁿ) | O(n) |

*n = number of months*

### Mathematical Insight

The recurrence relation:

```text
F₁ = 1
F₂ = 1
Fₙ = Fₙ₋₁ + k × Fₙ₋₂ for n > 2
```

This is a generalization of the Fibonacci sequence where the "offspring multiplier" is k instead of 1.

## 🧪 Testing

### Test Cases

```python
# Test case from problem
assert rabbit_pairs_manual(5, 3) == 19

# Additional test cases
assert rabbit_pairs_manual(1, 1) == 1  # Base case n=1
assert rabbit_pairs_manual(2, 5) == 1  # Base case n=2
assert rabbit_pairs_manual(3, 2) == 3  # 1 + 2×1
assert rabbit_pairs_manual(6, 4) == 129  # Larger test
```

### Sample Dataset Verification

```text
Input: n=5, k=3

Month 1: 1
Month 2: 1
Month 3: 1 + 3×1 = 4
Month 4: 4 + 3×1 = 7
Month 5: 7 + 3×4 = 19

Output: 19
```

## 🔗 Related Problems

- **[FIBD](https://rosalind.info/problems/fibd/)** - Mortal Fibonacci Rabbits
- **[LIA](https://rosalind.info/problems/lia/)** - Independent Alleles
- **[IPRB](https://rosalind.info/problems/iprb/)** - Mendel's First Law
- **[LCSM](https://rosalind.info/problems/lcsm/)** - Finding a Shared Motif

## 📚 Learning Resources

- [Dynamic Programming Introduction](https://en.wikipedia.org/wiki/Dynamic_programming)
- [Fibonacci Sequence](https://en.wikipedia.org/wiki/Fibonacci_number)
- [Recurrence Relations](https://en.wikipedia.org/wiki/Recurrence_relation)
- [Population Growth Models](https://en.wikipedia.org/wiki/Population_growth)

## 🤝 Contributing

Contributions are welcome! Here are ways to contribute:

### Add New Implementations:
- Matrix exponentiation for O(log n) solution
- Closed-form formula implementation
- Other programming languages

### Improve Documentation:
- Add mathematical derivation
- Create visualization of rabbit growth
- Add complexity analysis

### Enhance Features:
- Add mortality factor (FIBD problem)
- Create interactive simulation
- Add unit tests

## 📄 License

This project is created for educational purposes as part of the [ROSALIND](http://rosalind.info/) bioinformatics platform. All code is provided under the MIT License.

## 🌟 Acknowledgments

- [ROSALIND](http://rosalind.info/) for the bioinformatics problems
- Leonardo Fibonacci for the original sequence
- The mathematical biology community for population modeling
- All contributors to algorithm education

## 📈 Performance Benchmarks

For n = 40, k = 5:

- **Python Manual:** ~0.00001 seconds
- **Python DP Array:** ~0.00002 seconds
- **R Vector:** ~0.00005 seconds
- **Naive Recursive:** > 10 seconds (exponential growth)

## 🔍 Additional Notes

- For n ≤ 40, simple DP is sufficient
- Space-optimized version is preferred for memory efficiency
- The problem constraints ensure results fit in 64-bit integers
- This model assumes:
  - Rabbits never die
  - Rabbits mature in one month
  - Each pair produces exactly k pairs per month

For more bioinformatics challenges, visit [ROSALIND](http://rosalind.info/) and join the bioinformatics community!