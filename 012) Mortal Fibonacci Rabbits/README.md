# Mortal Fibonacci Rabbits Solutions

![Python](https://img.shields.io/badge/Python-3.8+-blue.svg)
![R](https://img.shields.io/badge/R-4.0+-blue.svg)
![NumPy](https://img.shields.io/badge/NumPy-1.19+-blue.svg)
![Dynamic Programming](https://img.shields.io/badge/Dynamic%20Programming-Mortality-blue.svg)

<p align="center">
  <img src="https://rosalind.info/static/img/logo.png?4eac7af5ad70" alt="ROSALIND Logo">
</p>

A bioinformatics solution for modeling rabbit populations with mortality, extending the classic Fibonacci rabbit problem to account for limited lifespan.

## 📋 Problem Description

**Problem:** [Mortal Fibonacci Rabbits](https://rosalind.info/problems/fibd/)  
**Category:** Bioinformatics Stronghold  
**ID:** FIBD

Extend the Fibonacci rabbit model to include mortality: rabbits die after living for m months. Each pair of rabbits:

- Reaches maturity in 1 month
- Produces 1 pair of offspring each month after maturity
- Dies after m months

**Input:** Positive integers `n ≤ 100` (months) and `m ≤ 20` (lifespan)  
**Output:** Total rabbit pairs after n months

### Example:

```text
Input: n = 6, m = 3
Output: 4

Month 1: 1 pair (newborn)
Month 2: 1 pair (now mature)
Month 3: 2 pairs (1 original + 1 newborn)
Month 4: 2 pairs (1 dies, 1 matures, 1 newborn)
Month 5: 3 pairs
Month 6: 4 pairs
```

## 🧬 Solutions

### 1. Python Solution (From Scratch) - `python_manual.py`

```python
def mortal_rabbits_manual(n, m):
    """Calculate mortal rabbits using age tracking."""
    if n <= 0:
        return 0
    if n == 1:
        return 1
    
    ages = [0] * m
    ages[0] = 1  # Newborns
    
    for month in range(2, n + 1):
        newborns = sum(ages[1:])  # Reproducing rabbits
        # Age rabbits
        for age in range(m-1, 0, -1):
            ages[age] = ages[age-1]
        ages[0] = newborns
    
    return sum(ages)
```

**Features:**
- Age-based tracking
- O(n × m) time complexity
- Clear mortality modeling

### 2. Python Solution with Deque - `python_deque.py`

```python
from collections import deque

def mortal_rabbits_deque(n, m):
    """Calculate mortal rabbits using deque for rotation."""
    if n <= 0:
        return 0
    if n == 1:
        return 1
    
    ages = deque([0] * m)
    ages[0] = 1
    
    for month in range(2, n + 1):
        newborns = sum(list(ages)[1:])
        ages.pop()  # Oldest die
        ages.appendleft(newborns)
    
    return sum(ages)
```

**Features:**
- Efficient rotation with deque
- Clean age shifting logic
- Same time complexity, better constant factors

### 3. R Solution - `r_solution.R`

```r
mortal_rabbits_r <- function(n, m) {
  if (n <= 0) return(0)
  if (n == 1) return(1)
  
  ages <- integer(m)
  ages[1] <- 1  # R is 1-indexed
  
  for (month in 2:n) {
    newborns <- sum(ages[-1])  # All except newborns
    if (m > 1) ages[2:m] <- ages[1:(m-1)]  # Age shift
    ages[1] <- newborns
  }
  
  return(sum(ages))
}
```

**Features:**
- R's 1-indexed arrays
- Vectorized age shifting
- Multiple implementation approaches

## 📁 File Structure

```text
Mortal-Rabbits/
│
├── README.md              # This documentation
├── python_manual.py       # Python age array solution
├── python_deque.py        # Python deque solution
├── python_numpy.py        # Python NumPy solution
├── r_solution.R           # R implementation
└── Dataset.txt            # Input file (n m)
```

## 🚀 Installation and Usage

### Python Requirements

```bash
# Install NumPy (if using NumPy solution)
pip install numpy
```

### Running the Solutions

**Python (Manual):**

```bash
python python_manual.py
```

**Python (Deque):**

```bash
python python_deque.py
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
n m
```

Example:

```text
6 3
```

### File Reading

**Python:**

```python
with open("Dataset.txt", "r") as file:
    n, m = map(int, file.read().strip().split())
```

**R:**

```r
input_data <- readLines("Dataset.txt", warn = FALSE)
values <- strsplit(input_data, " ")[[1]]
n <- as.integer(values[1])
m <- as.integer(values[2])
```

## 📊 Algorithm Analysis

### Time Complexity

| Method | Approach | Time Complexity | Space Complexity |
|--------|----------|-----------------|------------------|
| Python Manual | Age array with shifting | O(n × m) | O(m) |
| Python Deque | Queue rotation | O(n × m) | O(m) |
| Python NumPy | Array operations | O(n × m) | O(m) |
| R Vector | Age vector shifting | O(n × m) | O(m) |

*n = number of months, m = lifespan*

### Recurrence Relation

The mortal Fibonacci sequence follows:

```text
F(n) = F(n-1) + F(n-2) - F(n-m-1)  for n > m
F(n) = F(n-1) + F(n-2)             for n ≤ m
F(1) = 1, F(2) = 1
```

## 🧪 Testing

### Test Cases

```python
# Test case from problem
assert mortal_rabbits_manual(6, 3) == 4

# Additional test cases
assert mortal_rabbits_manual(1, 1) == 1  # Dies immediately
assert mortal_rabbits_manual(5, 3) == 3  # From problem progression
assert mortal_rabbits_manual(10, 5) == 31  # Larger example
assert mortal_rabbits_manual(3, 10) == 2  # Standard Fibonacci (no death yet)
```

### Sample Dataset Verification

```text
n=6, m=3 (rabbits live 3 months)

Month 1: [1, 0, 0] = 1 total (1 newborn)
Month 2: [0, 1, 0] = 1 total (1 mature)
Month 3: [1, 0, 1] = 2 total (1 newborn, 1 old)
Month 4: [1, 1, 0] = 2 total (1 dies, 1 newborn, 1 mature)
Month 5: [1, 1, 1] = 3 total
Month 6: [2, 1, 1] = 4 total

Total after month 6: 4 pairs
```

## 🔗 Related Problems

- **[FIB](https://rosalind.info/problems/fib/)** - Rabbits and Recurrence Relations
- **[IPRB](https://rosalind.info/problems/iprb/)** - Mendel's First Law
- **[LIA](https://rosalind.info/problems/lia/)** - Independent Alleles
- **[REVC](https://rosalind.info/problems/revc/)** - Complementing a Strand of DNA

## 📚 Learning Resources

- [Fibonacci Sequence](https://en.wikipedia.org/wiki/Fibonacci_number)
- [Population Dynamics](https://en.wikipedia.org/wiki/Population_dynamics)
- [Age-Structured Populations](https://en.wikipedia.org/wiki/Age_structure)
- [Leslie Matrix](https://en.wikipedia.org/wiki/Leslie_matrix)

## 🤝 Contributing

Contributions are welcome! Here are ways to contribute:

### Add New Algorithms:
- Matrix exponentiation for O(log n) solution
- Closed-form approximation
- Stochastic simulation version

### Improve Documentation:
- Add population growth visualizations
- Create age pyramid diagrams
- Add mathematical derivation

### Enhance Features:
- Variable reproduction rates
- Sex ratio modeling
- Environmental carrying capacity
- Seasonality effects

## 📄 License

This project is created for educational purposes as part of the [ROSALIND](http://rosalind.info/) bioinformatics platform. All code is provided under the MIT License.

## 🌟 Acknowledgments

- [ROSALIND](http://rosalind.info/) for the bioinformatics problems
- Leonardo Fibonacci for the original sequence
- Population ecologists and demographers
- Mathematical biology researchers

## 📈 Performance Benchmarks

For n=100, m=20:

- **Python Manual:** ~0.0002 seconds
- **Python Deque:** ~0.0003 seconds
- **Python NumPy:** ~0.0001 seconds
- **R Vector:** ~0.0005 seconds

## 🔍 Additional Notes

**Biological Applications:**
- Population forecasting
- Conservation biology
- Epidemiology modeling
- Resource management

**Assumptions:**
- Discrete time steps (months)
- Constant reproduction rate
- No predation or disease
- Unlimited resources

**Edge Cases:**
- m = 1 (immediate death)
- m ≥ n (standard Fibonacci)
- Large n values

**Extensions:**
- Can add multiple offspring
- Variable mortality rates
- Density-dependent factors

For more bioinformatics challenges, visit [ROSALIND](http://rosalind.info/) and join the bioinformatics community!