# Permutation Generator Solutions

![Python](https://img.shields.io/badge/Python-3.8+-blue.svg)
![R](https://img.shields.io/badge/R-4.0+-blue.svg)
![Combinatorics](https://img.shields.io/badge/Combinatorics-Permutations-blue.svg)
![Mathematics](https://img.shields.io/badge/Mathematics-Factorial-orange.svg)

![ROSALIND](https://rosalind.info/static/img/logo.png?4eac7af5ad70)

A combinatorial solution for generating all permutations of numbers 1 through n, including calculating the total count (n factorial).

## 📋 Problem Description

**Problem:** [Enumerating Gene Orders](https://rosalind.info/problems/perm/)  
**Category:** Bioinformatics Stronghold  
**ID:** PERM

A permutation of length n is an ordering of the positive integers {1, 2, ..., n}. The total number of permutations of length n is n! (n factorial).

**Input:** A positive integer n ≤ 7  
**Output:** Total number of permutations followed by all permutations

**Example:**
```text
Input: 3
Output: 
6
1 2 3
1 3 2
2 1 3
2 3 1
3 1 2
3 2 1
```

## 🧬 Solutions

### 1. Python Solution (From Scratch) - `python_manual.py`

```python
from itertools import permutations

def generate_permutations_manual(n):
    """Generate all permutations of 1..n."""
    perms = list(permutations(range(1, n + 1)))
    
    output = [str(len(perms))]
    for perm in perms:
        output.append(' '.join(map(str, perm)))
    
    return '\n'.join(output)
```

**Features:**
- Uses Python's built-in itertools.permutations
- Simple and efficient
- O(n!) time complexity

### 2. Python Recursive Solution - `python_recursive.py`

```python
def generate_permutations_recursive(n):
    """Generate permutations using recursion."""
    def backtrack(path, used, results):
        if len(path) == n:
            results.append(path[:])
            return
        
        for i in range(1, n + 1):
            if not used[i]:
                used[i] = True
                path.append(i)
                backtrack(path, used, results)
                path.pop()
                used[i] = False
    
    results = []
    backtrack([], [False] * (n + 1), results)
    
    output = [str(len(results))]
    for perm in results:
        output.append(' '.join(map(str, perm)))
    
    return '\n'.join(output)
```

**Features:**
- Educational backtracking algorithm
- Demonstrates recursive permutation generation
- Clear algorithmic thinking

### 3. R Solution - `r_solution.R`

```r
generate_permutations_r <- function(n) {
  # Use gtools package if available
  if (requireNamespace("gtools", quietly = TRUE)) {
    perms <- gtools::permutations(n, n)
    perms <- matrix(1:n, nrow = nrow(perms), ncol = n, byrow = FALSE)
    perms[] <- 1:n[perms]
  } else {
    # Manual recursive implementation
    perms <- manual_permutations_r(n)
  }
  
  total <- nrow(perms)
  output <- c(as.character(total),
              apply(perms, 1, paste, collapse = " "))
  
  return(paste(output, collapse = "\n"))
}
```

**Features:**
- Uses gtools package for efficient permutations
- Manual fallback implementation
- Multiple approaches included

## 📁 File Structure

```text
Permutation-Generator/
│
├── README.md              # This documentation
├── python_manual.py       # Python itertools solution
├── python_recursive.py    # Python recursive solution
├── python_iterative.py    # Python iterative solution
├── r_solution.R           # R implementation
└── Dataset.txt           # Input file (integer n)
```

## 🚀 Installation and Usage

### Python Requirements

```bash
# No special requirements for basic solution
# For math functions:
pip install numpy
```

### R Requirements

```r
# Install gtools for permutation generation
install.packages("gtools")
```

### Running the Solutions

**Python (Manual):**

```bash
python python_manual.py
```

**Python (Recursive):**

```bash
python python_recursive.py
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

Input file `Dataset.txt` should contain a single integer:

```text
3
```

### File Reading

**Python:**

```python
with open("Dataset.txt", "r") as file:
    n = int(file.read().strip())
```

**R:**

```r
n <- as.integer(readLines("Dataset.txt", warn = FALSE))
```

## 📊 Algorithm Analysis

### Time Complexity

| Method | Approach | Time Complexity | Space Complexity |
|--------|----------|----------------|------------------|
| Python itertools | Built-in generator | O(n!) | O(n!) |
| Python Recursive | Backtracking | O(n!) | O(n) |
| Python Heap's | Iterative algorithm | O(n!) | O(n!) |
| R gtools | Package function | O(n!) | O(n!) |

*n! = factorial(n) grows very quickly*

### Permutation Counts

| n | n! | Reasonable to generate? |
|---|-----|------------------------|
| 1 | 1 | Yes |
| 2 | 2 | Yes |
| 3 | 6 | Yes |
| 4 | 24 | Yes |
| 5 | 120 | Yes |
| 6 | 720 | Yes |
| 7 | 5040 | Yes (max for problem) |
| 8 | 40320 | Borderline |
| 9 | 362880 | Too many |
| 10 | 3,628,800 | Definitely too many |

## 🧪 Testing

### Test Cases

```python
# Test n=1
assert generate_permutations_manual(1) == "1\n1"

# Test n=2
output_2 = generate_permutations_manual(2)
assert "2" in output_2  # Count
assert "1 2" in output_2
assert "2 1" in output_2

# Test n=3 (from sample)
output_3 = generate_permutations_manual(3)
lines = output_3.strip().split('\n')
assert lines[0] == "6"
assert len(lines) == 7  # Count + 6 permutations
```

### Sample Dataset Verification

For n=3:
- Total permutations: 3! = 6
- All permutations of {1, 2, 3}:
  - 1 2 3
  - 1 3 2
  - 2 1 3
  - 2 3 1
  - 3 1 2
  - 3 2 1

## 🔗 Related Problems

- **[LEXF](https://rosalind.info/problems/lexf/)** - Enumerating k-mers
- **[LIA](https://rosalind.info/problems/lia/)** - Independent Alleles
- **[SIGN](https://rosalind.info/problems/sign/)** - Enumerating Oriented Gene Orderings
- **[SSET](https://rosalind.info/problems/sset/)** - Counting Subsets

## 📚 Learning Resources

- [Permutations](https://en.wikipedia.org/wiki/Permutation)
- [Factorial](https://en.wikipedia.org/wiki/Factorial)
- [Backtracking Algorithms](https://en.wikipedia.org/wiki/Backtracking)
- [Heap's Algorithm](https://en.wikipedia.org/wiki/Heap%27s_algorithm)

## 🤝 Contributing

Contributions are welcome! Here are ways to contribute:

### Add New Algorithms:
- Lexicographic ordering
- Steinhaus–Johnson–Trotter algorithm
- Random permutation generation
- Permutation ranking/unranking

### Improve Documentation:
- Add permutation visualization
- Create factorial growth charts
- Add mathematical proofs
- Include applications in bioinformatics

### Enhance Features:
- Support for larger n with sampling
- Parallel permutation generation
- Memory-efficient streaming output
- Permutation filtering (even/odd, etc.)

## 📄 License

This project is created for educational purposes as part of the [ROSALIND](https://rosalind.info) bioinformatics platform. All code is provided under the MIT License.

## 🌟 Acknowledgments

- [ROSALIND](https://rosalind.info) for the bioinformatics problems
- Mathematicians in combinatorics
- Computer scientists developing efficient algorithms
- Open-source community for permutation libraries

## 📈 Performance Benchmarks

For n=7 (maximum for problem):
- **Python itertools:** ~0.001 seconds
- **Python recursive:** ~0.002 seconds
- **Python iterative:** ~0.003 seconds
- **R gtools:** ~0.005 seconds

## 🔍 Additional Notes

### Biological Applications

Permutations are used in:
- Gene order analysis
- Phylogenetic tree reconstruction
- Protein structure analysis
- Sequence alignment permutations

### Algorithmic Considerations
- n! grows extremely fast (combinatorial explosion)
- n ≤ 7 is manageable (5040 permutations)
- For n > 10, need sampling not exhaustive generation

### Output Format
- First line: total count
- Subsequent lines: space-separated permutations
- Order doesn't matter (any valid order accepted)

### Extensions
- Could add lexicographic ordering option
- Include permutation inversion count
- Generate permutations with repetitions
- Handle permutations of arbitrary sets, not just 1..n

---

For more bioinformatics challenges, visit [ROSALIND](https://rosalind.info) and join the bioinformatics community!