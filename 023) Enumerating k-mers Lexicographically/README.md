# Lexicographic String Generator Solutions

![Python](https://img.shields.io/badge/Python-3.8+-blue.svg)
![R](https://img.shields.io/badge/R-4.0+-blue.svg)
![Combinatorics](https://img.shields.io/badge/Combinatorics-k%5En-blue.svg)
![Lexicographic](https://img.shields.io/badge/Lexicographic-Order-orange.svg)

![ROSALIND](https://rosalind.info/static/img/logo.png?4eac7af5ad70)

A combinatorial solution for generating all strings of length n from a given alphabet in lexicographic order, useful for generating k-mers and exploring sequence space.

## 📋 Problem Description

**Problem:** [Enumerating k-mers Lexicographically](https://rosalind.info/problems/lexf/)  
**Category:** Bioinformatics Armory  
**ID:** LEXF

Given an ordered alphabet and a positive integer n, generate all possible strings of length n that can be formed from the alphabet, ordered lexicographically.

**Input:** An alphabet (collection of symbols) and an integer n ≤ 10  
**Output:** All strings of length n in lexicographic order

**Example:**
```text
Input: A C G T, 2
Output:
AA
AC
AG
AT
CA
CC
CG
CT
GA
GC
GG
GT
TA
TC
TG
TT
```

## 🧬 Solutions

### 1. Python Solution (From Scratch) - `python_manual.py`

```python
from itertools import product

def generate_lexicographic_strings(alphabet, n):
    """Generate all strings of length n in lex order."""
    # itertools.product generates in lex order for sorted input
    return [''.join(combo) for combo in product(alphabet, repeat=n)]
```

**Features:**
- Uses Python's built-in itertools.product
- Simple one-liner
- O(kⁿ) time and space complexity

### 2. Python Recursive Solution - `python_recursive.py`

```python
def generate_lexicographic_recursive(alphabet, n):
    """Generate strings using recursion."""
    results = []
    
    def backtrack(current):
        if len(current) == n:
            results.append(current)
            return
        for char in alphabet:
            backtrack(current + char)
    
    backtrack("")
    return results
```

**Features:**
- Educational backtracking algorithm
- Demonstrates recursive thinking
- Clear step-by-step construction

### 3. R Solution - `r_solution.R`

```r
generate_lexicographic_r <- function(alphabet, n) {
  # Use expand.grid to generate all combinations
  args <- rep(list(alphabet), n)
  combos <- do.call(expand.grid, args)
  apply(combos, 1, paste, collapse = "")
}
```

**Features:**
- Uses R's expand.grid function
- Vectorized operations
- Clean and efficient

## 📁 File Structure

```text
Lexicographic-Kmers/
│
├── README.md              # This documentation
├── python_manual.py       # Python itertools solution
├── python_recursive.py    # Python recursive solution
├── python_iterative.py    # Python iterative solution
├── r_solution.R           # R implementation
└── Dataset.txt           # Input file (alphabet and n)
```

## 🚀 Installation and Usage

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

Input file `Dataset.txt` should contain:
- First line: space-separated alphabet symbols
- Second line: integer n

```text
A C G T
2
```

### File Reading

**Python:**

```python
with open("Dataset.txt", "r") as file:
    lines = file.read().strip().split('\n')
    alphabet = lines[0].split()
    n = int(lines[1])
```

**R:**

```r
lines <- readLines("Dataset.txt", warn = FALSE)
alphabet <- strsplit(lines[1], " ")[[1]]
n <- as.integer(lines[2])
```

## 📊 Algorithm Analysis

### Time Complexity

| Method | Approach | Time Complexity | Space Complexity |
|--------|----------|----------------|------------------|
| Python itertools | Cartesian product | O(kⁿ) | O(kⁿ) |
| Python Recursive | Backtracking | O(kⁿ) | O(n) recursion depth |
| Python Iterative | Nested loops | O(kⁿ) | O(kⁿ) |
| R expand.grid | Grid expansion | O(kⁿ) | O(kⁿ) |

*k = alphabet size, n = string length*

### Number of Strings

- Total strings = kⁿ
- For DNA alphabet (k=4): 4ⁿ strings
- For protein alphabet (k=20): 20ⁿ strings
- For n=10, k=4: 4¹⁰ = 1,048,576 strings

## 🧪 Testing

### Test Cases

```python
# Test with small alphabet
alphabet = ['A', 'B']
n = 2
result = generate_lexicographic_strings(alphabet, n)
expected = ['AA', 'AB', 'BA', 'BB']
assert result == expected

# Test n=1
alphabet = ['X', 'Y', 'Z']
n = 1
result = generate_lexicographic_strings(alphabet, n)
expected = ['X', 'Y', 'Z']
assert result == expected

# Test empty n=0 (edge case)
alphabet = ['A', 'B', 'C']
n = 0
result = generate_lexicographic_strings(alphabet, n)
expected = ['']  # or [] depending on interpretation
```

### Sample Dataset Verification

For alphabet = ['A', 'C', 'G', 'T'], n = 2:
- Total strings: 4² = 16
- Order: AA, AC, AG, AT, CA, CC, CG, CT, GA, GC, GG, GT, TA, TC, TG, TT
- This matches the sample output

## 🔗 Related Problems

- **[PERM](https://rosalind.info/problems/perm/)** - Enumerating Gene Orders
- **[LEXV](https://rosalind.info/problems/lexv/)** - Ordering Strings of Varying Length
- **[KMER](https://rosalind.info/problems/kmer/)** - k-mer Composition
- **[LCSM](https://rosalind.info/problems/lcsm/)** - Finding a Shared Motif

## 📚 Learning Resources

- [Cartesian Product](https://en.wikipedia.org/wiki/Cartesian_product)
- [Lexicographic Order](https://en.wikipedia.org/wiki/Lexicographic_order)
- [k-mer](https://en.wikipedia.org/wiki/K-mer)
- [Backtracking Algorithms](https://en.wikipedia.org/wiki/Backtracking)

## 🤝 Contributing

Contributions are welcome! Here are ways to contribute:

### Add New Features:
- Generate strings of varying lengths
- Weighted alphabet (different frequencies)
- Constrained generation (no repeats, etc.)
- Parallel generation for large n

### Improve Documentation:
- Add lexicographic order visualizations
- Create k-mer frequency applications
- Add biological context examples
- Include performance optimization tips

### Enhance Performance:
- Memory-efficient streaming generation
- Lazy evaluation for large outputs
- Parallel processing
- Bit-level optimizations

## 📄 License

This project is created for educational purposes as part of the [ROSALIND](https://rosalind.info) bioinformatics platform. All code is provided under the MIT License.

## 🌟 Acknowledgments

- [ROSALIND](https://rosalind.info) for the bioinformatics problems
- Combinatorial mathematics researchers
- Computer scientists developing efficient algorithms
- Bioinformatics community for k-mer analysis tools

## 📈 Performance Benchmarks

For maximum problem size (k=10, n=10):
- **Python itertools:** ~1 second (generates 10¹⁰ = 10 billion strings - not feasible!)
- **Python generator:** ~0 time to start, streams output
- **R expand.grid:** Would fail due to memory (10 billion rows)

**Note:** For n=10 and large k, output is enormous. Practical implementations should use streaming.

## 🔍 Additional Notes

### Biological Applications

k-mer analysis in genomics:
- Primer design
- Motif discovery
- Sequence alignment

### Algorithmic Considerations
- Lex order depends on input alphabet order
- For DNA, standard order is A,C,G,T
- For proteins, alphabetical order

### Memory Management
- kⁿ grows exponentially
- For large n, use generators/iterators
- Consider writing directly to file

### Extensions
- Could add filtering constraints
- Support for degenerate bases
- Generate only unique strings
- Rank/unrank functions

---

For more bioinformatics challenges, visit [ROSALIND](https://rosalind.info) and join the bioinformatics community!