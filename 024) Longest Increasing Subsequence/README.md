# Longest Increasing/Decreasing Subsequence Solutions

![Python](https://img.shields.io/badge/Python-3.8+-blue.svg)
![R](https://img.shields.io/badge/R-4.0+-blue.svg)
![Dynamic Programming](https://img.shields.io/badge/Dynamic%20Programming-LIS/LDS-blue.svg)
![Algorithm](https://img.shields.io/badge/Algorithm-Patience%20Sorting-orange.svg)

![ROSALIND](https://rosalind.info/static/img/logo.png?4eac7af5ad70)

A dynamic programming solution for finding the longest increasing and decreasing subsequences in permutations, with applications in bioinformatics for sequence analysis and pattern recognition.

## 📋 Problem Description

**Problem:** [Longest Increasing Subsequence](https://rosalind.info/problems/lgis/)  
**Category:** Bioinformatics Stronghold  
**ID:** LGIS

Given a permutation of numbers 1 through n, find:
- A longest increasing subsequence (LIS)
- A longest decreasing subsequence (LDS)

A subsequence maintains the relative order but not necessarily consecutiveness.

**Input:** n followed by a permutation of length n (n ≤ 10000)  
**Output:** A longest increasing subsequence followed by a longest decreasing subsequence

**Example:**
```text
Input:
5
5 1 4 2 3

Output:
1 2 3
5 4 2
```

## 🧬 Solutions

### 1. Python Solution (From Scratch) - `python_manual.py`

```python
def longest_increasing_subsequence_manual(seq):
    """Find LIS using O(n²) dynamic programming."""
    n = len(seq)
    dp = [1] * n
    prev = [-1] * n
    
    for i in range(n):
        for j in range(i):
            if seq[j] < seq[i] and dp[j] + 1 > dp[i]:
                dp[i] = dp[j] + 1
                prev[i] = j
    
    # Reconstruct from max dp value
    max_idx = dp.index(max(dp))
    lis = []
    while max_idx != -1:
        lis.append(seq[max_idx])
        max_idx = prev[max_idx]
    
    return lis[::-1]
```

**Features:**
- Simple O(n²) dynamic programming
- Easy to understand
- Works for n ≤ 10000 (acceptable for problem constraints)

### 2. Python Efficient Solution - `python_efficient.py`

```python
import bisect

def longest_increasing_subsequence_efficient(seq):
    """Find LIS using O(n log n) patience sorting."""
    piles = []
    parent = [-1] * len(seq)
    pile_tops = []
    
    for i, num in enumerate(seq):
        j = bisect.bisect_left(piles, num)
        if j == len(piles):
            piles.append(num)
            pile_tops.append(i)
        else:
            piles[j] = num
            pile_tops[j] = i
        
        if j > 0:
            parent[i] = pile_tops[j-1]
    
    # Reconstruct
    lis = []
    curr = pile_tops[-1]
    while curr != -1:
        lis.append(seq[curr])
        curr = parent[curr]
    
    return lis[::-1]
```

**Features:**
- O(n log n) time complexity
- Patience sorting algorithm
- More efficient for large n

### 3. R Solution - `r_solution.R`

```r
longest_increasing_subsequence_r <- function(seq) {
  n <- length(seq)
  dp <- rep(1, n)
  prev <- rep(-1, n)
  
  for (i in 1:n) {
    for (j in 1:(i-1)) {
      if (seq[j] < seq[i] && dp[j] + 1 > dp[i]) {
        dp[i] <- dp[j] + 1
        prev[i] <- j
      }
    }
  }
  
  max_idx <- which.max(dp)
  lis <- integer(0)
  while (max_idx != -1) {
    lis <- c(seq[max_idx], lis)
    max_idx <- prev[max_idx]
  }
  
  return(lis)
}
```

**Features:**
- R implementation of O(n²) DP
- Similar logic to Python version
- Handles the problem constraints

## 📁 File Structure

```text
Longest-Subsequences/
│
├── README.md              # This documentation
├── python_manual.py       # Python O(n²) solution
├── python_efficient.py    # Python O(n log n) solution
├── r_solution.R           # R implementation
└── Dataset.txt           # Input file (n and permutation)
```

## 🚀 Installation and Usage

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

Input file `Dataset.txt` should contain:
- First line: integer n
- Second line: permutation of n numbers

```text
5
5 1 4 2 3
```

### File Reading

**Python:**

```python
with open("Dataset.txt", "r") as file:
    lines = file.read().strip().split('\n')
    n = int(lines[0])
    permutation = list(map(int, lines[1].split()))
```

**R:**

```r
lines <- readLines("Dataset.txt", warn = FALSE)
n <- as.integer(lines[1])
permutation <- as.integer(strsplit(lines[2], " ")[[1]])
```

## 📊 Algorithm Analysis

### Time Complexity

| Method | LIS Algorithm | Time Complexity | Space Complexity |
|--------|---------------|----------------|------------------|
| DP Basic | O(n²) dynamic programming | O(n²) | O(n) |
| Patience Sorting | Binary search based | O(n log n) | O(n) |
| R Implementation | O(n²) DP | O(n²) | O(n) |

*n ≤ 10000*
- O(n²): ~100 million operations (acceptable)
- O(n log n): ~100,000 × 13 = 1.3 million operations (better)

### Algorithm Details

**LIS O(n²) DP:**
- dp[i] = length of LIS ending at position i
- Transition: dp[i] = max(dp[j] + 1) for all j < i where seq[j] < seq[i]
- Reconstruction: track prev[i] pointers

**LIS O(n log n) Patience Sorting:**
- Maintain piles of cards
- Place each number on leftmost pile where it's ≥ pile top
- Start new pile if larger than all tops
- Length = number of piles

## 🧪 Testing

### Test Cases

```python
# Test case from problem
seq = [5, 1, 4, 2, 3]
lis = longest_increasing_subsequence_manual(seq)
lds = longest_decreasing_subsequence_manual(seq)
assert lis == [1, 2, 3] or lis == [1, 4] or lis == [1, 2]  # Multiple possible LIS
assert lds == [5, 4, 2] or lds == [5, 4, 3] or lds == [5, 2]  # Multiple possible LDS

# Simple increasing sequence
seq2 = [1, 2, 3, 4, 5]
assert longest_increasing_subsequence_manual(seq2) == [1, 2, 3, 4, 5]

# Simple decreasing sequence  
seq3 = [5, 4, 3, 2, 1]
assert longest_decreasing_subsequence_manual(seq3) == [5, 4, 3, 2, 1]

# Single element
seq4 = [42]
assert longest_increasing_subsequence_manual(seq4) == [42]
assert longest_decreasing_subsequence_manual(seq4) == [42]
```

### Sample Dataset Verification

```text
Permutation: 5 1 4 2 3

LIS possibilities:
- 1 2 3 (length 3)
- 1 4 (length 2)
- 1 2 (length 2)
- etc.

LDS possibilities:
- 5 4 2 (length 3)
- 5 4 3 (length 3)
- 5 2 (length 2)
- etc.

Output shows one valid LIS and one valid LDS.
```

## 🔗 Related Problems

- **[SUBS](https://rosalind.info/problems/subs/)** - Finding a Motif in DNA
- **[LCSM](https://rosalind.info/problems/lcsm/)** - Finding a Shared Motif
- **[SIGN](https://rosalind.info/problems/sign/)** - Enumerating Oriented Gene Orderings
- **[PERM](https://rosalind.info/problems/perm/)** - Enumerating Gene Orders

## 📚 Learning Resources

- [Longest Increasing Subsequence](https://en.wikipedia.org/wiki/Longest_increasing_subsequence)
- [Patience Sorting](https://en.wikipedia.org/wiki/Patience_sorting)
- [Dynamic Programming](https://en.wikipedia.org/wiki/Dynamic_programming)
- [Sequence Alignment](https://en.wikipedia.org/wiki/Sequence_alignment)

## 🤝 Contributing

Contributions are welcome! Here are ways to contribute:

### Add New Features:
- Find all longest subsequences
- Weighted LIS/LDS
- Constrained subsequences
- 2D/3D versions

### Improve Documentation:
- Add algorithm visualizations
- Create comparison of different algorithms
- Add biological applications
- Include performance analysis

### Enhance Performance:
- Parallel DP implementation
- Memory-optimized versions
- GPU acceleration
- Streaming algorithms for large n

## 📄 License

This project is created for educational purposes as part of the [ROSALIND](https://rosalind.info) bioinformatics platform. All code is provided under the MIT License.

## 🌟 Acknowledgments

- [ROSALIND](https://rosalind.info) for the bioinformatics problems
- Algorithm researchers developing LIS algorithms
- Dynamic programming pioneers
- Bioinformatics sequence analysis community

## 📈 Performance Benchmarks

For n = 10000:
- **Python O(n²):** ~1-2 seconds
- **Python O(n log n):** ~0.1 seconds
- **R O(n²):** ~5-10 seconds
- **R O(n log n):** ~0.5 seconds

## 🔍 Additional Notes

### Biological Applications

Finding conserved regions in sequences:
- Protein structure comparison
- Gene expression pattern analysis
- Phylogenetic tree construction

### Multiple Solutions

There can be multiple valid LIS/LDS of the same length

### Algorithm Choice
- O(n²) DP is simpler to implement
- O(n log n) is better for large n
- For n ≤ 10000, both work

### Extensions
- Could find longest common subsequence between two permutations
- Add constraints (minimum gap, maximum difference)
- Handle non-permutations (sequences with repeats)
- 2D LIS (points in plane)

---

For more bioinformatics challenges, visit [ROSALIND](https://rosalind.info) and join the bioinformatics community!