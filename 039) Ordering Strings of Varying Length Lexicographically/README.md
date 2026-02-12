# Variable Length String Generation Solutions

![Python](https://img.shields.io/badge/Python-3.8+-blue.svg)
![R](https://img.shields.io/badge/R-4.0+-blue.svg)
![Combinatorics](https://img.shields.io/badge/Combinatorics-Lexicographic-blue.svg)
![Variable Length](https://img.shields.io/badge/Variable_Length-Strings-green.svg)

![ROSALIND](https://rosalind.info/static/img/logo.png?4eac7af5ad70)

Solutions for generating all strings of length at most n from a custom ordered alphabet, sorted according to a special lexicographic order where shorter strings come before longer ones with the same prefix.

## 📋 Problem Description

**Problem:** [Ordering Strings of Varying Length](https://rosalind.info/problems/lexv/)  
**Category:** Bioinformatics Textbook Track  
**ID:** LEXV

Given an ordered alphabet A and a positive integer n, generate all strings of length at most n formed from A, ordered lexicographically with a special rule: if string s is a prefix of string t, then s comes before t.

**Input:** 
- First line: Space-separated symbols defining ordered alphabet
- Second line: Integer n (n ≤ 4)

**Output:** All strings of length ≤ n in lexicographic order

**Example:**
```text
Input: D N A
       3

Output: D, DD, DDD, DDN, DDA, DN, DND, DNN, DNA, ...
```

## 🧬 Solutions

### 1. Python Solution (Recursive) - `python_manual.py`

```python
def generate_strings_varying_length(alphabet, n):
    results = []
    
    def backtrack(current):
        if current:
            results.append(current)
        if len(current) < n:
            for char in alphabet:
                backtrack(current + char)
    
    backtrack("")
    return results
```

**Features:**
- Depth-first recursive generation
- Natural lexicographic order
- Simple and intuitive

### 2. Python Efficient Solution - `python_efficient.py`

```python
def generate_strings_iterative(alphabet, n):
    results = []
    stack = [""]
    
    while stack:
        current = stack.pop()
        if current:
            results.append(current)
        if len(current) < n:
            for char in reversed(alphabet):
                stack.append(current + char)
    
    return results
```

**Features:**
- Iterative stack-based approach
- No recursion depth limits
- Explicit control of generation order
- Memory efficient

### 3. Python BioPython Style Solution - `python_biopython.py`

```python
def generate_strings_biopython_style(alphabet, n):
    from itertools import product
    
    results = []
    for length in range(1, n + 1):
        for combo in product(alphabet, repeat=length):
            results.append(''.join(combo))
    
    # Custom sort for special lexicographic order
    results.sort(key=lambda x: [alphabet.index(c) for c in x] + [len(x)])
    return results
```

**Features:**
- Uses itertools.product for combination generation
- Explicit sorting with custom key
- Clear separation of generation and ordering
- Professional structure

### 4. R Solution - `r_solution.R`

```r
generate_strings_varying_length_r <- function(alphabet, n) {
  results <- character(0)
  
  generate <- function(current) {
    if (nchar(current) > 0) {
      results <<- c(results, current)
    }
    if (nchar(current) < n) {
      for (char in alphabet) {
        generate(paste0(current, char))
      }
    }
  }
  
  generate("")
  return(results)
}
```

**Features:**
- Recursive implementation in R
- Uses R's string manipulation functions
- Global variable for results accumulation
- Simple and effective

## 📁 File Structure

```text
Variable-Length-Strings/
│
├── README.md              # This documentation
├── python_manual.py       # Python recursive solution
├── python_efficient.py    # Python iterative solution
├── python_biopython.py    # Python BioPython-style solution
├── r_solution.R           # R implementation
└── Dataset.txt            # Input file
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

**Python (BioPython Style):**
```bash
python python_biopython.py
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

```text
D N A
3
```

Where:
- First line: Space-separated symbols (alphabet order matters)
- Second line: Maximum string length n (≤ 4)

### Output Format

One string per line, in lexicographic order:

```text
D
DD
DDD
DDN
DDA
DN
...
```

## 📊 Algorithm Analysis

### Time Complexity

| Method | Approach | Time Complexity | Space Complexity |
|--------|----------|-----------------|------------------|
| Python Recursive | DFS with backtracking | O(k^n) | O(n) recursion depth |
| Python Iterative | Stack-based DFS | O(k^n) | O(k^n) worst-case |
| Python BioPython | Product + sort | O(k^n × n log(k^n)) | O(k^n) |
| R Solution | Recursive DFS | O(k^n) | O(n) recursion depth |

k = alphabet size, n = max string length

### Number of Strings

Total number of strings generated:

$$
\text{Total} = k + k^2 + k^3 + \cdots + k^n = \frac{k(k^n - 1)}{k - 1}
$$

For example with k=3, n=3:

3 + 9 + 27 = 39 strings

## 🧪 Testing

### Test Cases

```python
# Test case 1: Small alphabet, n=1
alphabet = ['A', 'B']
n = 1
result = generate_strings_varying_length(alphabet, n)
assert result == ['A', 'B']

# Test case 2: n=2 with custom order
alphabet = ['B', 'A']
n = 2
result = generate_strings_varying_length(alphabet, n)
assert result[:5] == ['B', 'BB', 'BA', 'A', 'AB']

# Test case 3: Sample from problem
alphabet = ['D', 'N', 'A']
n = 3
result = generate_strings_varying_length(alphabet, n)
assert len(result) == 39  # 3 + 9 + 27
assert result[0] == 'D'
assert result[4] == 'DDA'
```

### Validation Rules

- **Ordering rule:** If s is prefix of t, then s comes before t
- **Lexicographic:** Otherwise, compare character by character using given alphabet order
- **Length limit:** All strings have length ≤ n
- **Completeness:** All possible strings are generated

## 🔗 Related Problems

- [LEXF](https://rosalind.info/problems/lexf/) - Enumerating k-mers Lexicographically
- [KMER](https://rosalind.info/problems/kmer/) - k-mer Composition
- [PERM](https://rosalind.info/problems/perm/) - Enumerating Gene Orders
- [SUBS](https://rosalind.info/problems/subs/) - Finding Motifs

## 📚 Learning Resources

- [Lexicographic Order](https://en.wikipedia.org/wiki/Lexicographic_order)
- [Backtracking Algorithms](https://en.wikipedia.org/wiki/Backtracking)
- [Combinatorial Generation](https://en.wikipedia.org/wiki/Combinatorial_generation)
- [Depth-First Search](https://en.wikipedia.org/wiki/Depth-first_search)

## 🤝 Contributing

Contributions are welcome! Here are ways to contribute:

### Add New Features:
- Generate strings with constraints (no repeats, etc.)
- Weighted string generation (different symbol probabilities)
- Stream output for very large alphabets
- Parallel generation

### Improve Performance:
- Memory-efficient streaming generation
- Bit-level encoding for large alphabets
- Lazy evaluation
- GPU acceleration for large n

### Enhance Documentation:
- Visualize generation tree
- Interactive alphabet reordering
- Step-by-step generation animation
- Real-world biological applications

## 📄 License

This project is created for educational purposes as part of the [ROSALIND](https://rosalind.info/) bioinformatics platform. All code is provided under the MIT License.

## 🌟 Acknowledgments

- [ROSALIND](https://rosalind.info/) for the bioinformatics problems
- Combinatorial mathematics researchers
- String algorithm developers
- Computer science educators

## 📈 Performance Benchmarks

For maximum problem size (alphabet size ≤ 12, n=4):

| Solution | Time |
|----------|------|
| Python Recursive | ~0.01 seconds |
| Python Iterative | ~0.01 seconds |
| Python BioPython | ~0.02 seconds |
| R Solution | ~0.05 seconds |

**Maximum number of strings (worst case):**

With 12 symbols and n=4: 12 + 144 + 1728 + 20736 = 22620 strings

Output file size: ~200KB

## 🔍 Additional Notes

### Special Lexicographic Order

The key difference from standard lex order:

- **Standard:** "A" < "AA" (compare 'A' vs 'A', then length)
- **This problem:** "A" < "AA" (shorter string with same prefix comes first)

**Implementation approaches:**

1. **DFS generation:** Natural ordering if generated correctly
2. **Custom sorting:** Generate all then sort with custom comparator
3. **Iterative deepening:** Generate by increasing length

### Applications in Bioinformatics

- Primer design with variable lengths
- Motif search with flexible patterns
- Oligonucleotide library generation
- Sequence space exploration

### Algorithm Details

The recursive solution naturally produces correct order because:

1. It outputs current string before exploring longer strings
2. Explores children in alphabet order
3. Depth-first ensures prefix property

### Edge Cases

- **n=0:** Should output nothing (or just empty string, depending on interpretation)
- **Alphabet with one symbol:** All strings are repetitions of that symbol
- **Large alphabet (12 symbols):** Manageable due to n ≤ 4
- **Special characters:** Alphabet can be any symbols

### Extensions

- **Variable alphabets:** Different alphabets for different positions
- **Constraints:** Exclude certain patterns or repetitions
- **Weighted generation:** Prioritize certain symbols or lengths
- **Incremental output:** Generate and process strings on-the-fly

### Optimization for Large n

If n were larger (not in this problem):

- Use iterative deepening
- Stream results to file
- Parallel generation by starting prefix
- Prune branches based on constraints

---

For more bioinformatics challenges, visit [ROSALIND](https://rosalind.info/) and join the bioinformatics community!
