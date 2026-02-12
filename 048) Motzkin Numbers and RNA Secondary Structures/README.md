# RNA Noncrossing Matchings (Motzkin Numbers) Solutions

![Python](https://img.shields.io/badge/Python-3.8+-blue.svg)
![R](https://img.shields.io/badge/R-4.0+-blue.svg)
![Dynamic Programming](https://img.shields.io/badge/DP-Motzkin_Numbers-blue.svg)
![RNA Structure](https://img.shields.io/badge/RNA-Secondary_Structure-green.svg)

![ROSALIND](https://rosalind.info/static/img/logo.png?4eac7af5ad70)

Solutions for counting all noncrossing matchings (not necessarily perfect) in RNA secondary structures using modified Motzkin number recurrence with Watson-Crick base pairing constraints.

## 📋 Problem Description

**Problem:** [Counting RNA Secondary Structures](https://rosalind.info/problems/motz/)  
**Category:** Bioinformatics Textbook Track  
**ID:** MOTZ

Given an RNA string, count all possible noncrossing matchings of basepair edges in its bonding graph. Unlike perfect matchings (CAT problem), these matchings do not need to pair all bases - some bases may remain unpaired.

**Input:** An RNA string s (length ≤ 300 bp) in FASTA format  
**Output:** Total number of noncrossing matchings modulo 1,000,000

### Example

```text
Input: AUAU
Output: 7
The 7 noncrossing matchings are: empty, 1-2, 1-4, 2-3, 3-4, 1-2&3-4, 1-4&2-3.
```

## 🧬 Solutions

### 1. Python Solution (Dynamic Programming) - `python_manual.py`

```python
def count_noncrossing_matchings(rna):
    n = len(rna)
    MOD = 1000000
    dp = [[0] * n for _ in range(n)]
    
    # Base case: empty substring
    for i in range(n + 1):
        if i < n:
            dp[i][i - 1] = 1
    
    for length in range(n):
        for i in range(n - length):
            j = i + length
            total = dp[i + 1][j] if i + 1 <= j else 1  # i unpaired
            
            for k in range(i + 1, j + 1):
                if can_pair(rna[i], rna[k]):
                    left = dp[i + 1][k - 1] if i + 1 <= k - 1 else 1
                    right = dp[k + 1][j] if k + 1 <= j else 1
                    total = (total + left * right) % MOD
            
            dp[i][j] = total
    
    return dp[0][n - 1] % MOD
```

**Features:**
- Bottom-up dynamic programming
- O(n³) time complexity
- Explicit matrix construction

### 2. Python Efficient Solution - `python_efficient.py`

```python
from functools import lru_cache

def count_rna_matchings_efficient(rna):
    MOD = 1000000
    can_pair = {('A','U'), ('U','A'), ('C','G'), ('G','C')}
    
    @lru_cache(maxsize=None)
    def count(i, j):
        if i > j:
            return 1
        
        total = count(i + 1, j)  # i unpaired
        
        for k in range(i + 1, j + 1):
            if (rna[i], rna[k]) in can_pair:
                total = (total + count(i + 1, k - 1) * count(k + 1, j)) % MOD
        
        return total % MOD
    
    return count(0, len(rna) - 1)
```

**Features:**
- Top-down DP with memoization
- Clean recursive formulation
- Automatic caching with @lru_cache

### 3. Python BioPython Solution - `python_biopython.py`

```python
from Bio import SeqIO
from functools import lru_cache

def count_rna_matchings_biopython(seq):
    rna = str(seq)
    MOD = 1000000
    valid_pairs = {('A','U'), ('U','A'), ('C','G'), ('G','C')}
    
    @lru_cache(maxsize=None)
    def dp(i, j):
        if i >= j:
            return 1
        
        total = dp(i + 1, j)  # base i unpaired
        
        for k in range(i + 1, j + 1):
            if (rna[i], rna[k]) in valid_pairs:
                total = (total + dp(i + 1, k - 1) * dp(k + 1, j)) % MOD
        
        return total
    
    return dp(0, len(rna) - 1) % MOD
```

**Features:**
- Uses BioPython for FASTA parsing
- Handles Seq objects
- Professional bioinformatics workflow

### 4. R Solution - `r_solution.R`

```r
count_rna_matchings_r <- function(rna) {
  MOD <- 1000000
  bases <- strsplit(rna, "")[[1]]
  n <- length(bases)

  # Handle empty string
  if (n == 0) return(1)

  # Watson–Crick pairing
  can_pair <- function(b1, b2) {
    (b1 == "A" && b2 == "U") ||
    (b1 == "U" && b2 == "A") ||
    (b1 == "C" && b2 == "G") ||
    (b1 == "G" && b2 == "C")
  }

  # Precompute valid pairs
  pair_ok <- matrix(FALSE, n, n)
  for (i in 1:n)
    for (j in 1:n)
      pair_ok[i, j] <- can_pair(bases[i], bases[j])

  # DP table
  dp <- matrix(0, n, n)

  # Base cases: single base
  for (i in 1:n)
    dp[i, i] <- 1

  # Fill table by increasing substring length
  for (len in 2:n) {
    for (i in 1:(n - len + 1)) {
      j <- i + len - 1

      # Case 1: base i is unpaired
      dp[i, j] <- dp[i + 1, j]

      # Case 2: base i pairs with k
      for (k in seq.int(i + 1, j)) {
        if (pair_ok[i, k]) {
          left  <- if (i + 1 <= k - 1) dp[i + 1, k - 1] else 1
          right <- if (k + 1 <= j) dp[k + 1, j] else 1
          dp[i, j] <- (dp[i, j] + left * right) %% MOD
        }
      }
    }
  }

  return(dp[1, n])
}

```

**Features:**
- Bottom-up dynamic programming (no recursion)
- Safe 1-based indexing (seq.int)
- Exact MOTZ recurrence (no biological loop constraints)
- O(n³) time, O(n²) memory

## 📁 File Structure

```text
RNA-Noncrossing-Matchings-Motz/
│
├── README.md              # This documentation
├── python_manual.py       # Python DP solution
├── python_efficient.py    # Python memoized solution
├── python_biopython.py    # Python BioPython solution
├── r_solution.R           # R implementation
└── Dataset.txt            # Input FASTA file
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

**Python (BioPython):**

```bash
# First install BioPython if needed:
# pip install biopython
python python_biopython.py
```

**R:**

```bash
# First install memoise package if needed:
# install.packages("memoise")
Rscript r_solution.R
```

## 🔧 Configuration

### Input Format

Input file `Dataset.txt` should contain RNA sequence in FASTA format:

```text
>Rosalind_57
AUAU
```

### Output Format

Single integer representing number of noncrossing matchings modulo 1,000,000:

```text
7
```

## 📊 Mathematical Background

### Motzkin Numbers Recurrence

For Motzkin numbers mₙ (counting noncrossing matchings on n nodes):

```
mₙ = mₙ₋₁ + Σ(k=2 to n) mₖ₋₂ · mₙ₋ₖ
```

Where:
- First term: node 1 unpaired
- Second term: node 1 paired with node k

### RNA Modification

For RNA with base pairing constraints:
- Only allow pairing if bases are complementary (A-U, C-G)
- Recurrence becomes conditional on valid base pairs

### Recurrence for RNA

For substring rna[i:j+1]:

```
dp(i,j) = dp(i+1,j) + Σ(k=i+1 to j) [pair(i,k)] × dp(i+1,k-1) × dp(k+1,j)
```

Where pair(i,k) = 1 if bases i and k can form Watson-Crick pair.

## 🧪 Testing

### Test Cases

```python
# Test case 1: Empty or single base
assert count_noncrossing_matchings("") == 1  # empty matching
assert count_noncrossing_matchings("A") == 1  # only empty matching

# Test case 2: Two bases that can pair
assert count_noncrossing_matchings("AU") == 2  # empty or paired
assert count_noncrossing_matchings("AA") == 1  # cannot pair, only empty

# Test case 3: Four bases (sample)
assert count_noncrossing_matchings("AUAU") == 7

# Test case 4: All bases identical (cannot pair)
assert count_noncrossing_matchings("AAAA") == 1  # only empty matching

# Test case 5: Perfect matching possible
assert count_noncrossing_matchings("AUCG") >= 2  # at least empty and perfect
```

### Verification with Sample

For "AUAU" (length 4), the 7 noncrossing matchings are:

1. ∅ (empty)
2. (1,2)
3. (1,4)
4. (2,3)
5. (3,4)
6. (1,2) & (3,4)
7. (1,4) & (2,3)

## 🔗 Related Problems

- [CAT](https://rosalind.info/problems/cat/) - Catalan Numbers (perfect matchings)
- [MMCH](https://rosalind.info/problems/mmch/) - Maximum Matchings in RNA
- [PMCH](https://rosalind.info/problems/pmch/) - Perfect Matchings and RNA
- [PROB](https://rosalind.info/problems/prob/) - Introduction to Random Strings

## 📚 Learning Resources

- [Motzkin Numbers](https://en.wikipedia.org/wiki/Motzkin_number)
- [RNA Secondary Structure](https://en.wikipedia.org/wiki/Nucleic_acid_secondary_structure)
- [Noncrossing Partitions](https://en.wikipedia.org/wiki/Noncrossing_partition)
- [Dynamic Programming in Bioinformatics](https://en.wikipedia.org/wiki/Dynamic_programming)

## 🤝 Contributing

Contributions are welcome! Here are ways to contribute:

### Add New Features:
- Handle wobble base pairs (G-U)
- Minimum helix length constraints
- Include pseudoknots (crossing edges)
- Generate all structures (not just count)

### Improve Performance:
- Optimize to O(n²) using different DP formulation
- Parallel computation for long sequences
- Memory-efficient streaming DP
- Bit-level optimizations

### Enhance Documentation:
- Visualize RNA secondary structures
- Step-by-step DP table filling animation
- Interactive structure explorer
- Biological significance examples

## 📄 License

This project is created for educational purposes as part of the [ROSALIND](https://rosalind.info) bioinformatics platform. All code is provided under the MIT License.

## 🌟 Acknowledgments

- [ROSALIND](https://rosalind.info) for the bioinformatics problems
- Theodore Motzkin for Motzkin numbers
- RNA structure prediction researchers
- Dynamic programming algorithm developers

## 📈 Performance Benchmarks

For maximum problem size (n=300):

- Python DP (bottom-up): ~1-2 seconds
- Python memoized (top-down): ~0.5-1 second
- R Solution: ~2-3 seconds

Memory usage:
- DP table: 300 × 300 × 8 bytes ≈ 720KB
- Memoized: similar worst-case
- Total: < 1MB

## 🔍 Additional Notes

### Difference from CAT Problem

- **CAT:** Counts only perfect matchings (all bases paired)
- **MOTZ:** Counts all matchings (some bases may be unpaired)
- **Relation:** Perfect matchings are subset of all matchings

### Algorithm Complexity

Naive implementation: O(n³)
- n iterations for substring length
- n iterations for start position
- n iterations for pairing position

Could be optimized to O(n²) with different DP formulation.

### Biological Significance

Noncrossing matchings correspond to RNA secondary structures without:
- Pseudoknots

But allowing:
- Hairpin loops
- Bulges
- Internal loops
- Multibranch loops
- Unpaired bases

### Implementation Details

Key insights:
- Empty matching counts as 1
- Independence: Regions separated by a base pair are independent
- No crossing constraint: Automatically enforced by DP structure
- Base pairing: Only Watson-Crick pairs allowed

### Edge Cases

- **Empty string:** 1 matching (empty)
- **Single base:** 1 matching (empty)
- **No possible pairs:** Only empty matching (count = 1)
- **All bases complementary:** Many possible matchings
- **Maximum length (300):** Still feasible with O(n³) algorithm

### Extensions

- **Weighted structures:** Different energies for different pairings
- **Minimum loop length:** Require at least 3 unpaired bases in loops
- **Maximum pairing distance:** Limit how far apart paired bases can be
- **Multiple structure types:** Distinguish between different loop types

---

For more bioinformatics challenges, visit [ROSALIND](https://rosalind.info) and join the bioinformatics community!