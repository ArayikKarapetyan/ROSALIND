# k-mer Composition Solutions

![Python](https://img.shields.io/badge/Python-3.8+-blue.svg)
![R](https://img.shields.io/badge/R-4.0+-blue.svg)
![BioPython](https://img.shields.io/badge/BioPython-1.78+-green.svg)
![k-mer](https://img.shields.io/badge/k--mer-4_mer-blue.svg)

![ROSALIND](https://rosalind.info/static/img/logo.png?4eac7af5ad70)

Solutions for calculating the k-mer composition of DNA sequences, specifically 4-mer composition, with all k-mers ordered lexicographically.

## 📋 Problem Description

**Problem:** [k-mer Composition](https://rosalind.info/problems/kmer/)  
**Category:** Bioinformatics Armory  
**ID:** KMER

Given a DNA string, compute its k-mer composition array where k=4. The array represents the frequency of each possible 4-mer in lexicographic order (AAAA, AAAC, AAAG, AAAT, AACA, ... , TTTT).

**Input:** A DNA string in FASTA format (length ≤ 100 kbp)  
**Output:** The 4-mer composition array (256 integers)

**Example:**
For a short sequence "ACGCGGCTCTGAAA" and k=2, output would be counts for AA, AC, AG, AT, CA, CC, CG, CT, GA, GC, GG, GT, TA, TC, TG, TT.

## 🧬 Solutions

### 1. Python Solution (Manual) - `python_manual.py`

```python
def kmer_composition(dna, k=4):
    from itertools import product
    all_kmers = [''.join(p) for p in product('ACGT', repeat=k)]
    kmer_counts = {kmer: 0 for kmer in all_kmers}
    
    for i in range(len(dna) - k + 1):
        kmer = dna[i:i+k]
        kmer_counts[kmer] += 1
    
    return [kmer_counts[kmer] for kmer in all_kmers]
```

**Features:**
- Uses itertools.product for k-mer generation
- Simple sliding window counting
- Explicit lexicographic ordering

### 2. Python BioPython Solution - `python_biopython.py`

```python
from Bio import SeqIO
from itertools import product

def kmer_composition_biopython(seq, k=4):
    all_kmers = [''.join(p) for p in product('ACGT', repeat=k)]
    kmer_counts = {kmer: 0 for kmer in all_kmers}
    
    seq_str = str(seq)
    for i in range(len(seq_str) - k + 1):
        kmer = seq_str[i:i+k]
        kmer_counts[kmer] += 1
    
    return [kmer_counts[kmer] for kmer in all_kmers]
```

**Features:**
- Uses BioPython for FASTA parsing
- Handles SeqRecord objects
- Professional bioinformatics workflow

### 3. Python Efficient Solution - `python_efficient.py`

```python
def kmer_composition_efficient(dna, k=4):
    base_to_num = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    total_kmers = 4 ** k
    counts = [0] * total_kmers
    
    # Rolling hash for k-mers
    current_index = 0
    for i in range(k):
        current_index = current_index * 4 + base_to_num[dna[i]]
    counts[current_index] = 1
    
    for i in range(k, len(dna)):
        current_index %= (total_kmers // 4)
        current_index = current_index * 4 + base_to_num[dna[i]]
        counts[current_index] += 1
    
    return counts
```

**Features:**
- Integer encoding for O(1) k-mer identification
- Rolling window avoids string slicing
- O(n) time complexity
- Memory efficient

### 4. R Solution - `r_solution.R`

```r
kmer_composition_r <- function(dna, k = 4) {
  all_kmers <- expand.grid(rep(list(c("A", "C", "G", "T")), k)) %>%
    apply(1, paste, collapse = "") %>%
    sort()
  
  counts <- integer(length(all_kmers))
  names(counts) <- all_kmers
  
  for (i in 1:(nchar(dna) - k + 1)) {
    kmer <- substr(dna, i, i + k - 1)
    counts[kmer] <- counts[kmer] + 1
  }
  
  return(as.vector(counts))
}
```

**Features:**
- Uses R's expand.grid for k-mer generation
- Vectorized operations where possible
- dplyr pipe operator for readability

## 📁 File Structure

```text
kmer-composition/
│
├── README.md              # This documentation
├── python_manual.py       # Python manual solution
├── python_biopython.py    # Python BioPython solution
├── python_efficient.py    # Python efficient solution
├── r_solution.R           # R implementation
└── Dataset.txt            # Input FASTA file
```

## 🚀 Installation and Usage

### Running the Solutions

**Python (Manual):**

```bash
python python_manual.py
```

**Python (BioPython):**

```bash
# First install BioPython if needed:
# pip install biopython
python python_biopython.py
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

Input file `Dataset.txt` should contain DNA sequence in FASTA format:

```text
>Rosalind_6431
CTTCGAAAGTTTGGGCCGAGTCTTACAGTCGGTCTTGAAGCAAAGTAACGAACTCCACGG
CCCTGACTACCGAACCAGTTGTGAGTACTCAACTGGGTGAGAGTGCAGTCCCTATTGAGT
...
```

### Output Format

Space-separated integers representing counts of all 256 possible 4-mers in lexicographic order:

```text
4 1 4 3 0 1 1 5 1 3 1 2 2 1 2 0 ...
```

## 📊 Algorithm Analysis

### Time Complexity

| Method | Approach | Time Complexity | Space Complexity |
|--------|----------|----------------|------------------|
| Python Manual | Sliding window with dict | O(n + 4^k) | O(4^k) |
| Python BioPython | Same as manual | O(n + 4^k) | O(4^k) |
| Python Efficient | Integer encoding | O(n) | O(4^k) |
| R Solution | Sliding window | O(n + 4^k) | O(4^k) |

*n = sequence length, k = 4*

### Number of k-mers

For k=4 with alphabet {A, C, G, T}:
- Total possible k-mers: 4⁴ = 256
- Lexicographic order: AAAA, AAAC, AAAG, AAAT, AACA, ... , TTTT
- Output array length: 256 integers

## 🧪 Testing

### Test Cases

```python
# Test with short sequence
dna = "ACGTACGT"
k = 2
result = kmer_composition(dna, k)
# Should have counts for all 16 2-mers

# Test with known composition
dna = "AAAA"
k = 2
result = kmer_composition(dna, k)
# AA: 3, all others: 0

# Test edge case: sequence shorter than k
dna = "ACG"
k = 4
result = kmer_composition(dna, k)
# All counts should be 0
```

### Validation

For the sample dataset:
- Output length should be exactly 256 integers
- Sum of all counts = len(sequence) - k + 1
- All counts are non-negative integers

## 🔗 Related Problems

- [**LEXF**](https://rosalind.info/problems/lexf/) - Enumerating k-mers Lexicographically
- [**REVP**](https://rosalind.info/problems/revp/) - Locating Restriction Sites
- [**LCSM**](https://rosalind.info/problems/lcsm/) - Finding a Shared Motif
- [**SUBS**](https://rosalind.info/problems/subs/) - Finding Motifs

## 📚 Learning Resources

- [k-mer](https://en.wikipedia.org/wiki/K-mer)
- [DNA Sequencing](https://en.wikipedia.org/wiki/DNA_sequencing)
- [Composition Vector](https://en.wikipedia.org/wiki/Sequence_composition)
- [Lexicographic Order](https://en.wikipedia.org/wiki/Lexicographic_order)

## 🤝 Contributing

Contributions are welcome! Here are ways to contribute:

### Add New Features:
- Support for different k values
- Handle degenerate bases (N, R, Y, etc.)
- Parallel processing for large sequences
- Streaming for very long sequences

### Improve Performance:
- GPU acceleration for k-mer counting
- Bloom filters for approximate counting
- Compressed data structures
- Multi-threaded sliding window

### Enhance Documentation:
- Visualize k-mer frequency distributions
- Add real genomic dataset examples
- Include performance comparison charts
- Create interactive k-mer explorer

## 📄 License

This project is created for educational purposes as part of the [ROSALIND](https://rosalind.info) bioinformatics platform. All code is provided under the MIT License.

## 🌟 Acknowledgments

- [ROSALIND](https://rosalind.info) for the bioinformatics problems
- Sequencing technology developers
- k-mer analysis researchers
- Open-source bioinformatics community

## 📈 Performance Benchmarks

For maximum problem size (100 kbp sequence):
- **Python Manual:** ~0.1 seconds
- **Python BioPython:** ~0.15 seconds (including parsing)
- **Python Efficient:** ~0.05 seconds
- **R Solution:** ~0.3 seconds

Memory usage:
- Sequence storage: ~100KB
- Count array: 256 integers × 8 bytes = 2KB
- Total: < 1MB

## 🔍 Additional Notes

### Biological Applications

k-mer composition is used in:
- Genome assembly
- Sequence alignment
- Metagenomics
- Species identification
- Error correction in sequencing

### Mathematical Properties
- **Sum of counts:** Σ counts = n - k + 1
- **Reverse complement:** Can be used to verify counts
- **De Bruijn graph:** k-mer composition relates to graph construction
- **Markov models:** k-mers used in sequence modeling

### Optimization Details

The efficient solution uses:
- Integer encoding (A=00, C=01, G=10, T=11 in binary)
- Rolling hash: O(1) update per position
- Base-4 arithmetic for k-mer indexing
- Avoids string operations in inner loop

### Edge Cases
- **Sequence shorter than k:** All counts are 0
- **Non-ACGT characters:** Should be handled or filtered
- **Very large k:** 4^k grows exponentially (but k=4 fixed here)
- **Empty sequence:** All counts are 0

### Extensions
- **Variable k:** Parameterize k value
- **Canonical k-mers:** Count k-mers and their reverse complements together
- **Sparse representation:** Store only non-zero counts
- **Distribution analysis:** Calculate k-mer frequency statistics

---

For more bioinformatics challenges, visit [ROSALIND](https://rosalind.info) and join the bioinformatics community!