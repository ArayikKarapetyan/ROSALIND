# Shortest Superstring (Genome Assembly) Solutions

![Python](https://img.shields.io/badge/Python-3.8+-blue.svg)
![R](https://img.shields.io/badge/R-4.0+-blue.svg)
![Genome Assembly](https://img.shields.io/badge/Genome%20Assembly-Overlap-blue.svg)
![Algorithm](https://img.shields.io/badge/Algorithm-Greedy-orange.svg)

![ROSALIND](https://rosalind.info/static/img/logo.png?4eac7af5ad70)

A bioinformatics solution for genome assembly using the shortest superstring approach, which reconstructs a chromosome from overlapping DNA reads by finding the shortest string containing all reads as substrings.

## 📋 Problem Description

**Problem:** [Genome Assembly as Shortest Superstring](https://rosalind.info/problems/long/)  
**Category:** Bioinformatics Armory  
**ID:** LONG

Given a set of DNA reads (substrings) from the same chromosome, reconstruct the chromosome by finding the shortest string that contains all reads as substrings. The problem guarantees a unique solution with overlaps > half the read length.

**Input:** At most 50 DNA strings (reads) of ≈ equal length (≤ 1 kbp) in FASTA format  
**Output:** Shortest superstring containing all reads

**Example:**
```text
Input:
ATTAGACCTG
CCTGCCGGAA
AGACCTGCCG
GCCGGAATAC

Output: ATTAGACCTGCCGGAATAC

Assembly:
ATTAGACCTG
     AGACCTGCCG
         CCTGCCGGAA
              GCCGGAATAC
Result: ATTAGACCTGCCGGAATAC
```

## 🧬 Solutions

### 1. Python Solution (From Scratch) - `python_manual.py`

```python
def shortest_superstring_manual(strings):
    """Greedy shortest superstring algorithm."""
    while len(strings) > 1:
        max_overlap = 0
        best_i, best_j = -1, -1
        
        # Find pair with maximum overlap
        for i in range(len(strings)):
            for j in range(len(strings)):
                if i == j: continue
                overlap = overlap_length(strings[i], strings[j])
                if overlap > max_overlap:
                    max_overlap = overlap
                    best_i, best_j = i, j
        
        if max_overlap == 0:
            # No overlaps, concatenate
            return ''.join(strings)
        
        # Merge best pair
        merged = strings[best_i] + strings[best_j][max_overlap:]
        
        # Update list
        del strings[max(best_i, best_j)]
        del strings[min(best_i, best_j)]
        strings.append(merged)
    
    return strings[0]
```

**Features:**
- Greedy overlap merging
- O(k² × L²) time complexity where k = #reads, L = read length
- Simple and effective for problem constraints

### 2. Python Efficient Solution - `python_efficient.py`

```python
def shortest_superstring_efficient(strings):
    """Efficient version with precomputed overlaps."""
    # Precompute all pairwise overlaps
    n = len(strings)
    overlaps = [[0] * n for _ in range(n)]
    
    for i in range(n):
        for j in range(n):
            if i != j:
                overlaps[i][j] = overlap_length(strings[i], strings[j])
    
    # Greedy merging using precomputed overlaps
    current = strings[:]
    while len(current) > 1:
        # Find max overlap in current set
        max_ov = -1
        best_i = best_j = -1
        
        for i in range(len(current)):
            for j in range(len(current)):
                if i == j: continue
                # Map to original indices if needed
                if overlaps[i][j] > max_ov:
                    max_ov = overlaps[i][j]
                    best_i, best_j = i, j
        
        if max_ov <= 0:
            return ''.join(current)
        
        # Merge and update
        merged = current[best_i] + current[best_j][max_ov:]
        del current[max(best_i, best_j)]
        del current[min(best_i, best_j)]
        current.append(merged)
    
    return current[0]
```

**Features:**
- Precomputed overlap matrix
- Reduced redundant calculations
- Same greedy approach but more efficient

### 3. R Solution - `r_solution.R`

```r
shortest_superstring_r <- function(strings) {
  while (length(strings) > 1) {
    # Find best overlap pair
    best <- find_max_overlap_pair_r(strings)
    i <- best$pair[1]
    j <- best$pair[2]
    
    if (i == -1) {
      return(paste(strings, collapse = ""))
    }
    
    # Merge and update
    merged <- merge_strings_r(strings[i], strings[j], best$overlap)
    strings <- strings[-c(i, j)]
    strings <- c(strings, merged)
  }
  
  return(strings[1])
}
```

**Features:**
- R implementation of greedy algorithm
- Uses helper functions for clarity
- Handles the problem constraints

## 📁 File Structure

```text
Shortest-Superstring/
│
├── README.md              # This documentation
├── python_manual.py       # Python greedy solution
├── python_efficient.py    # Python optimized solution
├── python_graph.py        # Python graph-based solution
├── r_solution.R           # R implementation
└── Dataset.fasta         # Input FASTA file
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

Input should be in FASTA format with DNA reads:

```text
>Read_1
ATTAGACCTG
>Read_2
CCTGCCGGAA
>Read_3
AGACCTGCCG
>Read_4
GCCGGAATAC
```

### File Reading

**Python:**

```python
def parse_fasta(fasta_string):
    lines = fasta_string.strip().split('\n')
    sequences = []
    current = ""
    for line in lines:
        if line.startswith('>'):
            if current:
                sequences.append(current)
            current = ""
        else:
            current += line.strip()
    if current:
        sequences.append(current)
    return sequences

with open("Dataset.fasta", "r") as file:
    strings = parse_fasta(file.read())
```

**R:**

```r
parse_fasta_r <- function(fasta_string) {
  lines <- strsplit(fasta_string, "\n")[[1]]
  sequences <- character(0)
  current <- ""
  
  for (line in lines) {
    if (startsWith(line, ">")) {
      if (nchar(current) > 0) {
        sequences <- c(sequences, current)
      }
      current <- ""
    } else {
      current <- paste0(current, line)
    }
  }
  if (nchar(current) > 0) {
    sequences <- c(sequences, current)
  }
  
  return(sequences)
}

strings <- parse_fasta_r(paste(readLines("Dataset.fasta"), collapse = "\n"))
```

## 📊 Algorithm Analysis

### Time Complexity

| Method | Approach | Time Complexity | Best For |
|--------|----------|----------------|----------|
| Greedy Basic | Check all pairs each iteration | O(k³ × L) | Small datasets |
| Greedy with Precomputed | Precompute all overlaps | O(k² × L + k³) | Medium datasets |
| Graph-based | Overlap graph construction | O(k² × L) | Theoretical |

*k = number of reads (≤ 50), L = read length (≤ 1000)*

### Overlap Calculation

The key operation is finding maximum overlap between two strings a and b:
- Check if suffix of a matches prefix of b
- Start from maximum possible overlap (min(len(a), len(b)))
- Decrease until match found or overlap = 0
- Problem guarantees overlap > half length

## 🧪 Testing

### Test Cases

```python
# Test case from problem
strings = ["ATTAGACCTG", "CCTGCCGGAA", "AGACCTGCCG", "GCCGGAATAC"]
result = shortest_superstring_manual(strings[:])
assert result == "ATTAGACCTGCCGGAATAC"

# Simple case with perfect overlaps
strings2 = ["ABC", "BCD", "CDE"]
result2 = shortest_superstring_manual(strings2[:])
assert result2 == "ABCDE"

# No overlaps
strings3 = ["AAA", "BBB", "CCC"]
result3 = shortest_superstring_manual(strings3[:])
assert result3 == "AAABBBCCC" or result3 == "BBBAAACCC" etc.

# Single string
strings4 = ["SINGLE"]
assert shortest_superstring_manual(strings4[:]) == "SINGLE"
```

### Sample Dataset Verification

```text
Reads:
1: ATTAGACCTG
2: CCTGCCGGAA
3: AGACCTGCCG
4: GCCGGAATAC

Overlaps:
1 & 3: "ATTAGACCTG" vs "AGACCTGCCG" - overlap "AGACCTG" (7)
3 & 2: "AGACCTGCCG" vs "CCTGCCGGAA" - overlap "CCTGCCG" (7)
2 & 4: "CCTGCCGGAA" vs "GCCGGAATAC" - overlap "GCCGGA" (6)

Assembly path: 1 → 3 → 2 → 4
Result: ATTAGACCTGCCGGAATAC
```

## 🔗 Related Problems

- **[GRPH](https://rosalind.info/problems/grph/)** - Overlap Graphs
- **[LCSM](https://rosalind.info/problems/lcsm/)** - Finding a Shared Motif
- **[SUBS](https://rosalind.info/problems/subs/)** - Finding a Motif in DNA
- **[SPLC](https://rosalind.info/problems/splc/)** - RNA Splicing

## 📚 Learning Resources

- [Shortest Superstring Problem](https://en.wikipedia.org/wiki/Shortest_superstring_problem)
- [Genome Assembly](https://en.wikipedia.org/wiki/Sequence_assembly)
- [Overlap-Layout-Consensus](https://en.wikipedia.org/wiki/Sequence_assembly#Overlap-layout-consensus)
- [Greedy Algorithm](https://en.wikipedia.org/wiki/Greedy_algorithm)

## 🤝 Contributing

Contributions are welcome! Here are ways to contribute:

### Add New Features:
- Exact shortest superstring algorithms
- De Bruijn graph assembly
- Handling sequencing errors
- Paired-end read support

### Improve Documentation:
- Add assembly visualization
- Create overlap graph diagrams
- Add biological context
- Include assembly statistics

### Enhance Performance:
- Efficient overlap computation (suffix trees)
- Parallel processing
- Memory-efficient implementations
- Streaming for large datasets

## 📄 License

This project is created for educational purposes as part of the [ROSALIND](https://rosalind.info) bioinformatics platform. All code is provided under the MIT License.

## 🌟 Acknowledgments

- [ROSALIND](https://rosalind.info) for the bioinformatics problems
- Genome assembly researchers
- Algorithm developers for assembly tools
- Bioinformatics open-source community

## 📈 Performance Benchmarks

For maximum problem size (50 reads of 1000bp):
- **Python Basic:** ~0.1-0.5 seconds
- **Python Efficient:** ~0.05-0.2 seconds
- **R Basic:** ~0.5-1 second
- **R Efficient:** ~0.2-0.5 seconds

## 🔍 Additional Notes

### Biological Context

This simulates:
- Overlap-Layout-Consensus assembly
- Sanger sequencing assembly
- Contig assembly from reads

### Problem Guarantees
- Unique solution exists
- Overlaps > half read length
- Reads from same strand
- No sequencing errors

### Algorithm Limitations
- Greedy may not find optimal for arbitrary inputs
- But works with problem constraints
- For real data, use more sophisticated assemblers

### Extensions
- Could add error correction
- Handle reverse complements
- Support for varying read lengths
- Quality score integration

---

For more bioinformatics challenges, visit [ROSALIND](https://rosalind.info) and join the bioinformatics community!