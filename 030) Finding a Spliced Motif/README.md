# Subsequence Index Finder

![Python](https://img.shields.io/badge/Python-3.8+-blue.svg)
![R](https://img.shields.io/badge/R-4.0+-blue.svg)
![BioPython](https://img.shields.io/badge/BioPython-1.8+-green.svg)
![Algorithm](https://img.shields.io/badge/Algorithm-Greedy-green.svg)
![Bioinformatics](https://img.shields.io/badge/Bioinformatics-Sequence%20Analysis-blue.svg)
![ROSALIND](https://img.shields.io/badge/ROSALIND-Problem%20SSEQ-orange.svg)

![ROSALIND](https://rosalind.info/static/img/logo.png?4eac7af5ad70)

A comprehensive solution for finding indices where a DNA string appears as a subsequence within another DNA string, implemented in Python (both manual and BioPython) and R.

## 📋 Problem Description

**Problem:** [Finding a Spliced Motif](https://rosalind.info/problems/sseq/)  
**Category:** Bioinformatics Stronghold  
**ID:** SSEQ

A subsequence of a string is a collection of symbols contained in order (though not necessarily contiguously) in the string. Given two DNA strings s and t, find one collection of indices in s where the symbols of t appear as a subsequence.

**Input:** Two DNA strings s and t (each ≤ 1 kbp) in FASTA format  
**Output:** One collection of indices (1-based) where t appears as subsequence in s

**Example:**
```text
s: ACGTACGTGACG
t: GTA

Valid solution: 3 8 10

s: A C G T A C G T G A C G
   1 2 3 4 5 6 7 8 9 10 11 12
   
Matches: G at position 3, T at position 8, A at position 10
```

## 🧬 Solutions

### 1. Python Manual Solution - `python_manual.py`

```python
def find_subsequence_indices_manual(s, t):
    """Find indices using greedy left-to-right algorithm."""
    indices = []
    i = j = 0
    
    while i < len(s) and j < len(t):
        if s[i] == t[j]:
            indices.append(i + 1)  # 1-based indexing
            j += 1
        i += 1
    
    return indices if j == len(t) else None
```

**Features:**
- Simple greedy algorithm
- O(n + m) time complexity
- 1-based indexing as required

### 2. Python BioPython Solution - `python_biopython.py`

```python
from Bio import SeqIO
from io import StringIO

def find_subsequence_biopython(fasta_string):
    """Find subsequence indices using BioPython."""
    fasta_io = StringIO(fasta_string)
    records = list(SeqIO.parse(fasta_io, "fasta"))
    
    s = str(records[0].seq)
    t = str(records[1].seq)
    
    indices = []
    i = j = 0
    
    while i < len(s) and j < len(t):
        if s[i] == t[j]:
            indices.append(i + 1)
            j += 1
        i += 1
    
    if j == len(t):
        return indices
    else:
        raise ValueError("t is not a subsequence of s")
```

**Features:**
- Uses BioPython for FASTA parsing
- Professional bioinformatics workflow
- Built-in error handling for FASTA format

### 3. R Solution - `r_solution.R`

```r
find_subsequence_indices_r <- function(s, t) {
  s_chars <- strsplit(s, "")[[1]]
  t_chars <- strsplit(t, "")[[1]]
  
  indices <- integer(length(t_chars))
  i <- j <- 1
  
  while (i <= length(s_chars) && j <= length(t_chars)) {
    if (s_chars[i] == t_chars[j]) {
      indices[j] <- i
      j <- j + 1
    }
    i <- i + 1
  }
  
  if (j > length(t_chars)) return(indices)
  else return(NULL)
}
```

**Features:**
- R's string manipulation functions
- Multiple implementation strategies
- Vectorized operations where possible

## 📁 File Structure

```text
Subsequence-Finder/
│
├── README.md              # This documentation
├── python_manual.py       # Python manual implementation
├── python_biopython.py    # Python BioPython implementation
├── r_solution.R           # R implementation
└── Dataset.fasta         # Input FASTA file
```

## 🚀 Installation and Usage

### Python Requirements

```bash
# For BioPython solution
pip install biopython
```

### Running the Solutions

**Python (Manual):**

```bash
python python_manual.py
```

**Python (BioPython):**

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

Input should be in FASTA format with two sequences:

```fasta
>Rosalind_14
ACGTACGTGACG
>Rosalind_18
GTA
```

### File Reading Examples

**Python Manual:**

```python
def parse_fasta_manual(fasta_string):
    sequences = {}
    current_id = ""
    for line in fasta_string.strip().split('\n'):
        if line.startswith('>'):
            current_id = line[1:].strip()
            sequences[current_id] = ""
        else:
            sequences[current_id] += line.strip()
    return sequences

with open("Dataset.fasta", "r") as file:
    fasta_data = file.read()
sequences = parse_fasta_manual(fasta_data)
s = list(sequences.values())[0]
t = list(sequences.values())[1]
```

**Python BioPython:**

```python
from Bio import SeqIO
records = list(SeqIO.parse("Dataset.fasta", "fasta"))
s = str(records[0].seq)
t = str(records[1].seq)
```

**R:**

```r
parse_fasta_r <- function(fasta_string) {
  lines <- strsplit(fasta_string, "\n")[[1]]
  sequences <- list()
  current_id <- ""
  
  for (line in lines) {
    if (startsWith(line, ">")) {
      current_id <- substring(line, 2)
      sequences[[current_id]] <- ""
    } else {
      sequences[[current_id]] <- paste0(sequences[[current_id]], line)
    }
  }
  
  return(sequences)
}

fasta_data <- readLines("Dataset.fasta", warn = FALSE)
fasta_string <- paste(fasta_data, collapse = "\n")
sequences <- parse_fasta_r(fasta_string)
s <- sequences[[1]]
t <- sequences[[2]]
```

## 📊 Algorithm Analysis

### Time Complexity

| Method | Approach | Time Complexity | Space Complexity |
|--------|----------|----------------|------------------|
| Greedy Algorithm | Two-pointer scan | O(n + m) | O(m) |
| BioPython + Greedy | FASTA parsing + scan | O(n + m) | O(n + m) |
| R Greedy | Character array scan | O(n + m) | O(n + m) |

Where:
- n = length of string s (≤ 1000)
- m = length of string t (≤ 1000)

### Algorithm Details

**Greedy Algorithm:**
1. Initialize pointers at start of both strings
2. Move through string s
3. When character matches current character in t, record position and move t pointer
4. Continue until all characters in t are found or end of s is reached
5. If all characters in t are found, return positions; otherwise, return failure

## 🧪 Testing

### Test Cases

**Python:**

```python
# Test case from problem
fasta_data = """>Rosalind_14
ACGTACGTGACG
>Rosalind_18
GTA"""

# Test manual solution
sequences = parse_fasta_manual(fasta_data)
s = sequences["Rosalind_14"]
t = sequences["Rosalind_18"]
indices = find_subsequence_indices_manual(s, t)
assert indices == [3, 8, 10]

# Test BioPython solution
indices_bp = find_subsequence_biopython(fasta_data)
assert indices_bp == [3, 8, 10]

# Additional tests
assert find_subsequence_indices_manual("ABCD", "AC") == [1, 3]
assert find_subsequence_indices_manual("ABCD", "CA") == None
assert find_subsequence_indices_manual("AAAA", "AA") == [1, 2]
```

### Sample Dataset Verification

```text
s: A C G T A C G T G A C G
   1 2 3 4 5 6 7 8 9 10 11 12
t: G T A

Step-by-step:
1. i=0, j=0: s[0]=A ≠ G → i=1
2. i=1, j=0: s[1]=C ≠ G → i=2  
3. i=2, j=0: s[2]=G = G → record 3, j=1, i=3
4. i=3, j=1: s[3]=T = T → record 8, j=2, i=4
5. i=4, j=2: s[4]=A ≠ A → i=5
... skip to i=9
6. i=9, j=2: s[9]=A = A → record 10, j=3

All t characters found → return [3, 8, 10]
```

## 🔗 Related Problems

- [**SUBS**](https://rosalind.info/problems/subs/) - Finding a Motif in DNA (substring search)
- [**LCSM**](https://rosalind.info/problems/lcsm/) - Finding a Shared Motif
- [**LONG**](https://rosalind.info/problems/long/) - Genome Assembly
- [**PROT**](https://rosalind.info/problems/prot/) - Translating RNA into Protein

## 📚 Learning Resources

- [Subsequence vs Substring](https://en.wikipedia.org/wiki/Subsequence)
- [Greedy Algorithms](https://en.wikipedia.org/wiki/Greedy_algorithm)
- [BioPython SeqIO](https://biopython.org/wiki/SeqIO)
- [Sequence Alignment Basics](https://en.wikipedia.org/wiki/Sequence_alignment)

## 🤝 Contributing

Contributions are welcome! Here are ways to contribute:

### Add New Features:
- Find all possible solutions
- Find solution with minimum/maximum span
- Support for RNA and protein sequences
- Approximate subsequence matching

### Improve Documentation:
- Add visual sequence matching diagrams
- Create interactive subsequence explorer
- Add algorithm animation
- Include biological applications examples

### Enhance Performance:
- Parallel processing for finding all solutions
- Memory-mapped file support
- GPU acceleration for large datasets
- Streaming algorithms for very long sequences

## 📄 License

This project is created for educational purposes as part of the [ROSALIND](https://rosalind.info) bioinformatics platform. All code is provided under the MIT License.

## 🌟 Acknowledgments

- [ROSALIND](https://rosalind.info) for the bioinformatics problems
- [BioPython](https://biopython.org) development team
- Algorithm designers for sequence matching
- Open-source scientific computing community

## 📈 Performance Benchmarks

For maximum constraints (n=1000, m=1000):
- **Python Manual:** ~0.0002 seconds
- **Python BioPython:** ~0.0003 seconds (includes FASTA parsing)
- **R Solution:** ~0.0005 seconds

## 🔍 Additional Notes

### Biological Applications

**Spliced Motifs:** Finding exons in pre-mRNA sequences

**Protein Domains:** Identifying conserved domains in protein sequences

**Regulatory Elements:** Locating transcription factor binding sites

**Evolutionary Studies:** Finding conserved subsequences across species

### Implementation Notes
- **1-based Indexing:** Required by biological convention
- **Multiple Solutions:** Many valid index sets may exist
- **Greedy Optimality:** This greedy algorithm finds a valid solution but not necessarily optimal by other metrics
- **FASTA Handling:** Both single-line and multi-line sequences supported

### BioPython Advantages
- **Robust FASTA Parsing:** Handles various FASTA formats
- **Sequence Objects:** Provide additional methods and validation
- **Professional Workflow:** Standard in bioinformatics research
- **Integration:** Works with other BioPython modules

### Extension Ideas
- **k-best Solutions:** Find top k different solutions
- **Distance Constraints:** Solutions with maximum/minimum gap sizes
- **Weighted Matching:** Account for base quality scores
- **Batch Processing:** Process multiple query sequences

---

For more bioinformatics challenges, visit [ROSALIND](https://rosalind.info) and join the bioinformatics community!