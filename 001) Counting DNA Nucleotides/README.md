# DNA Nucleotide Counter

![Python](https://img.shields.io/badge/Python-3.8+-blue.svg)
![R](https://img.shields.io/badge/R-4.0+-blue.svg)
![ROSALIND](https://img.shields.io/badge/ROSALIND-Bioinformatics-green.svg)
<p align="center">
  <img src="https://rosalind.info/static/img/logo.png?4eac7af5ad70" alt="ROSALIND Logo">
</p>

A bioinformatics solution for counting nucleotide occurrences (A, C, G, T) in DNA sequences, implemented in both Python and R.

## 📋 Problem Description

**Problem:** [Counting DNA Nucleotides](https://rosalind.info/problems/dna/)  
**Category:** Bioinformatics Stronghold  
**ID:** DNA

Given a DNA string `s` (maximum length 1000 bp), return four integers separated by spaces counting the respective number of times that the symbols 'A', 'C', 'G', and 'T' occur in `s`.

### Example:

```text
Input: AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC
Output: 20 12 17 21
```

## 🧬 Solutions

### Python Solution (`simple_solution.py`)

```python
def count_nucleotides_dna(dna_string):
    """Count nucleotides using dictionary approach"""
    counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    for nucleotide in dna_string:
        if nucleotide in counts:
            counts[nucleotide] += 1
    return counts['A'], counts['C'], counts['G'], counts['T']
```

**Features:**
- Simple dictionary-based counting
- Alternative concise version using `str.count()`
- Robust handling of input files

### R Solution (`Solution_on_R.R`)

```r
count_nucleotides_r <- function(dna_string) {
    nucleotides <- strsplit(dna_string, "")[[1]]
    a_count <- sum(nucleotides == "A")
    c_count <- sum(nucleotides == "C")
    g_count <- sum(nucleotides == "G")
    t_count <- sum(nucleotides == "T")
    return(c(A = a_count, C = c_count, G = g_count, T = t_count))
}
```

**Features:**
- Three different implementations:
  - Basic vector comparison with `sum()`
  - Using `table()` for frequency counting
  - Factor-based approach with predefined levels
- Comprehensive error checking

## 📁 File Structure

```text
Counting DNA Nucleotides/
│
├── README.md              # This file
├── Python_solution.py     # Python implementation
├── R_Solution.R           # R implementation
└── Dataset.txt            # Example DNA sequence
```

## 🚀 Usage

### Python

```bash
# Make sure you're in the correct directory
cd /path/to/Counting\ DNA\ Nucleotides

# Run the Python script
python Python_solution.py
```

### R

```r
# From R console or RStudio
source("R_Solution.R")

# Or from command line
Rscript R_Solution.R
```

## 🔧 Path Configuration

### Current Structure (Recommended)

The scripts assume all files are in the same directory. Place your `Dataset.txt` in the same folder as the scripts.

### Custom Path Configuration

**For Python:**  
Edit line 24 in `Python_solution.py`:

```python
# Change this line to your dataset path
path = "your/custom/path/Dataset.txt"
```

**For R:**  
Edit line 40 in `R_Solution.R`:

```r
# Change this line to your dataset path
input_file <- "your/custom/path/Dataset.txt"
```

## 📊 Algorithm Explanation

Both implementations follow the same logical steps:

1. **Input Reading:** Read the DNA sequence from file, handling multi-line sequences
2. **Sequence Cleaning:** Remove whitespace and newline characters
3. **Counting:** Iterate through each character and count occurrences
4. **Output:** Format results as space-separated integers

**Time Complexity:** O(n) where n is the length of the DNA string  
**Space Complexity:** O(1) (constant space for counting)

## 📈 Performance Notes

- **Python:** The dictionary approach is O(n) and memory-efficient
- **R:** The `table()` method is fastest for large sequences, but all methods are O(n)
- Both implementations handle sequences up to 1000 bp efficiently

## 🤝 Contributing

Feel free to:
- Add implementations in other languages
- Optimize existing solutions
- Add unit tests
- Improve documentation

## 📚 Related Problems

- **[RNA](https://rosalind.info/problems/rna/)** - Transcribing DNA into RNA
- **[REVC](https://rosalind.info/problems/revc/)** - Complementing a Strand of DNA
- **[GC](https://rosalind.info/problems/gc/)** - Computing GC Content

## 📄 License

This project is for educational purposes as part of the [ROSALIND](http://rosalind.info/) bioinformatics platform.

## 🎯 Learning Outcomes

- String manipulation in Python and R
- Basic bioinformatics algorithms
- File I/O operations
- Cross-language implementation comparison

For more bioinformatics problems, visit [ROSALIND](http://rosalind.info/)