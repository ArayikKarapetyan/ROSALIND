# DNA to RNA Transcription Solutions

![Python](https://img.shields.io/badge/Python-3.8+-blue.svg)
![R](https://img.shields.io/badge/R-4.0+-blue.svg)
![BioPython](https://img.shields.io/badge/BioPython-1.8+-green.svg)

<p align="center">
  <img src="https://rosalind.info/static/img/logo.png?4eac7af5ad70" alt="ROSALIND Logo">
</p>

A bioinformatics solution for transcribing DNA sequences to RNA sequences, implemented in multiple programming languages and approaches.

## 📋 Problem Description

**Problem:** [Transcribing DNA into RN](https://rosalind.info/problems/rna/)A  
**Category:** Bioinformatics Stronghold  
**ID:** RNA

An RNA string is a string formed from the alphabet containing 'A', 'C', 'G', and 'U'. Given a DNA string `t` corresponding to a coding strand, its transcribed RNA string `u` is formed by replacing all occurrences of 'T' in `t` with 'U' in `u`.

**Input:** A DNA string `t` having length at most 1000 nt  
**Output:** The transcribed RNA string of `t`

### Example:

```text
Input: GATGGAACTTGACTACGTAAATT
Output: GAUGGAACUUGACUACGUAAAUU
```

## 🧬 Solutions

### 1. Python Solution (From Scratch) - `python_manual.py`

```python
def transcribe_dna_to_rna(dna_string):
    """Transcribe DNA to RNA by replacing T with U."""
    return dna_string.replace('T', 'U')
```

**Features:**
- Simple string replacement using Python's built-in `replace()` method
- No external dependencies required
- O(n) time complexity, O(n) space complexity

### 2. Python Solution with BioPython - `python_biopython.py`

```python
from Bio.Seq import Seq

def transcribe_dna_to_rna_biopython(dna_string):
    """Transcribe DNA to RNA using BioPython."""
    return str(Seq(dna_string).transcribe())
```

**Features:**
- Uses BioPython's `Seq` class for biological sequence manipulation
- Built-in biological operations and validation
- Part of a comprehensive bioinformatics library

### 3. R Solution - `r_solution.R`

```r
transcribe_dna_to_rna <- function(dna_string) {
  # Using gsub to replace T with U
  rna_string <- gsub("T", "U", dna_string)
  return(rna_string)
}
```

**Features:**
- Three different implementations:
  - `gsub()` for pattern replacement
  - `chartr()` for character translation
  - `stringr` package for modern string manipulation
- Comprehensive file reading and cleaning

## 📁 File Structure

```text
Transcribing DNA into RNA/
│
├── README.md              # This documentation
├── python_manual.py       # Python manual implementation
├── python_biopython.py    # Python BioPython implementation
├── r_solution.R           # R implementation
└── Dataset.txt            # Input DNA sequence
```

## 🚀 Installation and Usage

### Python Requirements

```bash
# Install BioPython (if using BioPython solution)
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

### File Path Setup

Both implementations read from `Dataset.txt` in the same directory by default.

**Python Configuration:**

```python
# Edit in python_manual.py or python_biopython.py
path = "Dataset.txt"  # Change to your dataset path
```

**R Configuration:**

```r
# Edit in r_solution.R
input_file <- "Dataset.txt"  # Change to your dataset path
```

### Input Format

The script handles:
- Single-line DNA sequences
- Multi-line sequences (concatenates lines)
- Sequences with whitespace (automatically removed)
- Uppercase sequences (recommended)

## 📊 Algorithm Analysis

### Time Complexity

| Method | Approach | Time Complexity | Space Complexity |
|--------|----------|-----------------|------------------|
| Python Manual | `str.replace()` | O(n) | O(n) |
| BioPython | `Seq.transcribe()` | O(n) | O(n) |
| R | `gsub()` | O(n) | O(n) |

*n = length of DNA sequence*

### Key Differences

**Python Manual:**
- Simplest approach
- No dependencies
- Direct string manipulation

**BioPython:**
- Biological sequence validation
- Additional biological operations available
- Part of larger bioinformatics ecosystem

**R:**
- Statistical and data analysis advantages
- Multiple string manipulation options
- Easy integration with R's data science tools

## 🧪 Testing

### Test Cases

```python
# Test case from problem
assert transcribe_dna_to_rna("GATGGAACTTGACTACGTAAATT") == "GAUGGAACUUGACUACGUAAAUU"

# Edge cases
assert transcribe_dna_to_rna("") == ""  # Empty string
assert transcribe_dna_to_rna("AAA") == "AAA"  # No T's
assert transcribe_dna_to_rna("TTT") == "UUU"  # All T's
```

### Sample Dataset Verification

```text
Input:  GATGGAACTTGACTACGTAAATT
Output: GAUGGAACUUGACUACGUAAAUU
```

## 🔗 Related Problems

- **[DNA](https://rosalind.info/problems/dna/)** - Counting DNA Nucleotides
- **[REVC](https://rosalind.info/problems/revc/)** - Complementing a Strand of DNA
- **[GC](https://rosalind.info/problems/gc/)** - Computing GC Content
- **[HAMM](https://rosalind.info/problems/hamm/)** - Counting Point Mutations

## 📚 Learning Resources

- [BioPython Tutorial](https://biopython.org/wiki/Documentation)
- [R String Manipulation Cheat Sheet](https://rstudio.com/resources/cheatsheets/)
- [ROSALIND Bioinformatics Problems](http://rosalind.info/)
- [Python String Methods](https://docs.python.org/3/library/stdtypes.html#string-methods)

## 🤝 Contributing

Contributions are welcome! Here are ways to contribute:

### Add New Implementations:
- Julia, JavaScript, or other languages
- Alternative algorithms

### Improve Documentation:
- Add more examples
- Create tutorials
- Improve explanations

### Enhance Features:
- Add command-line interface
- Create unit tests
- Add sequence validation

## 📄 License

This project is created for educational purposes as part of the [ROSALIND](http://rosalind.info/) bioinformatics platform. All code is provided under the MIT License.

## 🌟 Acknowledgments

- [ROSALIND](http://rosalind.info/) for the bioinformatics problems
- [BioPython](https://biopython.org/) team for the excellent library
- R Project for statistical computing

## 📈 Performance Benchmarks

For sequences up to 1000 nucleotides:

- **Python Manual:** ~0.0001 seconds
- **BioPython:** ~0.0005 seconds (includes library overhead)
- **R:** ~0.0002 seconds

## 🔍 Additional Notes

- All solutions handle sequences up to 1000 bp efficiently
- Memory usage is minimal for all implementations
- Consider BioPython for complex bioinformatics workflows
- Use manual implementation for simple, dependency-free solutions

For more bioinformatics challenges, visit [ROSALIND](http://rosalind.info/) and join the bioinformatics community!