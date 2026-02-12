# Protein Motif Finder Solutions

![Python](https://img.shields.io/badge/Python-3.8+-blue.svg)
![R](https://img.shields.io/badge/R-4.0+-blue.svg)
![Bioinformatics](https://img.shields.io/badge/Bioinformatics-Motif%20Finding-blue.svg)
![Web API](https://img.shields.io/badge/Web%20API-UniProt-green.svg)

![ROSALIND](https://rosalind.info/static/img/logo.png?4eac7af5ad70)

A bioinformatics solution for finding N-glycosylation motifs in protein sequences by fetching data from UniProt database and performing pattern matching.

## 📋 Problem Description

**Problem:** [Finding a Protein Motif](https://rosalind.info/problems/mprt/)  
**Category:** Bioinformatics Armory  
**ID:** MPRT

N-glycosylation is an important protein modification. The N-glycosylation motif is represented as N{P}[ST]{P}, where:

- **N** = Asparagine
- **{P}** = any amino acid except Proline
- **[ST]** = Serine or Threonine
- **{P}** = any amino acid except Proline

Given UniProt protein IDs, fetch their sequences and find all positions where this motif occurs.

**Input:** At most 15 UniProt Protein Database access IDs  
**Output:** For each protein with the motif, output ID followed by motif positions (1-based)

### Example

```text
Input:
A2Z669
B5ZC00
P07204_TRBM_HUMAN
P20840_SAG1_YEAST

Output:
B5ZC00
85 118 142 306 395
P07204_TRBM_HUMAN
47 115 116 382 409
P20840_SAG1_YEAST
79 109 135 248 306 348 364 402 485 501 614
```

## 🧬 Solutions

### 1. Python Solution (From Scratch) - `python_manual.py`

```python
import re
import requests

def find_n_glycosylation_motif(sequence):
    """Find N-glycosylation motif positions."""
    positions = []
    # Pattern: N followed by not P, then S or T, then not P
    pattern = re.compile(r'(?=(N[^P][ST][^P]))')
    
    for match in pattern.finditer(sequence):
        positions.append(match.start() + 1)
    
    return positions
```

**Features:**
- Web API fetching with requests
- Regex pattern matching with lookahead for overlaps
- 1-based position reporting

### 2. Python Solution with BioPython - `python_biopython.py`

```python
from Bio import SeqIO
import requests
from io import StringIO

def fetch_and_find_motif(protein_id):
    """Fetch sequence and find motifs using BioPython."""
    url = f"http://www.uniprot.org/uniprot/{protein_id}.fasta"
    response = requests.get(url)
    
    if response.status_code == 200:
        fasta_data = StringIO(response.text)
        record = SeqIO.read(fasta_data, "fasta")
        sequence = str(record.seq)
        
        return find_n_glycosylation_motif(sequence)
    return []
```

**Features:**
- BioPython for FASTA parsing
- Better sequence handling
- Professional bioinformatics approach

### 3. R Solution - `r_solution.R`

```r
find_n_glycosylation_motif_r <- function(sequence) {
  # Pattern: N[^P][ST][^P] with lookahead for overlaps
  matches <- str_locate_all(sequence, regex("(?=(N[^P][ST][^P]))"))[[1]]
  if (nrow(matches) > 0) matches[, 1] else integer(0)
}
```

**Features:**
- httr for web requests
- stringr for regex with lookahead
- Vectorized operations in R

## 📁 File Structure

```text
Protein-Motif-Finder/
│
├── README.md              # This documentation
├── python_manual.py       # Python manual solution
├── python_biopython.py    # Python BioPython solution
├── r_solution.R           # R implementation
└── Dataset.txt            # Input file (protein IDs)
```

## 🚀 Installation and Usage

### Python Requirements

```bash
pip install requests biopython
```

### R Requirements

```r
install.packages(c("httr", "stringr"))
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

Input file `Dataset.txt` should contain UniProt IDs, one per line:

```text
A2Z669
B5ZC00
P07204_TRBM_HUMAN
P20840_SAG1_YEAST
```

### File Reading

**Python:**

```python
with open("Dataset.txt", "r") as file:
    protein_ids = [line.strip() for line in file if line.strip()]
```

**R:**

```r
protein_ids <- readLines("Dataset.txt", warn = FALSE)
protein_ids <- protein_ids[protein_ids != ""]  # Remove empty lines
```

## 📊 Algorithm Analysis

### Time Complexity

| Step | Operation | Time Complexity |
|------|-----------|-----------------|
| Fetch Sequence | HTTP GET request | O(1) per protein |
| Pattern Search | String scanning | O(L) where L = sequence length |
| Total | For m proteins | O(m × L) |

*L = protein sequence length (≤ few thousand), m ≤ 15*

### Pattern Details

**N-glycosylation motif:** N{P}[ST]{P}

- **Position 1:** N (Asparagine) - required
- **Position 2:** Any amino acid except P (Proline)
- **Position 3:** S (Serine) or T (Threonine)
- **Position 4:** Any amino acid except P (Proline)

**Regex:** `N[^P][ST][^P]` with lookahead for overlapping matches

## 🧪 Testing

### Test Cases

```python
# Test motif finding
test_sequence = "MNSTPNSTNSTP"
# Positions: 1 (MNST), 5 (PNST - no, starts with P), 9 (NSTP - no, ends with P)
expected = [1]
assert find_n_glycosylation_motif(test_sequence) == expected

# Test with overlapping
test_sequence2 = "NNSTNST"
# Positions: 1 (NNST), 4 (NST - too short)
# Need 4-length window: positions 1 (NNST)
expected2 = [1]
assert find_n_glycosylation_motif(test_sequence2) == expected2
```

### Sample Dataset Verification

For protein B5ZC00:

- Motif found at positions: 85, 118, 142, 306, 395
- Each position satisfies: N, not P, S/T, not P

## 🔗 Related Problems

- [SUBS](https://rosalind.info/problems/subs/) - Finding a Motif in DNA
- [LCSM](https://rosalind.info/problems/lcsm/) - Finding a Shared Motif
- [CONS](https://rosalind.info/problems/cons/) - Consensus and Profile
- [ORF](https://rosalind.info/problems/orf/) - Open Reading Frames

## 📚 Learning Resources

- [N-glycosylation](https://en.wikipedia.org/wiki/Glycosylation#N-linked_glycosylation)
- [UniProt Database](https://www.uniprot.org/)
- [Protein Motifs](https://en.wikipedia.org/wiki/Sequence_motif)
- [Regular Expressions](https://en.wikipedia.org/wiki/Regular_expression)

## 🤝 Contributing

Contributions are welcome! Here are ways to contribute:

### Add New Features:
- Support for other protein motifs
- Batch processing with progress bar
- Caching fetched sequences
- Command-line interface with options

### Improve Documentation:
- Add visual motif explanation
- Create UniProt API usage guide
- Add protein structure visualization

### Enhance Performance:
- Parallel fetching of sequences
- Local sequence database
- Optimized pattern matching algorithms
- Memory-efficient sequence processing

## 📄 License

This project is created for educational purposes as part of the [ROSALIND](https://rosalind.info) bioinformatics platform. All code is provided under the MIT License.

## ⚠️ Important Notes

- **Internet Connection Required:** Solutions fetch data from UniProt
- **Rate Limiting:** Be respectful of UniProt servers
- **Error Handling:** Network issues or invalid IDs may occur
- **Overlapping Matches:** Motifs can overlap (use lookahead in regex)

## 🌟 Acknowledgments

- [ROSALIND](https://rosalind.info) for the bioinformatics problems
- UniProt Consortium for protein data
- Protein biochemistry researchers
- Web API developers and maintainers

## 📈 Performance Considerations

- **Network:** Most time spent fetching sequences
- **Pattern Matching:** Efficient regex engines used
- **Memory:** Sequences typically < 10k amino acids
- **Parallelism:** Can fetch multiple sequences concurrently

## 🔍 Additional Notes

### Biological Significance

N-glycosylation affects:
- Protein folding and stability
- Cellular localization
- Protein-protein interactions
- Immune recognition

### Technical Details:
- 1-based indexing for biological sequences
- Overlapping motifs allowed
- Empty output for proteins without motif
- HTTP timeout handling recommended

### Extensions:
- Could add other post-translational modification motifs
- Support for PROSITE pattern syntax
- Integration with local BLAST databases

---

For more bioinformatics challenges, visit [ROSALIND](https://rosalind.info) and join the bioinformatics community!