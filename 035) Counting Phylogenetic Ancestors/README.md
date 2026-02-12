# Binary Tree Internal Nodes Solutions

![Python](https://img.shields.io/badge/Python-3.8+-blue.svg)
![R](https://img.shields.io/badge/R-4.0+-blue.svg)
![Graph Theory](https://img.shields.io/badge/Graph_Theory-Trees-brightgreen.svg)
![Combinatorics](https://img.shields.io/badge/Combinatorics-Leaves_&_Nodes-blue.svg)

![ROSALIND](https://rosalind.info/static/img/logo.png?4eac7af5ad70)

A mathematical solution for counting internal nodes in unrooted binary trees based on graph theory principles and the handshake lemma.

## 📋 Problem Description

**Problem:** [Counting Internal Nodes in Binary Trees](https://rosalind.info/problems/inod/)  
**Category:** Algorithmic Heights  
**ID:** INOD

Given a positive integer n representing the number of leaves in an unrooted binary tree, calculate the number of internal nodes. In an unrooted binary tree:
- All leaves have degree 1
- All internal nodes have degree 3
- The tree is connected and acyclic

**Input:** A positive integer n (3 ≤ n ≤ 10000)  
**Output:** Number of internal nodes

### Example

```text
Input: 4
Output: 2

Tree structure:
   ○
  / \
 ○   ○
/ \ / \
L L L L  (4 leaves)
```

## 🧬 Solutions

### 1. Python Solution (Mathematical) - `python_manual.py`

```python
def count_internal_nodes(n_leaves):
    """
    Calculate internal nodes using handshake lemma.
    
    For unrooted binary tree:
    - Leaves: degree 1
    - Internal nodes: degree 3
    - Tree property: edges = nodes - 1
    
    Solving n + 3i = 2(n + i - 1)
    gives i = n - 2
    """
    return n_leaves - 2
```

**Features:**
- Direct mathematical solution
- O(1) time complexity
- Based on graph theory principles

### 2. Python Efficient Solution - `python_efficient.py`

```python
def count_internal_nodes_formula(n):
    """
    Efficient calculation using derived formula.
    
    Derivation:
    n*1 + i*3 = 2*(n + i - 1)
    n + 3i = 2n + 2i - 2
    i = n - 2
    """
    return n - 2
```

**Features:**
- Single arithmetic operation
- Constant time and space
- Clear formula derivation

### 3. R Solution - `r_solution.R`

```r
count_internal_nodes_r <- function(n_leaves) {
  # Simple formula: internal nodes = leaves - 2
  return(n_leaves - 2)
}
```

**Features:**
- One-line implementation
- Direct formula application
- Efficient and readable

## 📁 File Structure

```text
Binary-Tree-Nodes/
│
├── README.md              # This documentation
├── python_manual.py       # Python mathematical solution
├── python_efficient.py    # Python formula solution
├── r_solution.R           # R implementation
└── Dataset.txt            # Input file (single integer)
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

Input file `Dataset.txt` contains a single integer:

```text
4
```

### Range Constraints

- Minimum n: 3 (smallest binary tree)
- Maximum n: 10000
- Always produces positive integer result

## 📊 Mathematical Derivation

### Key Properties

**Tree Properties:**
- Connected graph
- No cycles
- Edges = Nodes - 1

**Unrooted Binary Tree:**
- Leaves: degree = 1
- Internal nodes: degree = 3

### Formula Derivation

Let:
- n = number of leaves (degree 1)
- i = number of internal nodes (degree 3)

From handshake lemma:
```
Sum of degrees = 2 × edges
n × 1 + i × 3 = 2 × edges
```

From tree property:
```
edges = (n + i) - 1
```

Combine equations:
```
n + 3i = 2(n + i - 1)
n + 3i = 2n + 2i - 2
3i - 2i = 2n - n - 2
i = n - 2
```

## 🧪 Testing

### Test Cases

```python
# Test cases based on known binary tree structures

# Smallest binary tree (3 leaves)
assert count_internal_nodes(3) == 1
# Structure: single internal node connected to 3 leaves

# Sample from problem
assert count_internal_nodes(4) == 2

# Larger trees
assert count_internal_nodes(5) == 3
assert count_internal_nodes(10) == 8
assert count_internal_nodes(100) == 98
assert count_internal_nodes(10000) == 9998
```

### Validation

The formula works because:
- For n = 3: i = 1 (star tree with 3 leaves)
- For n = 4: i = 2 (balanced binary tree)
- Always: i ≥ 1 for n ≥ 3
- Relationship is linear: Δi = Δn

## 🔗 Related Problems

- [TREE](https://rosalind.info/problems/tree/) - Completing a Tree
- [NKEW](https://rosalind.info/problems/nkew/) - Newick Format with Edge Weights
- [NWCK](https://rosalind.info/problems/nwck/) - Distances in Trees
- [INOD](https://rosalind.info/problems/inod/) - Counting Phylogenetic Ancestors

## 📚 Learning Resources

- [Handshake Lemma](https://en.wikipedia.org/wiki/Handshaking_lemma)
- [Binary Trees](https://en.wikipedia.org/wiki/Binary_tree)
- [Tree (Graph Theory)](https://en.wikipedia.org/wiki/Tree_(graph_theory))
- [Phylogenetic Trees](https://en.wikipedia.org/wiki/Phylogenetic_tree)

## 🤝 Contributing

Contributions are welcome! Here are ways to contribute:

### Add New Features:
- Visualize tree structures for given n
- Generate all possible tree topologies
- Calculate rooted vs unrooted differences
- Handle weighted edges

### Extend Mathematical Analysis:
- Derive formulas for degree constraints
- Analyze asymptotic behavior
- Compare with rooted binary trees
- Explore multi-furcating trees

### Enhance Documentation:
- Add tree visualization examples
- Include interactive tree builder
- Create step-by-step derivation guide
- Add biological applications

## 📄 License

This project is created for educational purposes as part of the [ROSALIND](https://rosalind.info) bioinformatics platform. All code is provided under the MIT License.

## 🌟 Acknowledgments

- [ROSALIND](https://rosalind.info) for the bioinformatics problems
- Graph theory mathematicians
- Phylogenetics researchers
- Computer scientists studying tree structures

## 📈 Performance Analysis

### Complexity Analysis

| Metric | Python Manual | Python Efficient | R Solution |
|--------|---------------|------------------|------------|
| Time Complexity | O(1) | O(1) | O(1) |
| Space Complexity | O(1) | O(1) | O(1) |
| Operations | 1 subtraction | 1 subtraction | 1 subtraction |

### Benchmark Results

For n = 10000:
- All solutions: < 0.001 seconds
- Memory usage: minimal
- Scalability: O(1) for any n

## 🔍 Additional Notes

### Biological Context

Unrooted binary trees are used in:
- Phylogenetic reconstruction
- Evolutionary biology
- Species relationship modeling
- Sequence evolution studies

### Mathematical Insights

The formula i = n - 2 reveals:
- **Conservation law:** Number of internal nodes is strictly determined by leaves
- **Growth pattern:** Each additional leaf requires one additional internal node
- **Structural constraint:** Minimum 1 internal node for n ≥ 3

### Alternative Derivation

Using Euler's formula for planar graphs or counting edges directly yields same result.

### Edge Cases

- **n = 3:** Minimum valid tree (1 internal node)
- **n = 2:** Not valid for unrooted binary tree (would need degree 2 node)
- **n = 1:** Not a tree in graph theory sense
- **Large n:** Formula remains valid up to computational limits

### Applications Beyond Biology

- Network design
- Computer science data structures
- Hierarchical clustering
- Decision trees in machine learning

### Extension Questions

- What if internal nodes could have degree 2 or 3?
- How does this change for rooted binary trees?
- What about trees with nodes of varying degrees?
- How to count all possible tree topologies?

---

For more bioinformatics challenges, visit [ROSALIND](https://rosalind.info) and join the bioinformatics community!