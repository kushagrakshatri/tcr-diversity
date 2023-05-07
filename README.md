# TCR Diversity Analyzer

This is a Python-based TCR diversity analysis tool that computes the MST-based diversity metrics for T cell repertoires using Smith-Waterman local alignment and BLOSUM-62 substitution matrix. The code constructs an undirected graph with TCR nodes and edge weights corresponding to pairwise Smith-Waterman local alignment scores.

## Table of Contents

- [Requirements](#requirements)
- [Usage](#usage)
- [Example](#example)

## Requirements

- Python 3.x
- NumPy

## Usage

1. Clone the repository:

```
git clone https://github.com/kushagrakshatri/tcr-diversity-analyzer.git
cd tcr-diversity-analyzer
```

2. Install the required dependencies:

```
pip install numpy
```

3. Create a CSV file containing the TCR sequences and place it in the `inputs` directory. The file should have the following format:

```
TCR_sequence_1
TCR_sequence_2
...
TCR_sequence_n
```

4. Import the `tcrdiversity` module in your Python script and use the provided functions and classes to analyze TCR diversity:

```python
from tcrdiversity import sw_align, Graph, diversity_mst

# Compute Smith-Waterman local alignment
score, alignment = sw_align(seq1, seq2)

# Create a TCR repertoire graph
graph = Graph('repertoire_id')

# Compute MST-based diversity metric
diversity = diversity_mst(graph)
```

## Example

Given a CSV file named `sample.csv` in the `inputs` directory containing the following TCR sequences:

```
CASSGRPGPTGDTGELFF
SAWRPPSGASWNIQYF
CASSHSAPGEIYNEQFF
CATTGEVEAQYF
CASSLTSATRGTEAFF
CASSFISSGVQEAQYF
CASTPLRGGGWPQHF
CASEQGEGTQYF
CASSLGLAVHSYEQYF
CASSLQGLAGGSTDTQYF
```
The following script will compute the MST-based diversity metric for the sample TCR repertoire:

```python
from tcrdiversity import Graph, diversity_mst

# Create a TCR repertoire graph
graph = Graph('sample')

# Compute MST-based diversity metric
diversity = diversity_mst(graph)

print(f"Diversity metric: {diversity}")
```

