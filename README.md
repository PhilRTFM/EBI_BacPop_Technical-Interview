# EBI_BacPop_Technical-Interview

![Python](https://img.shields.io/badge/python-v3.8+-blue.svg)
![Biopython](https://img.shields.io/badge/BioPython-v1.81+-blue.svg)
![Matplotlib](https://img.shields.io/badge/Matplotlib-v3.5+-blue.svg)

This repository contains my submission for the **Streptococcus pneumoniae Distance Coding Challenge** during my application for the Pathogen Informatics and Modelling group at EMBL-EBI , which implements the **MinHash Algorithm** to compute genetic distances between isolates and visualize the relationships using a **Neighbor-Joining Tree**.

---

## Features

- Efficient processing of genome sequences using Biopython.
- Implementation of **MinHash Sketching** for scalable computation.
- Calculation of **Jaccard Distances** (exact and approximate).
- Construction and visualization of a **Neighbor-Joining Tree**.
- Modular and well-documented Python code.

---

## Prerequisites

Before you begin, ensure you have the following:

- Python 3.8 or later installed on your system.
- Required Python packages installed. You can install them via:
  ```bash
  pip install -r requirements.txt
  ```

### Genome Files
Download the genomes archive `s_pneumoniae_genomes.tar.gz` provided with the challenge and extract it into the repository folder.

---

## Installation

Clone the repository and navigate to the project directory:

```bash
git clone https://github.com/yourusername/s_pneumoniae_distance_challenge.git
cd s_pneumoniae_distance_challenge
```

---

## Usage

### Running the Code
1. Ensure the genome files (`R6.fa`, `TIGR4.fa`, `14412_3#82.contigs_velvet.fa`, `14412_3#84.contigs_velvet.fa`) are present in the working directory.
2. Run the main script:
   ```bash
   python main.py
   ```

### Outputs
1. **Pairwise Distance Matrix**: Saved as `pairwise_distance_matrix.txt`.
2. **MinHash Sketches**: Each genome's sketch saved as `<genome>_sketch.txt`.
3. **Neighbor-Joining Tree**:
   - Newick format: `phylogenetic_tree.nwk`
   - PNG visualization: `phylogenetic_tree.png`

---

## Workflow

1. **FASTA Parsing**  
   Load genomes using Biopython and process sequences into k-mers of length 14.

2. **14-mers and Jaccard Distance**  
   Compute the exact Jaccard distances using k-mer sets.

3. **MinHash Sketching**  
   - Hash each 14-mer (forward and reverse complements).
   - Create sketches by selecting the smallest 1000 hashes.

4. **Approximate Jaccard Distance**  
   Use MinHash sketches to compute approximate distances.

5. **Neighbor-Joining Tree**  
   Construct and visualize a phylogenetic tree using the distance matrix.

---

## Example Results

### Pairwise Distance Matrix
```plaintext
        R6.fa   TIGR4.fa        14412_3#82.contigs_velvet.fa    14412_3#84.contigs_velvet.fa
R6.fa   0.0000  0.9744  0.9887  0.9887
TIGR4.fa        0.9744  0.0000  0.9885  0.9885
14412_3#82.contigs_velvet.fa    0.9887  0.9885  0.0000  0.5608
14412_3#84.contigs_velvet.fa    0.9887  0.9885  0.5608  0.0000
```

### Phylogenetic Tree Visualization
The Neighbor-Joining tree highlights relationships between the isolates. Example:

![phylogenetic_tree](https://github.com/user-attachments/assets/5d48707e-6491-4103-9a9c-9acf2ffc116f)

---

## Discussion

1. **Comparison of Full and MinHash Distances**  
   The MinHash distances closely approximate the full Jaccard distances while being computationally efficient.

2. **Effect of Sketch Size**  
   Increasing the sketch size improves accuracy but increases computational overhead. A sketch size of 1000 provides a good balance.

---

## Contact
For any queries, please reach out:
- **GitHub**: [yourusername](https://github.com/PhilRTFM)
