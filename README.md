# EBI_BacPop_Technical-Interview

![Python](https://img.shields.io/badge/python-v3.8+-blue.svg)
![Biopython](https://img.shields.io/badge/BioPython-v1.81+-blue.svg)
![Matplotlib](https://img.shields.io/badge/Matplotlib-v3.5+-blue.svg)

This repository contains my submission for the **Streptococcus pneumoniae Distance Coding Challenge** concerning my application for the Pathogen Informatics and Modelling group at EMBL-EBI, which implements the **MinHash Algorithm** to compute genetic distances between isolates and visualize the relationships using a **Neighbor-Joining Tree** as in (Ondov, B.D., Treangen, T.J., Melsted, P. et al. Mash: fast genome and metagenome distance estimation using MinHash. Genome Biol 17, 132 (2016). https://doi.org/10.1186/s13059-016-0997-x).

Here is a graphical workflow :

![Challenge drawio (1)](https://github.com/user-attachments/assets/dac8c0d9-6af4-4803-844b-e4e4ccdd99d0)


---

## Features

- Efficient processing of genome sequences using **Biopython**.
- Implementation of **MinHash Sketching** for scalable computation.
- Calculation of **Jaccard Distances**.
- Construction and visualization of a **Neighbor-Joining Tree**.
- Modular Python code.

---

## Prerequisites

Before you begin, ensure you have the following:

- Python 3.8 or later installed on your system.
- Required Python packages installed. You can install them via:
  ```bash
  pip install -r requirements.txt
  ```

### Genome Files
When cloning the repo, you will have 4 fasta files associated : 
`R6.fa` and `TIGR4.fa` : reference sequences.
`14412_3#82.contigs_velvet.fa` and `14412_3#84.contigs_velvet.fa` : draft sequences, assembled from shotgun sequencing runs.

---

## Installation

Clone the repository and navigate to the project directory:

```bash
git clone https://github.com/PhilRTFM/EBI_BacPop_Technical-Interview.git
cd EBI_BacPop_Technical-Interview
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
2. **Neighbor-Joining Tree**:
   - Newick format: `phylogenetic_{method}_tree.nwk`
   - PNG visualization: `phylogenetic_{method}_tree.png`
Where {method} is `full` for k-mer Jaccard or `minhash` with sketches.
---

## Workflow

1. **FASTA Parsing**  
   - **Function(s)**: 
     - `load_fasta(file_path)`: Loads genome sequences from FASTA files and concatenates them into a single string.
     - `extract_kmers(sequence, k=14)`: Extracts all unique k-mers (of length 14) from the concatenated sequence.
   - **Role**:
     - Converts genome sequences into k-mer sets for further analysis.
     - Prepares the dataset for exact and approximate distance calculations.

2. **14-mers and Jaccard Distance**  
   - **Function(s)**: 
     - `full_distance(kmers_a, kmers_b)`: Computes the exact Jaccard distance between two sets of k-mers.
     - `construct_distance_matrix(labels, distances)`: Constructs a symmetric distance matrix from pairwise distances.
   - **Role**:
     - Computes exact pairwise Jaccard distances for all genome pairs using the full k-mer sets.
     - Produces a numerical representation of genetic similarity between genomes in matrix form.

3. **MinHash Sketching**  
   - **Function(s)**: 
     - `hash_kmer(kmer)`: Generates a hash value for each k-mer, considering both forward and reverse complement.
     - `minhash_sketch(kmers, sketch_size=1000)`: Generates a sketch of the `sketch_size` smallest hash values from a set of k-mers.
   - **Role**:
     - Provides a compact representation of k-mer sets for efficient distance computation.
     - Handles large datasets by reducing memory requirements while retaining essential information.

4. **Approximate Jaccard Distance**  
   - **Function(s)**: 
     - `minhash_distance(sketch_a, sketch_b)`: Computes the approximate Jaccard distance between two MinHash sketches.
     - `construct_distance_matrix(labels, distances)`: Constructs a distance matrix for the approximate distances.
   - **Role**:
     - Approximates the genetic distances between genomes using MinHash sketches, saving computation time and memory.

5. **Neighbor-Joining Tree**  
   - **Function(s)**: 
     - `doNeighbourJoining(distance_matrix, labels)`: Constructs a phylogenetic tree using the Neighbor-Joining algorithm.
     - `save_newick(tree, output_path)`: Saves the tree in Newick format.
     - `plot_tree_from_newick(newick_path, output_path)`: Visualizes the tree and saves it as an image.
   - **Role**:
     - Constructs phylogenetic trees based on both full and MinHash distance matrices.
     - Saves the trees in both Newick format and graphical PNG format for further analysis and visualization.

---

## Example Results

### Pairwise Distance Matrices

#### MinHash Jaccard Distance Matrix
The MinHash-based distances, computed using sketches of size 1000, approximate the true Jaccard distances. This method is computationally efficient and well-suited for large datasets:
```plaintext
[[0.0000 0.2624 0.9868 0.9868]
 [0.2624 0.0000 0.9858 0.9863]
 [0.9868 0.9858 0.0000 0.2847]
 [0.9868 0.9863 0.2847 0.0000]]
```

#### Full Jaccard Distance Matrix
The exact Jaccard distances, calculated using full k-mer sets, provide the baseline for evaluating the accuracy of the MinHash approximation:
```plaintext
[[0.0000 0.9744 0.9887 0.9887]
 [0.9744 0.0000 0.9885 0.9885]
 [0.9887 0.9885 0.0000 0.5608]
 [0.9887 0.9885 0.5608 0.0000]]
```

### Phylogenetic Tree
The Neighbor-Joining tree, built from the computed distance matrices, highlights the evolutionary relationships between isolates. Below is the Newick format output:

#### Full Jaccard Tree
```
(((R6.fa:0.1316,TIGR4.fa:0.1316):0.7129,14412_3#82.contigs_velvet.fa:0.7129):0.1425,14412_3#84.contigs_velvet.fa:0.1425);
```

#### MinHash Jaccard Tree
```
(((R6.fa:0.1261,TIGR4.fa:0.1261):0.7202,14412_3#82.contigs_velvet.fa:0.7202):0.1408,14412_3#84.contigs_velvet.fa:0.1408);
```

### Tree Visualization
Both trees, visualized side by side, demonstrate a similarity between the Full Jaccard and MinHash methods:

![téléchargement](https://github.com/user-attachments/assets/82ffc663-d5fa-4e31-85ea-461e1a487305)

---

## Discussion

### Benchmarking

1. **Method: Full Jaccard Distance**
```
Time for Full Jaccard Distance Matrix: 3.7777 seconds
```

2. **Method: MinHash Jaccard Distance**
```
Time for MinHash Sketches: 31.3916 seconds
Time for MinHash Distance Matrix: 0.0030 seconds
```

1. **Comparison of Full and MinHash Distances**:
   - The Full Jaccard distances serve as the ground truth for assessing the MinHash approximation.
   - Despite slight numerical differences, the MinHash-based tree retains the same topology as the Full Jaccard tree, confirming its accuracy for phylogenetic inference.

2. **Effect of Sketch Size**:
   - Increasing the sketch size reduces the approximation error, aligning MinHash distances closer to Full Jaccard distances.
   - However, larger sketches increase memory usage and runtime. A size of 1000 offers a balanced trade-off for typical datasets.
![téléchargement (1)](https://github.com/user-attachments/assets/5ed5bce1-14fa-4fbc-a84b-ba035afc02b9)


3. **Biological Relevance**:
   - Phylogenetic trees elucidate evolutionary relationships between isolates, aiding in understanding pathogen lineages.
   - Even if the MinHash sketches computation takes time, computational efficiency of MinHash may enables analysis of large-scale genomic datasets, critical in epidemiological studies.

---

## Contact

For queries, reach out to:
- **GitHub**: [PhilRTFM](https://github.com/PhilRTFM)
