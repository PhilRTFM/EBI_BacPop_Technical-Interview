import hashlib
import numpy as np
import csv
from Bio import SeqIO
import matplotlib.pyplot as plt
from Bio import Phylo
from io import StringIO
import pandas as pd

# Existing Functions (No changes here)
def write_to_csv(file_name, headers, rows):
    with open(file_name, mode="w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(headers)
        writer.writerows(rows)

def load_fasta(file_path):
    sequences = []
    for record in SeqIO.parse(file_path, "fasta"):
        sequences.append(str(record.seq))
    concatenated_sequence = "".join(sequences)
    return concatenated_sequence

def extract_kmers(sequence, k=14):
    return {sequence[i:i+k] for i in range(len(sequence) - k + 1)}

def hash_kmer(kmer):
    reverse_comp = kmer[::-1].translate(str.maketrans("ACGT", "TGCA"))
    return min(int(hashlib.md5(kmer.encode()).hexdigest(), 16),
               int(hashlib.md5(reverse_comp.encode()).hexdigest(), 16))

def minhash_sketch(kmers, sketch_size=1000, genome_name=None, save_sketch=False):
    """
    Generate a MinHash sketch of a specified size and optionally save it to a file.

    Args:
        kmers (set of str): Set of k-mers to hash.
        sketch_size (int): Number of lowest hash values to include in the sketch.
        genome_name (str): Name of the genome (used for file naming).
        save_sketch (bool): If True, save the sketch to a file.

    Returns:
        list: MinHash sketch (list of lowest hash values).
    """
    hashes = {hash_kmer(kmer) for kmer in kmers}
    sketch = sorted(hashes)[:sketch_size]

    if save_sketch and genome_name:
        file_name = f"{genome_name}_sketch.txt"
        with open(file_name, "w") as f:
            for hash_value in sketch:
                f.write(f"{hash_value}\n")
        print(f"Sketch saved to {file_name}")

    return sketch

def full_distance(kmers_a, kmers_b):
    """
    Calculate the full Jaccard distance between two sets of k-mers.

    Args:
        kmers_a (set): Set of k-mers from genome A.
        kmers_b (set): Set of k-mers from genome B.

    Returns:
        float: Jaccard distance between the two sets.
    """
    set_a, set_b = set(kmers_a), set(kmers_b)
    intersection = len(set_a & set_b)
    union = len(set_a | set_b)
    return 1 - (intersection / union)

def minhash_distance(sketch_a, sketch_b):
    set_a, set_b = set(sketch_a), set(sketch_b)
    intersection = len(set_a & set_b)
    union = len(set_a | set_b)
    return 1 - (intersection / union)

def calculateQ(d):
    r = d.shape[0]
    q = np.zeros((r, r))
    for i in range(r):
        for j in range(r):
            if i != j:
                sumI = np.sum(d[i, :])
                sumJ = np.sum(d[j, :])
                q[i, j] = (r - 2) * d[i, j] - sumI - sumJ
    return q

def findLowestPair(q):
    minVal = np.inf
    minIndex = (0, 0)
    for i in range(q.shape[0]):
        for j in range(i + 1, q.shape[0]):
            if q[i, j] < minVal:
                minVal = q[i, j]
                minIndex = (i, j)
    return minIndex

def calculateNewDistanceMatrix(f, g, d):
    r = d.shape[0]
    indices = [i for i in range(r) if i != f and i != g]
    new_d = np.zeros((r - 1, r - 1))
    for i, ii in enumerate(indices):
        for j, jj in enumerate(indices):
            new_d[i + 1, j + 1] = d[ii, jj]
    for i, ii in enumerate(indices):
        new_distance = (d[f, ii] + d[g, ii] - d[f, g]) / 2
        new_d[0, i + 1] = new_distance
        new_d[i + 1, 0] = new_distance
    new_d[0, 0] = 0
    return new_d

def doNeighbourJoining(d, labels):
    while len(labels) > 2:
        q = calculateQ(d)
        i, j = findLowestPair(q)
        r = d.shape[0]
        sumI = np.sum(d[i, :])
        sumJ = np.sum(d[j, :])
        distToNewNode = (d[i, j] + (sumI - sumJ) / (r - 2)) / 2
        new_label = f"({labels[i]}:{distToNewNode:.4f},{labels[j]}:{distToNewNode:.4f})"
        labels[i] = new_label
        d = calculateNewDistanceMatrix(i, j, d)
        labels.pop(j)
    tree = f"({labels[0]}:{d[0, 1]:.4f},{labels[1]}:{d[0, 1]:.4f});"
    return tree

def save_newick(tree, output_path="tree.nwk"):
    with open(output_path, "w") as f:
        f.write(tree)
    print(f"Newick tree saved to {output_path}")

def plot_tree_from_newick(newick_path, output_path="tree.png"):
    with open(newick_path, "r") as f:
        tree = Phylo.read(f, "newick")
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(1, 1, 1)
    Phylo.draw(tree, axes=ax, do_show=False)
    plt.savefig(output_path)
    print(f"Tree plotted and saved as {output_path}")
    plt.close(fig)

def construct_distance_matrix(labels, distances, display_graphical=False, output_path="pairwise_distance_matrix.txt"):
    """
    Construct a square distance matrix from pairwise distances, optionally display it graphically,
    and save it to a text file.

    Args:
        labels (list of str): Names of the sequences.
        distances (dict): Pairwise distances as {(label_a, label_b): distance}.
        display_graphical (bool): If True, print the distance matrix as a table.
        output_path (str): Path to save the distance matrix as a text file.

    Returns:
        np.ndarray: Symmetric distance matrix.
    """
    size = len(labels)
    matrix = np.zeros((size, size))
    label_to_index = {label: idx for idx, label in enumerate(labels)}

    # Populate the matrix
    for (label_a, label_b), distance in distances.items():
        i, j = label_to_index[label_a], label_to_index[label_b]
        matrix[i, j] = distance
        matrix[j, i] = distance

    # Prepare the graphical output
    header = "\t" + "\t".join(labels)
    rows = [f"{labels[i]}\t" + "\t".join(f"{matrix[i, j]:.4f}" for j in range(size)) for i in range(size)]
    graphical_output = header + "\n" + "\n".join(rows)

    # Display the graphical representation
    if display_graphical:
        print("Distance Matrix (Graphical Representation):")
        print(graphical_output)

    # Save the distance matrix to a text file
    with open(output_path, "w") as f:
        f.write(graphical_output)
    print(f"Distance matrix saved to {output_path}")

    return matrix
