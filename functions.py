import hashlib
import numpy as np
import csv
from Bio import SeqIO
import matplotlib.pyplot as plt
from Bio import Phylo
from io import StringIO


def write_to_csv(file_name, headers, rows):
    """Write the given rows to a CSV file with specified headers."""
    with open(file_name, mode="w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(headers)
        writer.writerows(rows)

def load_fasta(file_path):
    """Load a FASTA file and return the concatenated sequences as a string."""
    sequences = []
    for record in SeqIO.parse(file_path, "fasta"):
        sequences.append(str(record.seq))
    concatenated_sequence = "".join(sequences)
    return concatenated_sequence

def extract_kmers(sequence, k=14):
    """Extract k-mers from a given sequence."""
    return {sequence[i:i+k] for i in range(len(sequence) - k + 1)}

def hash_kmer(kmer):
    """Return the hash (integer) of a k-mer and its reverse complement."""
    reverse_comp = kmer[::-1].translate(str.maketrans("ACGT", "TGCA"))
    return min(int(hashlib.md5(kmer.encode()).hexdigest(), 16),
               int(hashlib.md5(reverse_comp.encode()).hexdigest(), 16))

def minhash_sketch(kmers, sketch_size=1000):
    """Generate a MinHash sketch of a specified size."""
    hashes = {hash_kmer(kmer) for kmer in kmers}
    return sorted(hashes)[:sketch_size]

def minhash_distance(sketch_a, sketch_b):
    """Calculate the Jaccard distance between two MinHash sketches."""
    set_a, set_b = set(sketch_a), set(sketch_b)
    intersection = len(set_a & set_b)
    union = len(set_a | set_b)
    return 1 - (intersection / union)

def calculateQ(d):
    """Calculate the Q-matrix for the NJ algorithm."""
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
    """Find the pair of indices with the lowest Q value."""
    minVal = np.inf
    minIndex = (0, 0)
    for i in range(q.shape[0]):
        for j in range(i + 1, q.shape[0]):
            if q[i, j] < minVal:
                minVal = q[i, j]
                minIndex = (i, j)
    return minIndex

def calculateNewDistanceMatrix(f, g, d):
    """Calculate a new distance matrix after joining two nodes."""
    r = d.shape[0]
    indices = [i for i in range(r) if i != f and i != g]

    # Create a new distance matrix of reduced size
    new_d = np.zeros((r - 1, r - 1))

    # Copy over the old data to the new matrix
    for i, ii in enumerate(indices):
        for j, jj in enumerate(indices):
            new_d[i + 1, j + 1] = d[ii, jj]

    # Calculate distances for the new row and column
    for i, ii in enumerate(indices):
        new_distance = (d[f, ii] + d[g, ii] - d[f, g]) / 2
        new_d[0, i + 1] = new_distance
        new_d[i + 1, 0] = new_distance

    # Distance between the merged node and itself is always 0
    new_d[0, 0] = 0
    return new_d

def doNeighbourJoining(d, labels):
    """Perform the NJ algorithm to build a phylogenetic tree."""
    while len(labels) > 2:
        q = calculateQ(d)
        i, j = findLowestPair(q)

        # Calculate branch lengths
        r = d.shape[0]
        sumI = np.sum(d[i, :])
        sumJ = np.sum(d[j, :])
        distToNewNode = (d[i, j] + (sumI - sumJ) / (r - 2)) / 2

        # Update labels and distance matrix
        new_label = f"({labels[i]}:{distToNewNode:.4f},{labels[j]}:{distToNewNode:.4f})"
        labels[i] = new_label
        d = calculateNewDistanceMatrix(i, j, d)
        labels.pop(j)

    # Final two labels
    tree = f"({labels[0]}:{d[0, 1]:.4f},{labels[1]}:{d[0, 1]:.4f});"
    return tree

def save_newick(tree, output_path="tree.nwk"):
    """Save the tree in Newick format."""
    with open(output_path, "w") as f:
        f.write(tree)
    print(f"Newick tree saved to {output_path}")

def plot_tree_from_newick(newick_path, output_path="tree.png"):
    """Plot a phylogenetic tree from a Newick file and save as PNG."""
    with open(newick_path, "r") as f:
        tree = Phylo.read(f, "newick")
    
    # Set up the plot
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(1, 1, 1)

    # Draw the tree
    Phylo.draw(tree, axes=ax, do_show=False)
    
    # Save the plot
    plt.savefig(output_path)
    print(f"Tree plotted and saved as {output_path}")
    plt.close(fig)

def save_distance_matrix(matrix, labels, output_path="distance_matrix.txt"):
    """Save the pairwise distance matrix to a .txt file."""
    with open(output_path, "w") as f:
        # Write the header row
        f.write("\t" + "\t".join(labels) + "\n")
        # Write each row
        for i, row in enumerate(matrix):
            row_str = "\t".join(f"{val:.4f}" for val in row)
            f.write(f"{labels[i]}\t{row_str}\n")
    print(f"Distance matrix saved to {output_path}")
