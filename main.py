import numpy as np
from functions import *

def main():
    # Specify the genome files provided for the challenge
    genome_files = ["R6.fa", "TIGR4.fa", "14412_3#82.contigs_velvet.fa", "14412_3#84.contigs_velvet.fa"]

    # Load sequences and extract k-mers
    print("Loading genome sequences...")
    sequences = {file: load_fasta(file) for file in genome_files}
    print("Extracting k-mers...")
    kmers = {file: extract_kmers(seq, k=14) for file, seq in sequences.items()}

    # Generate MinHash sketches
    print("Generating MinHash sketches...")
    sketches = {file: minhash_sketch(kmers[file], sketch_size=1000) for file in genome_files}

    # Calculate pairwise MinHash distances
    print("Calculating MinHash distances...")
    distance_matrix = np.array(
        [[0 if i == j else minhash_distance(sketches[genome_files[i]], sketches[genome_files[j]])
          for j in range(len(genome_files))]
         for i in range(len(genome_files))]
    )

    # Save distance matrix to a .txt file
    save_distance_matrix(distance_matrix, genome_files, "pairwise_distance_matrix.txt")

    # Generate phylogenetic tree using NJ algorithm
    print("\nGenerating phylogenetic tree using NJ algorithm...")
    tree = doNeighbourJoining(distance_matrix, genome_files)
    print(f"Generated Phylogenetic Tree (Newick):\n{tree}")

    # Save the tree in Newick format
    save_newick(tree, "phylogenetic_tree.nwk")

    # Plot the tree as a PNG file
    plot_tree_from_newick("phylogenetic_tree.nwk", "phylogenetic_tree.png")

if __name__ == "__main__":
    main()
