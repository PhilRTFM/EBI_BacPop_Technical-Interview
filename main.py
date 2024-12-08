import time
from functions import *

# Step 1: Inputs and Preprocessing
def preprocess_genome_files(genome_files, k):
    print("Loading genome files and extracting k-mers...")
    sequences = {file: load_fasta(file) for file in genome_files}
    kmers = {file: extract_kmers(sequences[file], k=k) for file in genome_files}
    print(f"Genome files processed: {list(kmers.keys())}")
    return kmers

# Step 2: Full Jaccard Distance Matrix
def compute_full_distance_matrix(genome_files, kmers):
    print("\nCalculating Full Jaccard Distance Matrix...")
    start_time = time.time()
    full_distances = {
        (genome_files[i], genome_files[j]): full_distance(kmers[genome_files[i]], kmers[genome_files[j]])
        for i in range(len(genome_files)) for j in range(i + 1, len(genome_files))
    }
    full_distance_matrix = construct_distance_matrix(genome_files, full_distances, display_graphical=False)
    print(f"Full Jaccard Distance Matrix:\n{full_distance_matrix}")
    print(f"Time for Full Distance Matrix: {time.time() - start_time:.4f} seconds")
    return full_distance_matrix

# Step 3: MinHash Distance Matrix
def compute_minhash_distance_matrix(genome_files, kmers, sketch_size):
    print("\nGenerating MinHash sketches...")
    start_time = time.time()
    minhash_sketches = {file: minhash_sketch(kmers[file], sketch_size=sketch_size) for file in genome_files}
    print(f"Time for MinHash Sketches: {time.time() - start_time:.4f} seconds")
    
    print("\nCalculating MinHash Distance Matrix...")
    start_time = time.time()
    minhash_distances = {
        (genome_files[i], genome_files[j]): minhash_distance(
            minhash_sketches[genome_files[i]], minhash_sketches[genome_files[j]]
        )
        for i in range(len(genome_files)) for j in range(i + 1, len(genome_files))
    }
    minhash_distance_matrix = construct_distance_matrix(genome_files, minhash_distances, display_graphical=False)
    print(f"MinHash Jaccard Distance Matrix:\n{minhash_distance_matrix}")
    print(f"Time for MinHash Distance Matrix: {time.time() - start_time:.4f} seconds")
    return minhash_distance_matrix

# Step 4: Neighbor-Joining Trees
def generate_trees(full_matrix, minhash_matrix, genome_files):
    print("\nGenerating Neighbor-Joining Trees...")
    
    # Full Jaccard Tree
    start_time = time.time()
    full_tree = doNeighbourJoining(full_matrix, genome_files)
    print(f"Time for Full Jaccard Tree: {time.time() - start_time:.4f} seconds")
    save_newick(full_tree, "phylogenetic_full_tree.nwk")
    
    # MinHash Jaccard Tree
    start_time = time.time()
    minhash_tree = doNeighbourJoining(minhash_matrix, genome_files)
    print(f"Time for MinHash Jaccard Tree: {time.time() - start_time:.4f} seconds")
    save_newick(minhash_tree, "phylogenetic_minhash_tree.nwk")
    
    return full_tree, minhash_tree

# Step 5: Plotting Trees
def plot_trees(full_tree_path, minhash_tree_path):
    fig, axes = plt.subplots(1, 2, figsize=(16, 8))
    
    # Full Jaccard Tree
    ax1 = axes[0]
    with open(full_tree_path, "r") as f:
        full_tree = Phylo.read(f, "newick")
    Phylo.draw(full_tree, do_show=False, axes=ax1)
    ax1.set_title("Full Jaccard Distance Tree")
    fig.savefig("phylogenetic_full_tree.png")
    
    # MinHash Jaccard Tree
    ax2 = axes[1]
    with open(minhash_tree_path, "r") as f:
        minhash_tree = Phylo.read(f, "newick")
    Phylo.draw(minhash_tree, do_show=False, axes=ax2)
    ax2.set_title("MinHash Jaccard Distance Tree")
    fig.savefig("phylogenetic_minhash_tree.png")
    
    plt.tight_layout()
    plt.show()

# Main Execution
if __name__ == "__main__":
    # Inputs
    genome_files = ["R6.fa", "TIGR4.fa", "14412_3#82.contigs_velvet.fa", "14412_3#84.contigs_velvet.fa"]
    k = 14
    sketch_size = 1000
    
    # Steps 1-5
    kmers = preprocess_genome_files(genome_files, k)
    full_matrix = compute_full_distance_matrix(genome_files, kmers)
    minhash_matrix = compute_minhash_distance_matrix(genome_files, kmers, sketch_size)
    full_tree, minhash_tree = generate_trees(full_matrix, minhash_matrix, genome_files)
    plot_trees("phylogenetic_full_tree.nwk", "phylogenetic_minhash_tree.nwk")
