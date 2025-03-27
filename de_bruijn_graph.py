from collections import defaultdict

# Function to construct the de Bruijn graph from k-mers
def de_bruijn_graph(kmers):
    graph = defaultdict(list)

    for kmer in kmers:
        prefix = kmer[:-1]  # Prefix of length k-1
        suffix = kmer[1:]   # Suffix of length k-1
        graph[prefix].append(suffix)

    return graph

# Function to print the de Bruijn graph
def print_de_bruijn_graph(graph):
    for prefix, suffixes in graph.items():
        for suffix in suffixes:
            print(f"{prefix} -> {suffix}")

# Example usage
kmers = ["ATT", "TTA", "TAC", "ACC", "CCA", "CAC"]

# Build the de Bruijn graph
graph = de_bruijn_graph(kmers)

# Print the graph
print("De Bruijn Graph:")
print_de_bruijn_graph(graph)
