from collections import defaultdict, deque

def build_overlap_graph(kmers):
    graph = defaultdict(list)

    for kmer in kmers:
        prefix = kmer[:-1]
        suffix = kmer[1:]


        for other in kmers:
            if kmer != other and kmer[1:] == other[:-1]:
                graph[kmer].append(other)

    return graph

# Function to find a Hamiltonian path using depth-first search (DFS)
def find_hamiltonian_path(graph, start_node, visited, path):
    visited.add(start_node)
    path.append(start_node)

    if len(visited) == len(graph):
        return True

    for neighbor in graph[start_node]:
        if neighbor not in visited:
            if find_hamiltonian_path(graph, neighbor, visited, path):
                return True

    path.pop()
    visited.remove(start_node)
    return False

# Function to reconstruct the string from the Hamiltonian path
def reconstruct_from_path(path):
    # Start with the first k-mer
    reconstructed_string = path[0]

    # Append only the last character of each subsequent k-mer
    for kmer in path[1:]:
        reconstructed_string += kmer[-1]

    return reconstructed_string
kmers = ["ATT", "TTA", "TAC", "ACC", "CCA", "CAC"]

graph = build_overlap_graph(kmers)

# Find the Hamiltonian path
start_node = kmers[0]  # Start from the first k-mer
path = []
visited = set()

if find_hamiltonian_path(graph, start_node, visited, path):
    print("Hamiltonian path found:", path)

    # Reconstruct the original sequence from the path
    reconstructed_string = reconstruct_from_path(path)
    print("Reconstructed string:", reconstructed_string)
else:
    print("No Hamiltonian path found.")
