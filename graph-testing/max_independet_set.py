import itertools
import networkx as nx
import mis_benchmark as mb

MAX_HP = 2
PRIMER_BPS = 12
MAX_SELF_COMP = 4
MAX_INTER_COMP = 9
MIN_HAM = 6

complement_map = {
    'G': 'C',
    'C': 'G',
    'T': 'A',
    'A': 'T'
}

#grqph modeling and finding max independet set using mis solvers

def hamming_distance(strand1, strand2):
    dist = 0
    for i in range(0, len(strand1)):
        if strand1[i] != strand2[i]:
            dist += 1
    return dist

def complement_strand(strand):
    complement = ""
    for base in strand:
        complement += complement_map[base]
    return complement

def calc_inter_complementarity(strand1, strand2, length):
    strand1_comp = complement_strand(strand1)
    for i in range(0, len(strand1_comp) - length):
        if strand1_comp[i:(i + length)] in strand2:
            return True
    return False


def build_primer_graph(primer_list, min_hamming_distance, max_inter_complement_size):
    G = nx.Graph()

    # Add vertices (primers) to the graph
    G.add_nodes_from(primer_list)

    # Check Hamming distance and complementarity between primers
    for primer1, primer2 in itertools.combinations(primer_list, 2):
        # Check Hamming distance
        hamming_dist = hamming_distance(primer1, primer2)
        if hamming_dist <= min_hamming_distance:
            G.add_edge(primer1, primer2, weight=hamming_dist)

        # Check complementarity
        if calc_inter_complementarity(primer1, primer2, max_inter_complement_size):
            G.add_edge(primer1, primer2, weight=max_inter_complement_size)


    return G
def run():
    with open('output_12_len_primer_internal_advance_generator.txt', 'r') as f:
        all_primers = f.readlines()

    # Remove newline characters from each item
    all_primers = [item.strip() for item in all_primers]

    # Build the graph
    G = build_primer_graph(all_primers, MIN_HAM, MAX_INTER_COMP)

    # Convert the graph to a format suitable for mis-benchmark
    nodes = list(G.nodes)
    edges = list(G.edges)
    graph = mb.Graph()
    for node in nodes:
        graph.add_node(node)
    for edge in edges:
        graph.add_edge(*edge)

    # Solve the MIS problem using the framework
    mis_solver = mb.algorithms.MIS(graph)
    independent_set = mis_solver.run()

    # Print the size of the independent set and the primers in the set
    print("Size of the independent set:", len(independent_set))
    print("Independent set primers:", independent_set)

if __name__ == "__main__":
    run()