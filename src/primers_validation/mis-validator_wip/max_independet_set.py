import itertools

import networkx as nx
import sys
import pickle

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
    G.add_nodes_from(primer_list)
    # for 1 length primer there is 190000 choose 2 =18049905000 iteration !!!
    #for 10 length 12 choose 2 = 2483000 choose 2 = 3.e+12 iteration !!!

    for primer1, primer2 in itertools.combinations(primer_list, 2):
        hamming_dist = hamming_distance(primer1, primer2)
        if hamming_dist <= min_hamming_distance:
            G.add_edge(primer1, primer2, weight=hamming_dist)

        if calc_inter_complementarity(primer1, primer2, max_inter_complement_size):
            G.add_edge(primer1, primer2, weight=max_inter_complement_size)

    return G

def save_graph_to_file(graph, filename):
    with open(filename, 'wb') as f:
        pickle.dump(graph, f)

def run():
    with open('output_12_len_primer_internal_advance_generator.txt', 'r') as f:
        all_primers = f.readlines()

    all_primers = [item.strip() for item in all_primers]
    G = build_primer_graph(all_primers, MIN_HAM, MAX_INTER_COMP)
    save_graph_to_file(G, 'input/primer_graph.pkl')

if __name__ == "__main__":
    run()
