import random
import logging

MAX_HP = 2
PRIMER_BPS = 14
MAX_SELF_COMP = 4
MAX_INTER_COMP = 10
MIN_HAM = 6

complement_map = {
    'G': 'C',
    'C': 'G',
    'T': 'A',
    'A': 'T'
}

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class TrieNode:
    def __init__(self):
        self.children = {}
        self.is_end_of_primer = False

class Trie:
    def __init__(self):
        self.root = TrieNode()

    def insert(self, primer):
        node = self.root
        for char in primer:
            if char not in node.children:
                node.children[char] = TrieNode()
            node = node.children[char]
        node.is_end_of_primer = True
        
        
    # hamming distance is the number of positions at which the corresponding symbols are different
    #in our case we need every primer in the final set to have hamming distance of at least 6(max_mismatches)
    # basically we are searching (in tree that built from the primer added to the final set)
    # for a primer with hamming distance less than max_mismatches
    # if we find such primer then we return False
    #  ( means that in the final set primer with hamming distance less than 6 is already present in the tree
    # and that why we cant add him to the final set 
    def search_with_hamming_distance(self, node, primer, lvl, mismatches, max_mismatches, histogram):
        #case enough mismatches (hamming distance to all other in the set ) already found
        if mismatches > max_mismatches:
            return False
        #case EndOf primer chaeck if the primer have hamming distance less than max_mismatches
        # to all other in the final set (cant add it to the final set!)
        if lvl == len(primer):
            return node.is_end_of_primer and mismatches < max_mismatches

        char = primer[lvl]
        for child_char, child_node in node.children.items():
            new_mismatches = mismatches + (1 if char != child_char else 0)
            if self.search_with_hamming_distance(child_node, primer, lvl + 1, new_mismatches, max_mismatches,histogram):
                histogram[lvl]+=1
                return True#if found primer with hamming distance less than max_mismatches
        return False

    def is_valid_primer(self, primer, max_mismatches,histogram):
        return not self.search_with_hamming_distance(self.root, primer, 0, 0, max_mismatches,histogram)

def complement_strand(strand):
    try:
        complement = ''.join(complement_map[base] for base in strand)
        return complement
    except KeyError as e:
        logging.error(f"Invalid character found: {e.args[0]} in primer: {strand}")
        return None

def check_if_contains_complement_patterns_in_primers_set(primer, patterns_complement_of_max_inter_comp_len_set):
    complement_primer = complement_strand(primer)
    for i in range(len(primer) - MAX_INTER_COMP + 1):
        complement_pattern = complement_primer[i:i + MAX_INTER_COMP]
        if complement_pattern in patterns_complement_of_max_inter_comp_len_set:
            return False
    return True

def max_homopolymer(strand):
    max_count = 1
    cur_count = 1
    for i in range(1, len(strand)):
        if strand[i] == strand[i - 1]:
            cur_count += 1
        else:
            max_count = max(max_count, cur_count)
            cur_count = 1
    return max(max_count, cur_count)


#huristic filter
def unique_kmers(primer, k=4):
    kmers = set()
    for i in range(len(primer) - k + 1):
        kmer = primer[i:i + k]
        if kmer in kmers:
            return False
        kmers.add(kmer)
    return True

def heuristic_filter(primer):
    if not unique_kmers(primer):
        return False
    return True

def update_complement_set(primer, patterns_complement_of_max_inter_comp_len_set):
    complement_primer = complement_strand(primer)
    for i in range(PRIMER_BPS - MAX_INTER_COMP + 1):
        complement_pattern = complement_primer[i: i + MAX_INTER_COMP]
        patterns_complement_of_max_inter_comp_len_set.add(complement_pattern)

def process_primers(all_primers):
    primer_set = set()
    patterns_complement_of_max_inter_comp_len_set = set()
    trie = Trie()
    histogram = [0]*PRIMER_BPS
    countHowManyPrimersFilteredByMaxInterComp = 0

    for index, primer in enumerate(all_primers):
        if index % 100000 == 0:
            #(index + 1)/len(all_primers) with .2f precision
            precent = "{:.2f}".format((index + 1)/len(all_primers))
            logging.info(f"Processing primer {index + 1}/{len(all_primers)} ({precent}% of the primers), valid primers: {len(primer_set)}")
        #todo add heuristic filter?
        if check_if_contains_complement_patterns_in_primers_set(primer, patterns_complement_of_max_inter_comp_len_set):
            if trie.is_valid_primer(primer, MIN_HAM,histogram):
                trie.insert(primer)
                update_complement_set(primer,patterns_complement_of_max_inter_comp_len_set)
                primer_set.add(primer)
        else:
            countHowManyPrimersFilteredByMaxInterComp += 1

    for i in range(PRIMER_BPS):
        logging.info(f" at level {i} : {histogram[i]} primers throuwn away")
    logging.info(f"Number of primers filtered by max inter comp: {countHowManyPrimersFilteredByMaxInterComp}")
    return primer_set

def main():
    logging.info("Starting primer processing")
    with open('output_14_len_primer_internal_advance_generator.txt', 'r') as f:
        all_primers = f.readlines()

    all_primers = [item.strip() for item in all_primers]
    #random.shuffle(all_primers)

    valid_primers = process_primers(all_primers)

    logging.info(f"Total valid primers: {len(valid_primers)}")
    with open('output_14_len_primer_final_primers3.txt', 'w') as f:
        for primer in valid_primers:
            f.write(f"{primer}\n")



    logging.info("Primer processing completed")

if __name__ == "__main__":
    main()
