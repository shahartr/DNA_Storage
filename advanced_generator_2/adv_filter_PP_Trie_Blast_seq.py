import random
import logging
import primer3
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Blast import NCBIWWW, NCBIXML
from primer3 import calcHairpin, calcHomodimer

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

    def search_with_hamming_distance(self, node, primer, index, mismatches, max_mismatches):
        if mismatches > max_mismatches:
            return False
        if index == len(primer):
            return node.is_end_of_primer and mismatches <= max_mismatches

        char = primer[index]
        for child_char, child_node in node.children.items():
            new_mismatches = mismatches + (1 if char != child_char else 0)
            if self.search_with_hamming_distance(child_node, primer, index + 1, new_mismatches, max_mismatches):
                return True
        return False

    def is_valid_primer(self, primer, max_mismatches):
        return not self.search_with_hamming_distance(self.root, primer, 0, 0, max_mismatches)

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
        patterns_complement_of_max_inter_comp_len_set.add(complement_pattern)
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

def unique_kmers(primer, k=4):
    kmers = set()
    for i in range(len(primer) - k + 1):
        kmer = primer[i:i + k]
        if kmer in kmers:
            return False
        kmers.add(kmer)
    return True

def heuristic_filter(primer):
    max_hp = max_homopolymer(primer)
    if max_hp > MAX_HP:
        return False
    if not unique_kmers(primer):
        return False
    return True

def calculate_tm(primer):
    return mt.Tm_NN(primer)

def filter_by_tm(primer, min_tm=30, max_tm=65):
    tm = calculate_tm(primer)
    return min_tm <= tm <= max_tm

def blast_sequence(sequence):
    result_handle = NCBIWWW.qblast("blastn", "nt", sequence)
    blast_record = NCBIXML.read(result_handle)
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < 0.01:
                return False
    return True

def simulate_multiplex_pcr_validation(primer_set):
    # Placeholder function for actual Multiplex PCR validation simulation
    # In practice, this would involve checking for primer-primer interactions, compatibility, and more
    valid_primers = set()
    for primer in primer_set:
        if check_secondary_structures(primer):
            valid_primers.add(primer)
    return valid_primers

def check_secondary_structures(primer):
    hairpin = calcHairpin(primer)
    homodimer = calcHomodimer(primer)
    return hairpin.structure_found == 0 and homodimer.structure_found == 0

def process_primers(all_primers):
    primer_set = set()
    patterns_complement_of_max_inter_comp_len_set = set()
    trie = Trie()

    for index, primer in enumerate(all_primers):
        if index % 100000 == 0:
            logging.info(f"Processing primer {index + 1}/{len(all_primers)}")
        if filter_by_tm(primer):
            if check_if_contains_complement_patterns_in_primers_set(primer, patterns_complement_of_max_inter_comp_len_set):
                if trie.is_valid_primer(primer, MIN_HAM - 1):
                    if blast_sequence(primer):
                        trie.insert(primer)
                        primer_set.add(primer)
                        if len(primer_set) % 500 == 0:
                            logging.info(f"Primer set size: {len(primer_set)}")

    return primer_set

def main():
    logging.info("Starting primer processing")
    with open('output_14_len_primer_internal_advance_generator.txt', 'r') as f:
        all_primers = f.readlines()

    all_primers = [item.strip() for item in all_primers]
    random.shuffle(all_primers)

    valid_primers = process_primers(all_primers)

    logging.info(f"Total valid primers: {len(valid_primers)}")
    with open('output_14_len_primer_final_primers2.txt', 'w') as f:
        for primer in valid_primers:
            f.write(f"{primer}\n")

    logging.info("Primer processing completed")

if __name__ == "__main__":
    main()
