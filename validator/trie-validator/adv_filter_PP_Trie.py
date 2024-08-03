import datetime
import random
import logging
import trie as trie_utils
from Bio.SeqUtils import MeltingTemp as mt
from primer3 import calc_hairpin,calc_homodimer

MAX_HP = 4
PRIMER_BPS = 15
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

def calculate_tm(primer):
    return mt.Tm_NN(primer)

def filter_by_tm(primer, min_tm=30, max_tm=65):
    tm = calculate_tm(primer)
    return min_tm <= tm <= max_tm

def check_secondary_structures(primer):
    hairpin = calc_hairpin(primer) #hair pin is the primer binding to itself (hairpin like structure)
    homodimer = calc_homodimer(primer)#homodimer is the primer binding to itself
    return hairpin.structure_found == 0 and homodimer.structure_found == 0

#similaroty shown as hamming distance that we already checked in the trie
def blast_filters(primer):
    if filter_by_tm(primer):#filter by melting temperature
        if check_secondary_structures(primer):#filter by secondary structures
            return True
    return False

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
    trie = trie_utils.Trie()
    histogram = [0]*(PRIMER_BPS)
    count_primers_filtered_by_maxintercompliment = 0
    count_primers_filtered_by_blast = 0
    start_time = datetime.datetime.now()
    for index, primer in enumerate(all_primers):#(index, primer) = (0, 'TGGGAGGAGGAGGA')

        if index % 100000 == 0:
            show_logs(all_primers, count_primers_filtered_by_maxintercompliment,
                      histogram, index, primer_set, trie,start_time,count_primers_filtered_by_blast)
        #todo add heuristic filter?
        #todo tm and second structure...b4
        if not blast_filters(primer):
            count_primers_filtered_by_blast += 1
            continue
        if check_if_contains_complement_patterns_in_primers_set(primer, patterns_complement_of_max_inter_comp_len_set):
            if trie.is_valid_primer(primer, MIN_HAM,histogram):
                trie.insert(primer)
                update_complement_set(primer,patterns_complement_of_max_inter_comp_len_set)
                primer_set.add(primer)
        else:
            count_primers_filtered_by_maxintercompliment += 1

    for i in range(PRIMER_BPS):
        logging.info(f" at level {i} : {histogram[i]} primers throuwn away")
    logging.info(f"Number of primers filtered by max inter comp: {count_primers_filtered_by_maxintercompliment}")
    logging.info(f"Number of primers filtered by blast: {count_primers_filtered_by_blast}")
    return primer_set


def show_logs(all_primers, count_primers_filtered_by_maxintercompliment,
              histogram, index, primer_set, trie,start_time,count_primers_filtered_by_blast):
    leaves = trie.count_leaves()
    nodes = trie.count_tree_nodes()
    percent = "{:.2f}%".format(((index + 1) / len(all_primers)) * 100)
    logging.info(f" Processing primer {index + 1}/{len(all_primers)} ({percent}% of the primers), valid primers: {len(primer_set)}")
    if index % 5000000 == 0 and index != 0:
        logging.info(f" {histogram} primers disqualify by level")
        logging.info(f" tree status: {leaves} leaves and {nodes} nodes")
        logging.info(f" filtered by max inter comp: {count_primers_filtered_by_maxintercompliment} primers ")
        logging.info(f" filtered by blast: {count_primers_filtered_by_blast} primers ")
        elapsed_time = datetime.datetime.now() - start_time
        time_per_primer = elapsed_time / index
        remaining_primers = len(all_primers) - index
        remaining_time = time_per_primer * remaining_primers
        logging.info(f"Approximate time left: {remaining_time}")
        #show logs

#15 length file should have 700 million primers to process
def main():
    logging.info("Starting primer processing")
    with open('output_15_len_primer_internal_advance_generator.txt', 'r') as f:
        all_primers = f.readlines()

    all_primers = [item.strip() for item in all_primers]
    #random.shuffle(all_primers)

    valid_primers = process_primers(all_primers)

    logging.info(f"Total valid primers: {len(valid_primers)}")
    with open('output_15_len_primer_final_primers_tmp.txt', 'w') as f:
        for primer in valid_primers:
            f.write(f"{primer}\n")



    logging.info("Primer processing completed")

if __name__ == "__main__":
    main()
