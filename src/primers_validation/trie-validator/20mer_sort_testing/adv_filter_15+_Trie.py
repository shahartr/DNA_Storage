import datetime
import logging
import os
import signal
import sys

from nupack import mfe
from primer3 import calc_tm
from tqdm import tqdm
import nupack
import trie as trie_utils

# Configuration and Constants
MAX_HP = 4
PRIMER_BPS = 20
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


def save_valid_sofar(primer_set):
    with open('valid_primers_mid_run.txt', 'w') as f:
        for primer in primer_set:
            f.write(f"{primer}\n")




def calculate_tm(primer):
    return calc_tm(primer)

def heuristic_hairpin_and_homodimer(primer, min_hairpin_len=3, max_hairpin_loop=4, min_homodimer_len=4):
    # Hairpin check
    primer_len = len(primer)
    for i in range(primer_len - min_hairpin_len):
        for j in range(i + min_hairpin_len + max_hairpin_loop, primer_len):
            if primer[i:i+min_hairpin_len] == primer[j-min_hairpin_len+1:j+1][::-1]:
                return False  # Potential hairpin found

    # Homodimer check
    rev_primer = primer[::-1]
    for i in range(primer_len - min_homodimer_len + 1):
        subseq = primer[i:i+min_homodimer_len]
        if subseq in rev_primer:
            return False  # Potential homodimer found

    return True  # No potential hairpins or homodimers found

def filter_by_tm(primer, min_tm=55, max_tm=60):
    tm = calculate_tm(primer)
    return min_tm <= tm <= max_tm


import nupack

def check_secondary_structures(primer):
    # Using NUPACK to check for secondary structures
    model = nupack.Model(material='dna', celsius=37) # 37 is standard
    result = nupack.mfe(strands=[primer], model=model)#mfe is the minimum free energy
    # The free energy is now directly accessible from the result
    min_dG = result[0].energy

    # Adjust the threshold as necessary
    return min_dG > -3.0


def blast_filters(primer):
    tm_pass = filter_by_tm(primer)
    structure_pass = check_secondary_structures(primer)
    return tm_pass, structure_pass


def complement_strand(strand):
    try:
        return ''.join(complement_map[base] for base in strand)
    except KeyError as e:
        logging.error(f"Invalid character found: {e.args[0]} in primer: {strand}")
        return None


def show_logs(index, count_primers_filtered_by_maxintercompliment,
              histogram, primer_set, trie, count_primers_disqualify_by_tm,
              count_primers_disqualify_by_secondary_structure,count_primers_disqualify_by_secondary_structure_huristic):
    leaves = trie.count_leaves()
    nodes = trie.count_tree_nodes()
    # percent = "{:.2f}%".format(((index + 1) / len(all_primers)) * 100)
    # logging.info(f" Processing primer {index + 1}/{len(all_primers)} ({percent}% of the primers), valid primers: {len(primer_set)}")
    logging.info(f" Processing primer {index + 1}  Valid primers: {len(primer_set)}")
    if index % 4000000 == 0 and index != 0:
        logging.info(f" {histogram} primers disqualify by level")
        logging.info(f" tree status: {leaves} leaves and {nodes} nodes")
        logging.info(f" filtered by max inter comp: {count_primers_filtered_by_maxintercompliment} primers ")
        logging.info(f" filtered by tm: {count_primers_disqualify_by_tm} primers ")
        logging.info(f" filtered by secondary structure: {count_primers_disqualify_by_secondary_structure} primers ")
       # logging.info(f" filtered by secondary structure heuristic: {count_primers_disqualify_by_secondary_structure_huristic} primers ")


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


def unique_kmers(primer, k=4):
    kmers = set()
    for i in range(len(primer) - k + 1):
        kmer = primer[i:i + k]
        if kmer in kmers:
            return False
        kmers.add(kmer)
    return True


def heuristic_filter(primer):
    return unique_kmers(primer)


def update_complement_set(primer, patterns_complement_of_max_inter_comp_len_set):
    complement_primer = complement_strand(primer)
    for i in range(PRIMER_BPS - MAX_INTER_COMP + 1):
        complement_pattern = complement_primer[i: i + MAX_INTER_COMP]
        patterns_complement_of_max_inter_comp_len_set.add(complement_pattern)


def check_melting_temp(primer):
    tm_pass = filter_by_tm(primer)
    return tm_pass


def process_file(file_path, primer_set, patterns_complement_of_max_inter_comp_len_set, trie, histogram,
                 count_primers_filtered_by_inter_complement, count_primers_disqualify_by_tm,
                 count_primers_disqualify_by_secondary_structure, count_primers_disqualify_by_secondary_structure_huristic):
    with open(file_path, 'r') as f:
        lines = f.readlines()

    total_primers = len(lines)
    logging.info(f"Processing file: {file_path} with {total_primers} primers")

    for index, primer in tqdm(enumerate(lines), total=total_primers, desc=f"Processing {file_path}"):
        if index % 2000000 == 0:
            show_logs(index, count_primers_filtered_by_inter_complement,
                      histogram, primer_set, trie, count_primers_disqualify_by_tm,
                      count_primers_disqualify_by_secondary_structure, count_primers_disqualify_by_secondary_structure_huristic)
        primer = primer.strip()
        # if not heuristic_hairpin_and_homodimer(primer):
        #     count_primers_disqualify_by_secondary_structure_huristic += 1
        #     continue
        if not check_melting_temp(primer):#tm check
            count_primers_disqualify_by_tm += 1
            continue
        # if not check_secondary_structures(primer):
        #     count_primers_disqualify_by_secondary_structure += 1
        #     continue
        if check_if_contains_complement_patterns_in_primers_set(primer, patterns_complement_of_max_inter_comp_len_set):
            if trie.is_valid_primer(primer, MIN_HAM, histogram):
                trie.insert(primer)
                update_complement_set(primer, patterns_complement_of_max_inter_comp_len_set)
                primer_set.add(primer)
        else:
            count_primers_filtered_by_inter_complement += 1

    for i in range(PRIMER_BPS):
        logging.info(f"At level {i} : {histogram[i]} primers thrown away")

    logging.info(f"Number of primers filtered by max inter comp: {count_primers_filtered_by_inter_complement}")
    logging.info(f"Number of primers filtered by tm: {count_primers_disqualify_by_tm}")
    # logging.info(
    #     f"Number of primers filtered by secondary structure: {count_primers_disqualify_by_secondary_structure}")



def main():
    logging.info("Starting primer processing")

    primer_set = set()
    patterns_complement_of_max_inter_comp_len_set = set()
    trie = trie_utils.Trie()
    histogram = [0] * (PRIMER_BPS)
    count_primers_filtered_by_maxintercompliment = 0
    count_primers_disqualify_by_tm = 0
    count_primers_disqualify_by_secondary_structure = 0
    count_primers_disqualify_by_secondary_structure_huristic = 0


    temp_folder = 'temp/final_output'

    def signal_handler(sig, frame):
        logging.info("Keyboard interrupt received. Saving valid primers...")
        save_valid_sofar(primer_set)
        sys.exit(0)
    signal.signal(signal.SIGINT, signal_handler)

    try:
        for i in range(5):
            file_path = os.path.join(temp_folder, f'final_sorted_part_{i}.txt')
            process_file(file_path, primer_set, patterns_complement_of_max_inter_comp_len_set, trie, histogram,
                         count_primers_filtered_by_maxintercompliment, count_primers_disqualify_by_tm,
                         count_primers_disqualify_by_secondary_structure,count_primers_disqualify_by_secondary_structure_huristic)

        logging.info(f"Total valid primers: {len(primer_set)}")
        with open('final_20_length.txt', 'w') as f:
            for primer in primer_set:
                f.write(f"{primer}\n")

        logging.info("Primer processing completed")
    except Exception as e:
        logging.error(f"An error occurred: {e}")
        #is KeyboardInterrupt exception
        if e.__class__.__name__ == 'KeyboardInterrupt':
            logging.info("Saving valid primers...")
            save_valid_sofar(primer_set)
        else:
            logging.error("Exiting due to error")
            #y to save n to exit
            save = input("Do you want to save the valid primers? (y/n): ")
            if save.lower() == 'y':
                logging.info("Saving valid primers...")
                save_valid_sofar(primer_set)


        sys.exit(1)


if __name__ == "__main__":
    main()
