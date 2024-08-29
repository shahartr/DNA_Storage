
import logging
import os
import signal
import sys

from Bio.Seq import Seq
from nupack import mfe
from tqdm import tqdm
import nupack
from primer3 import calc_tm
import primers_trie_tree as trie_utils

# Configuration and Constants
MAX_HP = 2
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
    with open('../../../../data/valid_primers.txt', 'w') as f:
        for primer in primer_set:
            f.write(f"{primer}\n")


def complement_strand(strand):
    try:
        return ''.join(complement_map[base] for base in strand)
    except KeyError as e:
        logging.error(f"Invalid character found: {e.args[0]} in primer: {strand}")
        return None


def show_logs(index, count_primers_filtered_by_maxintercompliment,
              histogram, primer_set, trie, count_primers_disqualify_by_tm,
              count_primers_disqualify_by_secondary_structure):
    leaves = trie.count_leaves()
    nodes = trie.count_tree_nodes()
    # percent = "{:.2f}%".format(((index + 1) / len(all_primers)) * 100)
    # logging.info(f" Processing primer {index + 1}/{len(all_primers)} ({percent}% of the primers), valid primers: {len(primer_set)}")
    logging.info(f" Processing primer {index + 1}  Valid primers: {len(primer_set)}")
    if index % 80000 == 0 and index != 0:
        logging.info(f" {histogram} primers disqualify by level")
        logging.info(f" tree status: {leaves} leaves and {nodes} nodes")
        logging.info(f" filtered by max inter comp: {count_primers_filtered_by_maxintercompliment} primers ")
        logging.info(f" filtered by tm: {count_primers_disqualify_by_tm} primers ")
        logging.info(f" filtered by secondary structure: {count_primers_disqualify_by_secondary_structure} primers ")
       # logging.info(f" filtered by secondary structure heuristic: {count_primers_disqualify_by_secondary_structure_huristic} primers ")


def check_if_contains_complement_patterns_in_primers_set(primer, patterns_complement_of_max_inter_comp_len_set):

    complement_primer = Seq(primer).reverse_complement()
    complement_primer = str(complement_primer)
    for length in range(1, MAX_INTER_COMP + 1):
        for i in range(PRIMER_BPS - length + 1):
            complement_pattern = complement_primer[i:i + length]
            if complement_pattern in patterns_complement_of_max_inter_comp_len_set:
                return False
    return True


def check_melting_temp(primer):
    tm_pass = calc_tm(primer)
    return tm_pass


def update_complement_set(primer, patterns_complement_of_max_inter_comp_len_set):
# complement_primer = complement_strand(primer)
    complement_primer = Seq(primer).reverse_complement()
    complement_primer = str(complement_primer)
    for i in range(PRIMER_BPS - MAX_INTER_COMP + 1):
        complement_pattern = complement_primer[i: i + MAX_INTER_COMP]
        patterns_complement_of_max_inter_comp_len_set.add(complement_pattern)

def check_inter_primer_complementarity(new_primer, final_set, max_complementarity=10):
    """Check for inter-primer homology with the final set of primers"""
    reverse_complement = str(Seq(new_primer).reverse_complement())

    for existing_primer in final_set:
        for i in range(PRIMER_BPS - max_complementarity):
            subseq = existing_primer[i:i + max_complementarity+1]
            if subseq in reverse_complement:
                return False  # Significant complementarity found

    return True  # No significant complementarity found


def hairpin_homodimer_check(primer):
    """Check for hairpin and homodimer formation using NUPACK mfe"""
    # Hairpin check
    results = mfe([primer],model=nupack.Model(material='dna', celsius=37))
    if results:
        energy = results[0].energy
        if energy < -3.0:
            return False
    return True




def process_file(file_path, primer_set, patterns_complement_of_max_inter_comp_len_set, trie, histogram,
                 count_primers_filtered_by_maxintercompliment, count_primers_disqualify_by_tm,
                 count_primers_disqualify_by_secondary_structure,count_primers_disqualify_by_secondary_structure_huristic):
    with open(file_path, 'r') as f:
        lines = f.readlines()

    total_primers = len(lines)
    logging.info(f"Processing file: {file_path} with {total_primers} primers")

    for index, primer in tqdm(enumerate(lines), total=total_primers, desc=f"Processing {file_path}"):
        if index % 20000 == 0:
            show_logs(index, count_primers_filtered_by_maxintercompliment,
                      histogram, primer_set, trie, count_primers_disqualify_by_tm,
                      count_primers_disqualify_by_secondary_structure)
        primer = primer.strip()
        if check_inter_primer_complementarity(primer, primer_set, MAX_INTER_COMP):
        #if check_if_contains_complement_patterns_in_primers_set(primer, patterns_complement_of_max_inter_comp_len_set):
            if trie.is_valid_primer(primer, MIN_HAM, histogram):
                if check_melting_temp(primer):
                    if hairpin_homodimer_check(primer):
                        trie.insert(primer)
                        #update_complement_set(primer, patterns_complement_of_max_inter_comp_len_set)
                        primer_set.add(primer)
                    else:
                        count_primers_disqualify_by_secondary_structure += 1

        else:
            count_primers_filtered_by_maxintercompliment += 1

    for i in range(PRIMER_BPS):
        logging.info(f"At level {i} : {histogram[i]} primers thrown away")

    logging.info(f"Number of primers filtered by max inter comp: {count_primers_filtered_by_maxintercompliment}")

def main():

    # i already generated the primers that meet criteria (homopolimer , gc content , self-compliment and saved them in a file) include this in the graph
    logging.info("Starting primer processing")

    primer_set = set()
    patterns_complement_of_max_inter_comp_len_set = set()
    trie = trie_utils.Trie()
    histogram = [0] * (PRIMER_BPS)
    count_primers_filtered_by_maxintercompliment = 0
    count_primers_disqualify_by_tm = 0
    count_primers_disqualify_by_secondary_structure = 0
    count_primers_disqualify_by_secondary_structure_huristic = 0

    def signal_handler(sig, frame):
        logging.info("Keyboard interrupt received. Saving valid primers...")
        save_valid_sofar(primer_set)
        sys.exit(0)
    signal.signal(signal.SIGINT, signal_handler)

    try:
        file_path = os.path.join('../../../../data/1billionFiles/output_20_sorted_V.txt')
        process_file(file_path, primer_set, patterns_complement_of_max_inter_comp_len_set, trie, histogram,
                     count_primers_filtered_by_maxintercompliment, count_primers_disqualify_by_tm,
                     count_primers_disqualify_by_secondary_structure,count_primers_disqualify_by_secondary_structure_huristic)

        logging.info(f"Total valid primers: {len(primer_set)}")
        with open('../../../../data/final_20_length000002.txt', 'w') as f:
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
