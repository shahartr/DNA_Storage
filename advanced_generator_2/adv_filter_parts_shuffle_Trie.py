import random
import logging
import sys
from concurrent.futures import ProcessPoolExecutor
import trie as trie_utils
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

def process_primers(all_primers,run_id):
    primer_set = set()
    patterns_complement_of_max_inter_comp_len_set = set()
    trie = trie_utils.Trie()
    histogram = [0]*PRIMER_BPS
    countHowManyPrimersFilteredByMaxInterComp = 0

    for index, primer in enumerate(all_primers):
        if index % 100000 == 0:
            percent = "{:.2f}%".format(((index + 1) / len(all_primers)) * 100)
            logging.info(f"Run {run_id}: Processing primer {index + 1}/{len(all_primers)} ({percent}% of the primers), valid primers: {len(primer_set)}")
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

def run_primer_processing(run_id):
    logging.basicConfig(filename=f'logs/run_{run_id}_log_PP.txt', level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    logging.info("Starting primer processing")
    with open('output_14_len_primer_internal_advance_generator.txt', 'r') as f:
        all_primers = f.readlines()

    all_primers = [item.strip() for item in all_primers]

    num_parts = 320
    part_size = len(all_primers) // num_parts

    shuffled_primers = []
    for j in range(0, len(all_primers), part_size):
        part = all_primers[j:j + part_size]
        random.shuffle(part)
        shuffled_primers.extend(part)

    valid_primers = process_primers(shuffled_primers ,run_id )

    logging.info(f"Total valid primers: {len(valid_primers)}")
    with open(f'logs/output_14_len_primer_final_primers3_run_{run_id}.txt', 'w') as f:
        for primer in valid_primers:
            f.write(f"{primer}\n")

    logging.info("Primer processing completed")

if __name__ == '__main__':
    with ProcessPoolExecutor(max_workers=8) as executor:
        executor.map(run_primer_processing, range(1, 9))
