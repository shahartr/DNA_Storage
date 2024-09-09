import logging
import os
import random
import datetime
import time

import tqdm
from Bio.Seq import Seq
from nupack import mfe
import nupack
from primer3 import calc_tm
from utils import primers_trie_tree as trie_utils

# Configuration and Constants
MAX_HP = 2
PRIMER_BPS = 20
MAX_SELF_COMP = 4
MAX_INTER_COMP = 10
MIN_HAM = 6
nupack_model = nupack.Model(material='dna', celsius=37, sodium=0.05, magnesium=0.0015)
MIN_GC = 0.45
MAX_GC = 0.55
MIN_GC_CLAMPS = 3
NUM_PRIMERS = 5_000_000

complement_map = {
    'G': 'C',
    'C': 'G',
    'T': 'A',
    'A': 'T'
}
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


def complement_strand(strand):
    complement = ""
    for base in strand:
        complement += complement_map[base]
    return complement


def cg_percent(strand):
    """ calculate the percent of CG in a strand of DNA"""
    return (strand.count('C') + strand.count('G')) / len(strand)


def max_homo_polymer(strand):
    """ determine the maximum homo-polymers in a strand"""
    if not strand:
        return 0
    current_base = strand[0]
    max_length = current_length = 1
    for i in range(1, len(strand)):
        if strand[i] == current_base:
            current_length += 1
        else:
            current_base = strand[i]
            max_length = max(max_length, current_length)
            current_length = 1

    return max(max_length, current_length)


def reverse_complement_strand(strand):
    """ return the reverse complement of a strand of DNA"""
    seq = Seq(strand)
    return str(seq.reverse_complement())


def contains_self_complement(strand, max_self_comp=MAX_SELF_COMP, primer_length=PRIMER_BPS):
    """ determine if a strand of DNA contains a self-complement
     -return true if there is a longer than 4-base pair complement between the two strands"""
    strand_comp = str(Seq(strand).reverse_complement())
    for i in range(primer_length - max_self_comp):  # 16 times patterns of 5 bp
        if strand_comp[i:i + max_self_comp + 1] in strand:
            return True
    return False


def get_next_random_primer(length=PRIMER_BPS):
    """ generate a random DNA primer of a given length"""
    return ''.join(random.choices('AGCT', k=length))


def save_primers(primer_library):
    # save primer into 4 files each file contains 1/4 of the primers
    with open(f'output_{PRIMER_BPS}_sorted_{NUM_PRIMERS}.txt', 'w') as f:
        for item in primer_library:
            f.write("%s\n" % item)


def sort_primers(primers_library):
    # sort DNA (A, C, G, T) in lexicographical order
    primers_library.sort()

def print_logs2(i, length, start, time_list, count_primers_list, primers_library):
    currentTime = datetime.datetime.now()
    count_primers_list.append(i)
    time_list.append(currentTime)
    percent = "{:.2f}".format((i / length) * 100)
    logging.info(
        f"{currentTime} Primer count: {i}, of {length} ({percent} %) \n primer library size: {len(primers_library)}")
    if i % 2500000 == 0 and i != 0:
        # calc approx time left
        time_diff = currentTime - start
        time_diff = time_diff.total_seconds()
        time_diff = time_diff / i
        time_diff = time_diff * (length - i)
        time_diff = datetime.timedelta(seconds=time_diff)
        minutes = time_diff.seconds // 60
        logging.info(f"Approx {time_diff.seconds} seconds left")

def run():
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    print("Start running time:", datetime.datetime.now())
    start = datetime.datetime.now()
    length = NUM_PRIMERS
    count_gc_disqualified = 0
    count_hp_disqualified = 0
    count_self_comp_disqualified = 0
    count_gc_clamp_disqualified = 0
    min_gc = MIN_GC
    max_gc = MAX_GC
    max_hp = MAX_HP
    max_self_comp = MAX_SELF_COMP
    primers_library = []
    count_primers_list = []
    time_list = []
    time_gc_content = 0
    time_homopolymer = 0
    time_self_complement = 0
    primer_len = PRIMER_BPS

    for i in tqdm.tqdm(range(length)):
        # if (i % 500000) == 0:
        #     print_logs2(i, length, start, time_list, count_primers_list, primers_library)
        primer = get_next_random_primer(primer_len)
        # Measure time for GC Content
        start_gc = time.time()
        percent = cg_percent(primer)
        time_gc_content += time.time() - start_gc
        if percent < MIN_GC or percent > MAX_GC:
            count_gc_disqualified += 1
            continue

        # Measure time for homo-polymer check
        start_hp = time.time()
        if max_homo_polymer(primer) > MAX_HP:
            time_homopolymer += time.time() - start_hp
            count_hp_disqualified += 1
            continue
        time_homopolymer += time.time() - start_hp

        # Measure time for self-complementarity
        start_self_comp = time.time()
        if contains_self_complement(primer, MAX_SELF_COMP):
            time_self_complement += time.time() - start_self_comp
            count_self_comp_disqualified += 1
            continue
        time_self_complement += time.time() - start_self_comp

        primers_library.append(primer)

    logging.info(f"Total time for GC content check: {time_gc_content:.2f} seconds")
    logging.info(f"Total time for homo-polymer check: {time_homopolymer:.2f} seconds")
    logging.info(f"Total time for self-complement check: {time_self_complement:.2f} seconds")

    logging.info(f"Number of primers disqualified by GC content: {count_gc_disqualified}")
    logging.info(f"Number of primers disqualified by homo-polymers: {count_hp_disqualified}")
    # logging.info(f"Number of primers disqualified by GC clamps: {count_gc_clamp_disqualified}")
    logging.info(f"Number of primers disqualified by self-complementarity: {count_self_comp_disqualified}")
    logging.info(f"Number of primers passed first phase: {len(primers_library)}")
    print("sorting primers...")
    sort_primers(primers_library)
    # save_primers(list(primers_library))
    return primers_library



def print_logs(index, count_primers_filtered_by_max_inter_compliment
               , primer_set, trie, count_primers_disqualify_by_tm,
               count_primers_disqualify_by_secondary_structure,
               count_primers_disqualify_by_hamming_distance):
    leaves = trie.count_leaves()
    nodes = trie.count_tree_nodes()
    # percent = "{:.2f}%".format(((index + 1) / len(all_primers)) * 100)
    # logging.info(f" Processing primer {index + 1}/{len(all_primers)} ({percent}% of the primers), valid primers: {len(primer_set)}")
    logging.info(f" Processing primer {index + 1}  Valid primers: {len(primer_set)}")
    if index % 80000 == 0 and index != 0:
        logging.info(f" tree status: {leaves} leaves and {nodes} nodes")
        logging.info(f" filtered by max inter comp: {count_primers_filtered_by_max_inter_compliment} primers ")
        logging.info(f" filtered by hamming distance: {count_primers_disqualify_by_hamming_distance} primers ")
        logging.info(f" filtered by tm: {count_primers_disqualify_by_tm} primers ")
        logging.info(f" filtered by secondary structure: {count_primers_disqualify_by_secondary_structure} primers ")


def check_if_contains_fixed_length_complement_patterns_in_primers_set(primer,
                                                                      patterns_of_all_primers_in_set):
    """Check if the primer contains a fixed length complement pattern in the set of primers"""
    primer_len = PRIMER_BPS
    inter_comp_len = MAX_INTER_COMP
    rev_complement_primer = Seq(primer).reverse_complement()
    rev_complement_primer = str(rev_complement_primer)
    for i in range(primer_len - inter_comp_len):
        rev_complement_pattern = rev_complement_primer[i:i + inter_comp_len + 1]
        if rev_complement_pattern in patterns_of_all_primers_in_set:
            return False
    return True  # later add all parts of the primer to the set if the primer is valid


def check_melting_temp(primer, min_tm=55, max_tm=60):
    """Check for melting temperature of the primer"""
    tm_pass = calc_tm(primer)
    if tm_pass < min_tm or tm_pass > max_tm:
        return True
    return False


def update_complement_set(primer, patterns_complement_of_max_inter_comp_len_set):
    """Update the set of patterns of max inter complement length"""
    primer_len = PRIMER_BPS
    inter_comp_len = MAX_INTER_COMP
    for i in range(primer_len - inter_comp_len):
        pattern = primer[i: i + inter_comp_len + 1]
        patterns_complement_of_max_inter_comp_len_set.add(pattern)


def check_inter_primer_complementarity(new_primer, final_set, max_complementarity=4):
    """ Check for inter-primer homology with the final set of primers """
    reverse_complement = str(Seq(new_primer).reverse_complement())
    primer_length = PRIMER_BPS
    for existing_primer in final_set:
        for i in range(primer_length - max_complementarity + 1):
            subseq = existing_primer[i:i + max_complementarity]
            if subseq in reverse_complement:
                return False  # Significant complementarity found

    return True  # No significant complementarity ...found later added to the set if valid


def hairpin_homo_dimer_check(primer, model):
    """Check for hairpin and homo-dimer formation using NUPACK mfe , values might be little different from paper"""
    results = nupack.mfe([primer], model=model)
    if results:
        energy = results[0].energy
        if energy < -2:
            return False
        return True
    return False

def process_file( priemrs, primer_set, patterns_complement_of_max_inter_comp_len_set, trie,
                 count_primers_disqualify_by_inter_compliment, count_primers_disqualify_by_tm,
                 count_primers_disqualify_by_secondary_structure,
                 count_primers_disqualify_by_hamming_distance):

    primer_len = PRIMER_BPS

    total_primers = len(priemrs)
    min_ham_distance = 6
    time_hamming_distance = 0
    time_inter_compliment = 0
    time_melting_temp = 0
    time_secondary_structure = 0
    model=nupack_model
    logging.info(f"{total_primers} primers to check")
    for index, primer in tqdm.tqdm(enumerate(priemrs), total=total_primers, desc=f"Processing Primers"):
        if index % 20000 == 0:
            print_logs(index, count_primers_disqualify_by_inter_compliment,
                       primer_set, trie, count_primers_disqualify_by_tm,
                       count_primers_disqualify_by_secondary_structure,
                       count_primers_disqualify_by_hamming_distance)
            #print(f"valid primers: {len(primer_set)}")
        primer = primer.strip()

        # Measure time
        start_tm = time.time()
        if check_melting_temp(primer):
            time_melting_temp += time.time() - start_tm
            start_ham = time.time()
            if trie.is_valid_primer_no_info(primer, min_ham_distance):
                time_hamming_distance += time.time() - start_ham
                start_inter = time.time()
                if check_if_contains_fixed_length_complement_patterns_in_primers_set(primer,
                                                                                     patterns_complement_of_max_inter_comp_len_set):
                    time_inter_compliment += time.time() - start_inter
                    start_secondary = time.time()
                    if hairpin_homo_dimer_check(primer, model):
                        time_secondary_structure += time.time() - start_secondary
                        trie.insert(primer)
                        update_complement_set(primer, patterns_complement_of_max_inter_comp_len_set)
                        primer_set.add(primer)
                    else:
                        count_primers_disqualify_by_secondary_structure += 1
                        time_secondary_structure += time.time() - start_secondary
                else:
                    count_primers_disqualify_by_inter_compliment += 1
                    time_inter_compliment += time.time() - start_inter
            else:
                count_primers_disqualify_by_hamming_distance += 1
                time_hamming_distance += time.time() - start_ham

        else:
            count_primers_disqualify_by_tm += 1
            time_melting_temp += time.time() - start_tm


    logging.info(f"Total time for Hamming distance check: {time_hamming_distance:.2f} seconds")
    logging.info(f"Total time for inter-primer complementarity check: {time_inter_compliment:.2f} seconds")
    logging.info(f"Total time for melting temperature check: {time_melting_temp:.2f} seconds")
    logging.info(f"Total time for secondary structure check: {time_secondary_structure:.2f} seconds")
    logging.info(f"Number of primers filtered by max inter comp: {count_primers_disqualify_by_inter_compliment}")
    logging.info(f"Number of primers disqualify by tm: {count_primers_disqualify_by_tm}")
    logging.info(f"Number of primers disqualify by secondary structure: {count_primers_disqualify_by_secondary_structure}")
    logging.info(f"Number of primers disqualify by hamming distance: {count_primers_disqualify_by_hamming_distance}")


    #show paper relation
    valid_primers_without_tm_and_ss = len(primer_set)+count_primers_disqualify_by_tm+count_primers_disqualify_by_secondary_structure
    logging.info(f"box 1: {valid_primers_without_tm_and_ss}")
    logging.info(f"box 2: {len(primer_set)}")


def run2(primers):
    logging.info("Starting primer processing")
    primer_set = set()
    patterns_complement_of_max_inter_comp_len_set = set()
    trie = trie_utils.Trie()
    count_primers_filtered_by_max_inter_compliment = 0
    count_primers_disqualify_by_tm = 0
    count_primers_disqualify_by_secondary_structure = 0
    count_primers_disqualify_by_hamming_distance = 0

    process_file(primers, primer_set, patterns_complement_of_max_inter_comp_len_set, trie,
                 count_primers_filtered_by_max_inter_compliment, count_primers_disqualify_by_tm,
                 count_primers_disqualify_by_secondary_structure, count_primers_disqualify_by_hamming_distance)

    logging.info(f"Total valid primers: {len(primer_set)}")
    with open('final_20_length25.txt', 'w') as f:
        for primer in primer_set:
            f.write(f"{primer}\n")

    logging.info("Primer processing completed")


def load_primers():
    with open('output_20_sorted_2500000.txt', 'r') as f:
        primers = f.readlines()
    return primers


if __name__ == "__main__":
    primers = load_primers()
    run2(primers)

