import datetime
import logging
from primer3 import calc_tm
from tqdm import tqdm
import nupack


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
    with open('data/valid_primers.txt', 'w') as f:
        for primer in primer_set:
            f.write(f"{primer}\n")

def calculate_tm(primer):
    return calc_tm(primer)

def filter_by_tm(primer, min_tm=55, max_tm=60):
    tm = calculate_tm(primer)
    return min_tm <= tm <= max_tm

def check_melting_temp(primer):
    tm_pass = filter_by_tm(primer)
    return tm_pass

def check_secondary_structures(primer):
    # Using NUPACK to check for secondary structures
    model = nupack.Model(material='dna', celsius=37) # 37 is standard
    result = nupack.mfe(strands=[primer], model=model)#mfe is the minimum free energy
    # The free energy is now directly accessible from the result
    min_dG = result[0].energy

#min_dG is the minimum free energy means the stability of the secondary structure
#both hairpin and homodimer are secondary structures
#the threshold is set to -3.0 to filter out the primers
    #
    return min_dG > -3.0

def blast_filters(primer, count_tm, count_structure):
    if filter_by_tm(primer):
        if check_secondary_structures(primer):
            return True, count_tm, count_structure
        else:
            count_structure += 1
    else:
        count_tm += 1
    return False, count_tm, count_structure


    return tm_pass, structure_pass


def complement_strand(strand):
    try:
        return ''.join(complement_map[base] for base in strand)
    except KeyError as e:
        logging.error(f"Invalid character found: {e.args[0]} in primer: {strand}")
        return None



def main():
    logging.info("Starting primer processing")

    final_primer_set = set()
    primers=[]

    with open('data/final_20_length.txt', 'r') as f:
        for line in f:
            primer = line.strip()
            primers.append(primer)

    count_tm=0
    count_structure=0

    for primer in tqdm(primers, desc='Processing primers'):
        if check_melting_temp(primer):
            if check_secondary_structures(primer):
                final_primer_set.add(primer)
            else:
                count_structure += 1
        else:
            count_tm += 1

    logging.info(f"Total valid primers: {len(final_primer_set)}")
    logging.info(f"Number of primers filtered by tm: {count_tm}")
    logging.info(f"Number of primers filtered by secondary structure: {count_structure}")
    with open('data/final_after_tm_2nd_structure_20_length.txt', 'w') as f:
        for primer in final_primer_set:
            f.write(f"{primer}\n")




if __name__ == "__main__":
    main()
