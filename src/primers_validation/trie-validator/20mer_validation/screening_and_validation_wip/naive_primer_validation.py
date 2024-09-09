import logging
from tqdm import tqdm


# Constants
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

def complement_strand(strand):
    return ''.join(complement_map[base] for base in strand)

def cg_percent(strand):
    return (strand.count('C') + strand.count('G')) / len(strand)

def max_homopolymer(strand):
    max_len = cur_len = 1
    for i in range(1, len(strand)):
        if strand[i] == strand[i - 1]:
            cur_len += 1
            if cur_len > max_len:
                max_len = cur_len
        else:
            cur_len = 1
    return max_len

def contains_complement(strand1, strand2, length):
    strand1_comp = complement_strand(strand1)
    for i in range(len(strand1_comp) - length + 1):
        if strand1_comp[i:i + length] in strand2:
            return True
    return False

def hamming_distance(strand1, strand2):
    return sum(b1 != b2 for b1, b2 in zip(strand1, strand2))



def main():
    primer_library = []
    with open('../../../../../data/1millionFiles/final_20_length.txt', 'r') as f:
        for line in f:
            primer = line.strip()
            primer_library.append(primer)
    count_filter_by_gc = 0
    count_filter_by_hp = 0
    count_filter_by_self_comp = 0
    count_filter_by_inter_comp = 0
    count_filter_by_ham = 0
    final_primers = []

    #for each primer in the library (with tqdm)
    for primer in tqdm(primer_library, desc='Processing primers'):
        # GC Content check
        if not (0.45 <= cg_percent(primer) <= 0.55):
            count_filter_by_gc += 1
            continue
        # Max homopolymer check
        if max_homopolymer(primer) > MAX_HP:
            count_filter_by_hp += 1
            continue

        # Self-complementarity check
        if contains_complement(primer, primer, MAX_SELF_COMP):
            count_filter_by_self_comp += 1
            continue

        # Check against existing primers
        valid = True
        for p in final_primers:
            if hamming_distance(p, primer) < MIN_HAM:
                count_filter_by_ham += 1
                valid = False
                break
            if contains_complement(p, primer, MAX_INTER_COMP):
                count_filter_by_inter_comp += 1
                valid = False
                break

        if valid:
            final_primers.append(primer)


    logging.info(f"Number of primers: {len(final_primers)}")
    logging.info(f"Number of primers filtered by GC content: {count_filter_by_gc}")
    logging.info(f"Number of primers filtered by max homopolymer: {count_filter_by_hp}")
    logging.info(f"Number of primers filtered by self-complementarity: {count_filter_by_self_comp}")
    logging.info(f"Number of primers filtered by inter-complementarity: {count_filter_by_inter_comp}")
    logging.info(f"Number of primers filtered by Hamming distance: {count_filter_by_ham}")



    with open('data/final_primers333.txt', 'w') as f:
        for primer in primer_library:
            f.write(f"{primer}\n")

if __name__ == "__main__":
    main()
