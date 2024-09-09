import logging
import random
import datetime
import tqdm
from Bio.Seq import Seq

primer_library = []
count_primers_list = []
time_list = []
MAX_HP = 2
PRIMER_BPS = 20
MAX_SELF_COMP = 4
MAX_INTER_COMP = 10
MIN_HAM = 6
MIN_GC = 0.45
MAX_GC = 0.55
MIN_GC_CLAMPS = 3
NUM_PRIMERS = 10_000_000_000_000


complement_map = {
    'G': 'C',
    'C': 'G',
    'T': 'A',
    'A': 'T'
}

# return the complement of a strand of DNA l
def complement_strand(strand):
    complement = ""
    for base in strand:
        complement += complement_map[base]
    return complement

def gc_clamp(primer, clamp_length=5, min_gc_count=3):
    """
    Validates if the primer has a proper GC clamp at the 3' end for PCR.

    Parameters:
    primer (str): The DNA sequence of the primer.
    clamp_length (int): The number of bases at the 3' end to consider. Default is 5.
    min_gc_count (int): The minimum number of G or C bases required in the last `clamp_length` bases. Default is 3.

    Returns:
    bool: True if the primer has a valid GC clamp, False otherwise.
    """
    # Extract the last `clamp_length` bases from the primer
    end_sequence = primer[-clamp_length:]

    # Count the number of G or C bases in the end sequence
    gc_count = sum(base in 'GC' for base in end_sequence)

    # Validate that the number of G or C bases is at least `min_gc_count`
    return gc_count >= min_gc_count
# calculate the percent of CG in a strand of DNA
def cg_percent(strand):
    return (strand.count('C') + strand.count('G')) / len(strand)

# determine the maximum homopolymers in a strand
def max_homopolymer(strand):
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
    seq = Seq(strand)
    return str(seq.reverse_complement())

def contains_self_complement(strand, max_self_comp=4):
    strand_comp = str(Seq(strand).reverse_complement())
    for i in range(PRIMER_BPS - max_self_comp):#16 times patterns of 5 bp
        if strand_comp[i:i + max_self_comp + 1] in strand:
            return True
    return False
    # return true if there is a longer than 4-base pair complement between the two strands


# sequentially count through possible strands of DNA, using the lexicographic order of the nucleotides
def next_primer(length=PRIMER_BPS):
    return ''.join(random.choices('AGCT', k=length))

def save_primers(primer_library):
    #save primer into 4 files each file contains 1/4 of the primers
    with open(f'output_{PRIMER_BPS}_sorted_V.txt', 'w') as f:
        for item in primer_library:
            f.write("%s\n" % item)

def sort(primer_library):
    #sort DNA (A, C, G, T) in lexicographical order
    primer_library.sort()

def run():
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    print("Start running time:", datetime.datetime.now())
    start = datetime.datetime.now()
    length = NUM_PRIMERS
    count_gc_disqualified = 0
    count_hp_disqualified = 0
    count_self_comp_disqualified = 0
    count_gc_clamp_disqualified = 0
    #tqdm
    for i in tqdm.tqdm(range(length)):
        if (i % 500000) == 0:
            currentTime = datetime.datetime.now()
            count_primers_list.append(i)
            time_list.append(currentTime)
            percent = "{:.2f}".format((i / length) * 100)
            logging.info(f"{currentTime} Primer count: {i}, of {length} ({percent} %) \n primer library size: {len(primer_library)}")
            if i % 2500000 == 0 and i != 0:
                #calc approx time left
                time_diff = currentTime - start
                time_diff = time_diff.total_seconds()
                time_diff = time_diff / i
                time_diff = time_diff * (length - i)
                time_diff = datetime.timedelta(seconds=time_diff)
                minutes = time_diff.seconds // 60
                logging.info(f"Approx {time_diff.seconds} seconds left")

        primer = next_primer()

        # GC Content between 45 and 55%
        percent = cg_percent(primer)
        if (percent < MIN_GC) or (percent > MAX_GC):
            count_gc_disqualified += 1
            continue
        # if not gc_clamp(primer):
        #     count_gc_clamp_disqualified += 1
        #     continue

        # no homopolymers greater than length 2
        if max_homopolymer(primer) > MAX_HP:
            count_hp_disqualified += 1
            continue

        # no self complementing greater than 4
        if contains_self_complement(primer, MAX_SELF_COMP):
            count_self_comp_disqualified += 1
            continue

        primer_library.append(primer)

    logging.info(f"Number of primers disqualified by GC content: {count_gc_disqualified}")
    logging.info(f"Number of primers disqualified by homopolymers: {count_hp_disqualified}")
    logging.info(f"Number of primers disqualified by GC clamps: {count_gc_clamp_disqualified}")
    logging.info(f"Number of primers disqualified by self-complementarity: {count_self_comp_disqualified}")
    logging.info(f"Number of primers generated: {len(primer_library)}")
    print("start sorting...")
    sort(primer_library)

    save_primers(list(primer_library))



if __name__ == "__main__":
    run()

