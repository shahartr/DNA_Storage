import itertools
import random
import sys
primer_library = []
primer_to_hamming_distance_primer = {}
all_primers = []
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

characters = ['A', 'C', 'G', 'T']

# return the complement of a strand of DNA
def complement_strand(strand):
    complement = ""
    for base in strand:
        complement += complement_map[base]
    return complement

# calculate the percent of CG in a strand of DNA
def cg_percent(strand):
    return (strand.count('C') + strand.count('G')) / len(strand)

# determine the maximum homopolymers in a strand
def max_homopolymer(strand):
    bp = strand[0]
    max = 1
    cur = 1
    for i in range(1, len(strand)):
        if bp == strand[i]:
            cur += 1
        else:
            bp = strand[i]
            if cur > max:
                max = cur
            cur = 1
    if cur > max:
        max = cur
    return max

# return true if there is a longer than length-base pair complement between the two strands
def contains_complement(strand1, strand2, length):
    strand1_comp = complement_strand(strand1)
    for i in range(0, len(strand1_comp) - length):
        if strand1_comp[i:(i + length)] in strand2:
            return True
    return False

# return hamming distance between two strings, assuming the two strands are same length
def hamming_distance(strand1, strand2):
    dist = 0
    for i in range(0, len(strand1) - 1):
        if strand1[i] != strand2[i]:
            dist += 1
    return dist

def get_all_primers():
    all_primers = [''.join(p) for p in itertools.product(characters, repeat=PRIMER_BPS)]
    return all_primers


def run():

    all_primers = get_all_primers()

    i = 0
    while(len(all_primers) > 0):
        if (i % 500000) == 0:
            sys.stdout.write("")

        i+=1

        primer = random(all_primers)
        all_primers.remove(primer)

        # GC Content between 45 and 55%
        percent = cg_percent(primer)
        if (percent < 0.45) or (percent > 0.55):
            continue

        # no homopolymers greater than length 2
        if max_homopolymer(primer) > MAX_HP:
            continue

        # no self complementing greater than 4
        if contains_complement(primer, primer, MAX_SELF_COMP + 1):
            continue

        valid = True
        for p in primer_library:
            # hamming distance of at least 6 from other primers
            if hamming_distance(p, primer) < MIN_HAM:
                primer_to_hamming_distance_primer[p].append(primer)
                valid = False
                break

            # no interstrand complements of more than 10
            if contains_complement(p, primer, MAX_INTER_COMP + 1):
                valid = False
                break

        if valid:
            primer_library.append(primer)
            primer_to_hamming_distance_primer[primer] = []

    print("")
    print(primer_library)
    print(len(primer_library))


if __name__ == "__main__":
    run()