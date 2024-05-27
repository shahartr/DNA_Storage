import itertools
import datetime

all_primers = []
primer_set = set()
patterns_of_min_ham_len_set = set()
patterns_complement_of_max_inter_comp_len_set = set()
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
    max_count = 1
    cur_count = 1
    for i in range(1, len(strand)):
        if bp == strand[i]:
            cur_count += 1
        else:
            bp = strand[i]
            if cur_count > max_count:
                max_count = cur_count
            cur_count = 1
    if cur_count > max_count:
        max_count = cur_count
    return max_count

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
    for i in range(0, len(strand1)):
        if strand1[i] != strand2[i]:
            dist += 1
    return dist

# sequentially count through possible strands of DNA, using the lexigraphic order of the nucleotides
def next_primer(strand):
    suffix = ""
    for i in range((len(strand) - 1), -1, -1):
        if strand[i] == 'A':
            suffix = 'C' + suffix
            break
        elif strand[i] == 'C':
            suffix = 'G' + suffix
            break
        elif strand[i] == 'G':
            suffix = 'T' + suffix
            break
        else:
            suffix = 'A' + suffix
    return strand[0:len(strand) - len(suffix)] + suffix

def get_first_primer():
    first_primer = ""
    for i in range(PRIMER_BPS):
        first_primer += 'A'

    return first_primer

def get_all_patterns_of_min_ham():
    # Use itertools.product to generate all combinations of 'A', 'C', 'G', 'T' of length MIN_HAM
    bases = ['A', 'C', 'G', 'T']
    patterns_of_4_len_set = {''.join(pattern): False for pattern in itertools.product(bases, repeat=MIN_HAM)}

    for i in range(PRIMER_BPS - MIN_HAM + 1):
        patterns_of_min_ham_len_set.update(patterns_of_4_len_set) # not used

def check_hamming_distance(primer):
    for existing_primer in primer_set:
        if hamming_distance(primer, existing_primer) < MIN_HAM:
            return False
    return True

def update_patterns_of_min_ham_len_set(primer):
    for i in range(PRIMER_BPS - MIN_HAM + 1):
        pattern = primer[i: i + MIN_HAM]
        patterns_of_min_ham_len_set.add(pattern)

def get_all_patterns_of_max_inter_comp():
    bases = ['A', 'C', 'G', 'T']
    patterns_of_MIC_set = {''.join(pattern): False for pattern in itertools.product(bases, repeat=MAX_INTER_COMP)}
    patterns_complement_of_max_inter_comp_len_set.update(patterns_of_MIC_set)

def check_if_contains_complement_patterns_in_primers_set(primer):
    complement_primer = complement_strand(primer)

    for i in range(PRIMER_BPS - MAX_INTER_COMP + 1):
        complement_pattern = complement_primer[i: i + MAX_INTER_COMP]
        if complement_pattern in patterns_complement_of_max_inter_comp_len_set:
            return False
        else:
            patterns_complement_of_max_inter_comp_len_set.add(complement_pattern)

    return True

def update_patterns_complement_of_max_inter_complement_len_set(primer):
    for i in range(PRIMER_BPS - MAX_INTER_COMP + 1):
        pattern = primer[i: i + MAX_INTER_COMP]
        patterns_of_min_ham_len_set.add(pattern)


def run():
    global primer_set
    with open('output_14_len_primer_internal_advance_generator.txt', 'r') as f:
        all_primers = f.readlines()

    # Remove newline characters from each item
    all_primers = [item.strip() for item in all_primers]
   # get_all_patterns_of_max_inter_comp()
   # get_all_patterns_of_min_ham()  # bad!

    i = 0
    for primer in all_primers:
        if (i % 250000) == 0:
            print("Time: ", datetime.datetime.now(), " after: ", i, " strings from: ", len(all_primers), "\nsum primers: ", len(primer_set))

        valid = check_hamming_distance(primer)

        if(valid):

            valid = check_if_contains_complement_patterns_in_primers_set(primer)

            if(valid):
                #update_patterns_of_min_ham_len_set(primer)
                #update_patterns_complement_of_max_inter_complement_len_set(primer)
                #just add all of his complimentry to the set
                primer_set.add(primer)
        i = i + 1

    print(len(primer_set))


if __name__ == "__main__":
    run()


# after all: 179 primers for PRIMER_BPS = 10, MAX_INTER_COMP = 7, MIN_HAM = 4.
# after all: 724 primers for PRIMER_BPS = 11, MAX_INTER_COMP = 8, MIN_HAM = 5.
# after all: 682 primers for PRIMER_BPS = 12, MAX_INTER_COMP = 9, MIN_HAM = 5.

# after changing some methods to fit i got 546 primers with settings 12,9,6
#and 2443 primer with settings 12,9,5



#default setting 50 %
#after:  16250000  primers from:  32804376
#sum primers:  2626

#73%   after:  25750000  primers from:  32804376
#sum primers:  3610


#100% = 4188 primers!