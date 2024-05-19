import sys
import datetime

all_primers = []
primer_library = []
patterns_of_min_ham_len_library_list = []
patterns_complement_of_max_inter_comp_len_library_set = {}
MAX_HP = 2
PRIMER_BPS = 13
MAX_SELF_COMP = 4
MAX_INTER_COMP = 9
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
    #todo: check the option of iteration tolls
    patterns_of_4_len_library_set = {}
    first_pattern = ""
    for i in range(MIN_HAM):
            first_pattern += 'A'

    patterns_of_4_len_library_set[first_pattern] = False
    prev_pattern = first_pattern

    for i in range(4 ** MIN_HAM):
        suffix = ""
        for i in range((len(prev_pattern) - 1), -1, -1):
            if prev_pattern[i] == 'A':
                suffix = 'C' + suffix
                break
            elif prev_pattern[i] == 'C':
                suffix = 'G' + suffix
                break
            elif prev_pattern[i] == 'G':
                suffix = 'T' + suffix
                break
            else:
                suffix = 'A' + suffix
        current_pattern = prev_pattern[0:len(prev_pattern) - len(suffix)] + suffix
        patterns_of_4_len_library_set[current_pattern] = False
        prev_pattern = current_pattern

    for i  in range(PRIMER_BPS - MIN_HAM + 1):
        patterns_of_min_ham_len_library_list.append(patterns_of_4_len_library_set.copy())

def check_hamming_distance(primer):
    for i in range(PRIMER_BPS - MIN_HAM + 1):
        pattern = primer[i: i + MIN_HAM]
        if(patterns_of_min_ham_len_library_list[i][pattern]):
            return False

    return True

def update_patterns_of_min_ham_len_library(primer):
    for i in range(PRIMER_BPS - MIN_HAM + 1):
        pattern = primer[i: i + MIN_HAM]
        patterns_of_min_ham_len_library_list[i][pattern] = True

def get_all_patterns_of_max_inter_comp():
    patterns_complement_of_7_len_library_set = {}
    first_pattern = ""
    for i in range(MAX_INTER_COMP):
        first_pattern += 'A'

    patterns_complement_of_7_len_library_set[first_pattern] = False
    prev_pattern = first_pattern

    for i in range(4 ** MAX_INTER_COMP):
        suffix = ""
        for i in range((len(prev_pattern) - 1), -1, -1):
            if prev_pattern[i] == 'A':
                suffix = 'C' + suffix
                break
            elif prev_pattern[i] == 'C':
                suffix = 'G' + suffix
                break
            elif prev_pattern[i] == 'G':
                suffix = 'T' + suffix
                break
            else:
                suffix = 'A' + suffix
        current_pattern = prev_pattern[0:len(prev_pattern) - len(suffix)] + suffix
        patterns_complement_of_7_len_library_set[current_pattern] = False
        prev_pattern = current_pattern


    patterns_complement_of_max_inter_comp_len_library_set = patterns_complement_of_7_len_library_set.copy()

def check_if_contains_complement_patterns_in_primers_library(primer):
    complement_primer = complement_strand(primer)

    for i in range(PRIMER_BPS - MAX_INTER_COMP + 1):
        complement_pattern = complement_primer[i: i + MAX_INTER_COMP]
        if(patterns_complement_of_max_inter_comp_len_library_list[i][complement_pattern]):
            return False

    return True

def update_patterns_complement_of_max_inter_complement_len_library(primer):
        for i in range(PRIMER_BPS - MAX_INTER_COMP + 1):
            pattern = primer[i: i + MAX_INTER_COMP]
            patterns_of_min_ham_len_library_list[i][pattern] = True


def run():
    with open('output_13_len_primer_internal_advance_generator.txt', 'r') as f:
        all_primers = f.readlines()

    # Remove newline characters from each item
    all_primers = [item.strip() for item in all_primers]

    get_all_patterns_of_min_ham()

    i = 0
    for primer in all_primers:
        if (i % 250000) == 0:
            print("Time: ", datetime.datetime.now(), " after: ", i, " strings from: ", len(all_primers), "\nsum primers: ", len(primer_library))

        valid = check_hamming_distance(primer)

        if(valid):
            get_all_patterns_of_max_inter_comp()
            valid = check_if_contains_complement_patterns_in_primers_library(primer)

            if(valid):
                update_patterns_of_min_ham_len_library(primer)
                update_patterns_complement_of_max_inter_complement_len_library(primer)
                primer_library.append(primer)
        i = i + 1


    print("")
    #print(primer_library)
    print(len(primer_library))


if __name__ == "__main__":
    run()

# after all: 179 primers for PRIMER_BPS = 10, MAX_INTER_COMP = 7, MIN_HAM = 4.
# after all: 724 primers for PRIMER_BPS = 11, MAX_INTER_COMP = 8, MIN_HAM = 5.
# after all: 682 primers for PRIMER_BPS = 12, MAX_INTER_COMP = 9, MIN_HAM = 5.


# Time:  2024-05-19 12:32:27.829146  after:  6000000  strings from:  17171840
# sum primers:  1140