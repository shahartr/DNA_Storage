import sys

all_primers = []
primer_library = []
patterns_of_6_len_library = {}
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

def get_all_patterns_of_6():
    first_pattern = ""
    for i in range(MIN_HAM):
            first_pattern += 'A'

    patterns_of_6_len_library[first_pattern] = False
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
        prev_pattern = prev_pattern[0:len(prev_pattern) - len(suffix)] + suffix
        patterns_of_6_len_library[prev_pattern] =  True

def run():
    with open('output_14_len_primer_internal_advance_generator.txt', 'r') as f:
        all_primers = f.readlines()

    # Remove newline characters from each item
    all_primers = [item.strip() for item in all_primers]

    get_all_patterns_of_6()

    for primer in all_primers:
        valid = True
        for p in primer_library:
            # hamming distance of at least 6 from other primers
            if hamming_distance(p, primer) < MIN_HAM:
                valid = False
                break

            # no interstrand complements of more than 10
            if contains_complement(p, primer, MAX_INTER_COMP + 1):
                valid = False
                break

        if valid:
            primer_library.append(primer)


    print("")
    #print(primer_library)
    print(len(primer_library))


if __name__ == "__main__":
    run()