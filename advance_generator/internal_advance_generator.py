import sys
import datetime
primer_library = []
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


def run():
    print("Start running time:" , datetime.datetime.now())
    primer = get_first_primer()

    for i in range(4 ** PRIMER_BPS):
        if (i % 500000) == 0:
            print("Time: " , datetime.datetime.now() , " after: " , i , " strings")

        primer = next_primer(primer)

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

        primer_library.append(primer)

    print("")
    #print(primer_library)
    print(len(primer_library))

    with open('output_14_len_primer_internal_advance_generator.txt', 'w') as f:
        for item in primer_library:
            f.write("%s\n" % item)


if __name__ == "__main__":
    run()

#after the internal checking we have 32804376 primers.