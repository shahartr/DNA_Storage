import sys
primer_library = []
MAX_HP = 2
PRIMER_BPS = 10 # 14
MAX_SELF_COMP = 4
MAX_INTER_COMP = 7 #10
MIN_HAM = 4 #6

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

first_primer = ""
for i in range(PRIMER_BPS):
    first_primer += 'A'

primer = first_primer
for i in range(4 ** PRIMER_BPS):
    if (i % 500000) == 0:
        sys.stdout.write(".")


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


if __name__ == "__main__":
    print("")
    print(primer_library)
    print(len(primer_library))
