import logging
import sys
import datetime
import matplotlib.pyplot as plt
primer_library = []
count_primers_list = []
time_list = []
MAX_HP = 4
PRIMER_BPS = 15
MAX_SELF_COMP = 4
MAX_INTER_COMP = 10
MIN_HAM = 6
MIN_GC = 0.4
MAX_GC = 0.6

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


def save_primers(primer_library):
    #save primer into 4 files each file contains 1/4 of the primers
    for i in range(4):
        with open(f'../mis-validator/output_{PRIMER_BPS}_len_primer_internal_advance_generator_{i}.txt', 'w') as f:
            for item in primer_library[i::4]:
                f.write("%s\n" % item)



def run():
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    print("Start running time:" , datetime.datetime.now())
    primer = get_first_primer()

    start = datetime.datetime.now()
    length= 4 ** PRIMER_BPS
    for i in range(length):
        if (i % 5000000) == 0:
            currentTime = datetime.datetime.now()
            count_primers_list.append(i)
            time_list.append(currentTime)
            percent = "{:.2f}".format((i / length) * 100)
            logging.info(f"{currentTime} Primer count: {i}, of {length} ({percent} %) \n primer library size: {len(primer_library)}")
            if i % 25000000 == 0 and i != 0:
                #calc approx time left
                time_diff = currentTime - start
                time_diff = time_diff.total_seconds()
                time_diff = time_diff / i
                time_diff = time_diff * (length - i)
                time_diff = datetime.timedelta(seconds=time_diff)
                logging.info(f"Approx days left: {time_diff.days} days")

        primer = next_primer(primer)

        # GC Content between 45 and 55%
        percent = cg_percent(primer)
        if (percent < MIN_GC) or (percent > MAX_GC):
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

    # plt.figure(figsize=(10, 6))
    plt.plot(count_primers_list, time_list, marker='o', linestyle='-', color='b')

    # Step 3: Format the plot
    plt.title('Counts Primers Over Time')
    plt.xlabel('Time')
    plt.ylabel('Count')
    plt.grid(True)

    # Rotate date labels for better readability
    plt.xticks(rotation=45)

    # Show the plot
    plt.tight_layout()
    plt.show()

    save_primers(primer_library)
    # with open('../mis-validator/output_20_len_primer_internal_advance_generator.txt', 'w') as f:
    #     for item in primer_library:
    #         f.write("%s\n" % item)


if __name__ == "__main__":
    run()
#after the internal checking we have 0 primers of 9 len with the default values of other parameters.s
#after the internal checking we have 189600 primers of 10 len with the default values of other parameters.
#after the internal checking we have 1294544 primers of 11 len with the default values of other parameters.
#after the internal checking we have 2483096 primers of 12 len with the default values of other parameters.
#after the internal checking we have 17171840 primers of 13 len with the default values of other parameters.
#after the internal checking we have 32804376 primers of 14 len with the default values of other parameters.
#after the internal checking we have 228915528 primers of 15 len with the default values of other parameters.




#length 15 (4^15) = 1073741824
#passed filters: 228915528
#should be 712635024
