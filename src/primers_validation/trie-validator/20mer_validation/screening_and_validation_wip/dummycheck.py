#check final list to validate all criteria for primers
# gc content, melting temp, secondary structure, homology with other primers ...
from Bio.Seq import Seq
from tqdm import tqdm

MAX_HP = 2
PRIMER_BPS = 20
MAX_SELF_COMP = 4
MAX_INTER_COMP = 10
MIN_HAM = 6  #30% of PRIMER_BPS
MIN_GC = 0.45
MAX_GC = 0.55
MIN_GC_CLAMPS = 3
GC_CLAPM_LENGTH = 5


def gc_clamp(primer):
    return True
    """Check if the primer has a GC clamp at the 3' end."""
    gc_clamp_count = sum(1 for base in primer[-GC_CLAPM_LENGTH:] if base in "GC")
    return gc_clamp_count >= MIN_GC_CLAMPS


def GC(primer):
    """Calculate the GC content of the primer."""
    return (primer.count("G") + primer.count("C")) / PRIMER_BPS


def gc_content(primer):
    """Check if the primer's GC content is within the specified range."""
    gc_ratio = GC(primer)
    return MIN_GC <= gc_ratio <= MAX_GC


def homopolymer(primer):
    """Check for homopolymers in the primer."""
    for base in "ACGT":
        if base * (MAX_HP + 1) in primer:
            return False
    return True


def self_complementarity(primer):
    """Check for self-complementarity (intra-primer hairpins)."""
    strand_rev_comp = Seq(primer).reverse_complement()
    strand_rev_comp = str(strand_rev_comp)
    # strand_comp = complement_strand(strand)
    #check for self complement
    for i in range(PRIMER_BPS - MAX_SELF_COMP):  #16 times patterns of 5 bp
        if strand_rev_comp[i:i + MAX_SELF_COMP + 1] in primer:
            return False
    return True


def inter_complementarity(primer, final_primers):
    """Check for inter-primer complementarity (dimerization)."""
    reverse_complement = str(Seq(primer).reverse_complement())

    for existing_primer in final_primers:
        for i in range(PRIMER_BPS - MAX_INTER_COMP):
            subseq = existing_primer[i:i + MAX_INTER_COMP + 1]
            if subseq in reverse_complement:
                return False  # Significant complementarity found

    return True  # No significant complementarity found


def hamming_distance(primer, final_primers):
    """Check if the primer has a Hamming distance greater than MIN_HAM with all other primers."""

    def hamming(s1, s2):
        return sum(c1 != c2 for c1, c2 in zip(s1, s2))

    for other_primer in final_primers:
        if hamming(primer, other_primer) < MIN_HAM:
            return False
    return True


def reverse_complement(seq):
    """Generate the reverse complement of a sequence."""
    complement = str.maketrans('ACGT', 'TGCA')
    return seq.translate(complement)[::-1]


def run():
    count_gc_clamp = 0
    count_gc_content = 0
    count_homopolymer = 0
    count_intra_complementarity = 0
    count_inter_complementarity = 0
    count_hamming_distance = 0

    with open('../../../../../data/valid_primers_so_far.txt', 'r') as f:
        primers = [line.strip() for line in f]
    final_primers = []
    for primer in tqdm(primers):
        if gc_clamp(primer):
            if gc_content(primer):
                if homopolymer(primer):
                    if self_complementarity(primer):
                        if inter_complementarity(primer, final_primers):
                            if hamming_distance(primer, final_primers):
                                final_primers.append(primer)
                            else:
                                count_hamming_distance += 1
                        else:
                            count_inter_complementarity += 1
                    else:
                        count_intra_complementarity += 1
                else:
                    count_homopolymer += 1
            else:
                count_gc_content += 1
        else:
            count_gc_clamp += 1
    print(f"Number of valid primers: {len(final_primers)}")
    print(f"Number of primers filtered by GC clamp: {count_gc_clamp}")
    print(f"Number of primers filtered by GC content: {count_gc_content}")
    print(f"Number of primers filtered by homopolymer: {count_homopolymer}")
    print(f"Number of primers filtered by intra-complementarity: {count_intra_complementarity}")
    print(f"Number of primers filtered by inter-complementarity: {count_inter_complementarity}")
    print(f"Number of primers filtered by hamming distance: {count_hamming_distance}")
    with open('final_after_all_filters_20_length.txt', 'w') as f:
        for primer in final_primers:
            f.write(f"{primer}\n")


if __name__ == '__main__':
    run()
