import sys
from Bio.Seq import Seq
PRIMER_BPS = 20

def check_intra_homology(primer, max_self_complementarity=4):
    """Check for intra-primer homology (self-complementarity)"""
    for i in range(PRIMER_BPS - max_self_complementarity + 1):
        subseq = primer[i:i + max_self_complementarity]
        reverse_complement = str(Seq(subseq).reverse_complement())
        if reverse_complement in primer:
            return False
    return True

def check_inter_homology(primer, primers, max_complementarity=10):
    """Check for inter-primer homology with an existing primer library"""
    reverse_complement = str(Seq(primer).reverse_complement())
    for existing_primer in primers:
        for i in range(PRIMER_BPS - max_complementarity + 1):
            subseq = existing_primer[i:i + max_complementarity]
            if subseq in reverse_complement:
                return False
    return True


def generate_primers(num_primers=5):
    """Generate a set of primers and filter them"""
    with open('final_20_length002.txt', 'r') as f:
        primers = [line.strip() for line in f]
    valid_primers = []

    for primer in primers:
        if check_intra_homology(primer): #and check_inter_homology(primer, valid_primers):
            valid_primers.append(primer)
        else:
            print(f"Primer {primer} failed homology check and was discarded")

    return valid_primers


def run():
    valid_primers = generate_primers()



if __name__ == '__main__':
    run()
# Generate and filter primers

