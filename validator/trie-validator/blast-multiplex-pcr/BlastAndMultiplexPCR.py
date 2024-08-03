import logging
import random
import primer3
import subprocess
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Blast import NCBIWWW, NCBIXML
from primer3 import calcHairpin, calcHomodimer
from Bio.Blast.Applications import NcbiblastnCommandline

#todo
#we need to get each primer "blast" score



def calculate_tm(primer):
    return mt.Tm_NN(primer)

def filter_by_tm(primer, min_tm=30, max_tm=65):
    tm = calculate_tm(primer)
    return min_tm <= tm <= max_tm


def local_blast_sequence(sequence):
    pass
    # Save the sequence to a file
    with open("input.fasta", "w") as f:
        f.write(f">{sequence}\n{sequence}\n")

    # blast+ database
    blast_db_path = "/Users/saarco/blastdb/custom_db"

    blastn_cline = f"blastn -query input.fasta -db {blast_db_path} -evalue 0.01 -outfmt 5 -out blast_results.xml"

    # Use subprocess to run the BLAST+ command
    process = subprocess.run(blastn_cline, shell=True, capture_output=True)

    if process.returncode != 0:
        raise RuntimeError(f"BLAST+ returned an error: {process.stderr.decode()}")

    # Parse the BLAST results
    with open("blast_results.xml") as result_handle:
        blast_record = NCBIXML.read(result_handle)
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < 0.01:
                    return False
    return True
def blast_sequence(sequence):
    result_handle = NCBIWWW.qblast("blastn", "nt", sequence)#add the sequence to the blast queue to nt database
    blast_record = NCBIXML.read(result_handle)
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < 0.01:
                return False
    return True
#test only the final primer set result from first file
#and now will run the blast sequencing and melt TM for the final primer set
#then do then multiplex PCR \

#load output_14_len_primer_final_primers2.txt file

def filter_TM_and_blast_sequencing(final_primers):
    valid_primers = []
    for primer in final_primers:
        if filter_by_tm(primer):
            if local_blast_sequence(primer):
                if len(valid_primers) % 20 == 0:
                    logging.info(f"Primer set size: {len(valid_primers)}")
                valid_primers.append(primer)
    return valid_primers

def simulate_multiplex_pcr_validation(primer_set):
    # Placeholder function for actual Multiplex PCR validation simulation
    # In practice, this would involve checking for primer-primer interactions, compatibility, and more
    valid_primers = set()

    for primer in primer_set:
        if check_secondary_structures(primer):
            valid_primers.add(primer)
    return valid_primers

def check_secondary_structures(primer):
    hairpin = calcHairpin(primer)
    homodimer = calcHomodimer(primer)#homodimer is the primer binding to itself
    return hairpin.structure_found == 0 and homodimer.structure_found == 0



if __name__ == '__main__':
    with open('final.txt', 'r') as f:
        final_primers = f.readlines()

    final_primers = [item.strip() for item in final_primers]
    #for each primer do filter TM and blast sequencing
    #valid_primers = filter_TM_and_blast_sequencing(final_primers)
    valid_primers = simulate_multiplex_pcr_validation(final_primers)



    logging.info(f"Total valid primers: {len(valid_primers)}")
    with open('output_14_len_primer_after_blast.txt', 'w') as f:
        for primer in valid_primers:
            f.write(f"{primer}\n")

    print(f"Total valid primers: {len(valid_primers)}")

    logging.info("Primer processing completed")