import logging
import random

from advanced_generator_2.adv_filter_PP_Trie_Blast_seq import process_primers, simulate_multiplex_pcr_validation, \
    filter_by_tm, blast_sequence


#test only the final primer set result from first file
#and now will run the blast sequencing and melt TM for the final primer set
#then do then multiplex PCR \

#load output_14_len_primer_final_primers2.txt file

def filter_TM_and_blast_sequencing(final_primers):
    valid_primers = []
    for primer in final_primers:
        if filter_by_tm(primer):
            if blast_sequence(primer):
                valid_primers.append(primer)
    return valid_primers



if __name__ == '__main__':
    with open('output_14_len_primer_final_primers2.txt', 'r') as f:
        final_primers = f.readlines()

    final_primers = [item.strip() for item in final_primers]
    #for each primer do filter TM and blast sequencing
    valid_primers = filter_TM_and_blast_sequencing(final_primers)
    valid_primers = simulate_multiplex_pcr_validation(valid_primers)


    logging.info(f"Total valid primers: {len(valid_primers)}")
    with open('output_14_len_primer_after_blast.txt', 'w') as f:
        for primer in valid_primers:
            f.write(f"{primer}\n")

    print(f"Total valid primers: {len(valid_primers)}")

    logging.info("Primer processing completed")