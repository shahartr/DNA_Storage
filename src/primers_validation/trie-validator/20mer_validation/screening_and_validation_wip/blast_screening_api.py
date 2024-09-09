import os
import logging
from tqdm import tqdm
import subprocess
# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def load_primers(filename):
    with open(filename, 'r') as f:
        primers = [line.strip() for line in f]
    return primers

def save_primers_to_fasta(primers, filename):
    with open(filename, 'w') as f:
        for i, primer in enumerate(primers):
            f.write(f">primer_{i}\n{primer}\n")

def parse_blast_results(filename, evalue_threshold=0.01):
    valid_primers = set()
    with open(filename) as f:
        for line in f:
            cols = line.strip().split("\t")
            primer_id = cols[0]
            evalue = float(cols[10])
            if evalue > evalue_threshold:  # E-value threshold for filtering
                valid_primers.add(primer_id)
    return valid_primers

def run_blast_screening(primer_file, output_file):
    # Run BLAST using NCBI's remote service
    command = f"blastn -query {primer_file} -db nt -remote -out {output_file} -outfmt 6"
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    if process.returncode != 0:
        raise Exception(f"BLAST failed: {stderr.decode('utf-8')}")
    return parse_blast_results(output_file)

def main():
    input_file = "data/final_after_tm_2nd_structure_20_length.txt"
    primer_file = "../../../../../data/primers_for_blast.fasta"
    output_file = "blast_results.out"

    logging.info("Loading primers from file...")
    primers = load_primers(input_file)

    logging.info("Saving primers to FASTA file for BLAST screening...")
    save_primers_to_fasta(primers, primer_file)

    logging.info("Running BLAST screening...")
    valid_primers = set()
    batch_size = 1000  # Number of primers to process in each batch
    total_batches = (len(primers) + batch_size - 1) // batch_size

    for i in tqdm(range(total_batches), desc="BLAST Screening Progress"):
        start_idx = i * batch_size
        end_idx = min(start_idx + batch_size, len(primers))
        batch_file = f"batch_{i}.fasta"
        batch_output_file = f"batch_{i}_results.out"

        save_primers_to_fasta(primers[start_idx:end_idx], batch_file)
        valid_primers_batch = run_blast_screening(batch_file, batch_output_file)
        valid_primers.update(valid_primers_batch)

        os.remove(batch_file)
        os.remove(batch_output_file)

    valid_primers_list = [primers[int(primer_id.split('_')[1])] for primer_id in valid_primers]

    logging.info(f"Total valid primers after BLAST screening: {len(valid_primers_list)}")
    logging.info(f"Number of primers filtered by BLAST screening: {len(primers) - len(valid_primers_list)}")

    with open('data/final_valid_primers.txt', 'w') as f:
        for primer in valid_primers_list:
            f.write(f"{primer}\n")

    logging.info("Primer BLAST screening completed and results saved.")

if __name__ == "__main__":
    main()
