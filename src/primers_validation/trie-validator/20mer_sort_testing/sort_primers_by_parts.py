import os
import heapq
from tqdm import tqdm
from tempfile import NamedTemporaryFile

PRIMER_BPS = 20
INPUT_FILES = 5
OUTPUT_FILES = 10
CHUNK_SIZE = 200000000  # Estimated based on available memory

# Define the order of DNA bases
base_order = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

def sort_and_save_chunk(primers, chunk_idx):
    primers.sort()
    temp_file = NamedTemporaryFile(delete=False, mode='w+t')
    for primer in primers:
        temp_file.write(f"{primer}\n")
    temp_file.flush()
    temp_file.seek(0)
    return temp_file.name

def merge_files(file_list, output_file):
    min_heap = []
    file_pointers = [open(file, 'r') for file in file_list]

    # Initialize heap
    for idx, fp in enumerate(file_pointers):
        line = fp.readline().strip()
        if line:
            heapq.heappush(min_heap, (line, idx))

    with open(output_file, 'w') as f_out:
        while min_heap:
            smallest, idx = heapq.heappop(min_heap)
            f_out.write(smallest + '\n')
            next_line = file_pointers[idx].readline().strip()
            if next_line:
                heapq.heappush(min_heap, (next_line, idx))

    for fp in file_pointers:
        fp.close()
    for file in file_list:
        os.remove(file)

def main():
    chunk_files = []
    for i in range(INPUT_FILES):
        with open(f'output_{PRIMER_BPS}_len_primer_internal_advance_generator_{i}.txt', 'r') as f:
            while True:
                current_file_primers = [f.readline().strip() for _ in tqdm(range(CHUNK_SIZE), desc=f'Reading chunk from file {i}')]
                if not current_file_primers[0]:
                    break
                print(f"Sorting chunk of primers from file {i}...")
                chunk_file = sort_and_save_chunk(current_file_primers, i)
                chunk_files.append(chunk_file)
                print("Chunk sorted and saved.")

    # Merging sorted chunks into final output files
    chunk_size = len(chunk_files) // OUTPUT_FILES
    for i in range(OUTPUT_FILES):
        output_file = f'output_{PRIMER_BPS}_len_primer_sorted_part_{i}.txt'
        merge_files(chunk_files[i*chunk_size:(i+1)*chunk_size], output_file)
        print(f"Final sorted primers saved in {output_file}")

if __name__ == '__main__':
    main()
