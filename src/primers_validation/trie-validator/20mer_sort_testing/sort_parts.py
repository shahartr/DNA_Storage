import heapq
from tqdm import tqdm

PRIMER_BPS = 20
INPUT_FILES = 4
OUTPUT_FILES = 4
CHUNK_SIZE = 200_000_000

# Define the order of DNA bases (used for lexicographical order)
base_order = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

def save_chunk(output_file, chunk):
    with open(output_file, 'w') as f_out:
        for primer in chunk:
            f_out.write(primer + '\n')

def merge_files(file_list, output_files, chunk_size):
    print(f"Merging {len(file_list)} files into {output_files}")

    # Open all files
    file_pointers = [open(file, 'r') for file in file_list]
    min_heap = []

    # Initialize heap with the first line from each file
    for idx, fp in enumerate(file_pointers):
        line = fp.readline().strip()
        if line:
            heapq.heappush(min_heap, (line, idx))

    chunk = []
    file_count = 0
    while min_heap:
        smallest, idx = heapq.heappop(min_heap)
        chunk.append(smallest)

        if len(chunk) >= chunk_size:
            output_file = output_files[file_count]
            print(f"Saving sorted chunk to {output_file}")
            save_chunk(output_file, chunk)
            chunk = []
            file_count += 1

        next_line = file_pointers[idx].readline().strip()
        if next_line:
            heapq.heappush(min_heap, (next_line, idx))

    # Save any remaining items in the chunk
    if chunk:
        output_file = output_files[file_count]
        print(f"Saving final sorted chunk to {output_file}")
        save_chunk(output_file, chunk)

    for fp in file_pointers:
        fp.close()
    print("Merging complete.")

def main():
    # Input sorted files
    input_files = [f'output_{PRIMER_BPS}_len_primer_sorted_{i}.txt' for i in range(INPUT_FILES)]

    # Output files for final sorted chunks
    output_files = [f'output_{PRIMER_BPS}_len_primer_final_sorted_part_{i}.txt' for i in range(OUTPUT_FILES)]

    merge_files(input_files, output_files, CHUNK_SIZE)

if __name__ == '__main__':
    main()
