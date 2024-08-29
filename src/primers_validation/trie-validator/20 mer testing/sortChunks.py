import heapq
import os
from tqdm import tqdm

PRIMER_BPS = 20
CHUNKS_PER_FILE = 197
NUM_FILES = 5
FINAL_OUTPUT_FILES = 5
CHUNK_SIZE = 200_000_000  # Adjust the chunk size if necessary
FINAL_FOLDER = 'final_output'  # Change this to your desired final output directory

def save_chunk(output_file, chunk):
    with open(output_file, 'a') as f_out:
        for primer in chunk:
            f_out.write(primer + '\n')

def merge_all_chunks(chunk_files, output_files, chunk_size):
    print(f"Merging {len(chunk_files)} chunk files into {len(output_files)} output files")

    # Open all chunk files
    file_pointers = [open(file, 'r') for file in chunk_files]
    min_heap = []

    # Initialize heap with the first line from each chunk file
    for idx, fp in enumerate(file_pointers):
        line = fp.readline().strip()
        if line:
            heapq.heappush(min_heap, (line, idx))

    chunk = []
    file_count = 0
    total_chunks = len(chunk_files)

    with tqdm(total=total_chunks, desc="Merging chunks") as pbar:
        while min_heap:
            smallest, idx = heapq.heappop(min_heap)
            chunk.append(smallest)

            if len(chunk) >= chunk_size:
                output_file = output_files[file_count % FINAL_OUTPUT_FILES]
                save_chunk(output_file, chunk)
                chunk = []
                file_count += 1

            next_line = file_pointers[idx].readline().strip()
            if next_line:
                heapq.heappush(min_heap, (next_line, idx))

            pbar.update(1)

    # Save any remaining items in the chunk
    if chunk:
        output_file = output_files[file_count % FINAL_OUTPUT_FILES]
        save_chunk(output_file, chunk)

    for fp in file_pointers:
        fp.close()

    print("Merging complete.")

def main():
    # Collect all chunk files
    chunk_files = []
    for i in range(NUM_FILES):
        for j in range(CHUNKS_PER_FILE):
            chunk_files.append(f'sorted_chunk_{i}.txt_chunk_{j}.txt')

    # Ensure the final output directory exists
    if not os.path.exists(FINAL_FOLDER):
        os.makedirs(FINAL_FOLDER)

    # Final output files
    output_files = [os.path.join(FINAL_FOLDER, f'final_sorted_part_{i}.txt') for i in range(FINAL_OUTPUT_FILES)]

    # Clear the content of final output files if they already exist
    for output_file in output_files:
        open(output_file, 'w').close()

    merge_all_chunks(chunk_files, output_files, CHUNK_SIZE)

if __name__ == '__main__':
    main()
