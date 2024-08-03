def count_primers_in_file(file_path):
    with open(file_path, 'r') as file:
        return sum(1 for line in file)

file_path = 'output_15_len_primer_internal_advance_generator_0.txt'
primer_count = count_primers_in_file(file_path)

print(f"The file {file_path} contains {primer_count} primers.")
print(f"Expected number of primers: 178158756")
print(f"Check result: {'Match' if primer_count == 178158756 else 'Mismatch'}")
