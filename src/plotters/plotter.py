import re
import pandas as pd
from IPython.display import display

# Initialize an empty list to store log data
log_data = []

# Regular expressions for extracting data
run_re = re.compile(r'run (\d+) is now processing primer (\d+)/(\d+) \(([\d\.]+)%% of the primers\), valid primers: (\d+)')
disqualify_re = re.compile(r'primers disqualify by level')
inter_comp_re = re.compile(r'filtered by max inter comp: (\d+) primers')

# Variables to hold summary statistics
disqualify_by_level = None
filtered_by_max_inter_comp = []

# Read the log file
with open('plot/logsshuffle.txt', 'r') as file:
    for line in file:
        run_match = run_re.search(line)
        if run_match:
            run_id, current_primer, total_primers, percentage, valid_primers = run_match.groups()
            log_data.append({
                'run_id': int(run_id),
                'current_primer': int(current_primer),
                'total_primers': int(total_primers),
                'percentage': float(percentage),
                'valid_primers': int(valid_primers)
            })

        if disqualify_re.search(line):
            disqualify_by_level = line.split('primers disqualify by level')[1].strip()

        inter_comp_match = inter_comp_re.search(line)
        if inter_comp_match:
            filtered_by_max_inter_comp.append(int(inter_comp_match.group(1)))

# Convert log data to a DataFrame
df = pd.DataFrame(log_data)

# Perform analysis
average_valid_primers_per_run = df.groupby('run_id')['valid_primers'].mean()
total_primers_disqualified_by_level = [int(val) for val in disqualify_by_level.split(', ')] if disqualify_by_level else []
total_filtered_by_max_inter_comp = sum(filtered_by_max_inter_comp)

# Print summary statistics
print("Average valid primers per run:")
print(average_valid_primers_per_run)
print("\nTotal primers disqualified by level:")
print(total_primers_disqualified_by_level)
print("\nTotal primers filtered by max inter comp:")
print(total_filtered_by_max_inter_comp)

# displaying the DataFrame
display(df)
