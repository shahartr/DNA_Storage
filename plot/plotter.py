import re
import matplotlib.pyplot as plt
import pandas as pd
from datetime import datetime

# Naive logs from file output_14_len_primer_final_primers.txt
naive_logs = open('output_14_len_primer_final_primers.txt').read()
# Trie logs from file output_14_len_primer_final_primers2.txt
trie_logs = open('output_14_len_primer_final_primers2.txt').read()

# Function to parse naive logs
def parse_naive_logs(log_str):
    log_entries = log_str.split('Time: ')[1:]
    log_data = []
    for entry in log_entries:
        time_match = re.search(r'\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}\.\d{6}', entry)
        after_match = re.search(r'after:  (\d+)', entry)
        primers_match = re.search(r'sum primers:  (\d+)', entry)
        if time_match and after_match and primers_match:
            timestamp = datetime.strptime(time_match.group(), '%Y-%m-%d %H:%M:%S.%f')
            after = int(after_match.group(1))
            primers = int(primers_match.group(1))
            log_data.append((timestamp, after, primers))
    return pd.DataFrame(log_data, columns=['Time', 'After', 'Sum Primers'])

# Function to parse trie logs
def parse_trie_logs(log_str):
    log_entries = log_str.split('\n')
    log_data = []
    current_timestamp = None
    current_after = None
    for entry in log_entries:
        time_match = re.search(r'\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2},\d{3}', entry)
        if 'Processing primer' in entry:
            after_match = re.search(r'Processing primer (\d+)/32804376', entry)
            if time_match and after_match:
                current_timestamp = datetime.strptime(time_match.group(), '%Y-%m-%d %H:%M:%S,%f')
                current_after = int(after_match.group(1))
        if 'Primer set size' in entry:
            primers_match = re.search(r'Primer set size: (\d+)', entry)
            if current_timestamp and current_after and primers_match:
                primers = int(primers_match.group(1))
                log_data.append((current_timestamp, current_after, primers))
    return pd.DataFrame(log_data, columns=['Time', 'After', 'Sum Primers'])

# Parse both logs
naive_df = parse_naive_logs(naive_logs)
trie_df = parse_trie_logs(trie_logs)

# Normalize the timestamps to start from zero
def normalize_time(df):
    start_time = df['Time'].iloc[0]
    df['Time'] = (df['Time'] - start_time).dt.total_seconds()
    return df

naive_df = normalize_time(naive_df)
trie_df = normalize_time(trie_df)

# Calculate processing time per batch
def calculate_batch_time(df):
    df['Batch Time'] = df['Time'].diff().fillna(0)
    return df

naive_df = calculate_batch_time(naive_df)
trie_df = calculate_batch_time(trie_df)

# Plot Strings Processed Over Time
plt.figure(figsize=(14, 7))
plt.plot(naive_df['Time'], naive_df['After'], label='Naive', marker='o')
plt.plot(trie_df['Time'], trie_df['After'], label='Trie', marker='o')
plt.xlabel('Time (seconds)')
plt.ylabel('Strings Processed')
plt.title('Strings Processed Over Time')
plt.legend()
plt.grid(True)
plt.show()

# Plot Sum Primers Over Strings Processed
plt.figure(figsize=(14, 7))
plt.plot(naive_df['After'], naive_df['Sum Primers'], label='Naive', marker='o')
plt.plot(trie_df['After'], trie_df['Sum Primers'], label='Trie', marker='o')
plt.xlabel('Strings Processed')
plt.ylabel('Sum Primers')
plt.title('Sum Primers Over Strings Processed')
plt.legend()
plt.grid(True)
plt.show()

# Plot Processing Time per Batch
plt.figure(figsize=(14, 7))
plt.bar(naive_df['After'], naive_df['Batch Time'], width=250000, alpha=0.5, label='Naive')
plt.bar(trie_df['After'], trie_df['Batch Time'], width=250000, alpha=0.5, label='Trie')
plt.xlabel('Strings Processed')
plt.ylabel('Processing Time (seconds)')
plt.title('Processing Time per Batch')
plt.legend()
plt.grid(True)
plt.show()
