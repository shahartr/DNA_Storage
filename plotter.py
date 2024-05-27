import pandas as pd
import matplotlib.pyplot as plt

# Read the log file
log_file = 'logs.txt'
with open(log_file, 'r') as file:
    log_data = file.readlines()

timestamps = []
strings_processed = []
primers_sum = []

# Process each pair of lines in the log data
for i in range(0, len(log_data), 2):
    if i + 1 < len(log_data):
        time_part = log_data[i].split('Time:  ')[1].split('  after: ')
        timestamp = time_part[0].strip()
        strings = int(time_part[1].split('  strings from: ')[0].strip())
        primers = int(log_data[i + 1].split('sum primers: ')[1].strip())

        timestamps.append(pd.to_datetime(timestamp))
        strings_processed.append(strings)
        primers_sum.append(primers)

# Creating DataFrame
data = pd.DataFrame({
    'Timestamp': timestamps,
    'StringsProcessed': strings_processed,
    'PrimersSum': primers_sum
})

# Calculating rate of change
data['StringsProcessedRate'] = data['StringsProcessed'].diff()
data['PrimersSumRate'] = data['PrimersSum'].diff()

# Plotting
fig, axs = plt.subplots(4, 1, figsize=(14, 20), sharex=True)

# Plot 1: Number of strings processed over time
axs[0].plot(data['Timestamp'], data['StringsProcessed'], label='Strings Processed', color='blue')
axs[0].set_title('Number of Strings Processed Over Time')
axs[0].set_ylabel('Strings Processed')
axs[0].legend()

# Plot 2: Sum of primers over time
axs[1].plot(data['Timestamp'], data['PrimersSum'], label='Sum of Primers', color='green')
axs[1].set_title('Sum of Primers Over Time')
axs[1].set_ylabel('Sum of Primers')
axs[1].legend()

# Plot 3: Rate of change in number of strings processed
axs[2].plot(data['Timestamp'], data['StringsProcessedRate'], label='Rate of Change in Strings Processed', color='red')
axs[2].set_title('Rate of Change in Number of Strings Processed')
axs[2].set_ylabel('Rate of Change')
axs[2].legend()

# Plot 4: Rate of change in sum of primers
axs[3].plot(data['Timestamp'], data['PrimersSumRate'], label='Rate of Change in Sum of Primers', color='purple')
axs[3].set_title('Rate of Change in Sum of Primers')
axs[3].set_ylabel('Rate of Change')
axs[3].set_xlabel('Time')
axs[3].legend()

plt.tight_layout()
plt.show()
