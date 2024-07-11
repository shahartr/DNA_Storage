import re
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

def parse_log_file(log_file_path):
    # Define the regex pattern
    pattern = re.compile(
        r'(?P<date>\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2},\d{3}) - INFO - Run (?P<run>\d+): Processing primer (?P<primer>\d+)/(?P<total>\d+) \((?P<percent>\d+\.\d{2})%% of the primers\), valid primers: (?P<valid>\d+)'
    )

    runs = {}

    # Read the log file and extract relevant data
    with open(log_file_path, 'r') as file:
        for line in file:
            match = pattern.match(line)
            if match:
                data = match.groupdict()
                run = int(data['run'])
                primer = int(data['primer'])
                total = int(data['total'])
                percent = float(data['percent'])
                valid = int(data['valid'])

                if run not in runs:
                    runs[run] = {
                        'primers': [],
                        'valid': []
                    }

                runs[run]['primers'].append(primer)
                runs[run]['valid'].append(valid)

    return runs

def plot_runs(runs):
    sns.set(style="whitegrid")
    plt.figure(figsize=(14, 10))

    for run, data in runs.items():
        primers = data['primers']
        valid = data['valid']

        # Create x-axis values as percentage of progress
        x = [100 * p / primers[-1] for p in primers]

        sns.lineplot(x=x, y=valid, label=f'Run {run}')

    plt.xlabel('Percentage of Primers Processed', fontsize=14)
    plt.ylabel('Valid Primers', fontsize=14)
    plt.title('Valid Primers vs. Percentage of Primers Processed', fontsize=16)
    plt.legend(title='Run Number')
    plt.grid(True)
    plt.show()

def plot_total_valid_primers(runs):
    run_numbers = []
    total_valid = []

    for run, data in runs.items():
        run_numbers.append(run)
        total_valid.append(data['valid'][-1])

    df = pd.DataFrame({'Run': run_numbers, 'Total Valid Primers': total_valid})
    df = df.sort_values('Run')

    plt.figure(figsize=(12, 8))
    sns.barplot(x='Run', y='Total Valid Primers', data=df, palette='viridis')

    plt.xlabel('Run Number', fontsize=14)
    plt.ylabel('Total Valid Primers', fontsize=14)
    plt.title('Total Valid Primers Per Run', fontsize=16)
    plt.xticks(rotation=45)
    plt.show()

def plot_valid_primers_distribution(runs):
    all_valid_counts = []

    for data in runs.values():
        all_valid_counts.extend(data['valid'])

    plt.figure(figsize=(12, 8))
    sns.histplot(all_valid_counts, bins=50, kde=True, color='purple')

    plt.xlabel('Valid Primers', fontsize=14)
    plt.ylabel('Frequency', fontsize=14)
    plt.title('Distribution of Valid Primers', fontsize=16)
    plt.show()

# Path to your log file
log_file_path = 'logs.txt'

# Parse the log file and plot the data
runs_data = parse_log_file(log_file_path)
plot_runs(runs_data)
plot_total_valid_primers(runs_data)
plot_valid_primers_distribution(runs_data)
