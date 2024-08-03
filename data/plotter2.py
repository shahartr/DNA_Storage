import re
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats

# Read the log file
with open('plot/log_file.txt', 'r') as file:
    log_data = file.read()

# Extract run information
run_info = re.findall(r'run (\d+) is now processing primer (\d+)/(\d+) \(([\d.]+)% of the primers\), valid primers: (\d+)', log_data)
run_df = pd.DataFrame(run_info, columns=['Run', 'Processed', 'Total', 'Percent', 'Valid'])
run_df = run_df.astype({'Run': int, 'Processed': int, 'Total': int, 'Percent': float, 'Valid': int})

# Extract disqualification levels
disq_levels = re.findall(r'\[([^\]]+)\] primers disqualify by level', log_data)
disq_df = pd.DataFrame([list(map(int, level.split(','))) for level in disq_levels])
disq_df.columns = [f'Level_{i}' for i in range(len(disq_df.columns))]
disq_df['Run'] = range(1, len(disq_df) + 1)

# Extract filtered by max inter comp
filtered_comp = re.findall(r'filtered by max inter comp: (\d+) primers', log_data)
filtered_df = pd.DataFrame({'Filtered': list(map(int, filtered_comp))})
filtered_df['Run'] = range(1, len(filtered_df) + 1)

# Set up the plot style
sns.set(style="whitegrid")

# Ensure there is data to plot
if not disq_df.empty:
    # 1. Distribution of Disqualifications by Level
    plt.figure(figsize=(12, 6))
    disq_melted = pd.melt(disq_df, id_vars=['Run'], var_name='Level', value_name='Count')
    sns.boxplot(x='Level', y='Count', data=disq_melted, palette="Set3")
    plt.title('Distribution of Disqualifications by Level')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig('disqualifications_by_level.png', dpi=300)
    plt.close()

if not run_df.empty:
    # 2. Progress of Valid Primers (Differences between runs)
    plt.figure(figsize=(12, 6))
    for run in run_df['Run'].unique():
        run_data = run_df[run_df['Run'] == run]
        if not run_data.empty:
            plt.plot(run_data['Processed'], run_data['Valid'], label=f'Run {run}', marker='o')
    plt.title('Progress of Valid Primers (Differences between runs)')
    plt.xlabel('Processed Primers')
    plt.ylabel('Number of Valid Primers')
    plt.legend()
    plt.tight_layout()
    plt.savefig('valid_primers_progress.png', dpi=300)
    plt.close()

if not filtered_df.empty:
    # 3. Filtered Primers by Max Inter Comp (Differences between runs)
    plt.figure(figsize=(12, 6))
    sns.barplot(x='Run', y='Filtered', data=filtered_df, palette="Set2")
    plt.title('Filtered Primers by Max Inter Comp (Differences between runs)')
    plt.ylabel('Number of Filtered Primers')
    plt.tight_layout()
    plt.savefig('filtered_primers.png', dpi=300)
    plt.close()

if not run_df.empty:
    # 4. Valid Primers Percentage (Differences between runs)
    plt.figure(figsize=(12, 6))
    for run in run_df['Run'].unique():
        run_data = run_df[run_df['Run'] == run]
        if not run_data.empty:
            plt.plot(run_data['Processed'], run_data['Percent'], label=f'Run {run}', marker='x')
    plt.title('Percentage of Valid Primers (Differences between runs)')
    plt.xlabel('Processed Primers')
    plt.ylabel('Percentage of Valid Primers')
    plt.legend()
    plt.xscale('log')
    max_percent = run_df['Percent'].max()
    if not pd.isna(max_percent) and max_percent > 0:
        plt.ylim(0, max_percent * 1.1)
    plt.tight_layout()
    plt.savefig('valid_primers_percentage.png', dpi=300)
    plt.close()

# Calculate and print differences between runs
if not run_df.empty:
    last_checkpoints = run_df.groupby('Run').last().reset_index()

    # Calculate summary statistics
    summary_stats = last_checkpoints[['Run', 'Valid']].describe()

    print("\nSummary statistics for valid primers at last checkpoint:")
    print(summary_stats)

    print("\nDifferences in valid primers at last checkpoint:")
    for i in range(1, len(last_checkpoints)):
        diff = last_checkpoints.iloc[i]['Valid'] - last_checkpoints.iloc[i-1]['Valid']
        print(f"Run {i+1} vs Run {i}: {diff}")

    if not filtered_df.empty:
        print("\nDifferences in filtered primers:")
        for i in range(1, len(filtered_df)):
            diff = filtered_df.iloc[i]['Filtered'] - filtered_df.iloc[i-1]['Filtered']
            print(f"Run {i+1} vs Run {i}: {diff}")

    # Perform t-tests between runs for valid primers
    print("\nT-tests for valid primers between runs:")
    for i in range(1, len(last_checkpoints)):
        run1 = last_checkpoints.iloc[i-1]['Valid']
        run2 = last_checkpoints.iloc[i]['Valid']
        t_stat, p_value = stats.ttest_ind([run1], [run2])
        print(f"Run {i+1} vs Run {i}: t-statistic = {t_stat:.3f}, p-value = {p_value:.3f}")
