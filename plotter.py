import matplotlib.pyplot as plt

# Data provided
lengths = [9, 10, 11, 12, 13, 14, 15]
primers_count = [0, 189600, 1294544, 2483096, 17171840, 32804376, 228915528]
primers_count_after_filtering = [0, 179, 724, 682, 0, 0, 0]

# Plotting the data in two separate plots
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12))

# First plot: Before Filtering
ax1.plot(lengths, primers_count, marker='o', linestyle='-', color='b', label='Before Filtering')
ax1.set_title('Number of Primers by Length (Before Filtering)')
ax1.set_xlabel('Primer Length')
ax1.set_ylabel('Number of Primers')
ax1.set_yscale('log')
ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
ax1.legend()

# Second plot: After Filtering
ax2.plot(lengths, primers_count_after_filtering, marker='o', linestyle='-', color='r', label='After Filtering')
ax2.set_title('Number of Primers by Length (After Filtering)')
ax2.set_xlabel('Primer Length')
ax2.set_ylabel('Number of Primers')
ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
ax2.legend()

# Show the plots
plt.tight_layout()
plt.show()
