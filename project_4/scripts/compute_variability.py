from Bio import AlignIO
import math
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Parameters
alignment_file = "results/msa/aligned_scaffolds.fasta"
output_raw = "results/variability/variability.txt"
output_windowed = "results/variability/variability_windowed.txt"
output_plot = "results/variability/variability_plot.png"
window_size = 100  # adjustable smoothing window

def shannon_entropy(column):
    """
    Computes Shannon entropy of a column (alignment position).
    Ignores gaps ("-").
    """
    counts = {}
    for base in column:
        if base != '-':
            counts[base] = counts.get(base, 0) + 1

    total = sum(counts.values())
    if total == 0:
        return 0.0

    entropy = 0.0
    for count in counts.values():
        p = count / total
        entropy -= p * math.log2(p)
    return entropy

# Load alignment
alignment = AlignIO.read(alignment_file, "fasta")
length = alignment.get_alignment_length()

# Compute entropy per position
entropies = []
for i in range(length):
    column = alignment[:, i]
    entropies.append(shannon_entropy(column))

# Save raw entropy
with open(output_raw, "w") as f:
    for val in entropies:
        f.write(f"{val:.5f}\n")

# Compute smoothed entropy using a sliding window
smoothed = [
    np.mean(entropies[i:i + window_size])
    for i in range(0, length - window_size + 1, window_size)
]

# Save smoothed values
with open(output_windowed, "w") as f:
    for i, val in enumerate(smoothed):
        pos = i * window_size
        f.write(f"{pos}\t{val:.5f}\n")

# Plot
sns.set(style="whitegrid")
plt.figure(figsize=(14, 5))
plt.plot(
    range(0, length - window_size + 1, window_size),
    smoothed,
    lw=1.2
)
plt.title("Sequence Variability Across Genome")
plt.xlabel("Alignment Position")
plt.ylabel("Shannon Entropy (Smoothed)")
plt.tight_layout()
plt.savefig(output_plot)
