import glob
import os
import re
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def get_vcf_header_lines(filepath):
    """Finds how many lines of metadata (##) to skip in a VCF/VCD file."""
    count = 0
    with open(filepath, "r") as f:
        for line in f:
            if line.startswith("##"):
                count += 1
            else:
                break
    return count


def extract_sample_number(sample_name):
    """Extracts the first continuous block of digits from the sample name."""
    match = re.search(r"\d+", sample_name)
    return int(match.group()) if match else float("inf")


# ==========================================
# 1. LOAD, MERGE, AND SORT VCF/VCD FILES
# ==========================================
vcf_files = glob.glob("/media/simonray/data24/data/uio_dropbox_sr/DResearch/beathe/mybaits_vs_swift/VCF-files_for_Simon/swift+NF-viralrecon amp/*.vcf")  # Change extension if needed

sample_dict = {}
for file_path in vcf_files:
    sample_name = os.path.splitext(os.path.basename(file_path))[0]
    skip_rows = get_vcf_header_lines(file_path)

    # Read VCF, extract POS column
    df = pd.read_csv(file_path, sep="\t", skiprows=skip_rows)
    sample_dict[sample_name] = df["POS"].dropna().astype(int).unique()

# Melt into long-form and pivot into a presence/absence matrix
data = []
for sample, positions in sample_dict.items():
    for pos in positions:
        data.append({"Sample": sample, "POS": pos})

df_long = pd.DataFrame(data)
final_matrix = pd.crosstab(index=df_long["POS"], columns=df_long["Sample"])
final_matrix = final_matrix.sort_index()  # Sorts rows (Genomic Positions)

# --- NEW: Sort columns (Samples) numerically by the first number found ---
sorted_columns = sorted(final_matrix.columns, key=extract_sample_number)
final_matrix = final_matrix[sorted_columns]

# Save sorted matrix to a file
final_matrix.to_csv("/media/simonray/data24/data/uio_dropbox_sr/DResearch/beathe/mybaits_vs_swift/VCF-files_for_Simon/swift+NF-viralrecon amp/merged_vcf_mutations.csv")

# 4. Plot the Seaborn Heatmap
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# 1. Set up figure size
plt.figure(figsize=(10, 10))  # Can be a bit shorter now without a wall of text on the Y-axis

# 2. Define colors
my_colors = ["#e0e0e0", "#1f77b4"]
cmap = sns.color_palette(my_colors)

# 3. Create heatmap (turn OFF default yticklabels)
ax = sns.heatmap(
    final_matrix,
    cmap=cmap,
    cbar=False,
    linewidths=0,  # Setting to 0 looks cleaner when not using explicit row lines
    yticklabels=False,
)

# 4. Generate the custom scale for the Y-axis
# Find the actual genomic range from your data
min_pos = 0
max_pos = final_matrix.index.max()
# Round max_pos up to the nearest 1000 for a clean scale boundary
max_genome_limit = int(np.ceil(max_pos / 1000.0) * 1000)

# Create tick positions every 1000 nt, and labels every 5000 nt
tick_positions = []
tick_labels = []

for pos in range(0, max_genome_limit + 1, 1000):
    # Map the genomic position to the closest matrix row index position
    # (This ensures the tick line points to the correct relative spot on the heatmap)
    idx_pos = np.searchsorted(final_matrix.index, pos)

    if pos % 5000 == 0:
        tick_positions.append(idx_pos)
        tick_labels.append(f"{pos:,}")
    else:
        tick_positions.append(idx_pos)
        tick_labels.append("")  # Major tick, but no text label

# 5. Apply the custom genomic ticks and labels to the Y-axis
ax.set_yticks(tick_positions)
ax.set_yticklabels(tick_labels, size=10)

# Optional: Add small tick tick-marks pointing outwards or inwards
ax.tick_params(axis="y", which="major", left=True, length=6, color="black")

# 6. Final Polish
plt.title(
    "Viral Mutation Map Across Samples", fontsize=16, fontweight="bold", pad=20
)
plt.xlabel("Sample Names", fontsize=12, fontweight="bold", labelpad=10)
plt.ylabel(
    "Genomic Position (Nucleotides)", fontsize=12, fontweight="bold", labelpad=10
)
plt.xticks(rotation=45, ha="right")

plt.tight_layout()
plt.savefig("/media/simonray/data24/data/uio_dropbox_sr/DResearch/beathe/mybaits_vs_swift/VCF-files_for_Simon/swift+NF-viralrecon amp/vcf_mutation_heatmap.png", dpi=300)
plt.show()