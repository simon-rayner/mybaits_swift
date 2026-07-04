import glob
import os
import re
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from adjustText import adjust_text

vcf_files = glob.glob("/media/simonray/data24/data/uio_dropbox_sr/DResearch/beathe/mybaits_vs_swift/VCF-files_for_Simon/mybaits_NF-viralrecon amp incl iVar trim/*.vcf")
ref_vcf_path = "/media/simonray/data24/data/uio_dropbox_sr/DResearch/beathe/mybaits_vs_swift/VCF-files_for_Simon/B_1_1_529.vcf"
fig_path = "/media/simonray/data24/data/uio_dropbox_sr/DResearch/beathe/mybaits_vs_swift/plots/mybaits_NF_viralrecon_amp_vcf_mutation_heatmap_wrects.png"
plt_title = "SNPs by Sample: MyBaits/ViralRecon"
highlight_txt_path = "/media/simonray/data24/data/uio_dropbox_sr/DResearch/beathe/mybaits_vs_swift/plots/mybaits_NF__highlights.txt"

vcf_files = glob.glob("/media/simonray/data24/data/uio_dropbox_sr/DResearch/beathe/mybaits_vs_swift/VCF-files_for_Simon/swift+NF-viralrecon amp/*.vcf")
ref_vcf_path = "/media/simonray/data24/data/uio_dropbox_sr/DResearch/beathe/mybaits_vs_swift/VCF-files_for_Simon/B_1_1_529.vcf"
fig_path = "/media/simonray/data24/data/uio_dropbox_sr/DResearch/beathe/mybaits_vs_swift/plots/swift_NF_viralrecon_amp_vcf_mutation_heatmap_wrects.png"
highlight_txt_path = "/media/simonray/data24/data/uio_dropbox_sr/DResearch/beathe/mybaits_vs_swift/plots/swift_NF__highlights.txt"  # Adjust path to your file location
plt_title = "SNPs by Sample: Swift/ViralRecon"



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
# 1. LOAD, MERGE, AND SORT USER VCF/VCD FILES
# ==========================================

sample_dict = {}
for file_path in vcf_files:
    sample_name = os.path.splitext(os.path.basename(file_path))[0]
    skip_rows = get_vcf_header_lines(file_path)

    df = pd.read_csv(file_path, sep="\t", skiprows=skip_rows)
    sample_dict[sample_name] = df["POS"].dropna().astype(int).unique()

data = []
for sample, positions in sample_dict.items():
    for pos in positions:
        data.append({"Sample": sample, "POS": pos})

df_long = pd.DataFrame(data)
final_matrix = pd.crosstab(index=df_long["POS"], columns=df_long["Sample"])
final_matrix = final_matrix.sort_index()

# Sort columns numerically
sorted_columns = sorted(final_matrix.columns, key=extract_sample_number)
final_matrix = final_matrix[sorted_columns]

# ==========================================
# 2. LOAD REFERENCE VCF FILE
# ==========================================
 # Adjust to your reference file path
ref_skip = get_vcf_header_lines(ref_vcf_path)
ref_df = pd.read_csv(ref_vcf_path, sep="\t", skiprows=ref_skip)

# Create a dictionary mapping the mutation position to its label (e.g., 23403 -> "D614G")
# Note: VCF standard puts custom annotations in the INFO or ID column. Adjust if needed.
ref_mutations = {}
for _, row in ref_df.iterrows():
    pos = int(row["POS"])
    # Fallback to the ID or position string if an explicit amino acid label isn't present
    label = (
        row["ID"]
        if pd.notna(row["ID"]) and row["ID"] != "."
        else f"Mutation_{pos}"
    )

    # Clean up standard INFO text if you stored labels there (e.g., AA=S:D614G)
    if "AA=" in str(row.get("INFO", "")):
        match = re.search(r"AA=([^;]+)", str(row["INFO"]))
        if match:
            label = match.group(1)

    ref_mutations[pos] = label

# ==========================================
# 3. SET UP A TWO-GRID SUBPLOT SYSTEM
# ==========================================
# Create a figure with two side-by-side plots:
# [Heatmap Matrix] [Reference Label Column]
fig, ax1 = plt.subplots(figsize=(10.5, 12))

# Define presence/absence colors
my_colors = ["#e0e0e0", "#1f77b4"]
cmap = sns.color_palette(my_colors)

# Provide the specific sample names or 0-indexed column positions to highlight
# Indices match your target sample columns (e.g., 3rd, 5th, 6th column...)
#highlight_indices = [3, 5, 6, 11, 14, 15, 16, 20, 21, 22]
observed_positions = final_matrix.index.values
highlight_indices = []
box_color = "#6E7D8C"  # Default to Slate Grey if file read fails

if os.path.exists(highlight_txt_path):
    with open(highlight_txt_path, "r") as f:
        lines = [line.strip() for line in f.readlines() if line.strip()]

        if lines:
            # If the first line starts with #, treat it as the hex code
            if lines[0].startswith("#"):
                box_color = lines[0]
                indices_line = lines[1] if len(lines) > 1 else ""
            else:
                indices_line = lines[0]

            # Parse the indices out of the remaining text
            if indices_line:
                highlight_indices = [
                    int(idx.strip()) for idx in re.split(r"[,\s]+", indices_line) if idx.strip().isdigit()
                ]
else:
    print(f"Warning: Highlight configuration file '{highlight_txt_path}' not found. Skipping box highlights.")

total_rows = len(final_matrix.index)


for idx in highlight_indices:
    if idx < len(sorted_columns):
        rect = plt.Rectangle(
            (idx, 0),               # Start at top left corner of column
            1.0,                    # Exactly 1 column wide
            total_rows,             # Span down the full row space
            edgecolor=box_color,    # #00AA4F Crisp green #
                                    # #E53935 (A slightly brighter, classic cherry red)
                                    # #C62828 (A deeper, slightly darker crimson for high contrast against light backgrounds)
            facecolor="none",
            linewidth=1.5,
            zorder=3                # Forces border overlay on top of background grids
        )
        ax1.add_patch(rect)
# -------------------------------------------------------------

# ==========================================
# 4. LEFT PANEL: THE HEATMAP MATRIX
# ==========================================
sns.heatmap(
    final_matrix,
    cmap=cmap,
    cbar=False,
    linewidths=0.5,
    linecolor="white",
    yticklabels=False,  # Checked: Managed by the custom left ruler loop below
    ax=ax1,
)

# Left Y-Axis: The Clean Genomic Ruler (1,000 nt ticks / 5,000 nt labels)

max_pos = observed_positions.max()

left_tick_positions = []
left_tick_labels = []

for pos in range(0, max_pos + 1000, 1000):
    if pos == 0:
        continue
    idx_pos = np.searchsorted(observed_positions, pos)
    if idx_pos < len(observed_positions):
        left_tick_positions.append(idx_pos)
        if pos % 5000 == 0:
            left_tick_labels.append(f"{pos:,}")
        else:
            left_tick_labels.append("")

ax1.set_yticks(left_tick_positions)
ax1.set_yticklabels(left_tick_labels, size=10)
ax1.tick_params(axis="y", which="major", left=True, length=6, color="black")
ax1.set_ylabel(
    "Genome location", fontsize=12, fontweight="bold"
)

# Heatmap X-Axis Labels
ax1.set_xlabel("Sorted Sample Names", fontsize=12, fontweight="bold", labelpad=10)
ax1.set_xticklabels(ax1.get_xticklabels(), rotation=45, ha="right")

# ==========================================
# 5. RIGHT SIDE: RAW GEOMETRY VECTOR ARROWS (NO EXTRA OUTLINES)
# ==========================================
num_samples = len(sorted_columns)

# Get reference positions and sort them top-to-bottom
ref_positions_in_data = [
    p for p in ref_mutations.keys() if p in observed_positions
]

# We will collect the coordinates for all arrows and plot them simultaneously
X_points = []
Y_points = []

for ref_pos in ref_positions_in_data:
    matrix_row_idx = np.where(observed_positions == ref_pos)[0][0]
    actual_y = matrix_row_idx + 0.5

    # Base of the arrow starts just outside the heatmap
    X_points.append(num_samples + 0.5)
    Y_points.append(actual_y)

# Convert to numpy arrays for the vector engine
X = np.array(X_points)
Y = np.array(Y_points)

# Vector directions: We want them pointing left (U = -1, V = 0)
U = np.ones_like(X) * -1.0
V = np.zeros_like(Y)

# Draw the arrows as raw data vectors. This guarantees absolutely no secondary axis borders.
ax1.quiver(
    X,
    Y,
    U,
    V,
    color="#d95f02",
    scale=1,
    scale_units="x",  # Locks the arrow length strictly to the X-axis grid units
    angles="xy",
    width=0.004,  # Adjusts line thickness of the arrow stem
    headwidth=4,  # Adjusts how flared the arrowhead is
    headlength=5,
    clip_on=False,
)

# Add a clean text indicator running parallel to the heatmap edge
total_rows = len(observed_positions)
ax1.text(
    x=num_samples + 0.8,
    y=total_rows / 2.0,
    s="Omicron: B_1_1_529",
    color="#d95f02",
    weight="bold",
    size=11,
    rotation=270,
    va="center",
    ha="left",
    clip_on=False,
)

# ==========================================
# 6. FINAL POLISH, RESIZE & GHOST AXIS CLEANUP
# ==========================================
ax1.set_xlabel("Sorted Sample Names", fontsize=12, fontweight="bold", labelpad=10)
ax1.set_xticklabels(ax1.get_xticklabels(), rotation=45, ha="right")

plt.title(
    plt_title,
    fontsize=16,
    fontweight="bold",
    pad=20,
)

# Clean up any residual rogue layout frames
for extra_ax in fig.axes:
    if extra_ax != ax1:
        extra_ax.set_axis_off()
        for spine in extra_ax.spines.values():
            spine.set_visible(False)

# NEW LAYOUT BOUNDARIES:
# Setting right=0.95 forces the heatmap matrix to stretch out and occupy
# 95% of the total canvas width, leaving a crisp 5% strip for the arrows.
plt.subplots_adjust(left=0.08, right=0.95, top=0.92, bottom=0.18)

plt.savefig(fig_path, dpi=300, bbox_inches="tight")
plt.show()