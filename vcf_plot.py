import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from vcf2pandas import vcf2pandas

# 1. Define your 5 VCF files and give them clean display names
vcf_files = {
    "MyBaits_VR_Amp-iVar": "/home/simonray/PycharmProjects/mybaits_swift/data/1_baits_NF-amp skip iVar.vcf",
    "MyBaits_NSC": "/home/simonray/PycharmProjects/mybaits_swift/data/1_mybaits_NSC.vcf",
    "Swift_VR": "/home/simonray/PycharmProjects/mybaits_swift/data/1_swift_NF amp.vcf",
    "MyBaits_VR_Meta": "/home/simonray/PycharmProjects/mybaits_swift/data/1_mybaits_NF-meta.vcf",
    "MyBaits_VR_Amp+iVar": "/home/simonray/PycharmProjects/mybaits_swift/data/1_NF-amp incl iVar.vcf",
    "Swift_NSC": "/home/simonray/PycharmProjects/mybaits_swift/data/1_swift_NSC.vcf"
}


all_dfs = []

# 2. Load each VCF file using vcf2pandas
for sample_name, file_path in vcf_files.items():
    if os.path.exists(file_path):
        print(f"Loading {sample_name}...")
        # vcf2pandas outputs standard VCF columns like '#CHROM' and 'POS'
        df = vcf2pandas(file_path)

        # Add a column so we can group/offset them on the y-axis
        df['Sample'] = sample_name

        # Keep only what we need for the genomic position plot to save memory
        all_dfs.append(df[['CHROM', 'POS', 'Sample']])
    else:
        print(f"Warning: {file_path} not found. Skipping.")

# Combine all data into one master dataframe
master_df = pd.concat(all_dfs, ignore_index=True)

# Clean up chromosome names if necessary (e.g., matching 'chr1' vs '1')
master_df['CHROM'] = master_df['CHROM'].astype(str)

# 3. Set up the Seaborn plot
plt.figure(figsize=(14, 6))
sns.set_theme(style="whitegrid")

# Using stripplot automatically handles the y-axis offset per sample!
# Adjust 'size' for marker thickness and 'alpha' for density transparency
sns.stripplot(
    data=master_df,
    x="POS",
    y="Sample",
    hue="CHROM",  # Colors the dots by chromosome
    palette="tab20",  # Handles multiple chromosomes gracefully
    jitter=False,  # Keeps them strictly on their y-axis line
    size=3,
    alpha=0.6,
    dodge=False
)

# 4. Polish the layout
plt.title("Genomic Variant Locations Across Samples", fontsize=16, fontweight='bold')
plt.xlabel("Genomic Position (bp)", fontsize=12)
plt.ylabel("Samples", fontsize=12)
plt.legend(title="Chromosome", bbox_to_anchor=(1.02, 1), loc='upper left')

# Optional: Format x-axis numbers to look like 100,000,000 instead of 1e8
plt.gca().xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: format(int(x), ',')))

plt.tight_layout()
plt.savefig("/home/simonray/PycharmProjects/mybaits_swift/plots/vcf_plot_sample1.tiff")
plt.show()