import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import glob
import os
from pathlib import Path
import re

data_home = "/media/simonray/data24/data/uio_dropbox_sr/DResearch/beathe/mybaits_vs_swift/data_for_publication/bed_files/"

def plot_bed_distributions(file_pattern):
    # 1. Load all BED files into a single DataFrame
    all_data = []

    # Finds all files matching the pattern (e.g., "*.bed")
    for file_path in glob.glob(file_pattern):
        file_name = os.path.basename(file_path)

        # Read BED file (skipping headers if present)
        df = pd.read_csv(file_path, sep='\t', header=None,
                         names=['chrom', 'start', 'end'], usecols=[0, 1, 2])

        # Calculate interval lengths
        df['length'] = df['end'] - df['start']
        df['sample'] = file_name  # Label for Seaborn hue
        all_data.append(df)

    # Combine everything
    combined_df = pd.concat(all_data, axis=0)

    # 2. Set up the Seaborn plot
    plt.figure(figsize=(10, 6))
    sns.set_style("whitegrid")

    # Create a KDE plot to compare distributions of lengths
    plot = sns.histplot(
        data=combined_df,
        x='length',
        hue='sample',
        element="step",
        kde=True,
        log_scale=True  # Genomic intervals often vary by orders of magnitude
    )

    plt.title('Distribution of Interval Lengths across BED Files')
    plt.xlabel('Interval Length (bp) - Log Scale')
    plt.ylabel('Frequency')
    plt.tight_layout()
    # plt.show()

    # Save the figure
    plt.savefig("baits_NSC.png", dpi=300, bbox_inches='tight')


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import glob
import os


def plot_normalized_viral_coverage(file_pattern, output_name="normalized_coverage.png"):
    all_samples = []

    # 1. Load and Normalize data
    for file_path in glob.glob(file_pattern):
        sample_name = os.path.basename(file_path).replace(".bed", "")

        # Expects: chrom, pos, depth
        df = pd.read_csv(file_path, sep='\t', header=None,
                         names=['chrom', 'start_pos', 'pos', 'depth'], usecols=[0, 1, 2, 3])

        # NORMALIZATION: Divide depth by the total sum of depth in this file
        # We multiply by 1e6 to get "Reads Per Million" (RPM) for readability
        total_depth = df['depth'].sum()
        df['normalized_depth'] = (df['depth'] / total_depth) * 1e6

        df['Sample'] = sample_name
        all_samples.append(df)

    # Combine into a single DataFrame
    combined_df = pd.concat(all_samples)

    # 2. Plotting
    plt.figure(figsize=(12, 6))
    sns.set_style("ticks")

    # Use 'normalized_depth' for the Y-axis
    sns.lineplot(
        data=combined_df.reset_index(),
        x='pos',
        y='normalized_depth',
        hue='Sample',
        linewidth=1.2,
        alpha=0.8
    )

    # Styling the plot
    title = "Genomic Coverage across Viral Genome:" + Path(output_name).stem
    plt.title(title, fontsize=14, fontweight='bold')
    plt.xlabel('Position (bp)', fontsize=12)
    plt.ylabel('Read Depth', fontsize=12)
    plt.xlim(0, combined_df['pos'].max())  # Set X-axis to genome length
    plt.grid(axis='y', linestyle='--', alpha=0.4)

    # Move legend outside the plot area
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

    # 3. Save and Show
    plt.tight_layout()
    plt.savefig(output_name, dpi=300, bbox_inches='tight')
    print(f"Plot saved as {output_name}")
    plt.show()
    plt.close()


def plot_normalized_viral_coverage_grid(file_pattern, output_name="normalized_coverage.png"):


    #file_pattern = "path/to/files/*.bed"  # Update this to your path
    all_samples = []
    file_paths = glob.glob(file_pattern)

    def get_leading_number(path):
        filename = os.path.basename(path)
        match = re.match(r'^(\d+)', filename)
        # Return the integer if found, otherwise return a huge number so unnumbered files go to the end
        return int(match.group(1)) if match else float('inf')

    # Sort the file paths numerically based on their leading numbers
    file_paths.sort(key=get_leading_number)
    sample_order = []

    for file_path in file_paths:




        sample_name = os.path.basename(file_path).replace(".bed", "")
        sample_order.append(sample_name)  # Saves the correctly ordered name

        # Safely read the 4 columns. If your file only has 3 columns, adjust names/usecols.
        df = pd.read_csv(
            file_path,
            sep='\t',
            header=None,
            names=['chrom', 'start_pos', 'pos', 'depth'],
            usecols=[0, 1, 2, 3]
        )

        # NORMALIZATION: Compute RPM
        total_depth = df['depth'].sum()
        if total_depth > 0:
            df['normalized_depth'] = (df['depth'] / total_depth) * 1e6
        else:
            df['normalized_depth'] = 0

        df['Sample'] = sample_name
        all_samples.append(df)

    # Combine into a single DataFrame
    combined_df = pd.concat(all_samples, ignore_index=True)

    # --- GRID PLOTTING ---
    sns.set_style("ticks")

    # Generate the facet grid automatically
    g = sns.relplot(
        data=combined_df,
        x='pos',
        y='normalized_depth',
        col='Sample',         # Separate plot for each sample
        col_order=sample_order,
        col_wrap=3,           # Adjust this to change how many columns per row
        kind='line',
        height=3.5,           # Height of each individual subplot
        aspect=1.5,           # Aspect ratio (width = height * aspect)
        linewidth=1.2,
        color='darkcyan',     # Uniform clean color across the grid
        facet_kws={'sharey': False},  # Change to True if you want a uniform Y-axis across all samples
        sort=False,           # Skips sorting optimization for faster rendering
        estimator=None        # Skips error-bar checks for immediate drawing
    )

    # Clean up axis labels and titles
    g.set_titles("{col_name}")


    g.set_axis_labels("Genomic Position", "Reads Per Million (RPM)")

    # Fit layout nicely so titles don't overlap
    plt.tight_layout()
    plt.savefig(output_name, dpi=300, bbox_inches='tight')
    print(f"Plot saved as {output_name}")
    plt.show()


print("mybaits_amplicon_NF")
plot_normalized_viral_coverage_grid(
    f"{Path(data_home, "mybaits_amplicon_NF")}/*.bed",
    f"{Path(data_home, "mybaits_amplicon_NF")}/mybaits_amplicon_NF_norm_coverage.png")

print("myBaits_NSC")
plot_normalized_viral_coverage_grid(
    f"{Path(data_home, "mybaits_NSC")}/*.bed",
    f"{Path(data_home, "mybaits_NSC")}/mybaits_NSC_norm_coverage.png")

print("Swift_NSC")
plot_normalized_viral_coverage_grid(
    f"{Path(data_home, "swift_NSC")}/*.bed",
    f"{Path(data_home, "swift_NSC")}/swift_nsc_norm_coverage.png")

print("swift_NF")
plot_normalized_viral_coverage_grid(
    f"{Path(data_home, "swift_NF")}/*.bed",
    f"{Path(data_home, "swift_NF")}/swift_NF_norm_coverage.png")

print("myBaits_NF_metagenomic")
plot_normalized_viral_coverage_grid(
    f"{Path(data_home, "mybaits_metagenome_NF")}/*.bed",
    f"{Path(data_home, "mybaits_metagenome_NF")}/mybaits_metagenome_NF_norm_coverage.png")
