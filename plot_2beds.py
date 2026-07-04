import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


# 1. Define a helper function to read standard BED/bedGraph coverage files
def load_bed_coverage(filepath, sample_name):
    # Standard bedGraph has columns: chrom, start, end, coverage
    # Adjust names if your BED format has different columns
    df = pd.read_csv(filepath, sep='\t', header=None,
                     names=['chrom', 'start', 'end', 'coverage'])

    # Calculate the midpoint of each genomic interval to use as the X coordinate
    df['position'] = (df['start'] + df['end']) / 2
    df['sample'] = sample_name
    return df


# 2. Load your two datasets
df1 = load_bed_coverage('/media/simonray/data24/data/uio_dropbox_sr/DResearch/beathe/mybaits_vs_swift/SV_Re_Files_in_TSD_p1672/all_beds/baits_NSC/1_baits_NSC_coverage.bed', 'MSC_MyBaits')
df2 = load_bed_coverage('/media/simonray/data24/data/uio_dropbox_sr/DResearch/beathe/mybaits_vs_swift/SV_Re_Files_in_TSD_p1672/all_beds/swift_NSC/1_swift_NSC_coverage.bed', 'NSC_Swift')
df3 = load_bed_coverage('/media/simonray/data24/data/uio_dropbox_sr/DResearch/beathe/mybaits_vs_swift/SV_Re_Files_in_TSD_p1672/all_beds/baits_NF_meta/1_baits_NF_meta_coverage.bed', 'VR_MyBaits_Meta')

# Combine them vertically into a single long-form DataFrame
combined_df = pd.concat([df1, df2, df3], ignore_index=True)

# Optional: Filter for a specific chromosome if your genome is large
# combined_df = combined_df[combined_df['chrom'] == 'chr1']

# 3. Create the offset dot plot using FacetGrid
# 'row="sample"' creates the vertical track offset you asked for
g = sns.FacetGrid(
    combined_df,
    row="sample",
    hue="sample",
    palette="Set1",
    aspect=4,  # Makes the tracks wide/stretched out like a genome browser
    height=3,
    sharex=True  # Crucial: ensures the genomic coordinates align perfectly
)

# Map a scatterplot (dot plot) onto each facet
# s=5 or smaller prevents massive overlapping if you have high-density data
g.map(sns.scatterplot, "position", "coverage", s=4, alpha=0.7, linewidth=0)

# 4. Clean up the aesthetics
g.set_titles(row_template="{row_name}")
g.set_axis_labels("Genomic Position (bp)", "Read Coverage")
plt.subplots_adjust(hspace=0.3)  # Adds a small gap between the two tracks
plt.savefig("/home/simonray/PycharmProjects/mybaits_swift/plots/bedPlot_MyBaits_vs_SWIFT_sample1.tiff")
plt.show()