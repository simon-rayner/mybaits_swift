import gzip
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

def load_bed(path):
    # Reads BED files, automatically accommodating 3 or 4-column variants
    df = pd.read_csv(path, sep='\t', comment='#', header=None)
    if df.shape[1] >= 4:
        df.columns = ['chrom', 'start', 'end', 'coverage'] + list(df.columns[4:])
        df['position'] = (df['start'] + df['end']) / 2
    else:
        df.columns = ['chrom', 'position', 'coverage']
    return df

def load_vcf_positions(path):
    positions = []
    open_func = gzip.open if path.endswith('.gz') else open
    mode = 'rt' if path.endswith('.gz') else 'r'
    with open_func(path, mode) as f:
        for line in f:
            if not line.startswith('#'):
                parts = line.split('\t')
                if len(parts) > 1:
                    positions.append(int(parts[1]))
    return positions


def plot_beds_and_vcfs(bed_file, vcf_file, ref_vcf_file, vcf3_file, plot_title, sample_name, vcf3_name, figfile):
    bed_df = load_bed(bed_file)

    vcf1_pos = load_vcf_positions(vcf_file)
    vcf2_ref = load_vcf_positions(ref_vcf_file)

    # Check if the third VCF is provided
    has_vcf3 = bool(vcf3_file and str(vcf3_file).strip())

    chrom_name = bed_df['chrom'].iloc[0]
    min_pos, max_pos = bed_df['position'].min(), bed_df['position'].max()

    # --- 2. Plotting Layout Adjustments ---
    sns.set_theme(style="white")

    # Adjust height ratio dynamically based on track count
    height_ratios = [1.2, 4] if has_vcf3 else [0.8, 4]
    fig, (ax_vcf, ax_cov) = plt.subplots(
        2, 1, figsize=(12, 5.5 if has_vcf3 else 5), sharex=True,
        gridspec_kw={'height_ratios': height_ratios}
    )

    # Base tracks setup
    tracks_y = [1, 2]
    track_labels = [sample_name, 'B.1.1.529']

    if has_vcf3:
        vcf3_pos = load_vcf_positions(vcf3_file)
        tracks_y.append(3)
        track_labels.append(vcf3_name)

    # Top Panel: Draw the matching number of horizontal track lines
    ax_vcf.hlines(y=tracks_y, xmin=min_pos, xmax=max_pos, colors='#e0e0e0', linestyle='-', linewidth=1)

    # Scatter plot data onto rows
    ax_vcf.scatter(vcf1_pos, [1] * len(vcf1_pos), color='#d32f2f', s=35, zorder=3, label=sample_name)
    ax_vcf.scatter(vcf2_ref, [2] * len(vcf2_ref), color='#2e7d32', s=35, zorder=3, label='B.1.1.529')

    if has_vcf3:
        ax_vcf.scatter(vcf3_pos, [3] * len(vcf3_pos), color='#7b1fa2', s=35, zorder=3, label=vcf3_name)

    # Dynamic styling for the top panel based on tracks loaded
    ax_vcf.set_yticks(tracks_y)
    ax_vcf.set_yticklabels(track_labels, fontsize=10, fontweight='bold')
    ax_vcf.set_ylim(0.4, 3.6 if has_vcf3 else 2.6)
    ax_vcf.set_title(plot_title, fontsize=14, pad=12, fontweight='bold')
    sns.despine(ax=ax_vcf, left=True, bottom=True)
    ax_vcf.tick_params(left=False, bottom=False)

    # Bottom Panel: Coverage Depth
    sns.lineplot(data=bed_df, x='position', y='coverage', color='#1976d2', linewidth=1.5, ax=ax_cov)
    ax_cov.fill_between(bed_df['position'], bed_df['coverage'], color='#1976d2', alpha=0.15)

    # Style the Bottom Panel
    ax_cov.set_ylabel('Coverage Depth', fontsize=11, fontweight='bold')
    ax_cov.set_xlabel('Genome Position (bp)', fontsize=11, fontweight='bold')
    ax_cov.set_xlim(min_pos, max_pos)
    ax_cov.set_ylim(0, bed_df['coverage'].max() * 1.05)
    sns.despine(ax=ax_cov)

    plt.tight_layout()
    plt.savefig(figfile, dpi=300)
    plt.show()

vcf_data_home = "/media/simonray/data24/data/uio_dropbox_sr/DResearch/beathe/mybaits_vs_swift/data_for_publication/vcf_files/"
meta_data_home = "/media/simonray/data24/data/uio_dropbox_sr/DResearch/beathe/mybaits_vs_swift/data_for_publication/vcf_files/mybaits_metagenome_NF"
bed_data_home =  "/media/simonray/data24/data/uio_dropbox_sr/DResearch/beathe/mybaits_vs_swift/data_for_publication/bed_files/"
plot_data_home = "/media/simonray/data24/data/uio_dropbox_sr/DResearch/beathe/mybaits_vs_swift/data_for_publication/plots/"

f"{Path(vcf_data_home)}/B_1_1_529.vcf"
vcf_ref_file = f"{Path(vcf_data_home)}/B_1_1_529.vcf"
# --- 1. Load your actual data here ---


print("mybaits_NF - sample 4")
mybaits_vcf_amp_sample4 = f"{Path(vcf_data_home, "mybaits_amplicon_NF")}/4_NF-amp incl iVar.vcf"
mybaits_vcf_meta_sample4 = f"{Path(vcf_data_home, "mybaits_metagenome_NF")}/4_mybaits_NF-meta.vcf"
mybaits_bed_sample4 = f"{Path(bed_data_home, "mybaits_amplicon_NF")}/4_baits_NF-amplicon_coverage.bed"
mybaits_plotfile_out_sample4 = f"{Path(plot_data_home)}/mybaits_sample4_coverage_map.png"
mybaits_plot_title_sample4 = f'Read Coverage - Sample 4 mybaits_amplicon_NF-viralrecon (B.1.617.2)'
mybaits_xtitle_sample4_amp = "sample 4: AMP"
mybaits_xtitle_sample4_meta = "sample 4:META"
plot_beds_and_vcfs(mybaits_bed_sample4, mybaits_vcf_amp_sample4, vcf_ref_file, mybaits_vcf_meta_sample4, mybaits_plot_title_sample4, mybaits_xtitle_sample4_amp, mybaits_xtitle_sample4_meta, mybaits_plotfile_out_sample4)

print("mybaits_NF - sample 18")
mybaits_vcf_amp_sample18 = f"{Path(vcf_data_home, "mybaits_amplicon_NF")}/18_NF-amp incl iVar.vcf"
mybaits_vcf_meta_sample18 = f"{Path(vcf_data_home, "mybaits_metagenome_NF")}/18_mybaits_NF-meta.vcf"
mybaits_bed_sample18 = f"{Path(bed_data_home, "mybaits_amplicon_NF")}/18_baits_NF-amplicon_coverage.bed"
mybaits_plotfile_out_sample18 = f"{Path(plot_data_home)}/mybaits_sample18_coverage_map.png"
mybaits_plot_title_sample18 = f'Read Coverage - Sample 18 mybaits_amplicon_NF-viralrecon (B.1.617.2)'
mybaits_xtitle_sample18_amp = "sample 18:AMP"
mybaits_xtitle_sample18_meta = "sample 18:META"
plot_beds_and_vcfs(mybaits_bed_sample18, mybaits_vcf_amp_sample18, vcf_ref_file, mybaits_vcf_meta_sample18, mybaits_plot_title_sample18, mybaits_xtitle_sample18_amp, mybaits_xtitle_sample18_meta, mybaits_plotfile_out_sample18)

print("swift_NF - sample 4")
swift_vcf_sample4 = f"{Path(vcf_data_home, "swift_NF")}/4_swift_NF amp.vcf"
swift_bed_sample4 = f"{Path(bed_data_home, "swift_NF")}/4_swift_NF_coverage.bed"
swift_plotfile_out_sample4 = f"{Path(plot_data_home)}/swift_sample4_coverage_map.png"
swift_plot_title_sample4 = f'Read Coverage - Sample 4 swift_NF-viralrecon (B.1.617.2)'
swift_xtitle_sample4 = "sample 4"
plot_beds_and_vcfs(swift_bed_sample4, swift_vcf_sample4, vcf_ref_file, "", swift_plot_title_sample4, swift_xtitle_sample4, "", swift_plotfile_out_sample4)

print("swift_NF - sample 18")
swift_vcf_sample18 = f"{Path(vcf_data_home, "swift_NF")}/18_swift_NF amp.vcf"
swift_bed_sample18 = f"{Path(bed_data_home, "swift_NF")}/18_swift_NF_coverage.bed"
swift_plotfile_out_sample18 = f"{Path(plot_data_home)}/swift_sample18_coverage_map.png"
swift_plot_title_sample18 = f'Read Coverage - Sample 18 swift_NF-viralrecon (B.1.617.2)'
swift_xtitle_sample18 = "sample 18"
plot_beds_and_vcfs(swift_bed_sample18, swift_vcf_sample18, vcf_ref_file, "", swift_plot_title_sample18, swift_xtitle_sample18, "", swift_plotfile_out_sample18)


