import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

df = pd.read_csv('/media/simonray/data24/data/uio_dropbox_sr/DResearch/beathe/mybaits_vs_swift/data/median_genome_coverage__vs__80x_penalty.tsv', sep='\t')
# Clean up any remaining whitespace from column names just in case
df.columns = df.columns.str.strip()

# 2. Melt down all variables, keeping your sample column as the identifier
# (Make sure 'SAMPLE' matches the exact case in your file, e.g., 'SAMPLE' or 'Sample')
sample_col_name = 'SAMPLE'
df_melted = df.melt(id_vars=[sample_col_name], var_name='Variable', value_name='Value')
# Assuming 'df' is your loaded DataFrame
for col in df.columns:
    if col == 'SAMPLE':  # Skip your identifier column
        continue

    parts = col.split('__')
    if len(parts) != 2:
        print(f"❌ VIOLATION: Column '{col}' splits into {len(parts)} pieces instead of 2!")
# 3. Cleanly split on the double underscore
df_melted[['Group', 'Metric']] = df_melted['Variable'].str.split('__', expand=True)

# 4. Pivot the metrics back into side-by-side columns
df_tidy = df_melted.pivot(
    index=[sample_col_name, 'Group'],
    columns='Metric',
    values='Value'
).reset_index()

df_tidy.columns.name = None

# 5. Generate the Seaborn plot
sns.set_theme(style="whitegrid")
plt.figure(figsize=(10, 6))

# Dynamically map the exact names from your columns
# Assumes Metric names are exactly 'GenomeCoverage' and '80x'
sns.scatterplot(
    data=df_tidy,
    x='med_genome_cov',
    y='80x',
    hue='Group',
    style='Group',
    s=120,
    alpha=0.85
)


plt.savefig('/media/simonray/data24/data/uio_dropbox_sr/DResearch/beathe/mybaits_vs_swift/data/genome_coverage_vs_penalty.png', dpi=300)

custom_palette0 = {
    'covidseq_swift_bad': '#D32F2F',  # Dark Red
    'viralrecon_swift_bad': '#FF5252',  # Bright Red

    'covidseq_mybaits': '#1976D2',  # Deep Blue
    'viralrecon_mybaits': '#64B5F6',  # Sky Blue

    'covidseq_swift': '#388E3C',  # Forest Green
    'viralrecon_swift': '#81C784'  # Mint Green
}

# 2. Define your custom symbols/markers
# 'o' = circle, 'X' = bold cross, 's' = square, '^' = triangle, 'D' = diamond
custom_markers = {
    'covidseq_swift_bad': 'X',
    'viralrecon_swift_bad': 'o',

    'covidseq_mybaits': 'o',  # Mybaits group 1 (Circle)
    'viralrecon_mybaits': 's',  # Mybaits group 2 (Square)

    'covidseq_swift': '^',  # Swift group 1 (Triangle)
    'viralrecon_swift': 'D'  # Swift group 2 (Diamond)
}

# 1. Grab distinct shades from pre-defined libraries
reds = sns.color_palette("Reds_r", 3)  # Extracts a nice range of pre-defined reds
blues = sns.color_palette("Blues_r", 3)  # Extracts a nice range of pre-defined blues
greens = sns.color_palette("Greens_r", 3)  # Extracts a nice range of pre-defined greens

# 2. Build your custom palette mapping using the pre-defined color indices
# (Index 0 is usually the darkest, Index 1 or 2 are lighter shades)
custom_palette = {
    'covidseq_swift_bad': reds[0],
    'viralrecon_swift_bad': reds[1],

    'covidseq_mybaits': blues[0],
    'viralrecon_mybaits': blues[1],

    'covidseq_swift': greens[0],
    'viralrecon_swift': greens[1]
}


# 2. Setup the split-axis figure layout
# sharex=True binds their X-axes together
# gridspec_kw adjusts the height ratio so the bottom plot gets more visual real estate
fig, (ax_top, ax_bottom) = plt.subplots(
    2, 1,
    sharex=True,
    figsize=(10, 7),
    gridspec_kw={'height_ratios': [1, 3]}
)
fig.subplots_adjust(hspace=0.08) # Narrow the gap between the two subplots
fig.subplots_adjust(right=0.70)
ax_top.legend(bbox_to_anchor=(1.02, 0.5), loc='upper left', title='Sample Groups')
# 3. Plot the EXACT SAME data on both subplots
sns.set_theme(style="whitegrid")

for ax in [ax_top, ax_bottom]:
    sns.scatterplot(
        data=df_tidy,
        x='med_genome_cov',
        y='80x',
        hue='Group',
        style='Group',
        palette=custom_palette,  # Overrides defaults with your red/blue/green shades
        markers=custom_markers,
        s=35,
        alpha=0.7,
        ax=ax,
        legend=(ax == ax_top) # Only generate the legend once (on the top plot)
    )

# 4. Zoom into the specific Y-axis regions
ax_top.set_ylim(230, 260)     # Top axis: focused exclusively on the extreme outlier
ax_bottom.set_ylim(0, 16)     # Bottom axis: focused on the main cluster (0 to 15)

# 5. Hide structural spines to make it look like a unified plot
ax_top.spines['bottom'].set_visible(False)
ax_bottom.spines['top'].set_visible(False)
ax_top.xaxis.tick_top()
ax_top.tick_params(labeltop=False)  # Don't put X-axis labels at the very top
ax_bottom.xaxis.tick_bottom()
ax_top.tick_params(axis='y', labelsize=8)
ax_bottom.tick_params(axis='both', labelsize=8)

# 6. Add the characteristic diagonal "break marks" to the axis lines
d = .5  # Proportion factor for how big the tick marks are
kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12,
              linestyle="none", color='k', mec='k', mew=1, clip_on=False)

# Draw marks on the bottom of the top plot and top of the bottom plot
ax_top.plot([0, 1], [0, 0], transform=ax_top.transAxes, **kwargs)
ax_bottom.plot([0, 1], [1, 1], transform=ax_bottom.transAxes, **kwargs)

# 7. Labels and Legend Formatting
ax_top.set_title('Genome Coverage vs 80x Base Penalty', fontsize=12, pad=15)
ax_bottom.set_xlabel('Genome Coverage', fontsize=10)

# Clear individual Y-labels so we don't duplicate text
ax_top.set_ylabel('')
ax_bottom.set_ylabel('')
# Create a single, centered Y-axis label for the entire figure
fig.text(0.04, 0.5, '80x Base Penalty', va='center', rotation='vertical', fontsize=10)

# Position the unified legend beautifully outside the graph
ax_top.legend(bbox_to_anchor=(1.05, 0.5),
              loc='upper left',
              title='Sample Groups',
              title_fontsize=11,  # Size of the legend title text
              fontsize=10  # Size of the group names inside the legend box
              )


plt.show()
plt.savefig('/media/simonray/data24/data/uio_dropbox_sr/DResearch/beathe/mybaits_vs_swift/data/genome_coverage_vs_penalty_split.png', dpi=300, bbox_inches='tight')