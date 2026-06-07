import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# 1. Create dummy data matching your description (Replace this with your DataFrame)
import numpy as np
np.random.seed(42)
data = {
    'Sample': [f'Sample_{i+1}' for i in range(20)],
    'Protocol_1': np.random.normal(loc=12, scale=2, size=20),
    'Protocol_2': np.random.normal(loc=18, scale=3, size=20)
}
#df = pd.DataFrame(data)
df = pd.read_csv('/media/simonray/data24/data/uio_dropbox_sr/DResearch/beathe/mybaits_vs_swift/data/10x_genome_coverage.tsv', sep='\t')
# 2. Extract column names
cols = df.columns


# 2. Reshape the DataFrame from wide format to long format (required for Seaborn)
df_long = df.melt(id_vars=['Sample'], value_vars=['viralrecon_swift', 'covidseq_swift'],
                  var_name='Protocol', value_name='Coverage_10x')

# 3. Plotting
sns.set_theme(style="whitegrid")

# Create a categorical line plot where 'units' ensures lines connect individual samples
plt.figure(figsize=(6, 6))
sns.lineplot(
    data=df_long,
    x='Protocol',
    y='Coverage_10x',
    hue='Sample',
    units='Sample',
    estimator=None,
    alpha=0.7,
    marker='o',
    legend=False  # Hide legend if 20 colors are too cluttered
)

# Overlay a boxplot or pointplot to show the summary statistics / trends


sns.pointplot(
    data=df_long,
    x='Protocol',
    y='Coverage_10x',
    color='black',
    linewidth=3.6,      # Replaces scale for line thickness
    markersize=7.2,     # Replaces scale for marker/dot size
    errorbar=None
)

plt.title('10x Coverage Comparison between Protocols')
plt.ylabel('10x Coverage Value')
plt.xlabel('Protocol')
plt.tight_layout()
plt.savefig('/media/simonray/data24/data/uio_dropbox_sr/DResearch/beathe/mybaits_vs_swift/data/paired_line_plot.tiff')