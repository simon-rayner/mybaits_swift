import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# 1. Dummy data structure including Ct values (Replace with your actual data)
np.random.seed(42)
data = {
    'Sample': [f'Sample_{i+1}' for i in range(20)],
    'Protocol_1': np.random.normal(loc=12, scale=2, size=20),  # 10x coverage
    'Ct_Value': np.random.normal(loc=22, scale=3, size=20),    # qPCR Ct value
    'Protocol_2': np.random.normal(loc=18, scale=3, size=20)   # 10x coverage
}
df = pd.DataFrame(data)
df = pd.read_csv('/media/simonray/data24/data/uio_dropbox_sr/DResearch/beathe/mybaits_vs_swift/data/10x_genome_coverage.tsv', sep='\t')
# 2. Reshape into long-format containing all three variables


# 2. Apply Z-score Standardization: (Value - Mean) / StdDev
# This centers everything around 0 and scales by variance
df_zscore = df.copy()

for col in ['viralrecon_swift',  'Ct_Value', 'covidseq_swift']:
    col_mean = df_zscore[col].mean()
    col_std = df_zscore[col].std()
    df_zscore[col] = (df_zscore[col] - col_mean) / col_std


df_zscore['Ct_Value'] = df_zscore['Ct_Value'] * -1
#df_zscore['Ct_Value'] = 40 - df_zscore['Ct_Value']

# 4. Melt into long format
df_long = df_zscore.melt(
    id_vars=['Sample'],
    value_vars=['viralrecon_swift',  'Ct_Value', 'covidseq_swift'],
    var_name='Measurement',
    value_name='Z_Score'
)

order = ['viralrecon_swift',  'Ct_Value', 'covidseq_swift']
df_long['Measurement'] = pd.Categorical(df_long['Measurement'], categories=order, ordered=True)

# 5. Plotting
sns.set_theme(style="whitegrid")
fig, ax = plt.subplots(figsize=(4, 6))


sns.lineplot(
    data=df_long,
    x='Measurement',
    y='Z_Score',
    hue='Sample',
    units='Sample',
    estimator=None,
    marker='o',
    markersize=6,
    linewidth=1.0,
    alpha=0.7,
    legend=False,
    ax=ax
)



# 2. Shrink the x-axis group labels (Protocol_1, Ct_Value, Protocol_2)
ax.tick_params(axis='x', labelsize=7)
ax.tick_params(axis='y', labelsize=7)

ax.set_title('Swift: 10x Coverage vs Ct Value (Z-Score)', fontsize=10, pad=12)
ax.set_ylabel('Relative Variance (Standard Deviations from Mean)', fontsize=8, color="darkblue")
ax.set_xlabel('Expt', fontsize=8, color="darkblue")

# Add a horizontal line at 0 to represent the average for each group
plt.axhline(0, color='gray', linestyle='--', alpha=0.5)

plt.tight_layout()


plt.savefig('/media/simonray/data24/data/uio_dropbox_sr/DResearch/beathe/mybaits_vs_swift/data/zscore_trajectory_plot.png', dpi=300)