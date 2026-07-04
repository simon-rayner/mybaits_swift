import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats

# 1. Load the Excel files
df_mybaits = pd.read_excel('/media/simonray/data24/data/uio_dropbox_sr/DResearch/beathe/mybaits_vs_swift/SV_Re_Files_in_TSD_p1672_appendices/Appendix_myBaits.xlsx')
df_swift = pd.read_excel('/media/simonray/data24/data/uio_dropbox_sr/DResearch/beathe/mybaits_vs_swift/SV_Re_Files_in_TSD_p1672_appendices/Appendix_swift.xlsx')

# 2. Join the dataframes on a specific column
# Replace 'ID' with the name of your common column
merged_df_raw = pd.merge(df_mybaits, df_swift, on='Sample No', how='inner')
merged_df = merged_df_raw.drop(index=1).reset_index(drop=True)
merged_df['log_2_x'] = np.log2(merged_df['Input reads_x'])
merged_df['log_2_y'] = np.log2(merged_df['Input reads_y'])
sns.set_theme(style="whitegrid", palette="muted")

# 2. Initialize the figure
plt.figure(figsize=(10, 6))

slope_x, intercept_x, r_value_x, p_value_x, std_err_x = stats.linregress(merged_df['Ct-value_x'], merged_df['log_2_x'])
r_squared_x = r_value_x**2
ax = sns.regplot(
        data=merged_df,
        x='Ct-value_x',
        y='log_2_x',
        line_kws={'color': 'red', 'label': f'y = {slope_x:.2f}x + {intercept_x:.2f}'},
        scatter_kws={'alpha': 0.6})

plt.text(0.05, 0.92, f'$R^2 = {r_squared_x:.3f}$',
         transform=ax.transAxes, fontsize=12, verticalalignment='top',
         bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
plt.title('MyBaits: Ct versus Input Reads', fontsize=14)
plt.legend(loc='upper left', bbox_to_anchor=(0.01, 0.85))
plt.show()

# 1. Set the overall aesthetic theme


# 3. Create the scatter plot
# We add alpha (transparency) and s (size) for a more modern look
#ax = sns.scatterplot(
#    data=merged_df,
#    x='Ct-value_x',
#    y='Input reads_x',
#    #hue='category_column', # Optional: adds color based on a third column
#    palette='viridis',     # High-contrast, colorblind-friendly palette
#    alpha=0.7,             # Slight transparency helps with density
#    edgecolor='w',         # White outlines make points "pop"
#    s=100                  # Larger marker size
#)

# 4. Refine titles and labels
#plt.title('Analysis of X vs Y', fontsize=16, fontweight='bold', pad=20)
#plt.xlabel('X Axis Label', fontsize=12)
#plt.ylabel('Y Axis Label', fontsize=12)

# 5. Clean up the legend (if using 'hue')
#plt.legend(title='Categories', bbox_to_anchor=(1.05, 1), loc='upper left')

# 6. Final layout adjustment
plt.tight_layout()
plt.savefig('my_baits.png')
plt.show()
plt.close()

print("here")

slope_y, intercept_y, r_value_y, p_value_y, std_err_y = stats.linregress(merged_df['Ct-value_y'], merged_df['log_2_y'])
r_squared_y = r_value_y**2
ax = sns.regplot(
        data=merged_df,
        x='Ct-value_y',
        y='log_2_y',
        line_kws={'color': 'red', 'label': f'y = {slope_y:.2f}x + {intercept_y:.2f}'},
        scatter_kws={'alpha': 0.6})

plt.text(0.05, 0.92, f'$R^2 = {r_squared_y:.3f}$',
         transform=ax.transAxes, fontsize=12, verticalalignment='top',
         bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
plt.title('Swift: Ct versus Input Reads', fontsize=14)
plt.legend(loc='upper left', bbox_to_anchor=(0.01, 0.85))
plt.show()



# 6. Final layout adjustment
plt.tight_layout()
plt.savefig('swift_plot.png')
plt.show()

print("here")
