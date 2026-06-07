import pandas as pd
from scipy.stats import wilcoxon
from itertools import combinations

# 1. Load the TSV file
# sep='\t' is the critical parameter to read tab-separated files
df = pd.read_csv('/media/simonray/data24/data/uio_dropbox_sr/DResearch/beathe/mybaits_vs_swift/data/mean_genome_coverage.tsv', sep='\t')

# 2. Extract column names
cols = df.columns

# 3. Create combinations and run the test
# combinations(cols, 2) generates every possible unique pair of columns
results = []
for col1, col2 in combinations(cols, 2):
    # Perform the Wilcoxon test
    clean_pair = df[[col1, col2]].dropna()
    if len(clean_pair) < 1:
        print(f"Skipping {col1} vs {col2}: Not enough data after dropping NaNs.")
        continue

    # Perform the Wilcoxon test on the cleaned data
    stat, p_value = wilcoxon(clean_pair[col1], clean_pair[col2])

    results.append({
        'Column 1': col1,
        'Column 2': col2,
        'Statistic': stat,
        'P-Value': p_value
    })

# 4. Display results in a new DataFrame
results_df = pd.DataFrame(results)
print(results_df)