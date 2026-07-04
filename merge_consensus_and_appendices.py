import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats

#/media/simonray/data24/data/uio_dropbox_sr/DResearch/beathe/mybaits_vs_swift/lineage_report__pangolin_anonymised.xlsx
# /media/simonray/data24/data/uio_dropbox_sr/DResearch/beathe/mybaits_vs_swift/SV_Re_Files_in_TSD_p1672_appendices/Appendix_myBaits_metrics.xlsx
# /media/simonray/data24/data/uio_dropbox_sr/DResearch/beathe/mybaits_vs_swift/SV_Re_Files_in_TSD_p1672_appendices/Appendix_Swift_metrics.xlsx
# The "One" DataFrame (Each store ID appears exactly once)


df_mybaits = pd.read_excel('/media/simonray/data24/data/uio_dropbox_sr/DResearch/beathe/mybaits_vs_swift/SV_Re_Files_in_TSD_p1672_appendices/Appendix_myBaits.xlsx')
df_swift = pd.read_excel('/media/simonray/data24/data/uio_dropbox_sr/DResearch/beathe/mybaits_vs_swift/SV_Re_Files_in_TSD_p1672_appendices/Appendix_swift.xlsx')
merged_df_raw = pd.merge(df_mybaits, df_swift, on='Sample No', how='inner')
df_merged_appendices = merged_df_raw.drop(index=1).reset_index(drop=True)
df_merged_appendices['log_2_x'] = np.log2(df_merged_appendices['Input reads_x'])
df_merged_appendices['log_2_y'] = np.log2(df_merged_appendices['Input reads_y'])

df_pangolin = pd.read_excel('/media/simonray/data24/data/uio_dropbox_sr/DResearch/beathe/mybaits_vs_swift/data/lineage_report__pangolin_anonymised.xlsx')


df_merged_final= pd.merge(df_merged_appendices, df_pangolin, left_on='Sample No',  # Column name in df_employees
                          right_on='Sample')         # Column name in df_store)

print(len(df_merged_final))

df_merged_final.to_excel("/media/simonray/data24/data/uio_dropbox_sr/DResearch/beathe/mybaits_vs_swift/data/lineage_report__pangolin_anonymised_merged.xlsx", index=False)


