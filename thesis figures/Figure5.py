import numpy as np
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
df = adata.obs.groupby(["Observed", "Predicted"]).size().unstack(fill_value=0)
norm_df = df / df.sum(axis=0)
plt.figure(figsize=(12, 12))
ax = plt.gca()
c = ax.pcolor(norm_df)
plt.colorbar(c, ax=ax)
ax.set_xticks(np.arange(0.5, len(norm_df.columns), 1))
ax.set_yticks(np.arange(0.5, len(norm_df.index), 1))
ax.set_xticklabels(norm_df.columns, rotation=90)
ax.set_yticklabels(norm_df.index)
plt.ylabel("final_level")
plt.xlabel("twoPassAnnotation_clean")



for i in range(len(norm_df.index)):
    for j in range(len(norm_df.columns)):
        text = ax.text(j + 0.5, i + 0.5, int(df.iloc[i, j]),
                       ha="center", va="center", color="w", fontsize=7)  

plt.show()
