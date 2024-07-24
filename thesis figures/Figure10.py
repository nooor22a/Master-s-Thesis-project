import anndata as ad
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
adata =sc.read_h5ad('adata.h5ad')
# Extract necessary information from adata
df = adata.obs.copy()

# Rename columns to match your provided example
df = df.rename(columns={'cell_type_pred': 'cell_type', 'Day_fixed': 'day_fixed', 'batch': 'batch'})

# Add a count column assuming each row represents one cell
df['count'] = 1

# Calculate total counts per batch and day
df['total'] = df.groupby(['batch', 'day_fixed'])['count'].transform('sum')

# Calculate the ratio
df['ratio'] = df['count'] / df['total']

# Split the data by batch
df_fetal = df[df['batch'] == 'fetal']
df_2d = df[df['batch'] == '2d']
df_3d = df[df['batch'] == '3d']

# Aggregate the days into weeks for 3D batch
df_3d['week_fixed'] = (df_3d['day_fixed'] // 7) * 7

#  total counts per batch and week for 3D
df_3d['total'] = df_3d.groupby(['batch', 'week_fixed'])['count'].transform('sum')

# Calculate the ratio for 3D with weeks
df_3d['ratio'] = df_3d['count'] / df_3d['total']

pivot_fetal = df_fetal.pivot_table(index='day_fixed', columns='cell_type', values='ratio', aggfunc='sum', fill_value=0)
pivot_2d = df_2d.pivot_table(index='day_fixed', columns='cell_type', values='ratio', aggfunc='sum', fill_value=0)
pivot_3d = df_3d.pivot_table(index='week_fixed', columns='cell_type', values='ratio', aggfunc='sum', fill_value=0)

cmap = plt.get_cmap('tab20')
cell_types = pd.Index(list(pivot_fetal.columns) + list(pivot_2d.columns) + list(pivot_3d.columns)).unique()
colors = {cell_type: cmap(i) for i, cell_type in enumerate(cell_types)}

# Plotting
fig, axes = plt.subplots(1, 3, figsize=(20, 8), sharey=True)

# Plot each batch in separate subplots with reduced width for bars
for ax, (data, batch) in zip(axes, [(pivot_fetal, 'fetal'), (pivot_2d, '2d'), (pivot_3d, '3d')]):
    data.plot(kind='bar', stacked=True, ax=ax, width=0.95, color=[colors[cell_type] for cell_type in data.columns])
    ax.set_title(batch)
    ax.set_xlabel('Days' if batch != '3d' else 'Weeks')
    ax.set_ylabel('Cell Ratio')
    ax.legend().set_visible(False)  # Hide individual legends

# shared legend
handles = [plt.Line2D([0], [0], color=colors[cell_type], lw=4) for cell_type in cell_types]
labels = cell_types
fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol=3)

plt.tight_layout(rect=[0, 0.03, 1, 0.95])  
plt.show()

