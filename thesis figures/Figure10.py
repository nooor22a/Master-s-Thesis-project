
##Figure 10 A
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


##Figure 10B
import anndata as ad
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

adata_new = adata_2d[adata_2d.obs['cell_type_pred'] == 'Fibroblast'].copy()

df = adata_new.obs.copy()

df = df.rename(columns={'cell_type_pred': 'cell_type_pred', 'Day_fixed': 'day_fixed', 'level21': 'level21', 'org': 'org'})

df['count'] = 1

# Pivot the data for plotting
pivot_cell_type_pred = df.pivot_table(index='day_fixed', columns='cell_type_pred', values='count', aggfunc='sum', fill_value=0)
pivot_level21 = df.pivot_table(index='day_fixed', columns='level21', values='count', aggfunc='sum', fill_value=0)
pivot_org = df.pivot_table(index='day_fixed', columns='org', values='count', aggfunc='sum', fill_value=0)

# Define color maps with enough distinct colors
cmap_tab20 = plt.get_cmap('tab20')
cmap_tab20b = plt.get_cmap('tab20b')
colors_list = [cmap_tab20(i) for i in range(20)] + [cmap_tab20b(i) for i in range(3)]

# Create color dictionaries
cell_type_pred_colors = {cell_type: colors_list[i % len(colors_list)] for i, cell_type in enumerate(pivot_cell_type_pred.columns)}
level21_colors = {cell_type: colors_list[i % len(colors_list)] for i, cell_type in enumerate(pivot_level21.columns)}
org_colors = {org: colors_list[i % len(colors_list)] for i, org in enumerate(pivot_org.columns)}

# Plotting
fig, axes = plt.subplots(1, 3, figsize=(24, 8), sharey=True)

# Plot cell_type_pred
pivot_cell_type_pred.plot(kind='bar', stacked=True, ax=axes[0], color=[cell_type_pred_colors[cell_type] for cell_type in pivot_cell_type_pred.columns], width=0.9)
axes[0].set_title('Cell Type scPoli Prediction Distribution Over Days')
axes[0].set_xlabel('Days')
axes[0].set_ylabel('Cell Count')
axes[0].legend(loc='upper right')

# Plot snapseed annotations
pivot_level21.plot(kind='bar', stacked=True, ax=axes[1], color=[level21_colors[cell_type] for cell_type in pivot_level21.columns], width=0.9)
axes[1].set_title('Snapseed Cell Type Distribution Over Days')
axes[1].set_xlabel('Days')
axes[1].set_ylabel('Cell Count')
axes[1].legend(loc='upper right')

# Plot org
pivot_org.plot(kind='bar', stacked=True, ax=axes[2], color=[org_colors[org] for org in pivot_org.columns], width=0.9)
axes[2].set_title('Previous Distribution Over Days')
axes[2].set_xlabel('Days')
axes[2].set_ylabel(' Cell Count')
axes[2].legend(loc='upper right')

plt.tight_layout()
plt.show()

