

df = adata.obs.copy()

df = df.rename(columns={'cell_type_pred_y': 'cell_type', 'Day_fixed': 'day_fixed', 'batch': 'batch', 'donor': 'donor'})


df['count'] = 1

# Calculate total counts per donor
df['total'] = df.groupby(['donor'])['count'].transform('sum')

# Calculate the ratio
df['ratio'] = df['count'] / df['total']

# Pivot the data for plotting
pivot_data = df.pivot_table(index='donor', columns='cell_type', values='ratio', aggfunc='sum', fill_value=0)

# Define the color palette based on the extracted colors
color_palette = {
    'Erythrocyte': (31/255, 119/255, 180/255),
    'Fibroblast': (174/255, 199/255, 232/255),
    'Glioblast': (255/255, 127/255, 14/255),
    'Immune': (255/255, 187/255, 120/255),
    'Neural crest': (44/255, 160/255, 44/255),
    'Neuroblast': (152/255, 223/255, 138/255),
    'Neuron': (214/255, 39/255, 40/255),
    'Neuronal IPC': (255/255, 152/255, 150/255),
    'Oligo': (148/255, 103/255, 189/255),
    'Placodes': (197/255, 176/255, 213/255),
    'Radial glia': (140/255, 86/255, 75/255),
    'Vascular': (196/255, 156/255, 148/255)
}


cell_types = pd.Index(list(pivot_data.columns)).unique()
colors = {cell_type: color_palette.get(cell_type, 'gray') for cell_type in cell_types}

# Plotting
fig, ax = plt.subplots(figsize=(12, 8))

pivot_data.plot(kind='bar', stacked=True, ax=ax, width=0.95, color=[colors[cell_type] for cell_type in pivot_data.columns])

ax.set_title('Cell Type Ratios by Donor')
ax.set_xlabel('Donor')
ax.set_ylabel('Cell Ratio')
ax.legend().set_visible(False)  # Hide individual legends

# Create a shared legend
handles = [plt.Line2D([0], [0], color=colors[cell_type], lw=4) for cell_type in cell_types]
labels = cell_types
fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol=3)

plt.tight_layout(rect=[0, 0.03, 1, 0.95])

# Save the plot as a PDF file
plt.savefig('cell_type_ratios_by_donor.pdf', dpi=300)

# Show the plot
plt.show()
