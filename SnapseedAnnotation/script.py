#!/usr/bin/env python
# coding: utf-8
# %%


import snapseed as snap


# %%


from snapseed.utils import read_yaml

# %%
import jax
print(jax.devices())
###check if jax recognizes cuda to use GPUs 
# %%


snap.annotate


# %%
import scanpy as sc
##save raw counts in a seperate layer, normalize the counts, and run leiden clustering, for our analysis we ran sc.tl.leiden(adata, resolution=1.0)
adata = sc.read_h5ad('leidenclustered.h5ad')



# %%

##read the marker genes
marker_genes = read_yaml("Data_S1_snapseed_markers.yaml")


# %%


marker_genes



# %%
snap_annot = snap.annotate_hierarchy(
    adata,
    marker_genes,
    group_name="leiden",
    layer="counts",
)


## printing snap_annot will show a list of the annotations and their confidence values, save as needed
##adding annotations
adata.obs = adata.obs.join(snap_annot["assignments"], on=f"leiden")
print(adata.obs["level_1"].value_counts())
print(adata.obs["level_2"].value_counts())
print(adata.obs["level_3"].value_counts())


adata.obs['level_1'] = adata.obs['level_1'].astype(str)
adata.obs['level_2'] = adata.obs['level_2'].astype(str)
adata.obs['level_3'] = adata.obs['level_3'].astype(str)
adata.obs['level_4'] = adata.obs['level_4'].astype(str)
adata.obs['level_5'] = adata.obs['level_5'].astype(str)
adata.obs['final_level'] = adata.obs['level_3'].astype(str)


adata.write('annotated.h5ad')
