##Harmonizing before Integration
import scanpy as sc
import numpy as np
import gc

ot = sc.read_h5ad('SnapseedAnnotated_OT.h5ad')
hnoca = sc.read_h5ad('SnapseedAnnotated_HNOCA.h5ad')
linnarson = sc.read_h5ad('SnapseedAnnotated_Linnarson.h5ad')


ot.obs_names = ot.obs_names.astype(str)
hnoca.obs_names = hnoca.obs_names.astype(str)
linnarson.obs_names = linnarson.obs_names.astype(str)

ot.obs_names_make_unique()
hnoca.obs_names_make_unique()
linnarson.obs_names_make_unique()

ot.var_names_make_unique()
hnoca.var_names_make_unique()
linnarson.var_names_make_unique()

##some cleaning, matching obs and metadata
ot.obs.rename(columns={'twoPassAnnotation_clean': 'original_labels'}, inplace=True)
hnoca.obs.drop(columns= "batch", inplace= True)
ot.obs['assay_sc'] = "10x 3' v3"
columns_to_drop = ['organism_original', 'disease', 'cell_type_original']
hnoca.obs = hnoca.obs.drop(columns=columns_to_drop)
linnarson.obs['robustID'] = linnarson.obs.index
ot.obs = ot.obs.drop(columns=['final_intermediate', 'leiden'])


##intersecting genes
common_genes = np.intersect1d(np.intersect1d(ot.var.index, linnarson.var.index), hnoca.var.index)
ot_harm = ot[:, np.isin(ot.var.index, common_genes)].copy()
linnarson_harm = linnarson[:, np.isin(linnarson.var.index, common_genes)].copy()
hnoca_harm = hnoca[:, np.isin(hnoca.var.index, common_genes)].copy()

##concat
adata = sc.concat(
    [ot_harm, hnoca_harm, linnarson_harm],
    label='batch',
    keys=['2d', '3d', 'fetal']
)

print('current adata')
adata
adata.X
adata.write("Object_concat_before_integration.h5ad")

##delete to save memory
del ot
del hnoca
del linnarson
del ot_harm
del hnoca_harm
del linnarson_harm
gc.collect()

# +
##Normalize the raw counts and add to adata.X
adata.layers['normalized'] = sc.pp.normalize_total(adata, target_sum=1e4, inplace=False, layer='counts')['X']


adata.layers['log_normalized'] = sc.pp.log1p(adata.layers['normalized'])


adata.X = adata.layers['log_normalized'].copy()
# -

if 'normalized' in adata.layers:
    del adata.layers['normalized']

##subset for Highlyvariable genes, 2000, can be edited
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="batch",subset= True)

##pca, D=50
sc.pp.pca(adata, n_comps=50)
