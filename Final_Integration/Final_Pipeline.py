import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scvi
import seaborn as sns
import torch
from rich import print
import time
import psutil
import scvi
import GPUtil
import torch
from scib_metrics.benchmark import Benchmarker

adata = sc.read_h5ad("HarmonizedObject.h5ad')

##Basic preprocessing:
                     
# +
##Normalize the raw counts and add to adata.X
adata.layers['normalized'] = sc.pp.normalize_total(adata, target_sum=1e4, inplace=False, layer='counts')['X']


adata.layers['log_normalized'] = sc.pp.log1p(adata.layers['normalized'])


adata.X = adata.layers['log_normalized'].copy()
# -

if 'normalized' in adata.layers:
    del adata.layers['normalized']

##subset for Highlyvariable genes
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="batch",subset= True)

##pca, D=50
sc.pp.pca(adata, n_comps=50)

custom_palette = ['#F8766D', '#D39200', '#93AA00', '#00BA38', '#00C19F',
                  '#00B9E3', '#619CFF', '#DB72FB', '#FF61C3']
##scvi integration

scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="batch")
model = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")



##script used to measure computational resources
total_cpu_usage_start = psutil.cpu_percent(interval=None)
total_mem_usage_start = psutil.virtual_memory().used
total_gpu_before = GPUtil.getGPUs() if torch.cuda.is_available() else None
start_time_scvi = time.time()

model.train()

elapsed_time_scvi = time.time() - start_time_scvi
cpu_end_scvi = psutil.cpu_percent(interval=None)
mem_end_scvi = psutil.virtual_memory().used
mem_usage_diff = mem_end_scvi - total_mem_usage_start
# +
print("for scvi training: ")

print(f"Elapsed time for model training: {elapsed_time_scvi} seconds")
print(f"CPU usage at start: {total_cpu_usage_start}%")
print(f"CPU usage after training: {cpu_end_scvi}%")
print(f"Memory usage at start: {total_mem_usage_start / (1024 ** 3):.2f} GB")
print(f"Memory usage after training: {mem_end_scvi / (1024 ** 3):.2f} GB")
print(f"Memory usage difference: {mem_usage_diff / (1024 ** 3):.2f} GB")

if total_gpu_before:
    total_gpu_after = GPUtil.getGPUs()
    for i, gpu in enumerate(total_gpu_after):
        print(f"GPU {i} usage before: {total_gpu_before[i].memoryUsed / 1024:.2f} GB")
        print(f"GPU {i} usage after: {gpu.memoryUsed / 1024:.2f} GB")
        print(f"GPU {i} usage difference: {(gpu.memoryUsed - total_gpu_before[i].memoryUsed) / 1024:.2f} GB")
        


# +
SCVI_LATENT_KEY = "X_scVI"
adata.obsm[SCVI_LATENT_KEY] = model.get_latent_representation()

adata.write("integratedwith_XPCA_SCVI.h5ad")

# +
sc.pp.neighbors(adata, use_rep=SCVI_LATENT_KEY)
sc.tl.leiden(adata)
sc.tl.umap(adata)
### code used to save umaps generated

sc.pl.umap(
    adata,
    color=["batch"],
    frameon=False,
    save='batch.pdf',
    palette = custom_palette
)

##scPoli integration
from scarches.models.scpoli import scPoli
import torch

# +
scpoli_model = scPoli(
    adata=adata,
    unknown_ct_names=["nan"],
    condition_keys="batch",
    cell_type_keys=["level_1", "level_2", "final_level"],
    embedding_dims=5
)

early_stopping_kwargs = {
    "early_stopping_metric": "val_prototype_loss",
    "mode": "min",
    "threshold": 0,
    "patience": 20,
    "reduce_lr": True,
    "lr_patience": 13,
    "lr_factor": 0.1,
}
total_cpu_usage_start = psutil.cpu_percent(interval=None)
total_mem_usage_start = psutil.virtual_memory().used
total_gpu_before = GPUtil.getGPUs() if torch.cuda.is_available() else None
start_time_scvi = time.time()

scpoli_model.train(
    unlabeled_prototype_training=False,
    n_epochs=7,
    pretraining_epochs=5,
    early_stopping_kwargs=early_stopping_kwargs,
    eta=10,
    alpha_epoch_anneal=100
)
elapsed_time_scvi = time.time() - start_time_scvi
cpu_end_scvi = psutil.cpu_percent(interval=None)
mem_end_scvi = psutil.virtual_memory().used
mem_usage_diff = mem_end_scvi - total_mem_usage_start


print("for scvi training: ")

print(f"Elapsed time for model training: {elapsed_time_scvi} seconds")
print(f"CPU usage at start: {total_cpu_usage_start}%")
print(f"CPU usage after training: {cpu_end_scvi}%")
print(f"Memory usage at start: {total_mem_usage_start / (1024 ** 3):.2f} GB")
print(f"Memory usage after training: {mem_end_scvi / (1024 ** 3):.2f} GB")
print(f"Memory usage difference: {mem_usage_diff / (1024 ** 3):.2f} GB")

if total_gpu_before:
    total_gpu_after = GPUtil.getGPUs()
    for i, gpu in enumerate(total_gpu_after):
        print(f"GPU {i} usage before: {total_gpu_before[i].memoryUsed / 1024:.2f} GB")
        print(f"GPU {i} usage after: {gpu.memoryUsed / 1024:.2f} GB")
        print(f"GPU {i} usage difference: {(gpu.memoryUsed - total_gpu_before[i].memoryUsed) / 1024:.2f} GB")




scpoli_model.save("modelscpoli")
scpoli_model.get_conditional_embeddings().write_h5ad(
    "sample_embedding.h5ad")

adata.obsm["X_scpoli"] = scpoli_model.get_latent(
    adata,
    mean=True
)

# -

###view
sc.pp.neighbors(adata, use_rep='X_scpoli')
sc.tl.leiden(adata)
sc.tl.umap(adata)

custom_palette = ['#F8766D', '#D39200', '#93AA00', '#00BA38', '#00C19F',
                  '#00B9E3', '#619CFF', '#DB72FB', '#FF61C3']


# +
sc.pl.umap(
    adata,
    color=["final_level"],
    frameon=False,
    save='x_scpoli_final_level.pdf',
    palette = custom_palette
)


sc.pl.umap(
    adata,
    color=["leiden"],
    frameon=False,
    save='sscpol_leiden.pdf',
    palette = custom_palette
)

sc.pl.umap(
    adata,
    color=["level_1"],
    frameon=False,
    save='sscpol_level1.pdf',
    palette = custom_palette
)
sc.pl.umap(
    adata,
    color=["level_2"],
    frameon=False,
    save='sscpol_level2.pdf',
    palette = custom_palette
)

sc.pl.umap(
    adata,
    color=["batch"],
    frameon=False,
    save='scpolibatch.pdf',
    palette = custom_palette
)



print("UMAP plots saved successfully.")

adata.write("final_scpoli.h5ad")
