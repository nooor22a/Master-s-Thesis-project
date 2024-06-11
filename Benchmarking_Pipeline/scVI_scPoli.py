import numpy as np
import scanpy as sc
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

##adata is preprocessed 
adata = sc.read_h5ad('Final_object.h5ad')

custom_palette = ['#F8766D', '#D39200', '#93AA00', '#00BA38', '#00C19F',
                  '#00B9E3', '#619CFF', '#DB72FB', '#FF61C3']
##HVG selected, here subset = True. If you wanna keep all the genes, use subset =False.Later in scvi model, use adata=adata[:, adata.var["highly_variable"]].copy()

sc.pp.highly_variable_genes(
    adata,
    flavor="seurat_v3",
    n_top_genes=3000,
    layer="counts",
    batch_key="unifiedSampleID",
    subset = True
)   
###make sure to have a raw counts layer of your data before normalizing
#### here: adata.layers['counts'] = adata.X.copy()


scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="unifiedSampleID")
model = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")


#### you can delete these
total_cpu_usage_start = psutil.cpu_percent(interval=None)
total_mem_usage_start = psutil.virtual_memory().used
total_gpu_before = GPUtil.getGPUs() if torch.cuda.is_available() else None
start_time_scvi = time.time()
##### Train SCVI model


model.train()



# +
####### Calculate resources
elapsed_time_scvi = time.time() - start_time_scvi
cpu_end_scvi = psutil.cpu_percent(interval=None)
mem_end_scvi = psutil.virtual_memory().used



SCVI_LATENT_KEY = "X_scVI"
adata.obsm[SCVI_LATENT_KEY] = model.get_latent_representation()
sc.pp.neighbors(adata, use_rep=SCVI_LATENT_KEY)
sc.tl.leiden(adata)
sc.tl.umap(adata)
# -

sc.pl.umap(
    adata,
    color=["final_level"],
    frameon=False,
    save='lineagetrue.pdf',
    palette = custom_palette
)

sc.pl.umap(
    adata,
    color=["leiden"],
    frameon=False,
    save='leiden.pdf',
    palette = custom_palette
)
sc.pl.umap(
    adata,
    color=["final_level", "leiden"],
    frameon=False,
    save='lineagetrue.pdf',
    palette = custom_palette
)

""
###scPoli
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


print("for scpoli training: ")

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

adata.obsm["X_scpoli"] = scpoli_model.get_latent(adata,mean=True) 

""
##calc metrics
from scib_metrics.benchmark import Benchmarker
bm = Benchmarker(
    adata,
    batch_key="unifiedSampleID",
    label_key="final_level",
    embedding_obsm_keys=['X_pca', 'X_scVI', 'X_scpoli'],
    n_jobs=-1,
)
bm.benchmark()

df2 = bm.get_results(min_max_scale=False)
print(df2)
df2.to_csv('resultsintegration.csv', index=False)

