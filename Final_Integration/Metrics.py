import scanpy as sc
from scib_metrics.benchmark import Benchmarker

# Load your .h5ad file
adata = sc.read_h5ad('finalobject.h5ad')

# Initialize the Benchmarker
bm = Benchmarker(
    adata,
    batch_key="batch",
    label_key="final_level",
    embedding_obsm_keys=['X_pca', 'X_scVI', "X_scpoli"],
    n_jobs=-1,
)

# Run the benchmark
bm.benchmark()

# Get and print the results
df2 = bm.get_results(min_max_scale=False)
print(df2)

# Save the results to a CSV file
df2.to_csv('resultsintegration.csv', index=False)
