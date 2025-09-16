import seaborn as sns

import matplotlib
import matplotlib.pyplot as plt

print(sns.__version__)

print(matplotlib.__version__)

print(sc.__version__)

import scanpy as sc
import scvi

adata = sc.read_h5ad('.Results/monocytes.h5ad')

# Basic filtering and normalization
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Set up the AnnData object for scvi-tools (specify batch key if needed)
scvi.model.SCVI.setup_anndata(adata, batch_key='MT_CellType')  # Replace with your batch column or remove if none

# Initialize the model
model = scvi.model.SCVI(adata)

# Train the model
model.train(max_epochs=400)  # Adjust epochs based on dataset size

# Extract the latent space
adata.obsm["X_scVI"] = model.get_latent_representation()

# Run the clustering algorithm (Leiden or Louvain)
sc.pp.neighbors(adata, use_rep="X_scVI")  # Neighbors on latent space
sc.tl.leiden(adata, resolution=0.5)       # Adjust resolution for cluster granularity

# Compute UMAP for visualization
sc.tl.umap(adata)
sc.pl.umap(adata, color=['leiden'])  # Use 'leiden' or the column storing your cluster labels

sc.pl.umap(adata, color=['condition2'])

