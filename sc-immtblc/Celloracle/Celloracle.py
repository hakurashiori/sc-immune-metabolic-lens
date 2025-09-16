import celloracle as co
print(co.__version__)

import matplotlib.pyplot as plt

%matplotlib inline

# Instantiate Oracle object
oracle = co.Oracle()

import scanpy as sc

import numpy as np


# Load raw count matrix
adata = sc.read_h5ad("raw_count.h5ad")
raw_backup = adata.X.copy()
full_gene_names = adata.var_names.copy()

# Store raw counts in a separate layer and normalize/log-transform
adata.layers["raw_count"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4, layer="raw_count")
sc.pp.log1p(adata, layer="raw_count")

# Select highly variable genes
sc.pp.highly_variable_genes(adata, n_top_genes=4000, flavor="seurat", layer="raw_count")
adata = adata[:, adata.var.highly_variable].copy()

# Map variable gene indices to raw backup
var_idx = full_gene_names.get_indexer(adata.var_names)

adata.X = raw_backup[:, var_idx].copy()


for g in ["FOXO1", "SERPB1", "TCF7L2", "ATF4", "RUX2"]:
# Check presence of selected genes
    print(f"{g} in adata.var_names:", g in adata.var_names)

# Check data in anndata
# Check metadata structure
print("Metadata columns :", list(adata.obs.columns))
print("Dimensional reduction: ", list(adata.obsm.keys()))

# Instantiate Oracle object
oracle = co.Oracle()

oracle.import_anndata_as_raw_count(adata=adata,
                                   cluster_column_name="condition2",
                                   embedding_name="X_umap")

# Load TF info which was made from mouse cell atlas dataset.
# Load human promoter-based reference GRN
base_GRN = co.data.load_human_promoter_base_GRN()

# Check data
base_GRN.head()

# You can load TF info dataframe with the following code.
oracle.import_TF_data(TF_info_matrix=base_GRN)


# Perform PCA
oracle.perform_PCA()

# Select important PCs
# Plot cumulative explained variance
plt.plot(np.cumsum(oracle.pca.explained_variance_ratio_)[:100])
# Plot cumulative explained variance
n_comps = np.where(np.diff(np.diff(np.cumsum(oracle.pca.explained_variance_ratio_))>0.002))[0][0]
plt.axvline(n_comps, c="k")
plt.show()
print(n_comps)
n_comps = min(n_comps, 50)

n_cell = oracle.adata.shape[0]
print(f"cell number is :{n_cell}")

k = int(0.025*n_cell)
print(f"Auto-selected k is :{k}")

# Perform KNN imputation
oracle.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True, b_sight=k*8,
                      b_maxl=k*4, n_jobs=4)


# Save Oracle object to file
oracle.to_hdf5("data.celloracle.oracle")

# Load file.
#oracle = co.load_hdf5("data.celloracle.oracle")

# Construct neighborhood graph
sc.pp.neighbors(oracle.adata, use_rep="X_pca")

# Perform graph layout embedding
sc.tl.draw_graph(oracle.adata)

# Check clustering data
# Perform graph layout embedding
# Visualize embedding by condition
# Visualize gene expression on embedding
sc.pl.draw_graph(oracle.adata, color="condition2", show=False)
plt.savefig("./results/output.png", dpi=300, bbox_inches="tight")
plt.show()

%%time
# Calculate GRN for each population in clustering unit.
# This step may take some time.(~30 minutes)
# Infer gene regulatory network per cluster
links = oracle.get_links(cluster_name_for_GRN_unit="condition2", alpha=10,
                         verbose_level=10)

links.links_dict.keys()

links.links_dict["Pre-vaccination non-responder"]

# Save Links object.
# Save Oracle object to file
links.to_hdf5(file_path="links.celloracle.links")

# Filter GRN links
links.filter_links(p=0.001, weight="coef_abs", threshold_number=2000)

plt.rcParams["figure.figsize"] = [9, 4.5]

# Plot degree distributions
links.plot_degree_distributions(plot_model=True)

# Calculate network scores
links.get_network_score()

links.merged_score.head()

# Save Links object.
# Save Oracle object to file
links.to_hdf5(file_path="links.celloracle.links")

# You can load files with the following command.
#links = co.load_hdf5(file_path="links.celloracle.links")

# Filter GRN links
links.filter_links()
oracle.get_cluster_specific_TFdict_from_Links(links_object=links)
# Fit GRN for simulation
oracle.fit_GRN_for_simulation(alpha=10,
                              use_cluster_specific_TFdict=True)

for g in ["E2F8", "FOXO1", "TCF7L2", "HNF4A", "USF1", "SP1"]:
# Check presence of selected genes
    print(f"{g} in adata.var_names:", g in adata.var_names)

# Check gene expression
goi = "FOXO1"
# Perform graph layout embedding
# Visualize embedding by condition
# Visualize gene expression on embedding
sc.pl.draw_graph(oracle.adata, color=[goi, oracle.cluster_column_name],
                 layer="imputed_count", use_raw=False, cmap="viridis", show=False)

#plt.savefig("./results/output.png", dpi=300, bbox_inches="tight")


plt.show()

# Plot gene expression in histogram
# Plot histogram of gene expression
sc.get.obs_df(oracle.adata, keys=[goi], layer="imputed_count").hist()
plt.show()

# Enter perturbation conditions to simulate signal propagation after the perturbation.
# Simulate perturbation by gene knockout
oracle.simulate_shift(perturb_condition={goi: 0.0},
                      n_propagation=3)

# Get transition probability
# Construct neighborhood graph
# Estimate transition probabilities
oracle.estimate_transition_prob(n_neighbors=200,
                                knn_random=True,
                                sampled_fraction=1)

# Calculate embedding
# Calculate embedding shift after perturbation
oracle.calculate_embedding_shift(sigma_corr=0.05)

print(oracle.adata.obsm.keys())

# Perform graph layout embedding
#oracle.adata.obsm["X_umap"] = oracle.adata.obsm["X_draw_graph_fr"]

fig, ax = plt.subplots(1, 2,  figsize=[13, 6])

scale = 10
# Show quiver plot
# Plot quiver showing simulated vs randomized vector field
oracle.plot_quiver(scale=scale, ax=ax[0])
ax[0].set_title(f"Simulated cell identity shift vector: {goi} KO")

# Show quiver plot that was calculated with randomized graph.
# Plot quiver showing simulated vs randomized vector field
oracle.plot_quiver_random(scale=scale, ax=ax[1])
ax[1].set_title(f"Randomized simulation vector")

fig.savefig("./results/output.png", dpi=300, bbox_inches="tight")

plt.show()

# n_grid = 40 is a good starting value.
n_grid = 40
# Construct neighborhood graph
# Calculate probability mass for simulation
oracle.calculate_p_mass(smooth=0.8, n_grid=n_grid, n_neighbors=200)

# Search for best min_mass.
# Suggest thresholds for mass filtering
oracle.suggest_mass_thresholds(n_suggestion=12)

min_mass = 1.5
# Apply mass filter
oracle.calculate_mass_filter(min_mass=min_mass, plot=True)

fig, ax = plt.subplots(1, 2,  figsize=[13, 6])

scale_simulation = 10
# Show quiver plot
# Plot simulated vector field on grid
oracle.plot_simulation_flow_on_grid(scale=scale_simulation, ax=ax[0])
ax[0].set_title(f"Simulated cell identity shift vector: {goi} KO")

# Show quiver plot that was calculated with randomized graph.
oracle.plot_simulation_flow_random_on_grid(scale=scale_simulation, ax=ax[1])
ax[1].set_title(f"Randomized simulation vector")

fig.savefig("./results/output.png", dpi=300, bbox_inches="tight")

plt.show()

# Plot vector field with cell cluster
fig, ax = plt.subplots(figsize=[8, 8])

# Plot clusters with simulation overlay
oracle.plot_cluster_whole(ax=ax, s=10)
# Plot simulated vector field on grid
oracle.plot_simulation_flow_on_grid(scale=scale_simulation, ax=ax, show_background=False)

fig.savefig("./results/output.png", dpi=300, bbox_inches="tight")

plt.show()

# Perform pseudotime analysis
from celloracle.applications import Pseudotime_calculator
# Perform pseudotime analysis
pt = Pseudotime_calculator(oracle_object=oracle)
print("Clustering name: ", pt.cluster_column_name)
print("Cluster list", pt.cluster_list)

# Check data
pt.plot_cluster(fontsize=8)



clusters_in_lineage = ['Healthy control', 'Pre-vaccination non-responder']

# Make a dictionary
lineage_dictionary = {"Lineage": clusters_in_lineage}

# Input lineage information into pseudotime object
# Define lineage dictionary
pt.set_lineage(lineage_dictionary=lineage_dictionary)

# Visualize lineage information
pt.plot_lineages()

# Estimated root cell name for each lineage
root_cells = {"Lineage": "H_c5_TGAAAGAAGTCAATAG-1_2"}
#H_c7_TTCGGTCCACAGCGTC-1_2
# Define root cell for pseudotime
pt.set_root_cells(root_cells=root_cells)

# Check root cell and lineage
pt.plot_root_cells()

#import scanpy as sc

# Construct neighborhood graph
sc.pp.neighbors(pt.adata, n_neighbors=15, use_rep="X_umap")

# Compute diffusion map
sc.tl.diffmap(pt.adata, n_comps=15, neighbors_key=None, random_state=0, copy=False)

# Check diffusion map data.
"X_diffmap" in pt.adata.obsm


# Calculate pseudotime
pt.get_pseudotime_per_each_lineage()

# Check results
# Plot pseudotime trajectory
pt.plot_pseudotime(cmap="rainbow")

# Check result
pt.adata.obs[["Pseudotime"]].head()

# Visualize pseudotime
fig, ax = plt.subplots(figsize=[6,6])

sc.pl.embedding(adata=pt.adata, basis=oracle.embedding_name, ax=ax, cmap="rainbow",
                color=["Pseudotime"])

fig.savefig("./results/output.png", dpi=300, bbox_inches="tight")

plt.show()


oracle.adata.obs["Pseudotime"] = pt.adata.obs["Pseudotime"]

# Gradient analysis along pseudotime
from celloracle.applications import Gradient_calculator

# Instantiate Gradient calculator object
# Gradient analysis along pseudotime
gradient = Gradient_calculator(oracle_object=oracle, pseudotime_key="Pseudotime")

# Construct neighborhood graph
# Calculate probability mass for simulation
gradient.calculate_p_mass(smooth=0.8, n_grid=n_grid, n_neighbors=200)
# Apply mass filter
gradient.calculate_mass_filter(min_mass=min_mass, plot=True)

gradient.transfer_data_into_grid(args={"method": "polynomial", "n_poly":3}, plot=True)


# Calculate gradients
gradient.calculate_gradient()

# Show results
scale_dev = 40
# Visualize gradient results
gradient.visualize_results(scale=scale_dev, s=5)

plt.tight_layout()  # 避免圖擠在一起

# 儲存圖片
plt.savefig("./results/output.png", dpi=300, bbox_inches="tight")

# （選擇性）顯示圖
plt.show()

# Visualize results
fig, ax = plt.subplots(figsize=[6, 6])
gradient.plot_dev_flow_on_grid(scale=scale_dev, ax=ax)

# Development module to compare reference vs perturbation
from celloracle.applications import Oracle_development_module

# Development module to compare reference vs perturbation
# Make Oracle_development_module to compare two vector field
# Development module to compare reference vs perturbation
dev = Oracle_development_module()

# Load development flow
dev.load_differentiation_reference_data(gradient_object=gradient)

# Load simulation result
dev.load_perturb_simulation_data(oracle_object=oracle)


# Calculate inner produc scores
# Calculate similarity between flows
dev.calculate_inner_product()
dev.calculate_digitized_ip(n_bins=10)

# Show perturbation scores
vm = 0.02

fig, ax = plt.subplots(1, 2, figsize=[12, 6])
# Plot inner product map
dev.plot_inner_product_on_grid(vm=0.02, s=50, ax=ax[0])
ax[0].set_title(f"PS")

dev.plot_inner_product_random_on_grid(vm=vm, s=50, ax=ax[1])
ax[1].set_title(f"PS calculated with Randomized simulation vector")


fig.savefig("./results/output.png", dpi=300, bbox_inches="tight")

plt.show()

# Show perturbation scores with perturbation simulation vector field
fig, ax = plt.subplots(figsize=[6, 6])
# Plot inner product map
dev.plot_inner_product_on_grid(vm=vm, s=60, ax=ax)
# Plot simulated vector field on grid
dev.plot_simulation_flow_on_grid(scale=scale_simulation, show_background=False, ax=ax)

fig.savefig("./results/output.png", dpi=300, bbox_inches="tight")

plt.show()

# Let's visualize the results
# Visualize combined developmental module layout
dev.visualize_development_module_layout_0(s=5,
                                          scale_for_simulation=scale_simulation,
                                          s_grid=50,
                                          scale_for_pseudotime=scale_dev,
                                          vm=vm)

plt.savefig("./results/output.png", dpi=300, bbox_inches="tight")

plt.show()



