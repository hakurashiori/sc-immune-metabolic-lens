import palantir
import scanpy as sc
import pandas as pd
import os

# Plotting
import matplotlib
import matplotlib.pyplot as plt

# warnings
import warnings
from numba.core.errors import NumbaDeprecationWarning

warnings.filterwarnings(action="ignore", category=NumbaDeprecationWarning)
warnings.filterwarnings(
    action="ignore", module="scanpy", message="No data for colormapping"
)

# Inline plotting
%matplotlib inline

ad = sc.read_h5ad('./data/input.h5ad')

print(ad)

sc.pp.normalize_per_cell(ad)
palantir.preprocess.log_transform(ad)
sc.pp.highly_variable_genes(ad, n_top_genes=2500, flavor="seurat")
sc.pp.pca(ad)

# Run diffusion maps
dm_res = palantir.utils.run_diffusion_maps(ad, n_components=5)
ms_data = palantir.utils.determine_multiscale_space(ad)

sc.pp.neighbors(ad)
sc.tl.umap(ad)

# Use scanpy functions to visualize umaps or FDL
sc.pl.embedding(ad, basis="umap", frameon=False)

imputed_X = palantir.utils.run_magic_imputation(ad)

sc.pl.embedding(
    ad,
    basis="umap",
    layer="MAGIC_imputed_data",
    color=["JUN", "CCL3", "GAB2", "TMEM176B", "JDP2"],
    frameon=False,
)
plt.show()



# Define terminal states based on conditions
terminal_states = pd.Series(
    ["Pre-vaccination non-responder", "Pre-vaccination responder", "Healthy control"],
    index=["ID0NR1_CTGTGCTCAGCATGAG_4", "ID0R6_ACTTACTTCTTGTCAT-1_5", "H8_GGTCTGGTCTCAAAGC_1"],
)

# Check metadata columns
print(ad.obs.columns)
cluster_counts = ad.obs['condition2'].value_counts()
print(cluster_counts)

# Use Healthy control as early cell reference
palantir.utils.early_cell(ad, "Healthy control", celltype_column="condition2")

cells = ["H6_ATTGGTGTCCCTCTTT_1_1"]
palantir.plot.highlight_cells_on_umap(ad, cells)
plt.show()

sc.pl.umap(ad, color='condition2')
palantir.plot.highlight_cells_on_umap(ad, terminal_states)
plt.show()

# Run Palantir again with Healthy control start cell
start_cell = "H6_ATTGGTGTCCCTCTTT_1_1"
pr_res = palantir.core.run_palantir(
    ad, early_cell=start_cell, num_waypoints=1200, terminal_states=terminal_states, use_early_cell_as_start=True
)

palantir.plot.plot_palantir_results(ad, s=3)
plt.savefig("./results/output.png", dpi=300, bbox_inches='tight', transparent=True)
plt.show()

cells = ["H11_TTCAATCCACGTTCGG_1","ID0NR1_CTGTGCTCAGCATGAG_4", "ID0R6_ACTTACTTCTTGTCAT-1_5"]
palantir.plot.plot_terminal_state_probs(ad, cells)

fig = plt.gcf()
axes = fig.get_axes()
for ax in axes:
    plt.sca(ax)
    plt.xticks(rotation=45, ha='right')

plt.tight_layout()
plt.savefig("./results/output.png", dpi=300, bbox_inches='tight', transparent=True)
plt.show()

palantir.plot.highlight_cells_on_umap(ad, cells)
plt.savefig("./results/output.png", dpi=300, bbox_inches='tight', transparent=True)
plt.show()

masks = palantir.presults.select_branch_cells(ad, q=.01, eps=.01)
palantir.plot.plot_branch_selection(ad)
plt.tight_layout()
plt.show()



palantir.plot.plot_trajectory(
    ad,
    "Healthy control",
    cell_color="palantir_entropy",
    n_arrows=5,
    color="red",
    scanpy_kwargs=dict(cmap="viridis"),
    arrowprops=dict(arrowstyle="-|>,head_length=.5,head_width=.5"),
)
plt.savefig("./results/output.png", dpi=300, bbox_inches='tight', transparent=True)

# Compute gene trends
gene_trends = palantir.presults.compute_gene_trends(ad, expression_key="MAGIC_imputed_data")
genes = ["TNF"]
palantir.plot.plot_gene_trends(ad, genes)
plt.legend(loc='upper left', bbox_to_anchor=(1.05, 1), borderaxespad=0.)
plt.show()

palantir.plot.plot_gene_trend_heatmaps(ad, genes)
plt.show()

palantir.plot.plot_trend(ad, "Pre-vaccination non-responder", "ATP2B1",
                         color="n_counts", position_layer="MAGIC_imputed_data")
plt.show()

# Cluster gene trends
more_genes = ad.var_names[:1000]
communities = palantir.presults.cluster_gene_trends(ad, "Pre-vaccination non-responder", more_genes)

palantir.plot.plot_gene_trend_clusters(ad, "Pre-vaccination non-responder")
plt.show()

# Save processed object
file_path = os.path.join("./results", "output.h5ad")
ad.write(file_path)
