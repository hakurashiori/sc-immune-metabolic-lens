# ============================================================
# hdWGCNA workflow (Seurat + hdWGCNA)
# - Metacell-based coexpression network
# - Module eigengenes + connectivity + hub genes
# - Optional: export module genes
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(patchwork)
  library(ggplot2)
})

# -----------------------
# 0) User parameters
# -----------------------
SEURAT_RDS_IN  <- "pbmc.rds"                
SEURAT_RDS_OUT <- "pbmc_hdwgcna.rds"

WGCNA_NAME     <- "NR_hdwgcna"

GROUP_BY       <- "condition2"             # meta.data 
GROUP_NAME     <- "Pre-vaccination non-responder"

REDUCTION_USE  <- "UMAP_harmony"           # metacell KNN reduction
METACELL_K     <- 25
MAX_SHARED     <- 10

GENE_SELECT    <- "fraction"               # "fraction" recommended
FRACTION       <- 0.05                     # gene expressed in >= 5% cells

NETWORK_TYPE   <- "signed"                 # "signed" / "unsigned" / "signed hybrid"
TOM_NAME       <- "NR"

MODULE_PREFIX  <- "NR-M"
EXCLUDE_MODULE <- c("grey")                

N_HUBS         <- 10

# (Optional) export some modules
EXPORT_MODULES <- c("NR-M2", "NR-M3")
OUTDIR_FIG     <- "output/figures"
OUTDIR_TAB     <- "output/tables"

dir.create(OUTDIR_FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTDIR_TAB, recursive = TRUE, showWarnings = FALSE)

# -----------------------
# 1) Load Seurat object
# -----------------------
seurat_obj <- readRDS(SEURAT_RDS_IN)

stopifnot(GROUP_BY %in% colnames(seurat_obj@meta.data))
if (!REDUCTION_USE %in% names(seurat_obj@reductions)) {
  stop(paste0("Reduction not found: ", REDUCTION_USE,
              " (available: ", paste(names(seurat_obj@reductions), collapse = ", "), ")"))
}

# -----------------------
# 2) Setup hdWGCNA gene set
# -----------------------
seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = GENE_SELECT,
  fraction    = FRACTION,
  wgcna_name  = WGCNA_NAME
)

# -----------------------
# 3) Construct metacells by groups
# -----------------------
seurat_obj <- MetacellsByGroups(
  seurat_obj   = seurat_obj,
  group.by     = c(GROUP_BY),
  reduction    = REDUCTION_USE,
  k            = METACELL_K,
  max_shared   = MAX_SHARED,
  ident.group  = GROUP_BY
)

# -----------------------
# 4) Normalize metacells
# -----------------------
seurat_obj <- NormalizeMetacells(seurat_obj)

# -----------------------
# 5) Set expression matrix for target group
# -----------------------
seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = GROUP_NAME,
  group.by   = GROUP_BY,
  assay      = "RNA",
  layer      = "data"
)

# -----------------------
# 6) Choose soft power (diagnostics)
# -----------------------
seurat_obj <- TestSoftPowers(seurat_obj, networkType = NETWORK_TYPE)

plot_list <- PlotSoftPowers(seurat_obj)
p_soft <- wrap_plots(plot_list, ncol = 2)
ggsave(file.path(OUTDIR_FIG, "soft_powers.png"), p_soft, width = 10, height = 6, dpi = 300)

power_table <- GetPowerTable(seurat_obj)
write.csv(power_table, file.path(OUTDIR_TAB, "soft_power_table.csv"), row.names = FALSE)

# -----------------------
# 7) Construct co-expression network + TOM
# -----------------------
seurat_obj <- ConstructNetwork(
  seurat_obj,
  tom_name = TOM_NAME
)

# Dendrogram
p_dendro <- PlotDendrogram(seurat_obj, main = paste0(GROUP_NAME, " hdWGCNA Dendrogram"))
ggsave(file.path(OUTDIR_FIG, "dendrogram.png"), p_dendro, width = 10, height = 6, dpi = 300)

# -----------------------
# 8) Module eigengenes (harmonized) for ALL cells
# -----------------------
seurat_obj <- ModuleEigengenes(
  seurat_obj,
  group.by.vars = GROUP_BY
)

# -----------------------
# 9) Connectivity (kME) within target group
# -----------------------
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by   = GROUP_BY,
  group_name = GROUP_NAME
)

# Rename modules to "NR-M*"
seurat_obj <- ResetModuleNames(seurat_obj, new_name = MODULE_PREFIX)

# Plot kME ranks
p_kme <- PlotKMEs(seurat_obj, ncol = 5)
ggsave(file.path(OUTDIR_FIG, "kme_ranks.png"), p_kme, width = 12, height = 8, dpi = 300)

# -----------------------
# 10) Export module table + hub genes
# -----------------------
modules <- GetModules(seurat_obj)
modules <- modules %>% filter(!(module %in% EXCLUDE_MODULE))
write.csv(modules, file.path(OUTDIR_TAB, "module_assignments.csv"), row.names = FALSE)

hub_df <- GetHubGenes(seurat_obj, n_hubs = N_HUBS)
write.csv(hub_df, file.path(OUTDIR_TAB, paste0("hub_genes_top", N_HUBS, ".csv")), row.names = FALSE)

# -----------------------
# 11) Export gene lists for selected modules (optional)
# -----------------------
if (length(EXPORT_MODULES) > 0) {
  for (m in EXPORT_MODULES) {
    genes_m <- modules$gene[modules$module == m]
    out_m <- file.path(OUTDIR_TAB, paste0(m, "_genes.csv"))
    write.csv(data.frame(gene = genes_m), out_m, row.names = FALSE)
  }
}

# -----------------------
# 12) Save updated Seurat object
# -----------------------
saveRDS(seurat_obj, SEURAT_RDS_OUT)

message("Done. Outputs in: ", OUTDIR_FIG, " and ", OUTDIR_TAB)
message("Saved Seurat object: ", SEURAT_RDS_OUT)
