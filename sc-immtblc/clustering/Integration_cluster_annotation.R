# ============================================
# Integration pipeline (Seurat + Harmony)
# ============================================

library(Seurat)
library(harmony)
library(ggplot2)

# -------------------------------
# 1. Load individual sample objects
# -------------------------------
cell.list <- list(
  readRDS("./data/sample1.rds"),
  readRDS("./data/sample2.rds"),
  readRDS("./data/sample3.rds")
  # ... add as needed
)

# -------------------------------
# 2. Preprocess each sample
# -------------------------------
for (i in 1:length(cell.list)) {
  cell.list[[i]] <- NormalizeData(cell.list[[i]], verbose = FALSE)
  cell.list[[i]] <- FindVariableFeatures(
    cell.list[[i]],
    selection.method = "vst",
    nfeatures = 2000,
    verbose = FALSE
  )
}

# -------------------------------
# 3. Find anchors and integrate
# -------------------------------
k.filter <- min(200, min(sapply(cell.list, ncol)))
cell.list.anchors <- FindIntegrationAnchors(
  object.list = cell.list,
  k.filter = k.filter,
  dims = 1:20
)

# Save anchors if needed
# saveRDS(cell.list.anchors, "./results/anchors.rds")

# Integrate datasets
pbmc <- IntegrateData(anchorset = cell.list.anchors, dims = 1:20)

# -------------------------------
# 4. Standard Seurat workflow
# -------------------------------
DefaultAssay(pbmc) <- "RNA"
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc))
ElbowPlot(pbmc)

pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)
pbmc <- RunTSNE(pbmc, dims = 1:10)

DimPlot(pbmc, reduction = "umap")
DimPlot(pbmc, reduction = "tsne")

# -------------------------------
# 5. Harmony batch correction
# -------------------------------
DefaultAssay(pbmc) <- "integrated"
pbmc <- RunHarmony(
  pbmc,
  group.by.vars = "condition1",  # replace with your metadata column
  assay.use = "integrated",
  plot_convergence = TRUE
)

pbmc <- RunUMAP(pbmc, reduction = "harmony", dims = 1:30, reduction.name = "UMAP_harmony")
pbmc <- FindNeighbors(pbmc, reduction = "harmony", dims = 1:30)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc$harmony.cluster <- pbmc$seurat_clusters

DimPlot(pbmc, reduction = "UMAP_harmony", group.by = "condition", raster = FALSE)


# ===============================
# 6. Cell Type Annotation Pipeline
# ===============================

library(Seurat)
library(dplyr)
library(magrittr)
library(SingleR)
library(celldex)
library(SignacX)
library(ggplot2)

# Set default assay and cluster identity
DefaultAssay(pbmc) <- "RNA"
Idents(pbmc) <- pbmc$harmony.cluster

# Find marker genes for each cluster
features <- FindAllMarkers(
  pbmc,
  only.pos = TRUE,
  min.pct = 0,
  logfc.threshold = 0.25,
  verbose = TRUE
)


# Dot plot for marker genes across known immune cell types
DotPlot(pbmc, features = c(
  "CCR7", "TCF7", "LEF1", "NOSIP", "MAL", "LTB", "IL7R", "AQP3", "LDHB",     # Naive T cells
  "CD8B", "CD8A", "GZMK", "GZMH", "CCL5", "KLRB1",                           # Effector & Memory T cells
  "TRAV1-2", "FOXP3", "IL32", "ANK3", "BCL2", "TSHZ2",                       # Tregs / Helper T
  "GNLY", "GZMB", "PRF1", "SPON2", "PTGDS",                                 # NK cells
  "TRDV2", "TRGV9", "TRDC", "KLRG1",                                        # γδ T cells
  "S100A8", "S100A9", "LYZ", "VCAN", "LST1", "CD14",                         # Classical Monocytes
  "CDKN1C", "IFITM3", "SERPINA1", "MS4A7", "FCGR3A",                         # Non-classical Monocytes
  "PLD4", "TCF4", "ITM2C",                                                  # pDCs
  "IRF8", "HLA-DRA", "HLA-DRB1", "CD79A", "IGKC", "MS4A1",                  # B cells / DCs
  "PPBP", "NRGN", "TUBB1", "PF4", "CAVIN2",                                 # Megakaryocytes
  "HBB", "HBA2", "HBA1", "HBG2", "HBG1"                                     # Erythroid cells
), group.by = "harmony.cluster") +
  NoLegend() + RotatedAxis() + coord_flip()

# Manually assign identities based on marker expression
cell_types <- c(
  "Classical Monocytes", "Naïve CD4+ T cells", "Effector memory CD8+ T cells",
  "Memory CD4+ T cells", "Natural killer cells", "B cells", "Naïve CD8+ T cells",
  "Natural killer T cells", "Classical Monocytes", "Helper CD4+ T cells",
  "Gamma delta T cells", "Non-classical Monocytes", "Megakaryocytes",
  "Regulatory T cells", "Regulatory T cells", "Classical Monocytes",
  "Natural killer T cells", "Dendritic cells", "Classical Monocytes",
  "Natural killer cells", "Non-classical Monocytes", "Erythroid cells",
  "Naïve CD4+ T cells", "Plasmacytoid Dendritic cells", "B cells",
  "undefined", "Natural killer cells", "undefined", "Erythroid cells", "undefined", "undefined"
)

names(cell_types) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, cell_types)
pbmc$CellType <- Idents(pbmc)

# Save final cell type labels
pbmc$CellType_sample <- paste(pbmc$CellType, pbmc$sample, sep = "_")

# Reorder levels for plotting
new_order <- c(
  "Naïve CD4+ T cells", "Memory CD4+ T cells", "Naïve CD8+ T cells",
  "Effector memory CD8+ T cells", "Regulatory T cells", "Helper CD4+ T cells",
  "Natural killer cells", "Gamma delta T cells", "Classical Monocytes",
  "Non-classical Monocytes", "Plasmacytoid Dendritic cells", "Dendritic cells",
  "B cells", "Megakaryocytes", "Erythroid cells"
)

pbmc <- subset(pbmc, idents = "undefined", invert = TRUE)
Idents(pbmc) <- factor(pbmc$CellType, levels = new_order)

# Custom colors
colors <- c(
  "Naïve CD4+ T cells" = "#66c2a5",
  "Memory CD4+ T cells" = "#fc8d62",
  "Naïve CD8+ T cells" = "#e78ac3",
  "Effector memory CD8+ T cells" = "#8da0cb",
  "Regulatory T cells" = "#a6d854",
  "Helper CD4+ T cells" = "#ffd92f",
  "Natural killer cells" = "#e5c494",
  "Gamma delta T cells" = "#b3b3b3",
  "Classical Monocytes" = "#fb8072",
  "Non-classical Monocytes" = "#80b1d3",
  "Plasmacytoid Dendritic cells" = "#fdb462",
  "Dendritic cells" = "#b3de69",
  "B cells" = "#fccde5",
  "Megakaryocytes" = "#d9d9d9",
  "Erythroid cells" = "#bc80bd"
)

# Plot annotated cell types
DimPlot(pbmc, reduction = "umap_harmony", group.by = "CellType", label = TRUE) +
  scale_color_manual(values = colors)

# Save processed object
saveRDS(pbmc, "pbmc_annotated.rds")