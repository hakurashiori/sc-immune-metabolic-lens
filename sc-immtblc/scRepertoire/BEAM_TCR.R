# Required libraries
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)
library(VennDiagram)
library(SeuratObject)

#### Set directory paths ####
data_dir <- "~/path_to_data_directory"  
output_dir <- "~/path_to_output_directory"  

#### Load TCR data ####
clonotypes <- read.csv(file.path(data_dir, "clonotypes.csv"), stringsAsFactors = FALSE)
filtered_contigs <- read.csv(file.path(data_dir, "filtered_contig_annotations.csv"), stringsAsFactors = FALSE)

# Merge TCR metadata
filtered_contigs <- filtered_contigs %>% rename(clonotype_id = raw_clonotype_id)
merged_data_orig <- merge(filtered_contigs, clonotypes, by = 'clonotype_id')

#### TCR diversity analysis ####
clonotype_frequency <- merged_data_orig %>%
  group_by(cdr3s_aa) %>%
  summarise(frequency = n()) %>%
  arrange(desc(frequency))

shannon_div <- diversity(clonotype_frequency$frequency, index = "shannon")
simpson_div <- diversity(clonotype_frequency$frequency, index = "simpson")

#### Load and process GEX data ####
gex <- Read10X(file.path(data_dir, "sample_filtered_feature_bc_matrix_from_BEAM"))
seurat_obj <- CreateSeuratObject(counts = gex)
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 5)
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:8)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:8)

#### (Optional) Run cell state inference (e.g., SignacFast) ####
# labels_fast <- SignacFast(seurat_obj)
# celltypes <- GenerateLabels(labels_fast, E = seurat_obj)
# seurat_obj <- AddMetaData(seurat_obj, metadata = celltypes$CellStates, col.name = "celltype")

#### Merge TCR and GEX data ####
tcr_metadata <- merged_data_orig %>%
  group_by(barcode) %>%
  summarise(
    clonotype = paste(unique(clonotype_id), collapse = ","),
    cdr3 = paste(unique(cdr3), collapse = ","),
    cdr3_nt = paste(unique(cdr3_nt), collapse = ","),
    v_gene = paste(unique(v_gene), collapse = ","),
    d_gene = paste(unique(d_gene), collapse = ","),
    j_gene = paste(unique(j_gene), collapse = ","),
    chain = paste(unique(chain), collapse = ","),
    productive = paste(unique(productive), collapse = ",")
  )

# Make sure barcodes match Seurat format
seurat_obj <- RenameCells(seurat_obj, new.names = paste0(Cells(seurat_obj), "-1"))
common_barcodes <- intersect(Cells(seurat_obj), tcr_metadata$barcode)
seurat_obj <- subset(seurat_obj, cells = common_barcodes)
tcr_metadata <- tcr_metadata[tcr_metadata$barcode %in% Cells(seurat_obj), ]

# Add TCR metadata
for (col in colnames(tcr_metadata)[-1]) {
  seurat_obj <- AddMetaData(seurat_obj, metadata = tcr_metadata[[col]], col.name = col)
}

#### (Optional) Subset to CD8 T cells if needed ####
# seurat_obj <- subset(seurat_obj, idents = c("T.CD8.naive", "T.CD8.cm", "T.CD8.em"))

#### Clone analysis: Add clonotype ID and visualize ####
seurat_obj$clonotype <- tcr_metadata$clonotype[match(Cells(seurat_obj), tcr_metadata$barcode)]
Idents(seurat_obj) <- seurat_obj$clonotype

# Highlight specific clonotypes
DimPlot(seurat_obj, cells.highlight = WhichCells(seurat_obj, idents = c("sample0_clonotype1", "sample0_clonotype2", "sample0_clonotype3")))

#### Summary donut plot of top clonotypes (function) ####
topN_clonos <- function(seu, col = "clonotype", top_n = 10) {
  v <- as.character(seu@meta.data[[col]])
  v <- v[!is.na(v) & v != ""]
  tab <- sort(table(v), decreasing = TRUE)
  top_names <- names(tab)[seq_len(min(top_n, length(tab)))]
  lab <- ifelse(names(tab) %in% top_names, names(tab), "Other")
  out <- data.frame(clonotype = lab, n = as.integer(tab)) %>%
    group_by(clonotype) %>%
    summarise(n = sum(n), .groups = "drop") %>%
    mutate(prop = n / sum(n),
           clonotype = fct_reorder(clonotype, -prop))
  return(out)
}

#### Donut plot example ####
df <- topN_clonos(seurat_obj, col = "clonotype", top_n = 10) %>%
  mutate(clonotype = gsub("^sample0_", "", clonotype))
ggplot(df, aes(x = 2, y = prop, fill = clonotype)) +
  geom_col(width = 0.5, color = "darkgrey") +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) +
  theme_void() +
  theme(legend.position = "right") +
  labs(fill = "Clonotype") +
  scale_fill_manual(values = c(
    "clonotype1" = "red",
    "clonotype2" = "#A6CEE3",
    "clonotype3" = "#B2DF8A",
    "clonotype4" = "#CAB2D6",
    "clonotype5" = "#FDBF6F",
    "clonotype6" = "#FF7F00",
    "clonotype7" = "#6A3D9A",
    "clonotype8" = "#B15928",
    "clonotype9" = "#1F78B4",
    "clonotype10" = "#33A02C",
    "Other" = "lightgrey"
  ))

#### Save output ####
saveRDS(seurat_obj, file.path(output_dir, "BEAM_tcr_integrated.rds"))
