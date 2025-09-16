# Pathway Enrichment Analysis 

# Load required libraries
library(Seurat)
library(dplyr)
library(magrittr)
library(fgsea)
library(msigdbr)
library(gprofiler2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AUCell)
library(GSEABase)
library(pheatmap)

# ------------------------
# Data preparation
# ------------------------
DefaultAssay(obj) <- "RNA"
Idents(obj) <- obj$condition2

# Differential expression for multiple comparisons
NR.DEG <- FindMarkers(obj, ident.1 = 'Pre-vaccination non-responder', ident.2 = 'Healthy control', logfc.threshold = 0, min.pct = 0)
R.DEG  <- FindMarkers(obj, ident.1 = 'Pre-vaccination responder', ident.2 = 'Healthy control', logfc.threshold = 0, min.pct = 0)
HC.DEG <- FindMarkers(obj, ident.1 = 'Vaccinated healthy control', ident.2 = 'Post-vaccination responder', logfc.threshold = 0, min.pct = 0)

# ------------------------
# 1. GSEA using fgsea
# ------------------------
# Load gene sets
geneSet <- gmtPathways("PATH/TO/KEGG_GLYCOLYSIS.gmt")
pathway <- geneSet[["KEGG_MEDICUS_REFERENCE_GLYCOLYSIS"]]

# Create ranked gene lists
rankedGenes_NR <- setNames(NR.DEG$avg_log2FC, rownames(NR.DEG))
rankedGenes_R  <- setNames(R.DEG$avg_log2FC, rownames(R.DEG))
rankedGenes_HC <- setNames(HC.DEG$avg_log2FC, rownames(HC.DEG))

# Run fgsea
fgseaRes <- fgsea(pathways = list(Glycolysis = pathway), stats = rankedGenes_R, minSize = 3, maxSize = 350000, nperm = 350000)
plotEnrichment(pathway, rankedGenes_HC) + labs(title = "KEGG Glycolysis Pathway")

# ------------------------
# 2. gprofiler2 enrichment
# ------------------------
deg <- FindMarkers(obj, ident.1 = 'Pre-vaccination non-responder', ident.2 = 'Pre-vaccination responder', logfc.threshold = 0.2, min.pct = 0)
signif_genes <- deg[deg$p_val_adj < 0.05, ]
up_genes <- rownames(signif_genes[signif_genes$avg_log2FC > 0, ])
down_genes <- rownames(signif_genes[signif_genes$avg_log2FC < 0, ])

go_res_up <- gost(query = up_genes, organism = "hsapiens", sources = c("GO:BP", "GO:MF", "GO:CC", "REAC", "KEGG"))
go_res_down <- gost(query = down_genes, organism = "hsapiens", sources = c("GO:BP", "GO:MF", "GO:CC", "REAC", "KEGG"))

# ------------------------
# 3. clusterProfiler enrichment
# ------------------------
ego_up <- enrichGO(gene = up_genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
ego_down <- enrichGO(gene = down_genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)

# ------------------------
# 4. AUCell scoring
# ------------------------
# Load gene sets
gmt_files <- list.files("PATH/TO/GMT_FILES/", pattern = "\\.gmt$", full.names = TRUE)
gene_sets <- lapply(gmt_files, getGmt)
selected_sets <- list(
  Glycolysis = gene_sets[["REACTOME_GLUCOSE_METABOLISM.v2023.2.Hs.gmt"]][["REACTOME_GLUCOSE_METABOLISM"]],
  LipidMetabolism = gene_sets[["REACTOME_METABOLISM_OF_LIPIDS.v2023.2.Hs.gmt"]][["REACTOME_METABOLISM_OF_LIPIDS"]]
)
geneSetList <- lapply(selected_sets, function(x) x@geneIds)

# Prepare matrix
expr_mat <- GetAssayData(obj, slot = "data")
cell_ranking <- AUCell_buildRankings(expr_mat, plotStats = TRUE)
aucell_scores <- AUCell_calcAUC(geneSetList, cell_ranking, aucMaxRank = ceiling(0.2 * nrow(cell_ranking)))
aucMat <- getAUC(aucell_scores)

# Heatmap
annot <- data.frame(Condition = obj$condition2)
rownames(annot) <- colnames(aucMat)
pheatmap(
  aucMat,
  annotation_col = annot,
  show_colnames = FALSE,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("navy", "white", "firebrick"))(50)
)
