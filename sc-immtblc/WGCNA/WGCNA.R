library(WGCNA)
library(Seurat)
library(ggplot2)
library(ggrepel)


genes_of_interest <- c(
  "MT-ND3", "MT-ATP8", "RAC1", "IFI30", "GAPDH", "TMSB10", "MINOS1", 
  "RPL10", "ADIPOR1", "RCBTB2", "MT-CO2", "HNRNPH3", "NINJ1", "RPL7A", 
  "RPL24", "ADGRE2", "TMED5", "C8orf59", "IRF2BP2", "STRN4"
)

pbmc <- FindVariableFeatures(pbmc, nfeatures = 2000)
variable_genes <- VariableFeatures(pbmc)
genes_for_network <- unique(c(variable_genes, genes_of_interest))

pbmc <- subset(pbmc, features = genes_for_network)
pbmc <- ScaleData(pbmc, features = genes_for_network)
expr_data <- as.data.frame(GetAssayData(pbmc, assay = "RNA", slot = "scale.data"))
datExpr <- t(expr_data)

# Soft-thresholding
powers <- c(1:10, seq(12, 20, 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
softPower <- 1

# Network Construction 
adjacency <- adjacency(datExpr, power = softPower)
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1 - TOM
geneTree <- hclust(as.dist(dissTOM), method = "average")
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = 30)
dynamicColors <- labels2colors(dynamicMods)

# Module-Trait Relationship
trait_data <- pbmc@meta.data$condition2
trait_data_df <- as.data.frame(trait_data)
trait_data_numeric <- as.data.frame(model.matrix(~ trait_data - 1, data = trait_data_df))
trait_data_filtered <- trait_data_numeric[, -((ncol(trait_data_numeric) - 2):ncol(trait_data_numeric))]

MEList <- moduleEigengenes(datExpr, colors = dynamicColors)
MEs <- MEList$eigengenes
moduleTraitCor <- cor(MEs, trait_data_filtered, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

textMatrix <- paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
par(mar = c(16, 14, 4, 2) + 0.1)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = c('HC','NR','R'),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1,1),
               main = "Module-Trait Relationships")

# Module Membership vs Gene Significance
selected_module <- "turquoise"
module_genes <- dynamicColors == selected_module
MM <- abs(cor(datExpr, MEs[[paste0("ME", selected_module)]], use = "pairwise.complete.obs"))
GS <- abs(cor(datExpr, trait_data_filtered, use = "pairwise.complete.obs"))
MM_module <- MM[module_genes]
GS_module <- GS[module_genes, 1] # example with first trait (e.g., NR)

plotData <- data.frame(ModuleMembership = MM_module, GeneSignificance = GS_module)
ggplot(plotData, aes(x = ModuleMembership, y = GeneSignificance)) +
  geom_point(color = selected_module, alpha = 0.6) +
  geom_smooth(method = "lm", color = "red", linetype = "dashed") +
  labs(title = paste0("MM vs GS: ", selected_module, " module"),
       x = "Module Membership", y = "Gene Significance") +
  theme_minimal()

# Hub Gene Identification
adjacency_matrix <- adjacency(datExpr, power = softPower)
connectivity <- colSums(adjacency_matrix[module_genes, module_genes])
gene_names <- rownames(expr_data)
GS_specific <- GS[, "trait_dataPre-vaccination responder"]
hub_genes <- gene_names[abs(MM) > 0.2 & abs(GS_specific) > 0.2]
print(hub_genes)

# Highlight specific genes 
highlight_genes <- c("RBFOX2", "ZSWIM6", "ARHGAP26", "RABGEF1", "IFI30", "MT-CO2")
highlight_idx <- which(gene_names %in% highlight_genes)
plotData <- data.frame(ModuleMembership = MM,
                       GeneSignificance = GS_specific,
                       GeneName = gene_names,
                       Label = ifelse(gene_names %in% highlight_genes, gene_names, NA))
ggplot(plotData, aes(x = ModuleMembership, y = GeneSignificance)) +
  geom_point(color = "grey70", alpha = 0.5) +
  geom_point(data = plotData[!is.na(plotData$Label), ], color = "red", size = 3) +
  geom_label_repel(aes(label = Label), size = 3, fill = "white") +
  geom_smooth(method = "lm", color = "green", linetype = "dashed", se = FALSE) +
  labs(title = "Hub Gene Identification: Highlighted Genes",
       x = "Module Membership", y = "Gene Significance") +
  theme_light()

