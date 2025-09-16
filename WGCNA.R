# https://fuzzyatelin.github.io/bioanth-stats/module-F21-Group1/module-F21-Group1.html
library(WGCNA)

genes_of_interest <- c(
  #'PLAU', 'EPHX2', 'GSTA4', 'UBE2E1', 'OGDH', 'ILKAP', 'PTDSS2', 'CTRC', 'NRF1', 'SLC2A9', 'SLC9A8', 'NDST1', 'ITPK1', 'ACACA',
  #'MT-ND3', 'MT-ATP8', 
  #'MIDN',
  #'KLF2', 'JUND', 'TMEM176B', 'DUSP6', 'TSC22D3', 'ATP2B1', 'LST1', 'AKIRIN2', 'ITM2B', 'IRF1', 'C20orf27', 'TYMP', 'SULT1A1', 'MAFB', 'IFITM2', 'MYADM', 'RAB11FIP1', 'UBE2L6', 'BST2', 'CD302', 'LGALS9', 'NDUFA3', 'TNFSF10',
  'MT-ND3', 'MT-ATP8','RAC1', 'IFI30', 'GAPDH', 'TMSB10', 'MINOS1', 
  'RPL10', 'ADIPOR1', 'RCBTB2', 'MT-CO2', 'HNRNPH3', 'NINJ1', 'RPL7A', 
  'RPL24', 'ADGRE2', 'TMED5', 'C8orf59', 'IRF2BP2', 'STRN4'
)  

seurat_obj = cls_4
seurat_obj = FindVariableFeatures(seurat_obj, nfeatures = 2000) # 4000
variable_genes <- VariableFeatures(seurat_obj)
genes_for_network <- unique(c(variable_genes, genes_of_interest))

# Scale the data (important for WGCNA)
# Subset the Seurat object to include only the genes of interest
seurat_obj <- subset(seurat_obj, features = genes_for_network)
seurat_obj <- ScaleData(seurat_obj, features = genes_for_network)

# Extract the scaled expression matrix for variable genes
expr_data <- as.data.frame(GetAssayData(seurat_obj, assay = "RNA", slot = "scale.data"))

# Check the dimensions of expr_data
dim(expr_data)

datExpr <- t(expr_data)
# Pick a soft-thresholding power
powers <- c(seq(1, 10, by = 1), seq(12, 20, by = 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2", type = "n",
     main = "Scale independence")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = 0.9, col = "red")
abline(h = 0.75, col = "red")

plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = 0.9, col = "red")



softPower <- 1


# Construct the adjacency matrix
adjacency <- adjacency(datExpr, power = softPower)

# Turn adjacency into a topological overlap matrix (TOM)
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1 - TOM

# Perform hierarchical clustering on the TOM dissimilarity
geneTree <- hclust(as.dist(dissTOM), method = "average")

# Module identification using dynamic tree cut
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = 30)
table(dynamicMods)

# Convert module labels to colors
dynamicColors <- labels2colors(dynamicMods)

# Plot the dendrogram with module colors
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


# Assume we have a metadata column in Seurat object named 'condition'
trait_data <- seurat_obj@meta.data$condition2
trait_data <- as.data.frame(trait_data)

# Calculate module eigengenes
MEList <- moduleEigengenes(datExpr, colors = dynamicColors)
MEs <- MEList$eigengenes

# Use the model.matrix function to create dummy variables
trait_data_numeric <- as.data.frame(model.matrix(~ trait_data - 1, data = trait_data))
# Drop the last three columns
trait_data_filtered <- trait_data_numeric[, -((ncol(trait_data_numeric) - 2):ncol(trait_data_numeric))]
# Correlate module eigengenes with trait data
moduleTraitCor <- cor(MEs, trait_data_filtered, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

# Display correlations and p-values in a heatmap
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
par(mar = c(16, 14, 4, 2) + 0.1) 
xlabel = c('Healthy control','Pre-vaccination non-responder','Pre-vaccination responder')
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = xlabel,
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1,1),
               main = "Module-trait relationships")

# Adjust module-trait relationship heatmap
par(cex.lab = 0.3)
par(mar = c(10, 4, 4, 2) + 0.1) 
labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = names(trait_data_filtered),
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = colorRampPalette(c("blue", "white", "red"))(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.7,         # Font size for values in cells
  cex.lab.x = 1.0,        # Font size for x-axis labels
  cex.lab.y = 1.0,        # Font size for y-axis labels
  main = "Module-Trait Relationships",
  zlim = c(-1, 1)         # Color scale range
)



#Select the modules
modules = c("turquoise")
modules = c("blue", "turquoise")
#Select the gene modules
genes = colnames(datExpr)

#if you want export specific colors, substitute the second modulecolors by above modules
inModule = is.finite(match(dynamicColors, modules))
modGenes_R = genes[inModule]
modGenes_R
#write.csv(modGenes_R, '~/Desktop/Thesis_Progress/Output_Data/modGenes_R.csv')
#write.csv(modGenes_NR, '~/Desktop/Thesis_Progress/Output_Data/modGenes_NR.csv')


### Plot correlation scatter of genes
MM <- abs(cor(datExpr, MEs$MEblue, use = "pairwise.complete.obs"))
trait_data <- as.numeric(as.character(trait_data))
# Remove the last 3 columns
#trait_data_numeric <- trait_data_numeric[, -((ncol(trait_data_numeric) - 2):ncol(trait_data_numeric))]
GS <- abs(cor(datExpr, trait_data_numeric, use = "pairwise.complete.obs"))


# Filter only the genes in the "blue" module
blueGenes <- dynamicColors == "blue"
MM_blue <- MM[blueGenes]
GS_blue <- GS[blueGenes]

plotData <- data.frame(
  ModuleMembership = MM_blue,
  GeneSignificance = GS_blue
)
plot(
  plotData$ModuleMembership,
  plotData$GeneSignificance,
  xlab = "Module Membership in blue module",
  ylab = "Gene significance",
  main = "Blue module membership vs. gene significance",
  pch = 21,          # Use open circles
  bg = "lightblue",       # Fill color for points
  col = "grey",      # Outline color for points
  cex = 0.6          # Point size
)

# Add linear regression line
fit <- lm(GS ~ MM)  # Fit a linear model
abline(fit, col = "red", lwd = 2)  # Add regression line
specific_genes <- c("MT-ATP8",'GAPDH', "MT-ND2")  # Replace with your target genes
#gene_names <- rownames(expr_data)
specific_indices <- which(gene_names %in% specific_genes)

text(
  MM[specific_indices], GS[specific_indices], 
  labels = gene_names[specific_indices], 
  pos = 3,  # Position of the labels (4 = right)
  col = "red", cex = 1.2  # Color and size of text
)

# Highlight specific points
points(
  MM[specific_indices], GS[specific_indices], 
  pch = 21, bg = "red", col = "black", cex = 1.5  # Larger size and red fill
)

MM <- abs(cor(datExpr, MEs$MEturquoise, use = "pairwise.complete.obs"))
trait_data <- as.numeric(as.character(trait_data))
# Remove the last 3 columns
#trait_data_numeric <- trait_data_numeric[, -((ncol(trait_data_numeric) - 2):ncol(trait_data_numeric))]
GS <- abs(cor(datExpr, trait_data_numeric, use = "pairwise.complete.obs"))

tqGenes <- dynamicColors == "turquoise"
MM_tq <- MM[tqGenes]
GS_tq <- GS[tqGenes]

plotData <- data.frame(
  ModuleMembership = MM_tq,
  GeneSignificance = GS_tq
)
plot(
  plotData$ModuleMembership,
  plotData$GeneSignificance,
  xlab = "Module Membership in turquoise module",
  ylab = "Gene significance",
  main = "Turquoise module membership vs. gene significance",
  pch = 21,          # Use open circles
  bg = "turquoise",       # Fill color for points
  col = "grey",      # Outline color for points
  cex = 0.6          # Point size
)

# Add linear regression line
fit <- lm(GS ~ MM)  # Fit a linear model
abline(fit, col = "red", lwd = 2)  # Add regression line

specific_genes <- c("ATP2B1",'KLF12')  # Replace with your target genes
#gene_names <- rownames(expr_data)
specific_indices <- which(gene_names %in% specific_genes)

text(
  MM[specific_indices], GS[specific_indices], 
  labels = gene_names[specific_indices], 
  pos = 3,  # Position of the labels (4 = right)
  col = "red", cex = 1.2  # Color and size of text
)

# Highlight specific points
points(
  MM[specific_indices], GS[specific_indices], 
  pch = 21, bg = "red", col = "black", cex = 1.5  # Larger size and red fill
)


###### Scatter plot with regression line by ggplot
ggplot(plotData, aes(x = ModuleMembership, y = GeneSignificance)) +
  geom_point(color = "blue", alpha = 0.7) +  # Add scatter points
  geom_smooth(method = "lm", color = "red", se = FALSE, linetype = "dashed") +  # Add regression line
  labs(
    title = "Module membership vs. gene significance",
    x = "Module Membership in blue module",
    y = "Gene significance"
  ) +
  theme_minimal()


### HUB GENES ------
module_genes <- dynamicColors == "turquoise"  # Replace "blue" with the module of interest
adjacency_matrix <- adjacency(datExpr, power = 1)  # Use the chosen soft-thresholding power
connectivity <- colSums(adjacency_matrix[module_genes, module_genes])
moduleEigengene <- MEs$MEturquoise  
MM <- cor(datExpr, moduleEigengene, use = "pairwise.complete.obs")
GS <- cor(datExpr, trait_data_filtered, use = "pairwise.complete.obs")
gene_names <- rownames(expr_data)

GS_specific <- GS[, which(colnames(GS) == "trait_dataPre-vaccination responder")]
hub_genes <- gene_names[abs(MM) > 0.2 & abs(GS_specific) > 0.2]
hub_genes
# NR: "CD99"    "TMSB10"  "EEF1A1"  "MINOS1"  "SEPT9"   "PCBP2"   "S100A10" "CRIP1"   "HBB"
# R: "RBFOX2","ZSWIM6","ARHGAP26","RABGEF1","SIPA1L1","MT-CO2","ADAMTSL4-AS1","IFI30","MTRNR2L12"   

specific_genes <- c("RBFOX2","ZSWIM6","ARHGAP26","RABGEF1","SIPA1L1","MT-CO2","ADAMTSL4-AS1","IFI30","MTRNR2L12")  
#gene_names <- rownames(expr_data)
specific_indices <- which(gene_names %in% specific_genes)

# Create a data frame for ggplot2
plotData <- data.frame(
  ModuleMembership = MM,
  GeneSignificance = GS_specific,
  GeneName = gene_names
)

# Highlight specific genes
plotData$Label <- ifelse(plotData$GeneName %in% specific_genes, plotData$GeneName, NA)


ggplot(plotData, aes(x = ModuleMembership, y = GeneSignificance)) +
  # Background points
  geom_point(color = "blue", alpha = 0.5, size = 1) +
  # Highlighted points
  geom_point(data = plotData[!is.na(plotData$Label), ], color = "red", size = 3) +
  # Repelling labels with a white background
  geom_label_repel(aes(label = Label), size = 3, color = "red", fill = "white", 
                   box.padding = 0.5, max.overlaps = 15, segment.color = "gray50") +
  # Regression line
  geom_smooth(method = "lm", color = "lightgreen", size = 1, linetype = "dashed", se = FALSE) +
  # Title and axis labels
  labs(
    title = "Hub Gene Identification: Potential Responders",
    subtitle = "Correlation between Module Membership and Gene Significance",
    x = "Module Membership",
    y = "Gene Significance"
  ) +
  # Customize theme
  theme_light(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "none"
  ) +
  # Adjust axis limits and gridlines
  scale_x_continuous(limits = c(-0.1, 0.6), breaks = seq(-0.5, 0.6, 0.2)) +
  scale_y_continuous(limits = c(-0.05, 0.5), breaks = seq(-0.25, 0.5, 0.1))


#  Plot by each patients -----
trait_data <- seurat_obj@meta.data$patients
trait_data <- as.data.frame(trait_data)

# Calculate module eigengenes
MEList <- moduleEigengenes(datExpr, colors = dynamicColors)
MEs <- MEList$eigengenes

# Use the model.matrix function to create dummy variables
trait_data_numeric <- as.data.frame(model.matrix(~ trait_data - 1, data = trait_data))

# Correlate module eigengenes with trait data
moduleTraitCor <- cor(MEs, trait_data_numeric, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

# Display correlations and p-values in a heatmap
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
par(mar = c(10, 12, 4, 2))
moduleTraitCor <- moduleTraitCor[, order]
par(mar = c(14, 16, 4, 2))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = order,
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1,1),
               main = "Module-trait relationships")


order = c('52BA23',
          'HD1_GEX',
          'healthy_control_1',
          'healthy_control_10',
          'healthy_control_11',
          'healthy_control_2',
          'healthy_control_4',
          'healthy_control_5',
          'healthy_control_7',
          'healthy_control_8',
          'Kolab',
          'P0_S1',
          'P0_S2',
          'P0_S3',
          'P0_S4',
          'P0_S5',
          'P0_S6',
          'S74',
          'S76',
          'S78',
          'V02_ID0',
          'V16_ID0',
          'V24_ID0',
          'V36_ID0',
          'V47_ID0',
          'V50_ID0',
          'V51_ID0',
          'V54_ID0',
          'V56_ID0',
          'V20_ID0',
          'V23_ID0',
          'V26_ID0',
          'V44_ID0',
          'V45_ID0',
          'V52_ID0',
          'V53_ID0',
          'V57_ID0')





MEheatmap <- moduleEigengenes(expr_data, dynamicColors)$eigengenes
heatmap(t(MEheatmap), Rowv = NA, scale = "row",
        col = colorRampPalette(c("blue", "white", "red"))(50))


plot(moduleMembership[, "blue"], geneSignificance[, "trait"],
     xlab = "Module Membership in Blue Module",
     ylab = "Gene Significance for Trait",
     main = "Module Membership vs. Gene Significance")

TOMplot(dissTOM, geneTree, dynamicColors, main = "TOM Heatmap")

plotEigengeneNetworks(MEs, "Eigengene Adjacency Heatmap")

exportNetworkToCytoscape(TOM, edgeFile = "edges.txt", nodeFile = "nodes.txt",
                         weighted = TRUE, threshold = 0.02, nodeNames = genes_for_network)

barplot(MEs$MEblue, main = "Blue Module Eigengene Expression")

hist(kIN$blue, breaks = 50, main = "Intra-modular Connectivity in Blue Module")
