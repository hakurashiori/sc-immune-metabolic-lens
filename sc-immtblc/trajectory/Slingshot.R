# ================================
# Slingshot trajectory analysis
# ================================

library(Seurat)
library(ggplot2)
library(SingleCellExperiment)
library(slingshot)
library(tradeSeq)
library(dplyr)
library(tidyr)
library(patchwork)
library(viridis)

# -------------------------------
# Preprocessing with Seurat
# -------------------------------
Idents(pbmc) <- pbmc$condition2
HC <- subset(pbmc, idents = 'Healthy control')

DefaultAssay(NR) <- 'RNA'
NR <- NormalizeData(NR)
NR <- ScaleData(NR, verbose = FALSE)
NR <- FindVariableFeatures(NR)
NR <- RunPCA(NR, verbose = FALSE)
ElbowPlot(NR)

NR <- RunUMAP(NR, dims = 1:20, reduction.name = "UMAP")
NR <- FindNeighbors(NR, dims = 1:20)
NR <- FindClusters(NR, resolution = 0.8)

DimPlot(NR, group.by = 'patientID_pre_NR', reduction = 'UMAP', label = TRUE)

# -------------------------------
# Convert to SingleCellExperiment
# -------------------------------
sce_2 <- as.SingleCellExperiment(NR)

# Use Seurat clustering results
sce_2$clusters <- NR$seurat_clusters
table(sce_2$clusters)

set.seed(1)  # For reproducibility

# -------------------------------
# Slingshot trajectory inference
# -------------------------------
sce_2 <- slingshot(
  sce_2,
  clusterLabels = sce_2$clusters,
  reducedDim = "UMAP",
  start.clus = "2"   # Modify as appropriate
)

# Extract lineages and curves
lineages <- slingLineages(sce_2)
curves <- slingCurves(sce_2)

print(lineages)
print(curves)

# Plot UMAP with trajectories
plot(reducedDim(sce_2, "UMAP"), col = sce_2$clusters, pch = 16, cex = 0.5)
lines(SlingshotDataSet(sce_2), col = "black", lwd = 2)

# Extract pseudotime
pseudotime_2 <- slingPseudotime(sce_2)
NR$Pseudotime <- pseudotime_2

# -------------------------------
# Visualization: Pseudotime
# -------------------------------
df_2 <- data.frame(
  UMAP1 = reducedDim(sce_2, "UMAP")[, 1],
  UMAP2 = reducedDim(sce_2, "UMAP")[, 2],
  Pseudotime = pseudotime_2[, 1]   # Example: first lineage
)

curve_df_2 <- as.data.frame(slingCurves(sce_2)[[1]]$s)

ggplot(df_2, aes(x = UMAP1, y = UMAP2, color = Pseudotime)) +
  geom_point(size = 1.5, alpha = 0.9) +
  geom_path(data = curve_df_2, aes(x = UMAP_1, y = UMAP_2), color = "black", size = 1) +
  scale_color_viridis(option = "viridis", direction = -1, end = 0.95, name = "Pseudotime") +
  coord_fixed() +
  theme_classic(base_size = 13)

# Alternative gradient color
ggplot(df_2, aes(x = UMAP1, y = UMAP2, color = Pseudotime)) +
  geom_point(size = 2, alpha = 0.9) +
  geom_path(data = curve_df_2, aes(x = UMAP_1, y = UMAP_2), color = "black", size = 1) +
  scale_color_gradientn(colors = c("#2B004F", "#800080", "#FF66CC", "#FFD700", "#FFFF00"),
                        name = "Pseudotime") +
  coord_fixed() +
  theme_classic(base_size = 13)

# -------------------------------
# tradeSeq analysis
# -------------------------------
counts_2 <- assay(sce_2, "counts")
pseudotime_2 <- slingPseudotime(sce_2)

# Filter lowly expressed genes
filt_counts_2 <- counts_2[rowSums(counts_2 > 5) > ncol(counts_2)/100, ]

# Fit GAM
sds_2 <- slingshot(Embeddings(NR, "UMAP"), clusterLabels = sce_2$seurat_clusters, start.clus = "2")
sce_GAM <- fitGAM(counts = as.matrix(filt_counts_2), sds = sds_2)

# Association test
assocRes_2 <- associationTest(sce_GAM, lineages = TRUE)
topGenes_2 <- rownames(assocRes_2)[order(assocRes_2$pvalue)][1:10]

# Visualize smoothed expression
plotSmoothers(sce_GAM, filt_counts_2, gene = topGenes_2[1]) + ggtitle(topGenes_2[1])

# Identify lineage-specific best fits
lineage_pval_cols <- grep("^pvalue_", colnames(assocRes_2), value = TRUE)
assocRes_2$best_lineage <- apply(assocRes_2[, lineage_pval_cols, drop = FALSE], 1, function(x) {
  valid_pvals <- x[!is.na(x)]
  if (length(valid_pvals) == 0) return(NA)
  names(x)[which.min(valid_pvals)]
})

# -------------------------------
# Differential expression at lineage endpoints
# -------------------------------
diff_end_association_2 <- diffEndTest(sce_GAM)
diff_end_association_2$feature_id <- rownames(diff_end_association_2)

significant_genes_2 <- diff_end_association_2 %>%
  filter(pvalue < 0.05) %>%
  arrange(desc(waldStat)) %>%
  dplyr::slice(1:30) %>%
  pull(feature_id)

# Example visualization
plotSmoothers(sce_GAM, as.matrix(filt_counts_2), gene = significant_genes_2[1]) +
  ggtitle(significant_genes_2[1])

# -------------------------------
# Save results (placeholder paths)
# -------------------------------
saveRDS(sce_2, "./results/sce_slingshot.rds")
saveRDS(sce_GAM, "./results/sce_GAM.rds")
saveRDS(sds_2, "./results/sds.rds")
saveRDS(filt_counts_2, "./results/filt_counts.rds")


# ================================
# Lineage-specific late-phase regression
# ================================

# NR: extract pseudotime trajectories
df_NR <- as.data.frame(t(pred_NR))
df_NR$pseudotime <- seq(0, 1, length.out = nrow(df_NR))
df_NR_lineage1 <- df_NR[grep("^lineage1_", rownames(df_NR)), ]
df_NR_lineage2 <- df_NR[grep("^lineage2_", rownames(df_NR)), ]
df_NR_lineage3 <- df_NR[grep("^lineage3_", rownames(df_NR)), ]

# R: extract pseudotime
df_R <- as.data.frame(t(pred_R))
df_R$pseudotime <- seq(0, 1, length.out = nrow(df_R))

# HC: extract pseudotime trajectories
df_HC <- as.data.frame(t(pred_HC))
df_HC$pseudotime <- seq(0, 1, length.out = nrow(df_HC))
df_HC_lineage1 <- df_HC[grep("^lineage1_", rownames(df_HC)), ]
df_HC_lineage2 <- df_HC[grep("^lineage2_", rownames(df_HC)), ]
df_HC_lineage3 <- df_HC[grep("^lineage3_", rownames(df_HC)), ]
df_HC_lineage4 <- df_HC[grep("^lineage4_", rownames(df_HC)), ]

# Regression: late phase NR
late_phase_NR <- df_NR_lineage3[df_NR_lineage3$pseudotime > 0.886, ]
lm_late_NR <- lm(S100A9 ~ pseudotime, data = late_phase_NR)
summary(lm_late_NR)

# Regression: late phase R
late_phase_R <- df_R[df_R$pseudotime > 0.656, ]
lm_late_R <- lm(S100A9 ~ pseudotime, data = late_phase_R)
summary(lm_late_R)

# Regression: late phase HC (across lineages)
late_phase_HC <- df_HC_lineage4[df_HC_lineage4$pseudotime > 0.914, ]
lm_late_HC <- lm(S100A9 ~ pseudotime, data = late_phase_HC)
summary(lm_late_HC)

# -------------------------------
# Align pseudotime across groups
# -------------------------------
late_phase_NR <- df_NR_lineage3[df_NR_lineage3$pseudotime > 0.75, ]
late_phase_R  <- df_R[df_R$pseudotime > 0.75, ]
late_phase_HC <- df_HC_lineage4[df_HC_lineage4$pseudotime > 0.75, ]

late_phase_NR$group <- "NR"
late_phase_R$group  <- "R"
late_phase_HC$group <- "HC"

df_all <- rbind(late_phase_NR, late_phase_R, late_phase_HC)
df_all$group <- factor(df_all$group, levels = c("NR", "R", "HC"))

# Linear model with interaction
model <- lm(S100A9 ~ pseudotime * group, data = df_all)
summary(model)

# Visualization: late-phase expression vs pseudotime
ggplot(df_all, aes(x = pseudotime, y = S100A9, color = group)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_minimal() +
  labs(title = "S100A9 Expression ~ Pseudotime across Conditions (Late Phase)")

ggplot(df_all, aes(x = pseudotime, y = S100A9, color = group)) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_minimal() +
  labs(title = "S100A9 expression at transition phase from pre- to post-vaccination") +
  annotate("text", x = 0.95, y = 75, label = "italic(P)['NR vs HC'] < 2e-16", parse = TRUE) +
  annotate("text", x = 0.95, y = 70, label = "italic(P)['NR vs R'] < 2e-16", parse = TRUE)

# ================================
# Quartile-level abundance changes
# ================================

df <- data.frame(
  Group = rep(c("HC", "R", "NR"), each = 4),
  Quartile = rep(c("Q1", "Q2", "Q3", "Q4"), times = 3),
  Pre = c(0.451, 0.223, 0.229, 0.0847,   # HC
          0.223, 0.400, 0.112, 0.0328,   # R
          0.703, 0.708, 0.683, 0.107),   # NR
  Post = c(0.549, 0.777, 0.771, 0.915,   # HC
           0.777, 0.600, 0.888, 0.967,   # R
           0.297, 0.292, 0.317, 0.893)   # NR
)

# Compute difference (Post - Pre)
df_diff <- df %>%
  mutate(diff = Post - Pre) %>%
  select(Group, Quartile, diff)

# Plot abundance changes along pseudotime quartiles
ggplot(df_diff, aes(x = Quartile, y = diff, group = Group, color = Group)) +
  geom_hline(yintercept = 0, color = "black", size = 1.2, linetype = "dashed") +
  geom_line(size = 1.5) +
  geom_point(size = 3) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Post-Vaccination Gain in Monocyte Abundance Along Pseudotime Trajectory",
    y = expression(Delta~"Proportion (Post – Pre)"),
    x = "Trajectory Stage (Quartile)"
  ) +
  scale_color_manual(
    values = c("HC" = "gray40", "R" = "#377EB8", "NR" = "#E41A1C"),
    breaks = c("HC", "R", "NR")
  ) +
  theme(
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold"),
    legend.title = element_blank()
  )
