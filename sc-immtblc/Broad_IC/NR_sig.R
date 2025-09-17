# NR signature projecting: correlation w/ metabolic pathways

library(Seurat)
library(dplyr)
library(ggplot2)
library(tibble)
library(tidyr)
library(stringr)
library(purrr)
library(pheatmap)
library(reshape2)

NR_genes <- c(
  # NR cluster DEGs
  "MT-ATP8", "MT-ATP6", "MT-CO1", "MT-CO2", "MT-CYB", 
  "MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", "MT-ND5", "MT-ND6", "MT-ND4L",
  # repeated in every individual DEG
  'UBE2E1', 'OGDH', 'ILKAP', 'PTDSS2', 'CTRC', 'GNAL', 'PLAU', 'EPHX2', 'GSTA4',
  'GAPDH',
  'LILRB1', "ITGAM", 'CD4',
  # NR co-expression module genes
  'PLAUR', 'HSPB1', 'RNASE2', 'NRG1', 'CMTM2', 'NR4A1', 'HMGB2',
  'FOXJ1', 'ITGA4', 'SMPD3', 'CD9', 'CD99', 'RPL13A', 'CH25H',
  'THBS1', 'NINJ1', 'TRPM1', 'ANXA1', 'TNFAIP6', 'RAC1', 'IL1B',
  'DUSP1', 'CXCR4', 'CXCR2', 'CXCL3', 'SERPINE1', 'CALR',
  'CCL3', 'CCL4', 'CXCL1', 'GPR183', 'CD74', 'CXCL2','CXCL10', 'CXCL8',
  'S100A12', 'S100A9', 'S100A8', 'CCL3L1'
)

pbmc <- AddModuleScore(pbmc, features = list(NR_genes), name = "NR_Sig")

# Load REACTOME gene sets
gmt_files <- list.files("/path/to/REACTOME_db/", pattern = "\\.gmt$", full.names = TRUE)
all_pathways <- map(gmt_files, gmtPathways) %>% flatten()
named_gene_sets <- all_pathways %>%
  keep(~ length(.x) > 5) %>%
  map(~ unique(.x))

# Score all metabolic pathways
for (i in seq_along(named_gene_sets)) {
  gene_set <- named_gene_sets[[i]]
  set_name <- names(named_gene_sets)[i]
  score_name <- paste0("MS_", str_replace_all(set_name, "[^A-Za-z0-9]", "_"))
  pbmc <- AddModuleScore(object = pbmc, features = list(gene_set), name = score_name, assay = "RNA")
  message("✓ Done: ", set_name)
}

# Calculate group-wise correlations (Spearman)
score_cols <- grep("^MS_.*1$", colnames(pbmc@meta.data), value = TRUE)
meta_df <- pbmc@meta.data[, c("IC_group", "NR_Sig1", score_cols)]

cor_result <- meta_df %>%
  pivot_longer(cols = all_of(score_cols), names_to = "Pathway", values_to = "Score") %>%
  group_by(IC_group, Pathway) %>%
  summarise(
    cor = cor(NR_Sig1, Score, method = "spearman", use = "complete.obs"),
    n = n(),
    .groups = "drop"
  )

# Visualization
ggplot(cor_result, aes(x = IC_group, y = Pathway, fill = cor)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limits = c(-1, 1), name = "Spearman\nCorr") +
  theme_minimal(base_size = 10) +
  labs(title = "Correlation between NR_Sig1 and Module Scores", x = "Condition", y = "Pathway") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Filter pathway modules for Pearson correlation analysis
all_module_scores <- grep("^MS_", colnames(pbmc@meta.data), value = TRUE)
remove_keywords <- c("GSE40666", "GSE39110", "GSE25846", "REACTOME_TCR_SIGNALING", "GLUTAMATE")
filtered_module_scores <- all_module_scores[!grepl(paste(remove_keywords, collapse = "|"), all_module_scores)]

cor_results <- list()
for (group in unique(pbmc$IC_group)) {
  pbmc_group <- subset(pbmc, subset = IC_group == group)
  df <- pbmc_group@meta.data[, c("NR_Sig1", filtered_module_scores)]
  cor_list <- lapply(filtered_module_scores, function(module) {
    test <- cor.test(df$NR_Sig1, df[[module]], method = "pearson")
    tibble(module = module, cor = as.numeric(test$estimate), p = test$p.value)
  })
  cor_df <- bind_rows(cor_list)
  cor_df$IC_group <- group
  cor_results[[group]] <- cor_df
}

cor_all_df <- bind_rows(cor_results)
cor_all_df <- cor_all_df %>%
  mutate(signif = case_when(
    p <= 0.001 ~ "***",
    p <= 0.01  ~ "**",
    p <= 0.05  ~ "*",
    TRUE       ~ ""
  ))

cor_all_df$module <- factor(cor_all_df$module, levels = unique(cor_all_df$module))

# Plot correlation heatmap
ggplot(cor_all_df, aes(x = IC_group, y = module, fill = cor)) +
  geom_tile(color = "white") +
  geom_text(aes(label = signif), size = 3) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1, 1), name = "Pearson\nCorrelation") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title = element_blank())

# Standardize labels
cor_all_df$IC_group <- factor(cor_all_df$IC_group, levels = c("HC", "HIV", "HGSOC", "AD", "TNBC", "MD", "LC", "CLL", "KTR"))
cor_all_df$module_clean <- cor_all_df$module %>%
  str_replace("^MS_", "") %>%
  str_replace("^REACTOME_", "") %>%
  str_replace("^KEGG_", "") %>%
  str_replace("1$", "")
cor_all_df$module_clean <- factor(cor_all_df$module_clean, levels = unique(cor_all_df$module_clean))

# Cleaned heatmap
ggplot(cor_all_df, aes(x = IC_group, y = module_clean, fill = cor)) +
  geom_tile(color = "white") +
  geom_text(aes(label = signif), size = 3) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1, 1), name = "Pearson\nCorrelation") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title = element_blank())

# pheatmap matrix
my_palette <- colorRampPalette(c("#39FF14", "#FFFFFF", "#C800C8"))(100)
mat <- acast(cor_all_df, module_clean ~ IC_group, value.var = "cor")
signif_mat <- acast(cor_all_df, module_clean ~ IC_group, value.var = "signif")

pheatmap(mat,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = my_palette,
         display_numbers = signif_mat,
         fontsize = 10,
         border_color = NA)

# Correlation with response rate
plot_data <- data.frame(
  IC_group = c("HIV", "AD", "MD", "TNBC", "HGSOC", "LC", "CLL", "KTR"),
  NR_Sig1 = c(-0.193, 0.1205, 0.443, 0.554, 0.692, 0.763, 1.623, 2.017),
  response_rate = c(0.97, 0.715, 0.93, 0.58, 0.48, 0.6, 0.4, 0.39) # collect from ref citated and specified in articles
)

group_colors <- c(
  "HIV" = "#f7e21d", "AD" = "#47a4e6", "MD" = "#90cdb0",
  "TNBC" = "#141414", "HGSOC" = "#5be441", "LC" = "#ea670e",
  "CLL" = "#db5b15", "KTR" = "#ac5df7"
)

ggplot(plot_data, aes(x = response_rate, y = NR_Sig1, color = IC_group)) +
  geom_point(size = 3.5) +
  geom_smooth(method = "lm", se = FALSE, color = "gray40", linetype = "dashed", linewidth = 0.6) +
  geom_text(aes(label = IC_group), color = "black", nudge_y = 0.15, size = 4) +
  scale_color_manual(values = group_colors) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    legend.position = "none"
  ) +
  labs(x = "Seroconversion Rate", y = "Mean NR_Sig1 Score")
