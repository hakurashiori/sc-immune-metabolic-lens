library(Seurat)
library(dplyr)
library(Hmisc)
library(pheatmap)
library(ggplot2)

compute_coherence_score <- function(seurat_obj,
                                    gene_set_A,
                                    gene_set_B,
                                    cluster_name = NULL,
                                    reduction_name = "umap",
                                    cluster_label = NULL,
                                    heatmap_title = NULL) {
  
  # 1: Optional subset
  if (!is.null(cluster_name)) {
    seurat_obj <- subset(seurat_obj, idents = cluster_name)
  }
  
  # 2: Filter valid genes
  all_genes <- rownames(seurat_obj)
  genes_A <- gene_set_A[gene_set_A %in% all_genes]
  genes_B <- gene_set_B[gene_set_B %in% all_genes]
  
  if (length(genes_A) == 0 || length(genes_B) == 0) {
    stop("One or both gene sets not present in dataset.")
  }
  
  # 3: Fetch expression
  expr_A <- FetchData(seurat_obj, vars = genes_A)
  expr_B <- FetchData(seurat_obj, vars = genes_B)
  
  # 4: Compute correlation
  cor_mat <- cor(expr_A, expr_B, method = "spearman")
  cor_test <- Hmisc::rcorr(as.matrix(expr_A), as.matrix(expr_B), type = "spearman")
  pvals <- p.adjust(cor_test$P, method = "fdr")
  pmat <- matrix(pvals, nrow = nrow(cor_mat), ncol = ncol(cor_mat),
                 dimnames = list(rownames(cor_mat), colnames(cor_mat)))
  stars <- ifelse(pmat < 0.001, "***",
                  ifelse(pmat < 0.01, "**",
                         ifelse(pmat < 0.05, "*", "")))
  stars[is.na(stars)] <- ""
  
  # 5: Compute metrics
  sign_consistency <- mean(cor_mat > 0)
  
  compute_entropy <- function(mat) {
    signs <- ifelse(mat > 0, 1, -1)
    probs <- table(signs) / length(mat)
    -sum(probs * log2(probs + 1e-10))
  }
  
  entropy <- compute_entropy(cor_mat)
  ratio <- sign_consistency / entropy
  
  # 6: Heatmap
  if (is.null(heatmap_title)) {
    heatmap_title <- paste0("Correlation: ", cluster_name)
  }
  pheatmap(cor_mat,
           color = colorRampPalette(c("blue", "white", "red"))(50),
           display_numbers = stars,
           fontsize_number = 8,
           main = heatmap_title,
           cluster_rows = TRUE,
           cluster_cols = TRUE)
  
  # 7: Summary barplot (if multiple clusters used)
  df <- data.frame(
    cluster = ifelse(is.null(cluster_label), "Cluster", cluster_label),
    consistency = sign_consistency,
    entropy = entropy,
    ratio = ratio
  )
  
  print(
    ggplot(df, aes(x = cluster, y = ratio, fill = cluster)) +
      geom_bar(stat = "identity", width = 0.6, color = "black") +
      geom_text(aes(label = round(ratio, 2)), vjust = -0.8, size = 5) +
      scale_fill_manual(values = "#E41A1C") +
      labs(x = "Cluster", y = "Consistency / Entropy Ratio") +
      theme_minimal(base_size = 14) +
      theme(
        legend.position = "none",
        panel.grid = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 14, face = "bold"),
        axis.line = element_line(size = 0.6, color = "black")
      ) +
      ylim(0, ratio * 1.2) +
      ggtitle("Functional Coherence Score")
  )
  
  return(list(
    correlation_matrix = cor_mat,
    p_matrix = pmat,
    significance_stars = stars,
    sign_consistency = sign_consistency,
    entropy = entropy,
    coherence_ratio = ratio
  ))
}
