library(CellChat)
library(patchwork)
library(ComplexHeatmap)
library(future)

# Define identities
Idents(HC) <- HC$Annotation
Idents(NR) <- NR$Annotation
Idents(R)  <- R$Annotation

# Create CellChat objects
cellchat_HC <- createCellChat(object = HC, group.by = "ident", assay = "RNA")
cellchat_NR <- createCellChat(object = NR, group.by = "ident", assay = "RNA")
cellchat_R  <- createCellChat(object = R, group.by = "ident", assay = "RNA")

# Assign database
CellChatDB.use <- CellChatDB.human
cellchat_HC@DB <- CellChatDB.use
cellchat_NR@DB <- CellChatDB.use
cellchat_R@DB  <- CellChatDB.use

# Preprocessing function
preprocess_cellchat <- function(obj) {
  obj <- subsetData(obj)
  obj <- identifyOverExpressedGenes(obj)
  obj <- identifyOverExpressedInteractions(obj)
  plan("multisession", workers = 4)
  options(future.globals.maxSize = 8000 * 1024^2)
  obj <- computeCommunProb(obj, type = "triMean")
  obj <- filterCommunication(obj, min.cells = 10)
  obj <- computeCommunProbPathway(obj)
  obj <- aggregateNet(obj)
  return(obj)
}

# Apply preprocessing
cellchat_HC <- preprocess_cellchat(cellchat_HC)
cellchat_NR <- preprocess_cellchat(cellchat_NR)
cellchat_R  <- preprocess_cellchat(cellchat_R)

# Merge all into one object
object.list <- list(HC = cellchat_HC, NR = cellchat_NR, R = cellchat_R)
cellchat_merged <- mergeCellChat(object.list, add.names = names(object.list))

# Global interaction overview
p1 <- compareInteractions(cellchat_merged, show.legend = FALSE)
p2 <- compareInteractions(cellchat_merged, show.legend = FALSE, measure = "weight")
p1 + p2

par(mfrow = c(1, 2), xpd = TRUE)
netVisual_diffInteraction(cellchat_merged, weight.scale = TRUE)
netVisual_diffInteraction(cellchat_merged, weight.scale = TRUE, measure = "weight")

p3 <- netVisual_heatmap(cellchat_merged)
p4 <- netVisual_heatmap(cellchat_merged, measure = "weight")
p3 + p4

# Intra-group circle plots
weight.max <- getMaxWeight(object.list, attribute = c("idents", "count"))
par(mfrow = c(1, 3), xpd = TRUE)
for (i in seq_along(object.list)) {
  netVisual_circle(
    object.list[[i]]@net$count,
    weight.scale = TRUE,
    label.edge = FALSE,
    edge.weight.max = weight.max[2],
    edge.width.max = 12,
    title.name = paste("Interactions -", names(object.list)[i])
  )
}

# Specific pathway circle plots (example: SEMA4)
pathways.show <- "SEMA4"
weight.max <- getMaxWeight(object.list, slot.name = "netP", attribute = pathways.show)
par(mfrow = c(1, 3), xpd = TRUE)
for (i in seq_along(object.list)) {
  netVisual_aggregate(
    object.list[[i]],
    signaling = pathways.show,
    layout = "circle",
    edge.weight.max = weight.max[1],
    edge.width.max = 10,
    signaling.name = paste(pathways.show, names(object.list)[i])
  )
}

# Pathway-specific heatmaps (example: MHC-II)
pathways.show <- "MHC-II"
ht <- lapply(seq_along(object.list), function(i) {
  netVisual_heatmap(
    object.list[[i]],
    signaling = pathways.show,
    color.heatmap = "Reds",
    title.name = paste(pathways.show, names(object.list)[i])
  )
})
ComplexHeatmap::draw(ht[[1]] + ht[[2]] + ht[[3]], ht_gap = unit(0.5, "cm"))

# Centrality and signaling roles
for (i in seq_along(object.list)) {
  object.list[[i]] <- netAnalysis_computeCentrality(object.list[[i]])
}
num.link <- sapply(object.list, function(x) {
  rowSums(x@net$count) + colSums(x@net$count) - diag(x@net$count)
})
weight.MinMax <- range(num.link)
gg <- lapply(seq_along(object.list), function(i) {
  netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
})
wrap_plots(plots = gg)

# Signaling change in specific cell types
gg1 <- netAnalysis_signalingChanges_scatter(cellchat_merged, idents.use = "Classical Monocytes")
gg2 <- netAnalysis_signalingChanges_scatter(cellchat_merged, idents.use = "Naïve CD8+ T cells")
gg3 <- netAnalysis_signalingChanges_scatter(cellchat_merged, idents.use = "Effector memory CD8+ T cells")
wrap_plots(plots = list(gg1, gg2, gg3))

# Pathway ranking comparison
g1 <- rankNet(cellchat_merged, mode = "comparison", measure = "weight", stacked = TRUE, do.stat = TRUE)
g2 <- rankNet(cellchat_merged, mode = "comparison", measure = "weight", stacked = FALSE, do.stat = TRUE)
g1 + g2
