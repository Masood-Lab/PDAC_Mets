# Loading required packages
library(ISCHIA)
library(robustbase)
library(data.table)
library(ggplot2)
library(Seurat)
library(dplyr)
library(factoextra)
library(cluster)
library(showtext)
library(gridExtra)
library(pdftools)

# Set random seed for reproducibility
set.seed(123)

# Load data
pdac <- readRDS("/path/to/pdac_mets_rctd.rds")
assay_matrix <- pdac[["rctd_tier1"]]@data
norm_weights <- as.data.frame(t(assay_matrix))

# Elbow Method
k.values <- 1:20
wss_values <- sapply(k.values, function(k) kmeans(norm_weights, k, nstart = 10)$tot.withinss)

pdf("1_elbow_plot.pdf")
plot(k.values, wss_values, type = "b", pch = 19, frame = FALSE, 
     xlab = "Number of clusters K", ylab = "Total within-cluster sum of squares",
     main = "Elbow Method for Optimal K")
dev.off()

# Gap Statistic
gap_stat <- function(k) {
  km.res <- kmeans(norm_weights, k, nstart = 10)
  if (k == 1) return(NA)
  obs_disp <- sum(km.res$withinss)
  reference_disp <- mean(replicate(10, {
    km.null <- kmeans(matrix(rnorm(nrow(norm_weights) * ncol(norm_weights)), 
                             ncol = ncol(norm_weights)), k, nstart = 10)
    sum(km.null$withinss)
  }))
  log(reference_disp) - log(obs_disp)
}

gap_stat_values <- sapply(k.values, gap_stat)

pdf("2_gap_statistic_plot.pdf")
plot(k.values, gap_stat_values, type = "b", pch = 19, frame = FALSE, 
     xlab = "Number of Clusters (K)", ylab = "Gap Statistic",
     main = "Gap Statistic: Determining Optimal K")
dev.off()

# Calinski-Harabasz Index
calinski_harabasz_index <- function(data, labels) {
  num_clusters <- length(unique(labels))
  num_points <- nrow(data)
  centroids <- tapply(data, labels, FUN = colMeans)
  between_disp <- sum(sapply(1:num_clusters, function(i) {
    cluster_points <- data[labels == i, ]
    nrow(cluster_points) * sum((colMeans(cluster_points) - centroids[i, ]) ^ 2)
  }))
  within_disp <- sum(sapply(1:num_clusters, function(i) {
    cluster_points <- data[labels == i, ]
    sum((cluster_points - centroids[i, ]) ^ 2)
  }))
  (between_disp / (num_clusters - 1)) / (within_disp / (num_points - num_clusters))
}

ch_values <- sapply(k.values, function(k) {
  km.res <- kmeans(norm_weights, k, nstart = 10)
  calinski_harabasz_index(norm_weights, km.res$cluster)
})

pdf("3_calinski_harabasz_plot.pdf")
plot(k.values, ch_values, type = "b", pch = 19, frame = FALSE,
     xlab = "Number of Clusters (K)", ylab = "Calinski-Harabasz Index",
     main = "Calinski-Harabasz Index: Determining Optimal K")
dev.off()

# ISCHIA Analysis
pdf("4_composition_cluster_k_plot.pdf")
Composition.cluster.k(norm_weights, 20)
dev.off()

pdac <- Composition.cluster(pdac, norm_weights, 12)
pdac$cc_12 <- pdac$CompositionCluster_CC

# Spatial Dimension Plot
image_names <- c("IU_PDA_T1", "IU_PDA_T2", "IU_PDA_HM2", "IU_PDA_HM2_2", "IU_PDA_NP2", 
                 "IU_PDA_T3", "IU_PDA_HM3", "IU_PDA_T4", "IU_PDA_HM4", "IU_PDA_HM5", 
                 "IU_PDA_T6", "IU_PDA_HM6", "IU_PDA_LNM6", "IU_PDA_LNM7", "IU_PDA_T8", 
                 "IU_PDA_HM8", "IU_PDA_LNM8", "IU_PDA_T9", "IU_PDA_HM9", "IU_PDA_T10", 
                 "IU_PDA_HM10", "IU_PDA_LNM10", "IU_PDA_NP10", "IU_PDA_T11", "IU_PDA_HM11", 
                 "IU_PDA_NP11", "IU_PDA_T12", "IU_PDA_HM12", "IU_PDA_LNM12", "IU_PDA_HM13")

paletteMartin <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', 
                   '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff')

all_ccs <- unique(pdac$CompositionCluster_CC)
color_mapping <- setNames(paletteMartin[1:length(all_ccs)], all_ccs)

pdf("5_spatial_plots_K12.pdf", width = 10, height = 7)
for (image_name in image_names) {
  plot <- SpatialDimPlot(pdac, group.by = "CompositionCluster_CC", images = image_name) +
    scale_fill_manual(values = color_mapping) +
    theme_minimal() +
    ggtitle(image_name)
  print(plot)
}
dev.off()

# Enriched Cell Types
save_cc_plot <- function(cc) {
  plot <- Composition_cluster_enrichedCelltypes(pdac, cc, as.matrix(norm_weights))
  pdf_name <- paste0(cc, ".pdf")
  pdf(file = pdf_name)
  print(plot)
  dev.off()
}

ccs <- paste0("CC", 1:12)
for (cc in ccs) {
  save_cc_plot(cc)
}

pdf_files <- paste0("CC", 1:12, ".pdf")
pdf_combine(pdf_files, output = "6_enrichedCelltypes_CC_12.pdf")

# UMAP
pdac.umap <- Composition_cluster_umap(pdac, norm_weights)
pdf("7_umap_pie_chart.pdf")
print(pdac.umap$umap.deconv.gg)
dev.off()

# Add UMAP to Seurat object
emb.umap <- pdac.umap$umap.table
emb.umap$CompositionCluster_CC <- NULL
emb.umap$Slide <- NULL
emb.umap <- as.matrix(emb.umap)
colnames(emb.umap) <- c("UMAP1", "UMAP2")

pdac[['umap.ischia12']] <- CreateDimReducObject(embeddings = emb.umap, key = 'umap.ischia12_', assay = 'rctd_tier1')

pdf("8_seurat_ischia_umap_12.pdf")
DimPlot(pdac, reduction = "umap.ischia12", label = FALSE, group.by="cc_12")
dev.off()

# Bar plots
pdf("9_barplot_SampVsorig_12.pdf", height=12, width=20)
dittoBarPlot(pdac, "orig.ident", group.by = "cc_12")
dev.off()

pdf("10_barplot_origVsSamp_12.pdf", height=10, width=20)
dittoBarPlot(pdac, "cc_12", group.by = "orig.ident")
dev.off()

# Cell type co-occurrence
CC4.celltype.cooccur <- spatial.celltype.cooccurence(spatial.object=pdac, deconv.prob.mat=norm_weights, 
                                                     COI="CC4", prob.th= 0.05, 
                                                     Condition=unique(pdac$orig.ident))
pdf("11_celltype_cooccurrence_CC4.pdf")
plot.celltype.cooccurence(CC4.celltype.cooccur)
dev.off()
