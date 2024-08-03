# Load required libraries
library(spacexr)
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(progeny)

# Data Preprocessing

# Load the Spatial data
pdac <- readRDS("/path/to/pdac_most_updated.rds")

# Load the Single cell data
sc <- readRDS("/path/to/sc_rctd.rds")

# Clean up cell type names
sc$celltype_nicheDE <- gsub("T/NK cells", "T_NK cells", sc$celltype_nicheDE)
sc$celltype_nicheDE <- gsub("CAF-S1", "Fibroblasts", sc$celltype_nicheDE)

# Prepare data for RCTD
counts <- sc@assays$RNA@counts
Idents(sc) <- "celltype_nicheDE"
cluster <- as.factor(sc$celltype_nicheDE)
names(cluster) <- colnames(sc)
nUMI <- sc$nCount_RNA
names(nUMI) <- colnames(sc)

# Create the reference object
reference <- Reference(counts, cluster, nUMI)

# Prepare spatial transcriptomics data
counts <- pdac@assays$Spatial@counts

# Add coordinates for all images
image_names <- c("IU_PDA_HM9", "IU_PDA_HM10", "IU_PDA_HM11", "IU_PDA_HM12", "IU_PDA_HM13", 
                 "IU_PDA_HM2", "IU_PDA_HM3", "IU_PDA_HM4", "IU_PDA_HM5", "IU_PDA_HM6", 
                 "IU_PDA_HM8", "IU_PDA_LNM10", "IU_PDA_LNM12", "IU_PDA_LNM6", "IU_PDA_LNM7", 
                 "IU_PDA_LNM8", "IU_PDA_HM2_2", "IU_PDA_NP10", "IU_PDA_NP11", "IU_PDA_NP2", 
                 "IU_PDA_T1", "IU_PDA_T9", "IU_PDA_T10", "IU_PDA_T11", "IU_PDA_T12", 
                 "IU_PDA_T2", "IU_PDA_T3", "IU_PDA_T4", "IU_PDA_T6", "IU_PDA_T8")

coordinates_list <- lapply(image_names, function(image_name) {
  pos <- GetTissueCoordinates(pdac, image = image_name)
  colnames(pos) <- c('x','y')
  return(pos)
})

coords <- do.call(rbind, coordinates_list)
rownames(coords) <- gsub("^.+\\.", "", rownames(coords))

# Create SpatialRNA object
query <- SpatialRNA(coords, counts, colSums(counts))

# Run RCTD
RCTD <- create.RCTD(query, reference, max_cores = 8)
RCTD.full <- run.RCTD(RCTD, doublet_mode = "full")

# Add RCTD results to Seurat object
pdac <- AddMetaData(pdac, metadata = RCTD.full@results$results_df)

# Normalize weights
weights <- RCTD.full@results$weights
norm_weights <- normalize_weights(weights)

# Add RCTD results as a new assay
pdac[["rctd_full"]] <- CreateAssayObject(data = t(as.matrix(norm_weights)))
if (length(pdac@assays$rctd_full@key) == 0) {
    pdac@assays$rctd_full@key <- "rctd_full_"
}

# Plotting
DefaultAssay(pdac) <- "rctd_full"
cell_types <- c("B cells", "Fibroblasts", "Endothelial cells", "Epithelial cells", 
                "Hepatocytes", "Myeloid", "T_NK cells")

for (img in image_names) {
  plot <- SpatialFeaturePlot(pdac, features = cell_types, pt.size.factor = 1.6, 
                             ncol = 3, crop = TRUE, images = img)
  ggsave(paste0("plots/", img, "_rctd_full.png"), plot, width = 12, height = 8, dpi = 300)
}

# PROGENy analysis
pdac_progeny <- progeny(pdac, scale=FALSE, organism="Human", top=1000, 
                        perm=1, return_assay = TRUE, assay="Spatial")
pdac_progeny <- ScaleData(pdac_progeny, assay = "progeny")

progeny_scores_df <- as.data.frame(t(GetAssayData(pdac_progeny, slot = "scale.data", 
                                                  assay = "progeny"))) %>%
  rownames_to_column("Cell") %>%
  gather(Pathway, Activity, -Cell)

CellsClusters <- data.frame(Cell = names(Idents(pdac)), 
                            CellType = as.character(Idents(pdac)), 
                            stringsAsFactors = FALSE)

progeny_scores_df <- inner_join(progeny_scores_df, CellsClusters)

summarized_progeny_scores <- progeny_scores_df %>% 
  group_by(Pathway, CellType) %>%
  summarise(avg = mean(Activity), std = sd(Activity))

summarized_progeny_scores_df <- summarized_progeny_scores %>%
  select(-std) %>%   
  spread(Pathway, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE) 

# Heatmap of PROGENy scores
paletteLength <- 100
myColor <- colorRampPalette(c("#008080", "white","#FFA500"))(paletteLength)

progenyBreaks <- c(seq(min(summarized_progeny_scores_df), 0, 
                       length.out=ceiling(paletteLength/2) + 1),
                   seq(max(summarized_progeny_scores_df)/paletteLength, 
                       max(summarized_progeny_scores_df), 
                       length.out=floor(paletteLength/2)))

progeny_hmap <- pheatmap(t(summarized_progeny_scores_df[,-1]), fontsize=14, 
                         fontsize_row = 10, color=myColor, breaks = progenyBreaks, 
                         main = "PROGENy (1000)", angle_col = 45,
                         treeheight_col = 0,  border_color = NA)

# Save the heatmap
ggsave("progeny_heatmap.png", progeny_hmap, width = 12, height = 8, dpi = 300)
