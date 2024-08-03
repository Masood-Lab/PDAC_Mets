# Load required libraries
suppressPackageStartupMessages({
  library(escape)
  library(Seurat)
  library(dittoSeq)
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(patchwork)
  library(pheatmap)
})

# Read data
print("Reading pdac")
a <- readRDS("/path/to/Pdac_updated.rds")
Idents(a) <- a$cc_ischia_10
print(a)

# Define gene sets
gene.sets <- list(Angiogenesis = c('ANGPT1','ANGPT2','CDH5','CXCL5','CXCL8','CXCR2','FLT1','KDR','PDGFC','PGF','TEK','VEGFA','VEGFB','VEGFC','VWF'), Antitumor_cytokines = c('CCL3','IFNA2','IFNB1','IL21','TNF','TNFSF10'), Bcells = c('BLK','CD19','CD22','CD79A','CD79B','CR2','FCRL5','MS4A1','PAX5','STAP1','TNFRSF13B','TNFRSF13C','TNFRSF17'), Cancer_associated_fibroblasts = c('ACTA2','CD248','COL11A1','COL1A1','COL1A2','COL5A1','COL6A1','COL6A2','COL6A3','CXCL12','FAP','FBLN1','FGF2','FN1','LRP1','LUM','MFAP5','MMP2','MMP3','PDGFRA','PDGFRB'), Checkpoint_molecules = c('BTLA','CD274','CTLA4','HAVCR2','LAG3','PDCD1','PDCD1LG2','TIGIT','VSIR'), Co_stimulatory_ligands =  c('CD40LG','CD70','CD80','CD83','CD86','ICOSLG','TNFSF4','TNFSF9'), Co_stimulatory_receptors = c('CD27','CD28','CD40','ICOS','TNFRSF4','TNFRSF9'), Effector_cell_traffic = c('CCL3','CCL4','CCL5','CX3CL1','CX3CR1','CXCL10','CXCL11','CXCL9','CXCR3'), Effector_cells = c('CD8A','CD8B','EOMES','FASLG','GNLY','GZMA','GZMB','GZMK','IFNG','PRF1','TBX21','ZAP70'), Endothelium = c('CDH5','CLEC14A','ENG','FLT1','KDR','MMRN1','MMRN2','NOS3','VCAM1','VWF'), Granulocyte_traffic = c('CCL11','CCR3','CXCL1','CXCL2','CXCL5','CXCL8','CXCR1','CXCR2','KITLG'), Immune_Suppression_by_Myeloid_cells = c('ARG1','CYBB','IDO1','IL10','IL4I1','IL6','PTGS2'), M1_signature = c('CMKLR1','IL12A','IL12B','IL1B','IL23A','IRF5','NOS2','SOCS3','TNF'), Macrophage_and_DC_traffic = c('CCL2','CCL7','CCL8','CCR2','CSF1','CSF1R','XCL1','XCR1'), Matrix = c( 'COL11A1','COL1A1','COL1A2','COL3A1','COL4A1','COL5A1','ELN','FN1','LAMA3','LAMB3','LAMC2','LGALS7','LGALS9','TNC','VTN'), Matrix_remodeling   = c('ADAMTS4','ADAMTS5','CA9','LOX','MMP1','MMP11','MMP12','MMP2','MMP3','MMP7','MMP9','PLOD2'), MHCI = c('B2M','HLA_A','HLA_B','HLA_C','NLRC5','TAP1','TAP2','TAPBP'), MHCII = c('CIITA','HLA_DMA','HLA_DMB','HLA_DPA1','HLA_DPB1','HLA_DQA1','HLA_DQB1','HLA_DRA','HLA_DRB1'), Myeloid_cells_traffic = c('CCL15','CCL26','CSF1','CSF1R','CSF2','CSF2RA','CSF3','CSF3R','CXCL12','CXCL5','CXCL8','CXCR2','CXCR4','IL6','IL6R'), Neutrophil_signature = c('CD177','CTSG','CXCR1','CXCR2','ELANE','FCGR3B','FFAR2','MPO','PGLYRP1','PRTN3'), NK_cells  =c('CD160','CD226','CD244','EOMES','FGFBP2','GNLY','GZMB','GZMH','IFNG','KIR2DL4','KLRC2','KLRF1','KLRK1','NCR1','NCR3','NKG7','SH2D1B'), Protumor_cytokines = c('IL10','IL22','IL6','MIF','TGFB1','TGFB2','TGFB3'), Tcells = c('CD28','CD3D','CD3E','CD3G','CD5','ITK','TBX21','TRAC','TRAT1','TRBC1','TRBC2'), Th1_signature = c('CD40LG','IFNG','IL12RB2','IL2','IL21','STAT4','TBX21'), Th2_signature = c('CCR4','IL10','IL13','IL4','IL5'), Treg = c('CCR8','CTLA4','FOXP3','IKZF2','IKZF4','IL10','TNFRSF18'), Treg_and_Th2_traffic = c('CCL1','CCL17','CCL22','CCL28','CCR10','CCR4','CCR8'), Tumor_proliferation_rate = c('AURKA','AURKB','BUB1','CCNB1','CCND1','CCNE1','CDK2','CETN3','E2F1','ESCO2','MCM2','MCM6','MKI67','MYBL2','PLK1'), Tumor_associated_Macrophages = c('CD163','CD68','CSF1R','IL10','IL4I1','MRC1','MSR1','SIGLEC1'))

# Perform enrichment analysis
print("Starting Es")
DefaultAssay(a) <- "RNA"
Es <- enrichIt(obj = a, gene.sets = gene.sets, groups = 1000, cores = 8)
saveRDS(Es, "ES.rds")
print("Es Done!")

# Add enrichment scores to Seurat object
a[['scfea']] <- Es %>%
  pivot_longer(cols = everything(), names_to = "source", values_to = "score") %>%
  pivot_wider(id_cols = 'source', names_from = '.', values_from = 'score') %>%
  column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)

DefaultAssay(a) <- "scfea"
a <- ScaleData(a)
a@assays$scfea@data <- a@assays$scfea@scale.data
Idents(a) <- a$cc_ischia_10

# Prepare data for heatmap
df <- t(as.matrix(a@assays$scfea@data)) %>%
  as.data.frame() %>%
  mutate(cluster = Idents(a)) %>%
  pivot_longer(cols = -cluster, names_to = "source", values_to = "score") %>%
  group_by(cluster, source) %>%
  summarise(mean = mean(score), .groups = "drop")

top_acts_mat <- df %>%
  pivot_wider(id_cols = 'cluster', names_from = 'source', values_from = 'mean') %>%
  column_to_rownames('cluster') %>%
  as.matrix()

# Create heatmap
palette_length <- 100
my_color <- colorRampPalette(c("darkblue", "white", "red"))(palette_length)
my_breaks <- c(seq(-2, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 2, length.out=floor(palette_length/2)))

pdf("FGSEA_histology.pdf", width = 10, height = 6)
pheatmap(
  top_acts_mat,
  border_color = NA,
  color = my_color,
  breaks = my_breaks,
  cellwidth = 20,
  cellheight = 20,
  fontsize_row = 9,
  fontsize_col = 9,
  main = "Heatmap FGES Genesets",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  show_colnames = TRUE,
  show_rownames = TRUE,
  angle_col = 90
)
dev.off()

# Morpheus heatmap (if needed)
# Note: The 'morpheus' function is not a standard R function. 
# Make sure you have the correct package loaded for this to work.
rowAnnotations <- data.frame(
  annotation1 = 1:nrow(top_acts_mat),
  annotation2 = rownames(top_acts_mat)
)

morpheus(top_acts_mat, 
         colorScheme = list(scalingMode = "fixed", colors = heat.colors(3)), 
         rowAnnotations = rowAnnotations, 
         overrideRowDefaults = FALSE, 
         rows = list(list(field = 'annotation2', 
                          highlightMatchingValues = TRUE, 
                          display = list('color'))))

morpheus(top_acts_mat, 
         dendrogram='column', 
         colorScheme=list(scalingMode="fixed", colors=heat.colors(3)), 
         rowAnnotations=rowAnnotations, 
         tools=list(list(name='Hierarchical Clustering', 
                         params=list(group_rows_by=list('annotation1'), cluster='Rows'))), 
         rowGroupBy=list(list(field='annotation2')),
         rows=list(list(field='annotation2',display=list('color'))))
