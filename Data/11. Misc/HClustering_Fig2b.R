### This script performs the wilcoxon rank sum test and hierarchical clustering on the RCTD tier data to identify the significant abundant cell types between the clusters.
  
#### Load necessary packages
 
```{r}
library(Seurat)
library(compositions)
library(tidyverse)
library(clustree)
library(patchwork)
library(uwot)
library(scran)
library(cluster)
library(ggrastr)
library(cowplot)
# library(conflicted) # to be loaded in case of a conflict arises.
 
config <- config::get()
 
# source(here::here("pdac_nac", "visualization", "eda.R"))
```
 
### Load the seurat object and get the proportions data
 
```{r}
so <- readRDS(here::here(config$data_processed, "06-pdac_CC10_msig.rds"))
 
# Join with metadata if needed
metadata <- so@meta.data %>%
  select(orig.ident, patient_id, neoadjuvant_chemo, CompositionCluster_CC) %>%
  rownames_to_column("row_id")
 
 
# Get the proportions data
rctd_tier2 <- t(so@assays$rctd_tier2@data)
 
# Ensure the data is in the right format
rownames(rctd_tier2) <- make.unique(rownames(rctd_tier2))
 
# Log transformation of rctd_tier1
log_comps <- log10(rctd_tier2)
```
 
 
 
## Perform the summary statistics
We perform the summary statistics for the RCTD tier data. We perform hierarchical clustering and do the wilcoxon rank sum test to identify the differentially abundant cell types between the clusters.
 
 
#### Prepare the data for the summary statistics
 
```{r}
# Prepare data for summary statistics
cluster_summary_pat <- rctd_tier2 %>%
  as.data.frame() %>%
  rownames_to_column("row_id") %>%
  left_join(metadata, by = "row_id") %>%  # Join with meta_data using row_id as the key
  pivot_longer(-c(row_id, orig.ident, patient_id, neoadjuvant_chemo, CompositionCluster_CC), values_to = "ct_prop", names_to = "cell_type") %>%
  group_by(orig.ident, patient_id, neoadjuvant_chemo, CompositionCluster_CC, cell_type) %>%
  summarize(median_ct_prop = median(ct_prop, na.rm = TRUE))
```
 
 
```{r}
# Aggregate data for median ct prop
cluster_summary <- cluster_summary_pat %>%
  ungroup() %>%
  group_by(CompositionCluster_CC, cell_type) %>%
  summarize(patient_median_ct_prop = median(median_ct_prop, na.rm = TRUE))
 
# Prepare matrix for hierarchical clustering
cluster_summary_mat <- cluster_summary %>%
  pivot_wider(values_from = patient_median_ct_prop, names_from = cell_type, values_fill = list(patient_median_ct_prop = 0)) %>%
  column_to_rownames("CompositionCluster_CC") %>%
  as.matrix()
 
# conflicts_prefer(stats::"dist") # resolve conflicts between %*% functions
```
 
 
```{r}
# Perform hierarchical clustering
cluster_order <- hclust(dist(cluster_summary_mat))$labels[hclust(dist(cluster_summary_mat))$order] # use this if you want to order the clusters based on the hierarchical clustering
ct_order <- hclust(dist(t(cluster_summary_mat)))$labels[hclust(dist(t(cluster_summary_mat)))$order]
 
# Order Clusters in ascending order
# cluster_order1 <- c("CC1", "CC2", "CC3", "CC4", "CC5", "CC6", "CC7", "CC8", "CC9", "CC10")
 
# Wilcoxon test for characteristic cell types
run_wilcox_up <- function(prop_data) {
  prop_data_group <- prop_data[["CompositionCluster_CC"]] %>% unique() %>% set_names()
  map(prop_data_group, function(g) {
    test_data <- prop_data %>%
      mutate(test_group = ifelse(CompositionCluster_CC == g, "target", "rest")) %>%
      mutate(test_group = factor(test_group, levels = c("target", "rest")))
    wilcox.test(median_ct_prop ~ test_group, data = test_data, alternative = "greater") %>%
      broom::tidy()
  }) %>% enframe("CompositionCluster_CC") %>% unnest()
}
```
 
 
```{r}
 
wilcoxon_res <- cluster_summary_pat %>%
  ungroup() %>%
  group_by(cell_type) %>%
  nest() %>%
  mutate(wres = map(data, run_wilcox_up)) %>%
  dplyr::select(wres) %>%
  unnest() %>%
  ungroup() %>%
  mutate(p_corr = p.adjust(p.value)) %>%
  mutate(significant = ifelse(p_corr <= 0.15, "*", ""))
```
 
#### Save the summary statistics
 
```{r}
# give the path to save the summary statistics
file_path_cluster_summ <- here::here(config$data_interim,  "summary_of_clusters.txt")
file_path_wilcox_res <- here::here(config$data_interim, ß "wilcoxon_res_cells_clusters.txt")
 
# Save the summary statistics
write.table(cluster_summary_pat, file = file_path_cluster_summ, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(wilcoxon_res, file = file_path_wilcox_res, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
```
 
#### Plot the summary statistics
 
```{r}
# Plotting mean ct prop and barplots
mean_ct_prop_plt <- cluster_summary %>%
  left_join(wilcoxon_res, by = c("CompositionCluster_CC", "cell_type")) %>%
  mutate(cell_type = factor(cell_type, levels = ct_order), CompositionCluster_CC = factor(CompositionCluster_CC, levels = cluster_order)) %>%
  ungroup() %>%
  group_by(cell_type) %>%
  mutate(scaled_pat_median = (patient_median_ct_prop - mean(patient_median_ct_prop)) / sd(patient_median_ct_prop)) %>%
  ungroup() %>%
  ggplot(aes(x = cell_type, y = CompositionCluster_CC, fill = scaled_pat_median)) +
  geom_tile(color = "black") +
  geom_text(aes(label = significant)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12), legend.position = "bottom", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.text.y = element_text(size = 12)) +
  scale_fill_gradient2()
 
cluster_counts <- cluster_info %>%
  dplyr::select_at(c("row_id", "CompositionCluster_CC")) %>%
  group_by(CompositionCluster_CC) %>%
  summarize(nspots = length(CompositionCluster_CC)) %>%
  mutate(prop_spots = nspots / sum(nspots))
 
file_path_cluster_prop_summ <- here::here(config$data_interim, "cluster_prop_summary.csv")
 
write_csv(cluster_counts, file_path_cluster_prop_summ)
 
#barplots for cluster counts
barplts <- cluster_counts %>%
  mutate(CompositionCluster_CC = factor(CompositionCluster_CC, levels = cluster_order)) %>%
  ggplot(aes(y = CompositionCluster_CC, x = prop_spots)) +
  geom_bar(stat = "identity") +
  theme_classic() + ylab("") +
  theme(axis.text.y = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.text.x = element_text(size = 12))
 
cluster_summary_plt <- cowplot::plot_grid(mean_ct_prop_plt, barplts, align = "hv", axis = "tb") # use if barplots are needed to show the spot counts otherwise directly use mean_ct_prop_plt for the plot
```
 
#### plot the summary clusters
 
```{r}
pdf_path_summ_clust <- here::here(config$plots, "wilcox_summary_clusters.pdf")
 
pdf(pdf_path_summ_clust, width = 20, height = 10)
plot(cluster_summary_plt)
dev.off()
 
#box plot for median ct prop
pdf_path_boxplot <- here::here(config$plots,  "wilcox_bboxplot_median_ct_prop.pdf")
 
plt <- cluster_summary_pat %>%
  ggplot(aes(x = CompositionCluster_CC, y = median_ct_prop)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  facet_wrap(. ~ cell_type, ncol = 3, scales = "free_y")
 
pdf(pdf_path_boxplot, width = 20, height = 10)
plot(plt)
dev.off()
```
 
