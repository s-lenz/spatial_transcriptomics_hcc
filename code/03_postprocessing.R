# R version - 4.5.3
suppressPackageStartupMessages({
  library(Seurat) # 5.5.0
  library(SeuratObject) # 5.4.0
  library(dplyr) # 1.2.1
  library(tidyr) # 1.3.2
  library(patchwork) # 1.3.2
  library(ggplot2) # 4.0.3
  library(here) # 1.0.2
  library(data.table) # 1.18.4
  library(pheatmap) # 1.0.13
  library(Polychrome) # 1.6.1
  library(cowplot) # 1.2.0
  library(TCA) # 1.2.1
  library(Matrix) # 1.7-5
  library(stringr) # 1.6.0
  library(tibble) # 3.3.1
})

project_root <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
if (getwd() != project_root) {
  setwd(project_root)
}

source(here("code","utils.R"))
set.seed(123)

annotated_spatial_data <- readRDS(here("data", "merged_spatial_annotated.RDS"))
sample_names <- unique(annotated_spatial_data$orig.ident)

total_columns <- ncol(annotated_spatial_data@meta.data)
all_proportions <- annotated_spatial_data@meta.data[(total_columns-9):total_columns]
cell_types <- colnames(all_proportions)
all_proportions$sample_name <- annotated_spatial_data$orig.ident

# Plot proportions on spatial tissue slides
for (selected_sample in sample_names) {
  
  tissue_positions <- GetTissueCoordinates(annotated_spatial_data, 
                                           image = selected_sample) %>% 
    select(x,y)
  sample_proportions <- all_proportions %>% filter(sample_name == selected_sample) %>%
    select(all_of(cell_types))
  
  p <- plot_Proportions_RCTD(tissue_positions, 
                              sample_proportions, 
                              paste0("Sample ID: ", selected_sample))
  ggsave(here("results","plots",paste0(selected_sample,"_rctd_proportions.png")),p,
         width = 14, height = 15, dpi = 150, units = "in", device='png')
  
}

# Heatmaps with proportions by sample (per cell)
for (selected_sample in sample_names) {
  
  sample_proportions <- all_proportions %>% 
    filter(sample_name == selected_sample) %>% 
    select(all_of(cell_types))

  p <- plot_Heatmap_RCTD(sample_proportions, selected_sample)
  ggsave(here("results","plots",paste0(selected_sample,"_rctd_heatmap.png")),p,
         width = 7, height = 6, dpi = 150, units = "in", device='png')
}

# Spatial slides with annotated major cell types per spot
annotated_spatial_data@meta.data <- annotated_spatial_data@meta.data %>%
  mutate(main_cell_type = names(.[cell_types])[max.col(.[cell_types], ties.method = "first")])

p <- plotSpatialDim_by_celltypes(annotated_spatial_data, 
                            "orig.ident", 
                            cell_types, pt.size = 4)
p
ggsave(here("results","plots","cell_types_spatial.png"),p,
       width = 14, height = 6, dpi = 150, units = "in", device='png')

# Naive gene expression imputation based on proportion matrix
spatial_counts_sct <- as.matrix(LayerData(annotated_spatial_data, 
                                      assay = "SCT", 
                                      layer = "data"))

for (target_cell_type in cell_types) {
  print(target_cell_type)
  props_vector <- all_proportions[, target_cell_type]
  weight_diag <- Diagonal(x = props_vector)
  pseudo_expr_matrix <- spatial_counts_sct %*% weight_diag
  colnames(pseudo_expr_matrix) <- colnames(spatial_counts_sct)
  rownames(pseudo_expr_matrix) <- rownames(spatial_counts_sct)
  assay_name <- paste0("Pseudo_",
                       str_replace_all(target_cell_type, " ", "_"),
                       "_SCT")
  annotated_spatial_data[[assay_name]] <- CreateAssay5Object(data = pseudo_expr_matrix)
}

# Plot for FGF19 total expression vs expression by cell type (based on proportions)
p <- plotGEX_by_celltypes(annotated_spatial_data, 
                     "FGF19", 
                     "orig.ident", 
                     cell_types)
p
ggsave(here("results","plots","FGF19_GEX_by_celltype.png"),p,
       width = 18, height = 10, dpi = 150, units = "in", device='png')

# Heatmap illustrating cluster enrichment for estimated cell types
median_proportions <- annotated_spatial_data@meta.data %>% 
  group_by(integrated_clusters, orig.ident) %>%
  summarize(across(all_of(cell_types), median, 
                   .names = "{.col}"), 
            .groups = "drop")

p <- plotClusterHeatmap_by_celltype(median_proportions,
                               cell_types,
                               n_samples = length(unique(annotated_spatial_data$orig.ident)))

p
ggsave(here("results","plots","cluster_celltypes_heatmap.png"),p,
       width = 8, height = 3, dpi = 150, units = "in", device='png')


