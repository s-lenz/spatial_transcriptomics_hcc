# R version - 4.5.3
suppressPackageStartupMessages({
  library(Seurat) # 5.5.0
  library(SeuratObject) # 5.4.0
  library(dplyr) # 1.2.1
  library(patchwork) # 1.3.2
  library(ggplot2) # 4.0.3
  library(here) # 1.0.2
  library(data.table) # 1.18.4
  library(pheatmap) # 1.0.13
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
cell_types <- colnames(annotated_spatial_data@meta.data)[(total_columns-9):total_columns]

# Plot proportions on spatial tissue slides
for (sample_name in sample_names) {
  
  tissue_positions <- GetTissueCoordinates(annotated_spatial_data, 
                                           image = sample_name) %>% 
    select(x,y)
  sample_proportions <- annotated_spatial_data@meta.data %>% filter(orig.ident == sample_name) %>%
    select(all_of(cell_types))
  
  p <- plot_Proportions_RCTD(tissue_positions, 
                              sample_proportions, 
                              paste0("Sample ID: ", sample_name))
  ggsave(here("results","plots",paste0(sample_name,"_rctd_proportions.png")),p,
         width = 14, height = 15, dpi = 150, units = "in", device='png')
  
}