# R version - 4.5.3
suppressPackageStartupMessages({
  library(Seurat) # 5.5.0
  library(SeuratObject) # 5.4.0
  library(dplyr) # 1.2.1
  library(patchwork) # 1.3.2
  library(ggplot2) # 4.0.3
  library(here) # 1.0.2
  library(spacexr) # 2.2.1
  library(SpatialExperiment) # 1.40.0
  library(SummarizedExperiment) # 1.20.0
  library(data.table) # 1.18.4
})

project_root <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
if (getwd() != project_root) {
  setwd(project_root)
}

source(here("code","utils.R"))
set.seed(123)

# Load scRNAseq reference dataset
reference_data <- readRDS(here("data","GSE149614_final_ref.RDS"))
reference_data$custom_celltype <- gsub("/","-",reference_data$custom_celltype)
counts <- SeuratObject::LayerData(reference_data, 
                                  assay = "RNA", 
                                  layer = "counts")
cell_types <- as.factor(reference_data$custom_celltype)
names(cell_types) <- colnames(reference_data)
n_cell_types <- length(levels(cell_types))
rctd_reference <- Reference(counts, 
                            cell_types,
                            reference_data$nCount_RNA)

# Load spatial query dataset
spatial_data <- readRDS(here("data","merged_spatial_processed_harmony.RDS"))
sample_names <- unique(spatial_data$orig.ident)
all_proportions <- list()
for (sample_name in sample_names) {
  selected_sample_data <- subset(spatial_data, orig.ident == sample_name)
  spatial_counts <- SeuratObject::LayerData(selected_sample_data, 
                                            assay = "Spatial", 
                                            layer = "counts")
  tissue_positions <- GetTissueCoordinates(selected_sample_data, 
                                           image = Images(selected_sample_data)[1]) %>% 
    select(x,y)
  
  rctd_spatial <- SpatialRNA(tissue_positions,
                             counts = spatial_counts, 
                             nUMI = selected_sample_data$nCount_Spatial)
  
  rctd_obj <- create.RCTD(rctd_spatial, rctd_reference, max_cores = 8)
  rctd_obj <- run.RCTD(rctd_obj, doublet_mode = "full")
  gc()
  
  sample_proportions <- rctd_obj@results$weights 
  # Normalize weights to obtain proportions
  sample_proportions <- sweep(sample_proportions, 1, 
                              rowSums(sample_proportions), `/`) %>% 
    as.data.frame()
  
  all_proportions[[sample_name]] <- sample_proportions
}

rctd_proportions <- bind_rows(all_proportions)
saveRDS(rctd_proportions, here("results", "rctd_proportions.RDS"))

# all(rownames(spatial_data@meta.data) == rownames(rctd_proportions))
# TRUE

spatial_data <- AddMetaData(spatial_data, metadata = rctd_proportions)
saveRDS(spatial_data, file = here("data", "merged_spatial_annotated.RDS"))

