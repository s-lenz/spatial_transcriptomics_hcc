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

# Prepare reference dataset
# Load scRNAseq data
reference_data <- readRDS(here("data","GSE149614_final_ref.RDS"))
reference_data$custom_celltype <- gsub("/","-",reference_data$custom_celltype)
# Extract raw counts
counts <- SeuratObject::LayerData(reference_data, 
                                  assay = "RNA", 
                                  layer = "counts")
# Get factorized vector with cell types annotations
cell_types <- as.factor(reference_data$custom_celltype)
names(cell_types) <- colnames(reference_data)
n_cell_types <- length(levels(cell_types))
# Create RCTD reference
rctd_reference <- Reference(counts, 
                            cell_types,
                            reference_data$nCount_RNA)

# Load spatial query dataset
merged <- readRDS(here("data","merged_spatial_processed_harmony.RDS"))
sample_names <- unique(merged$orig.ident)
all_proportions <- list()
# For each sample, run RCTD algorithm and obtain estimated cell type proportions
for (sample_name in sample_names) {
  # Filter merged Seurat object
  selected_sample_data <- subset(merged, orig.ident == sample_name)
  # Extract raw counts
  spatial_counts <- SeuratObject::LayerData(selected_sample_data, 
                                            assay = "Spatial", 
                                            layer = "counts")
  # Get spot coordinates
  tissue_positions <- GetTissueCoordinates(selected_sample_data, 
                                           image = Images(selected_sample_data)[1]) %>% 
    select(x,y)
  # Create RCTD spatial data
  rctd_spatial <- SpatialRNA(tissue_positions,
                             counts = spatial_counts, 
                             nUMI = selected_sample_data$nCount_Spatial)
  # Create RCTD object from reference and spatial data
  # You can specify number of core to enable parallelization
  rctd_obj <- create.RCTD(rctd_spatial, rctd_reference, max_cores = 8)
  # Run RCTD deconvolution in full mode 
  # We are assuming that spot can contain any number of reference cell type
  rctd_obj <- run.RCTD(rctd_obj, doublet_mode = "full")
  gc()
  
  # Get estimated RCTD model weights
  sample_proportions <- rctd_obj@results$weights 
  # Normalize weights to obtain proportions (sum to 100%)
  sample_proportions <- sweep(sample_proportions, 1, 
                              rowSums(sample_proportions), `/`) %>% 
    as.data.frame()
  # Save sample proportions
  all_proportions[[sample_name]] <- sample_proportions
}

rctd_proportions <- bind_rows(all_proportions)
saveRDS(rctd_proportions, here("results", "rctd_proportions.RDS"))

# all(rownames(merged@meta.data) == rownames(rctd_proportions))
# TRUE

merged <- AddMetaData(merged, metadata = rctd_proportions)
saveRDS(merged, file = here("data", "merged_spatial_annotated.RDS"))

