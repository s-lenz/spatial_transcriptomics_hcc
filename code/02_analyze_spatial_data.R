suppressPackageStartupMessages({
  library(utils) # 4.5.3
  library(R.utils) # 2.13.0
  library(Seurat) # 5.5.0
  library(dplyr) # 1.2.1
  library(patchwork) # 1.3.2
  library(ggplot2) # 4.0.3
  library(spacexr)
})

source("code/utils.R")
set.seed(123)

## Load data files
# spatial
merged <- readRDS(file.path("data", "merged_spatial_processed_harmony.RDS"))

# reference
ref <- readRDS(file.path("data", "GSE149614_final_ref.RDS"))

table(ref$custom_celltype_clean)

cell_types <- ref$custom_celltype_clean
names(cell_types) <- colnames(ref)
table(cell_types)

ref_counts <- GetAssayData(ref, assay = "RNA", layer = "counts")

# make sure cells line up
ref_counts <- ref_counts[, names(cell_types)]
all(colnames(ref_counts) == names(cell_types))


# create RCTD reference object 
reference <- Reference(
  counts = ref_counts,
  cell_types = as.factor(cell_types)
)

## Loop over samples since RCTD runs per section (orig.ident)
sample_ids <- names(table(merged$orig.ident))
rctd_results <- list()   # one RCTD object per sample, named by orig.ident

# creat results folder
dir.create("results",         showWarnings = FALSE)
dir.create("results/figures", showWarnings = FALSE)

for (sample_id in sample_ids) {
  
  message("\n=== Running RCTD for sample: ", sample_id, " ===")
  
  merged_one <- subset(merged, subset = orig.ident == sample_id)
  
  count_layer    <- grep("^counts", Layers(merged_one[["Spatial"]]), value = TRUE)
  spatial_counts <- GetAssayData(merged_one, assay = "Spatial", layer = count_layer)
  
  coords      <- GetTissueCoordinates(merged_one)
  coords_rctd <- coords[, c("x", "y")]
  
  # CHECK: every count spot must have a coordinate
  spots_matched <- sum(colnames(spatial_counts) %in% rownames(coords_rctd))
  spots_total   <- ncol(spatial_counts)
  if (spots_matched != spots_total) {
    stop("Sample ", sample_id, ": only ", spots_matched, " of ", spots_total,
         " count spots have matching coordinates. Stopping the pipeline.")
  }
  message("  check passed: ", spots_matched, "/", spots_total, " spots have coordinates")
  
  common_spots   <- intersect(colnames(spatial_counts), rownames(coords_rctd))
  spatial_counts <- spatial_counts[, common_spots]
  coords_rctd    <- coords_rctd[common_spots, ]
  
  nUMI <- colSums(spatial_counts)
  puck <- SpatialRNA(coords = coords_rctd, counts = spatial_counts, nUMI = nUMI)
  
  myRCTD <- create.RCTD(spatialRNA = puck, reference = reference,
                        max_cores = 4, test_mode = FALSE)
  myRCTD <- run.RCTD(myRCTD, doublet_mode = "doublet")
  
  rctd_results[[sample_id]] <- myRCTD
  
  ## Analysis: pull results, write back into the object
  results_df <- myRCTD@results$results_df
  weights    <- myRCTD@results$weights_doublet
  
  # print a record to the log per sample
  cat("\n--- ", sample_id, " spot classes ---\n"); print(table(results_df$spot_class))
  cat("--- ", sample_id, " first types ---\n");   print(table(results_df$first_type))
  
  # create new metadata categories for primary and secondary cell types and weights
  merged_one$RCTD_spot_class   <- NA
  merged_one$RCTD_first_type   <- NA
  merged_one$RCTD_second_type  <- NA
  merged_one$RCTD_first_weight <- NA
  merged_one$RCTD_second_weight<- NA
  merged_one$RCTD_spot_class[rownames(results_df)]   <- as.character(results_df$spot_class)
  merged_one$RCTD_first_type[rownames(results_df)]   <- as.character(results_df$first_type)
  merged_one$RCTD_second_type[rownames(results_df)]  <- as.character(results_df$second_type)
  merged_one$RCTD_first_weight[rownames(weights)]    <- weights[, "first_type"]
  merged_one$RCTD_second_weight[rownames(weights)]   <- weights[, "second_type"]
  
  ## Plots: save to file instead of displaying
  p1 <- SpatialDimPlot(merged_one, group.by = "RCTD_first_type",  pt.size.factor = 3)
  p2 <- SpatialDimPlot(merged_one, group.by = "RCTD_spot_class",  pt.size.factor = 3)
  p3 <- SpatialDimPlot(merged_one, group.by = "RCTD_second_type", pt.size.factor = 3)
  
  ggsave(paste0("results/figures/", sample_id, "_first_type.png"),  p1, width = 8, height = 7, dpi = 150)
  ggsave(paste0("results/figures/", sample_id, "_spot_class.png"),  p2, width = 8, height = 7, dpi = 150)
  ggsave(paste0("results/figures/", sample_id, "_second_type.png"), p3, width = 8, height = 7, dpi = 150)
  
  ## Save objects
  saveRDS(myRCTD,     file = paste0("results/RCTD_",            sample_id, ".rds"))
  saveRDS(merged_one, file = paste0("results/merged_one_RCTD_", sample_id, ".rds"))
  
  message("  saved results, annotated object, and 3 figures for ", sample_id)
}

