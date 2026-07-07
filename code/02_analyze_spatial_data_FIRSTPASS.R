# This runs one sample at a time with no sample looping. Used as an initial 
# check for the data processing. Control sample looks good (mostly hepatocytes)

# R version - 4.5.3
suppressPackageStartupMessages({
  library(GEOquery) # 2.78.0
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

merged <- readRDS(file.path("data", "merged_spatial_processed_harmony.RDS"))

ref <- readRDS("/fs/ess/PAS0854/Sidney/classes/BMI5750/spatial_transcriptomics_hcc/data/GSE149614_final_ref.RDS")

table(ref$custom_celltype_clean)

cell_types <- ref$custom_celltype_clean
names(cell_types) <- colnames(ref)
table(cell_types)

ref_counts <- GetAssayData(ref, assay = "RNA", layer = "counts")

# make sure cells line up
ref_counts <- ref_counts[, names(cell_types)]
all(colnames(ref_counts) == names(cell_types))

# create RTCD reference object
reference <- Reference(
  counts = ref_counts,
  cell_types = as.factor(cell_types)
)


## Pull raw counts from visium
Assays(merged)
Layers(merged[["Spatial"]])

# run one sample first
table(merged$orig.ident)

sample_id <- names(table(merged$orig.ident))[1]
merged_one <- subset(merged, subset = orig.ident == sample_id)
Layers(merged_one[["Spatial"]])

# get raw counts
count_layer <- grep("^counts", Layers(merged_one[["Spatial"]]), value = TRUE)

spatial_counts <- GetAssayData(
  merged_one,
  assay = "Spatial",
  layer = count_layer
)

dim(spatial_counts)
head(colnames(spatial_counts))

# get coordinates
coords <- GetTissueCoordinates(merged_one)

dim(coords)
head(coords)
colnames(coords)
head(rownames(coords))

# compare to counts
head(colnames(spatial_counts))

sum(colnames(spatial_counts) %in% rownames(coords))
ncol(spatial_counts) # should equal sum

# ADD CHECK: true/false? else exit with warning message

coords_rctd <- coords[, c("x", "y")]

common_spots <- intersect(colnames(spatial_counts), rownames(coords_rctd))

spatial_counts <- spatial_counts[, common_spots]
coords_rctd <- coords_rctd[common_spots, ]

all(colnames(spatial_counts) == rownames(coords_rctd))

# create spatial RCTD object
# they call it puck
nUMI <- colSums(spatial_counts)

puck <- SpatialRNA(
  coords = coords_rctd,
  counts = spatial_counts,
  nUMI = nUMI
)

## Run RCTD
myRCTD <- create.RCTD(
  spatialRNA = puck,
  reference = reference,
  max_cores = 4,
  test_mode = FALSE
)

myRCTD <- run.RCTD(
  myRCTD,
  doublet_mode = "doublet" # first pass, other options are singlet and full
)

# look at results
results_df <- myRCTD@results$results_df
weights <- myRCTD@results$weights_doublet

head(results_df)
table(results_df$spot_class)
dim(weights)
head(weights)

# which cell types called overall?
table(results_df$first_type)
table(results_df$second_type)

summary(weights[, "first_type"])
summary(weights[, "second_type"])

# confidence?
table(results_df$spot_class, results_df$first_type)
table(results_df$spot_class, results_df$second_type)

# call types
merged_one$RCTD_spot_class <- NA
merged_one$RCTD_first_type <- NA
merged_one$RCTD_second_type <- NA
merged_one$RCTD_first_weight <- NA
merged_one$RCTD_second_weight <- NA

merged_one$RCTD_spot_class[rownames(results_df)] <- as.character(results_df$spot_class)
merged_one$RCTD_first_type[rownames(results_df)] <- as.character(results_df$first_type)
merged_one$RCTD_second_type[rownames(results_df)] <- as.character(results_df$second_type)

merged_one$RCTD_first_weight[rownames(weights)] <- weights[, "first_type"]
merged_one$RCTD_second_weight[rownames(weights)] <- weights[, "second_type"]

# plot
SpatialDimPlot(
  merged_one,
  group.by = "RCTD_first_type",
  pt.size.factor = 3 # not sure if this is too large but I like how it looks
)

SpatialDimPlot(
  merged_one,
  group.by = "RCTD_spot_class",
  pt.size.factor = 3
)

# plot secondary type
SpatialDimPlot(
  merged_one,
  group.by = "RCTD_second_type",
  pt.size.factor = 3
)

# save
saveRDS(myRCTD, file = paste0("results/RCTD_", sample_id, ".rds"))
saveRDS(merged_one, file = paste0("results/merged_one_RCTD_", sample_id, ".rds"))
