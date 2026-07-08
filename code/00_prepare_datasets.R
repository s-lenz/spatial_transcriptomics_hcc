# R version - 4.5.3
suppressPackageStartupMessages({
  library(GEOquery) # 2.78.0
  library(utils) # 4.5.3
  library(R.utils) # 2.13.0
  library(Seurat) # 5.5.0
  library(dplyr) # 1.2.1
  library(here) # 1.0.2
})

project_root <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
if (getwd() != project_root) {
  setwd(project_root)
}

source(here("code","utils.R"))
set.seed(123)

## Spatial dataset: GSE281759

gse <- "GSE281759"
save_dir <- here("data", gse)
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

# Get metadata that contains phenotypical data + other info about Visium Samples
geo_obj <- getGEO(gse)[[1]]
metadata <- pData(geo_obj)
# Save file with metadata so we can load it in later if we need
write.csv(metadata, file.path(save_dir, "GSE281759_metadata.csv"))

# Get Expression data files (raw counts) and tissue images (CytAssist files)
supp_files <- getGEOSuppFiles(gse,baseDir = "data")

# Unpack main data
untar(file.path(save_dir,"GSE281759_RAW.tar"), 
      exdir = file.path(save_dir,"GSE281759_RAW"))

# Unpack only spatial image files (keep GEX as gz files)
for (filename in list.files(file.path(save_dir,"GSE281759_RAW"),
                            full.names = TRUE)) {
  if (!grepl(pattern = "(matrix.mtx|barcodes.tsv|features.tsv)",
             x = filename)) {
    gunzip(filename)
  }
}

# Reorganize spatial data by sample (Seurat expected inputs)
old_spatial_folder <- file.path(save_dir, "GSE281759_RAW")
all_files <- list.files(old_spatial_folder, full.names = T)
sample_list <- unique(sapply(list.files(old_spatial_folder), 
                             FUN = extract_sample_names))
gsm_list <- unique(sapply(list.files(old_spatial_folder), 
                          FUN = extract_gsm_id))

for (idx in seq_along(sample_list)) {
  sample_name <- sample_list[idx]
  gsm_name <- gsm_list[idx]
  new_spatial_folder <- file.path(old_spatial_folder, sample_name)
  # Create sample target folder
  dir.create(file.path(new_spatial_folder,"spatial"), 
             recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(new_spatial_folder,"gex"), 
             recursive = TRUE, showWarnings = FALSE)
  
  # Make a list of samples files to move
  files_to_move <- all_files[grepl(sample_name,all_files)]
  new_names <- gsub(paste0(gsm_name,"_",sample_name,"_"),"",files_to_move)
  
  # Move GEX files
  selected_files <- grepl("(matrix.mtx|barcodes.tsv|features.tsv)",
                          x = files_to_move)
  files_to_move_gex <- files_to_move[selected_files]
  target_paths_gex <- file.path(new_spatial_folder, 
                                "gex",
                                basename(new_names[selected_files]))
  file.rename(from = files_to_move_gex, to = target_paths_gex)
  
  # Move Spatial files
  files_to_move_st <- files_to_move[!selected_files]
  target_paths_st <- file.path(new_spatial_folder, 
                               "spatial", 
                               basename(new_names[!selected_files]))
  file.rename(from = files_to_move_st, to = target_paths_st)
}

# Read 10X Visium objects to Seurat
objects.list <- list()
input_filepath <- file.path(save_dir, "GSE281759_RAW")
seurat_savedir <- file.path(save_dir, "seurat_objects")
dir.create(seurat_savedir, recursive = TRUE, showWarnings = FALSE)
for (sample_name in sample_list) {
  st_object <- Load10X_Visium(file.path(input_filepath, sample_name, "gex"),
                              file.path(input_filepath, sample_name, "spatial"),
                              sample_id = sample_name)
  sample_metadata <- metadata[metadata$description == sample_name,]
  st_object$Sex <- as.character(sample_metadata["Sex:ch1"])
  st_object$Tumour_Stage <- as.character(sample_metadata["tumour stage:ch1"])
  saveRDS(st_object, file.path(seurat_savedir,paste0(sample_name,".RDS")))
  objects.list[[sample_name]] <- st_object
}

# Merge and save Seurat object
merged <- merge(
  objects.list[[1]],
  y = objects.list[-1]
)
saveRDS(merged, here("data","merged_spatial_data.RDS"))
