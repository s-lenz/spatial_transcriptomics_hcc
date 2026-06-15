suppressPackageStartupMessages({
  library(Seurat)
})
set.seed(123)

# Functions to parse file name and extract sample information
# Extract sample name
extract_sample_names <- function(filepath) {
  return(strsplit(filepath, split="_")[[1]][2])
}
# Extract GSM ID
extract_gsm_id <- function(filepath) {
  return(strsplit(filepath, split="_")[[1]][1])
}

# Function to load 10X Visium data into Seurat object
Load10X_Visium <- function (
    data.dir,
    image.dir,
    assay = "Spatial",
    sample_id = "slice1",
    filter.matrix = TRUE,
    image.name = "tissue_lowres_image.png"
) {
  # read counts matrix
  counts <- Read10X(data.dir)
  # build a Seurat object
  seu <- CreateSeuratObject(counts, 
                            project = sample_id, 
                            assay = assay)
  # read the corresponding image and coordinate mapping
  img <- Read10X_Image(image.dir,
                       assay = assay, 
                       slice = sample_id,
                       image.name = image.name,
                       filter.matrix = filter.matrix)
  # align the image identifier with the object cells
  img <- img[Cells(seu)]
  # add the image to the corresponding Seurat instance
  seu[[sample_id]] <- img
  # return Seurat object
  return(seu)
}