suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
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

plot_Proportions_RCTD <- function(tissue_positions, proportions, title=NULL){
  proportions <- sample_proportions
  if (nrow(tissue_positions) != nrow(proportions)) {
    stop("Number of spatial spots doesn't match number of estimated proportions")
  }
  plotData <- cbind(sample_proportions, tissue_positions)
  allPlots <- list()
  for (ct in colnames(sample_proportions)) {
    plt <- ggplot(plotData, aes(x, y, col=.data[[ct]])) + 
      coord_equal() + theme_void() + 
      geom_point(size=0.3)
    allPlots[[ct]] <- plt
  }
  merged_plot <- allPlots |> wrap_plots(nrow=3) & theme(
    legend.key.width=unit(0.5, "lines"),
    legend.key.height=unit(1, "lines")) &
    scale_color_gradientn(colors=pals::jet())
  
  if (!is.null(title)) {
    merged_plot <- merged_plot + 
      plot_annotation(title,
                      theme = theme(plot.title = element_text(hjust = 0.5)))
  }
  
  return(merged_plot)
}