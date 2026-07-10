suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(stringr)
  library(pheatmap)
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

plot_Heatmap_RCTD <- function(proportions, title) {
  p <- pheatmap(as.matrix(proportions), 
                show_rownames=FALSE, show_colnames=TRUE, main=title)
                #cellwidth=12, treeheight_row=5, treeheight_col=5)
  return(p)
}

get_celltypes_colorpanel <- function(cell_types) {
  colors <- glasbey.colors(12)[c(2:4,6:12)]
  names(colors) <- cell_types
  return(colors)
}

get_celltypes_legend <- function(cell_types, color_pallete) {
  dummy_data <- data.frame(
    Cell_Type = cell_types,
    Values = c(1:length(cell_types))
  )
  
  custom_legend <- ggplot(dummy_data, aes(x = Values, y = Values, color = Cell_Type)) +
    scale_color_manual(values = color_pallete) +
    geom_point(size=4) + theme(legend.position = "bottom")
  pure_legend <- cowplot::get_legend(custom_legend)
  return(pure_legend)
}

plotSpatialDim_by_celltypes <- function(spatial_data, 
                                        sample_column, 
                                        cell_types,
                                        cell_type_column = "main_cell_type",
                                        n_rows = 3,
                                        pt.size = 3) {
  
  all_samples <- unique(spatial_data@meta.data[[sample_column]])
  
  celltype_colors <- get_celltypes_colorpanel(cell_types)
  
  allPlots <- list()
  for (sample_name in all_samples) {
    plt <- SpatialDimPlot(spatial_data, pt.size.factor = pt.size,
                          group.by = cell_type_column, 
                          images = sample_name, 
                          alpha = 1.5,
                          image.alpha = 0,
                          cols = celltype_colors) +
      NoLegend() + 
      labs(title = sample_name) + 
      theme(title=element_text(hjust = 0.5))
    allPlots[[sample_name]] <- plt
  }
  
  main_grid <- wrap_plots(allPlots, nrow=n_rows)
  pure_legend <- get_celltypes_legend(cell_types, celltype_colors)
  
  final_layout <- wrap_plots(
    main_grid,
    pure_legend,
    ncol = 1,
    heights = c(0.95, 0.05)
  )
  return(final_layout)
  
}

plotGEX_by_celltypes <- function(data, 
                                 gene_name,
                                 sample_column,
                                 cell_types, 
                                 n_rows = 3) {
  allPlots <- list()
  p <- VlnPlot(data, 
          features = c(gene_name), 
          group.by = sample_column,
          assay = "SCT") + xlab("") + ylab("") +
    labs(title = paste(gene_name," - Total")) + 
    NoLegend()
  allPlots[["total"]] <- p
  
  for (cell_type in cell_types) {
    assay_name <- paste0("Pseudo_",
                         str_replace_all(str_replace_all(cell_type, " ", "_"),
                                         "-","."), 
                         "_SCT")
    p <- VlnPlot(data, 
                 features = c(gene_name), 
                 group.by = sample_column,
                 assay = assay_name, 
                 layer="data") + xlab("") + ylab("") +
      labs(title = paste0(gene_name, " - ", cell_type)) + 
      NoLegend()
    allPlots[[assay_name]] <- p
  }
  return(wrap_plots(allPlots, nrow=n_rows))
  
}

plotClusterHeatmap_by_celltype <- function(median_proportions, 
                                           cell_types,
                                           n_samples){
  toPlot <- median_proportions %>%
    pivot_longer(
      cols = all_of(cell_types), 
      names_to = "cell_type", 
      values_to = "median_proportion"
    ) %>%
    complete(integrated_clusters, 
             orig.ident, cell_type, 
             fill = list(median_proportion = 0))
  
  toPlot <- toPlot %>%
    arrange(cell_type, orig.ident) %>%
    mutate(sample_ct_col = paste(cell_type, orig.ident, sep = "__"))
  
  wide_df <- toPlot %>%
    select(integrated_clusters, sample_ct_col, median_proportion) %>%
    pivot_wider(names_from = sample_ct_col, values_from = median_proportion) %>%
    arrange(as.numeric(as.character(integrated_clusters))) %>%
    column_to_rownames("integrated_clusters")
  
  plot_matrix <- as.matrix(wide_df)
  
  col_anno_df <- data.frame(sample_ct_col = colnames(plot_matrix)) %>%
    separate(sample_ct_col, into = c("cell_type", "sample"), sep = "__", remove = FALSE) %>%
    select(cell_type)
  rownames(col_anno_df) <- colnames(plot_matrix)
  col_anno_df$cell_type <- factor(
    col_anno_df$cell_type, 
    levels = unique(col_anno_df$cell_type) 
  )
  
  anno_colors <- list(cell_type = get_celltypes_colorpanel(cell_types))
  
  n_cell_types <- length(unique(col_anno_df$cell_type))
  col_gaps <- seq(n_samples, (n_samples * n_cell_types) - 1, by = n_samples)
  
  pheatmap(
    mat = plot_matrix,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    annotation_col = col_anno_df,
    show_colnames = FALSE, 
    gaps_col = col_gaps,   
    color = colorRampPalette(c("white", "indianred1" ,"red2", "darkred"))(100),
    annotation_colors = anno_colors,
    main = "Median Cell Type Proportions per Cluster"
  )
}

