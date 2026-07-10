# R version - 4.5.3
suppressPackageStartupMessages({
  library(Seurat) # 5.5.0
  library(dplyr) # 1.2.1
  library(patchwork) # 1.3.2
  library(ggplot2) # 4.0.3
  library(clustree) # 0.5.1
  library(here) # 1.0.2
})

project_root <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
if (getwd() != project_root) {
  setwd(project_root)
}

source(here("code","utils.R"))
set.seed(123)

# Read Seurat object with spatial data
merged <- readRDS(here("data","merged_spatial_data.RDS"))

# Add QC metrics
merged$log10_nCount_Spatial <- log(merged$nCount_Spatial, 10)
merged$log10_nFeature_Spatial <- log(merged$nFeature_Spatial, 10)
merged$percent.mt <- PercentageFeatureSet(merged, pattern = "^MT-")

s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

# Calculate cell cycle genes scores
## Make temporary joined object
merged <- JoinLayers(merged, assay = "Spatial")

# Log-normalize the joined RNA assay for CellCycleScoring
merged <- NormalizeData(merged, 
                        assay = "Spatial", 
                        verbose = FALSE)

# Calculate cell cycle scores
merged <- CellCycleScoring(object = merged,  
                          s.features = s.genes, 
                           g2m.features = g2m.genes, 
                           slot = "counts") # why counts and not data?

# Re-split layers for the multisample SCTransform
merged[["Spatial"]] <- split(merged[["Spatial"]], f = merged$orig.ident)
gc()

head(merged@meta.data)

# List of all sample images
names(merged@images)

# Spatial visualization
## Number of UMIs + Spatial distribution
plot1 <- VlnPlot(merged, 
                 features = "log10_nCount_Spatial", 
                 layer = "counts",
                 pt.size = 0) + 
  NoLegend()

selected_sample <- "A009N"
plot2 <- SpatialFeaturePlot(merged, 
                            features = "log10_nCount_Spatial", 
                            images = selected_sample) + 
  labs(title = selected_sample) + 
  theme(legend.position = "right", 
        plot.title = element_text(hjust = 0.5, face="bold"))
wrap_plots(plot1, plot2)

## Number of genes + MT content plots
VlnPlot(merged, 
        features = c("nFeature_Spatial","percent.mt"), 
        layer = "counts",
        pt.size = 0,) + 
  NoLegend()

## Number of genes + MT content plots
VlnPlot(merged, 
        features = c("S.Score","G2M.Score"), 
        layer = "counts",
        pt.size = 0,) + 
  NoLegend()

# Before QC: 36601 genes and 27224 spots
dim(merged)

merged <- subset(merged, log10_nCount_Spatial >= 3 & percent.mt <= 10)
# After QC: 36601 genes and 25096 spots
dim(merged)

# Apply SCTransform
## Make sure that sample are stored as layers
Layers(merged)
## Apply SCTransform on each sample layer
## If getting memory limit error - increase mem.maxVSize()
## usethis::edit_r_environ() to open .Renviron
## R_MAX_VSIZE=50Gb

merged <- SCTransform(merged, 
                      layer = "counts",
                      vst.flavor = "v2", 
                      assay = "Spatial",
                      vars.to.regress = c("nFeature_Spatial", 
                                          "percent.mt"),
                      verbose = FALSE)
gc()
merged <- RunPCA(merged, verbose = F) # Elbow plot?

# Unintegrated data processing
merged <- FindNeighbors(merged, dims = 1:30, reduction = "pca", verbose = F) %>%
  RunUMAP(dims = 1:30, 
          reduction = "pca", 
          reduction.name = "umap.unintegrated", 
          verbose = FALSE)

# Run each cluster resolution independently to enable custom names
cluster_resolutions <- seq(0.1, 1.0, 0.1)
for (res in cluster_resolutions) {
  merged <- FindClusters(merged,
                         resolution = res, 
                         cluster.name = paste0("unintegrated_clusters_",res), 
                         verbose = FALSE)
}
# Cluster tree plot to identify optimal
clustree(merged, prefix = "unintegrated_clusters_")

## Let's pick 0.2 cluster resolution as higher stability
merged$unintegrated_clusters <- merged$unintegrated_clusters_0.2
p1 <- DimPlot(merged, reduction = "umap.unintegrated", group.by = c("orig.ident"), label.size = 2)
p2 <- DimPlot(merged, reduction = "umap.unintegrated", group.by = c("unintegrated_clusters"), 
              label = T, label.box = T, label.size = 4) + NoLegend()
wrap_plots(p1, p2, ncol = 2, byrow = F)

# Integrated (Harmony) data processing
merged <- IntegrateLayers(object = merged, 
                          method = HarmonyIntegration, 
                          normalization.method = "SCT",
                          orig.reduction = "pca",
                          new.reduction = "integrated.harmony", 
                          verbose = FALSE)
merged <- FindNeighbors(merged, 
                        reduction = "integrated.harmony", 
                        dims = 1:30, 
                        verbose = FALSE) %>%
  RunUMAP(reduction = "integrated.harmony", 
          dims = 1:30, 
          reduction.name = "umap.harmony", 
          verbose = FALSE)

# Comparing sample mixing between unintegrated and integrated with Harmony versions
p1 <- DimPlot(merged, reduction = "umap.unintegrated", group.by = c("orig.ident"), label.size = 2)
p2 <- DimPlot(merged, reduction = "umap.harmony", group.by = c("orig.ident"), label.size = 2)
wrap_plots(p1, p2, ncol = 2, byrow = F)

# Find clusters on Harmony embedding
for (res in cluster_resolutions) {
  merged <- FindClusters(merged,
                         resolution = res, 
                         cluster.name = paste0("harmony_clusters_",res), 
                         verbose = FALSE)
}
# Cluster tree for harmony integrated clusters
clustree(merged, prefix = "harmony_clusters_")

# We will keep resolution 0.2 for harmony clusters
merged$integrated_clusters <- merged$harmony_clusters_0.2
p1 <- DimPlot(merged, reduction = "umap.harmony", group.by = c("orig.ident"), label.size = 2)
p2 <- DimPlot(merged, reduction = "umap.harmony", group.by = c("integrated_clusters"), 
        label = T, label.box = T, label.size = 4) + NoLegend()
wrap_plots(p1, p2, ncol = 2, byrow = F)

saveRDS(merged, file = here("data", "merged_spatial_processed_harmony.RDS"))