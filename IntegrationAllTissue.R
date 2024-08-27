library(Seurat)
library(ggplot2)
library(patchwork)
library(harmony)

setwd("/data/Blizard-AlazawiLab/rk")

# Load data
SeuObj <- readRDS("~/Seurat/SeuObj.rds")

#Perform analysis without integration
# Scaling
all.genes <- rownames(SeuObj)
SeuObj <- ScaleData(SeuObj, features = all.genes)

#PCA
SeuObj <- RunPCA(SeuObj)

# Determine dimensionality of the data
ElbowPlot(SeuObj, ndims = 50)

#Clustering and UMAP
SeuObj <- FindNeighbors(SeuObj, dims = 1:30, reduction = "pca")
SeuObj <- FindClusters(SeuObj, resolution = 0.5, cluster.name = "unintegrated_clusters")
SeuObj <- RunUMAP(SeuObj, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(SeuObj, reduction = "umap.unintegrated", group.by = 'Patient_ID', raster = F)
DimPlot(SeuObj, reduction = "umap.unintegrated", raster = F)

# Harmony integration
SeuObj <- IntegrateLayers(
 object = SeuObj, method = HarmonyIntegration,
 orig.reduction = "pca", new.reduction = "harmony",
 verbose = FALSE
)

# Re-join layers after integration
SeuObj[["RNA"]] <- JoinLayers(SeuObj[["RNA"]])

# Checking clusters of harmony integration
SeuObj <- FindNeighbors(SeuObj, reduction = "harmony", dims = 1:30)
SeuObj <- FindClusters(SeuObj, resolution = 0.5)
SeuObj <- RunUMAP(SeuObj, dims = 1:30, reduction = "harmony")

# Identify markers of each clusters
SeuObj@misc$markers <- FindAllMarkers(SeuObj)
saveRDS(SeuObj, file = "SeuObjI.rds")

# Ploting UMAP 
DimPlot(SeuObj, reduction = "umap", raster=FALSE, label = TRUE)
DimPlot(SeuObj, reduction = "umap", split.by = "Tissue", raster=FALSE)
DimPlot(SeuObj, reduction = "umap", group.by = "Stage", raster=FALSE)
DimPlot(SeuObj, reduction = "umap", group.by = "Batch", raster=FALSE)

# Visualisation
features <- c("CD3E", "CD19", "ALB")
FeaturePlot(SeuObj, features = features, reduction = "umap", raster = FALSE)+ RotatedAxis()
DotPlot(SeuObj, features = features)+ RotatedAxis()

DefaultAssay(SeuObj) <- 'ADT'
features <- c("Hu.CD3-UCHT1")
DotPlot(SeuObj, features = features)+ RotatedAxis()





