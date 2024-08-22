# Load packages
library(Seurat)
library(hdf5r)
library(ggplot2)
library(ggsci)
library(dplyr)

# Set working directory
setwd("~/Seurat")

# Load data
SeuObj <- readRDS("~/Seurat/SeuObj.rds")

# Isolate Liver tissue
Idents(SeuObj) <- "Tissue"
liverObj <- subset(x=SeuObj, idents = "LIVER")

# Perform analysis without integration
# Scaling
all.genes <- rownames(liverObj)
liverObj <- ScaleData(liverObj, features = all.genes)

#PCA
liverObj <- RunPCA(liverObj)

# Determine dimensionality of the data
ElbowPlot(liverObj, ndims = 50)

# Harmony integration
liverObj <- IntegrateLayers(
 object = liverObj, method = HarmonyIntegration,
 orig.reduction = "pca", new.reduction = "harmony",
 verbose = FALSE
)

# Save obj in case R studio crashes
saveRDS(liverObj, file = "liverObj.rds")

# Re-join layers after integration
liverObj[["RNA"]] <- JoinLayers(liverObj[["RNA"]])

# Checking clusters of harmony integration
liverObj <- FindNeighbors(liverObj, reduction = "harmony", dims = 1:30)
liverObj <- FindClusters(liverObj, resolution = 0.5)
liverObj <- RunUMAP(liverObj, dims = 1:30, reduction = "harmony")

# Identify markers of each clusters
liverObj@misc$markers <- FindAllMarkers(liverObj)
table(liverObj@misc$markers$cluster)
top10_markers <- as.data.frame(liverObj@misc$markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC))
top10_markers
write.csv(top10_markers, file="/data/home/hdx044/Files/liver/cluster_markers_top10_liver.csv")

#### Add metadata ####
# Add stage in meta data
# Create lookup table with stage and sample
lookup <- data.frame(Patient_ID = c('10113-1','10113-2','10202','10203','10205','10291-2','10380','10634','10738','10742','9680','9932', '9961', '9991', '9999'),
                     Stage = c('F1_2','F0','F1_2','F1_2','F1_2','F1_2','F1_2','F1_2','F1_2','H','F1_2','F1_2', 'F0', 'F1_2', 'F0'))

# Extract meta data
meta <- liverObj@meta.data

# Create new column in meta data for Stage
meta$Stage <- NA

# Run for all samples in a for loop to populate the whole column
for (x in 1:length(lookup$Patient_ID)) {
 meta$Stage[grep(lookup$Patient_ID[x], meta$Patient_ID)] <- lookup$Stage[x]
}

# Add meta data back to seurat object
liverObj <- AddMetaData(liverObj, metadata = meta)

# Add cluster.stage in metadata
liverObj$Cluster.Stage <- paste(liverObj$seurat_clusters , liverObj$Stage, sep = "_")


#### Plot UMAPs ####
# Ploting UMAP for clusters
ggp = DimPlot(liverObj, reduction = "umap", raster=FALSE, label = TRUE)+ 
      theme_classic()+ labs(title = "LIVER clusters")
ggsave(
 plot = ggp,
 filename = ("~/plots/liver/clusters.tiff"),
 height = 8,
 width = 8,
 dpi = 300,
 units = 'in'
)

#Ploting UMAP for clusters by stage
ggp = DimPlot(liverObj, reduction = "umap", group.by = "Stage", raster=FALSE)+
      theme_classic()+ labs(title = "LIVER clusters by Stage")
ggsave(
 plot = ggp,
 filename = ("~/plots/liver/clusters_by_stage.tiff"),
 height = 8,
 width = 8,
 dpi = 300,
 units = 'in'
)

# Ploting UMAP for clusters by Patient_ID
ggp = DimPlot(liverObj, reduction = "umap", group.by = "Patient_ID", raster=FALSE)+
      theme_classic()+ labs(title = "LIVER clusters by Patient_ID")
ggsave(
 plot = ggp,
 filename = ("~/plots/liver/clusters_by_patients.tiff"),
 height = 8,
 width = 8,
 dpi = 300,
 units = 'in'
)

#### Generate count plots ####
# Create tissue summary data frame
summary <- liverObj@meta.data %>%
 group_by(seurat_clusters) %>%
 count(Stage) %>%
 mutate(proportion = n / sum(n))

# Plot summary stacked bar chart
ggp <- ggplot(summary, aes(x = seurat_clusters, y = n, fill = Stage)) +
 geom_bar(stat = "identity", position = "stack") +
 guides(fill=guide_legend(title="Stage")) +
 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(
 plot = ggp,
 filename = ("~/plots/liver/cluster_count.tiff"),
 height = 5,
 width = 10,
 dpi = 300,
 units = 'in'
)

# Plot summary proportion bar chart
ggp <- ggplot(summary, aes(x = seurat_clusters, y = proportion, fill = Stage)) +
 geom_bar(stat = "identity") +
 guides(fill=guide_legend(title="Stage")) +
 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
 scale_y_continuous(expand = c(0,0))
ggsave(
 plot = ggp,
 filename = ("~/plots/liver/cluster_count_stacked.tiff"),
 height = 5,
 width = 10,
 dpi = 300,
 units = 'in'
)

#### MHC expression ####

features <- c('HLA-A','HLA-B','HLA-C','HLA-DRA', 'HLA-DP', 'HLA-DM', 'HLA-DOA', 'HLA-DOB','HLA-DQ')

# Split by Stage
Idents(object = liverObj) <- "Stage"
ggp = DotPlot(liverObj, features = features)+ RotatedAxis()+
      scale_fill_npg()+ labs(title = "Antigen presentation in Liver", x ='MHC genes', y = 'Stage')
ggsave(
 plot = ggp,
 filename = ("~/plots/liver/Antigen_presentation_markers_stage.tiff"),
 height = 5,
 width = 5,
 dpi = 300,
 units = 'in'
)

# Split by Patient_ID
Idents(object = liverObj) <- "Patient_ID"

# Arrange the patients according to stage
liverObj@active.ident <- factor (liverObj@active.ident,
                                 levels = c('10205','9991','9932', '9680','10113-1','10202','10203','10291-2','10380','10634','10738','10113-2','9999','9961','10742'))

ggp = DotPlot(liverObj, features = features)+ RotatedAxis()+
      scale_fill_npg()+ labs(title = "Antigen presentation in Liver", x ='MHC genes', y = 'Patient_ID')
ggsave(
 plot = ggp,
 filename = ("~/plots/liver/Antigen_presentation_markers_patients.tiff"),
 height = 5,
 width = 5,
 dpi = 300,
 units = 'in'
)

# change cluster name
liverObj$cluster <- paste0("C", liverObj$seurat_clusters)
Idents(liverObj) <- "cluster"

# Arrange the cluster according to no
liverObj@active.ident <- factor (liverObj@active.ident,
                                 levels = c('C0','C1','C2','C3','C4','C5','C6','C7','C8','C9','C10','C11','C12','C13','C14', 'C15', 'C16', 'C17', 'C18', 'C19', 'C20', 'C21', 'C22'))


features <- c('CD3E','CD4','CD8B','ALB')
ggp = DotPlot(liverObj, features = features)+ RotatedAxis()+
 scale_fill_npg()+ labs(title = "T cells and Hepatocytes in Liver tissue", x ='Genes', y = 'Cluster_ID')
ggsave(
 plot = ggp,
 filename = ("~/plots/liver/T cells and Hepatocytes in Liver tissue.tiff"),
 height = 5,
 width = 8,
 dpi = 300,
 units = 'in'
)
#### Perform DE analysis #### 

# within the same cell type across conditions
liverObj$stage.seurat_clusters <- paste(liverObj$Stage, liverObj$seurat_clusters, sep = "_")
Idents(liverObj) <- "stage.seurat_clusters"
liverObj.markers <- FindMarkers(liverObj, ident.1 = "F1_2_1", ident.2 = "F0_1")
head(liverObj.markers, n = 50)

# Find differentially expressed features between Cluster 1 and all other cells, only
# search for positive markers
liverObj.markers <- FindMarkers(liverObj, ident.1 = "F1_2_1", ident.2 = NULL, only.pos = TRUE)
head(liverObj.markers, n = 50)


saveRDS(liverObj, file = "liverObj.rds")