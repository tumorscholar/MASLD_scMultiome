# This script is for the loading and prepreocessing of CITEseq datasets using Seurat

# Set working directory
setwd("~/Seurat")

# Load packages
library(Seurat)
library(hdf5r)
library(ggplot2)

# Locate files 
dataDir <- '/data/Blizard-AlazawiLab/rk/cellranger'
files <- list.files(path = dataDir,
                    pattern = 'filtered_feature_bc_matrix.h5',
                    recursive = T,
                    full.names = T)

# Extract sample names from file names
pattern1 <- '/data/Blizard-AlazawiLab/rk/cellranger/'
samples <- gsub(pattern1,
                '',
                files)

pattern2 <- '/outs/filtered_feature_bc_matrix.h5'
samples <- gsub(pattern2,
                '',
                samples)

# Load all data
rawDat <- lapply(files, function(x){
 Read10X_h5(x)
}) 

# Specify list element names 
names(rawDat) <- samples

# Make seurat objects in for loop
datObjs <- list() 
for (x in samples) {
 
 # Create seurat object
 datObjs[[x]] <- CreateSeuratObject(counts = rawDat[[x]][['Gene Expression']], project = x)
 
 # Add antibody to new assay object
 datObjs[[x]][['ADT']] <- CreateAssayObject(counts = rawDat[[x]][['Antibody Capture']])
 
 datObjs[[x]]@meta.data$sample <- x
 
 # Plot mitochondria and hemoglobin genes to remove dead and doublet cells
 datObjs[[x]][["percent.mt"]] <- PercentageFeatureSet(datObjs[[x]], pattern = "^MT-")
 features <- c('HBA1','HBA2','HBB')
 datObjs[[x]][["percent.hb"]] <- PercentageFeatureSet(datObjs[[x]], features = features)
 plt <- VlnPlot(datObjs[[x]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.hb"), ncol = 4)
 ggsave(
  plot = plt,
  filename = paste0(x, '_qcPlot.tiff'),
  height = 8,
  width = 8,
  dpi = 300,
  units = 'in'
 )
 plt <- FeatureScatter(datObjs[[x]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')
 ggsave(
  plot = plt,
  filename = paste0(x, '_lmPlot.tiff'),
  height = 8,
  width = 8,
  dpi = 300,
  units = 'in'
 )
 # Filtering
 datObjs[[x]] <- subset(datObjs[[x]], subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)
 plt <-VlnPlot(datObjs[[x]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.hb"), ncol = 4)
 ggsave(
  plot = plt,
  filename = paste0(x, '_afterQCplot.tiff'),
  height = 8,
  width = 8,
  dpi = 300,
  units = 'in'
 )
}

saveRDS(datObjs, 'SeuObj.rds')

