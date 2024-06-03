# This script is for the loading and prepreocessing of CITEseq datasets using Seurat

# Set working directory
setwd('...')

# Load packages
library(Seurat)
library(hdf5r)

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
 datObjs[[x]] <- CreateSeuratObject(counts = rawDat[[x]][['Gene Expression']],
                                    project = x)
 
 # Add antibody to new assay object
 datObjs[[x]][['ADT']] <- CreateAssayObject(counts = rawDat[[x]][['Antibody Capture']])
 
 datObjs[[x]]@meta.data$sample <- x
 
}

