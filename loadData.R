# This script is for the loading and prepreocessing of scRNAseq datasets using Seurat

# Set working directory
setwd('...')

# Load packages
library(Seurat)

# Locate files 
dataDir <- '...'
files <- list.files(path = dataDir,
                    pattern = '*.h5',
                    recursive = T,
                    full.names = T)

# Extract sample names from file names
pattern <- '...'
samples <- gsub(pattern,
                '',
                files)

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
