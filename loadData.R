# This script is for the loading and prepreocessing of CITEseq datasets using Seurat

# Set working directory
setwd('...')

# Load packages
library(Seurat)
library(hdf5r)

# Locate files 
dataDir <- '/data/Blizard-AlazawiLab/rk/cellranger_counts_rk'
files <- list.files(path = dataDir,
                    pattern = 'filtered_feature_bc_matrix.h5',
                    recursive = T,
                    full.names = T)
files <- files[-c(10,16,18,19,28,30,32)]

# Extract sample names from file names
pattern1 <- '/data/Blizard-AlazawiLab/rk/cellranger_counts_rk/'
samples <- gsub(pattern1,
                '',
                files)

pattern2 <- '/outs/filtered_feature_bc_matrix.h5'
samples <- gsub(pattern2,
                '',
                samples)

pattern3 <- '/count/sample_filtered_feature_bc_matrix.h5'
samples <- gsub(pattern3,
                '',
                samples)

pattern4 <- 'new_samples/Analysis/Seurat_Demux/h5/'
samples <- gsub(pattern4,
                '',
                samples)

pattern5 <- '_sample_filtered_feature_bc_matrix.h5'
samples <- gsub(pattern5,
                '',
                samples)

pattern6 <- 'new_samples/Analysis/CellRanger_Total/'
samples <- gsub(pattern6,
                '',
                samples)

pattern7 <- 'GC-WL-10634-LIVER/outs/per_sample_outs/'
samples <- gsub(pattern7,
                '',
                samples)

pattern8 <- 'GC-WL-10738-LIVER/outs/per_sample_outs/'
samples <- gsub(pattern8,
                '',
                samples)

pattern9 <- 'GC-WL-10742-LIVER/outs/per_sample_outs/'
samples <- gsub(pattern9,
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
for (x in samples[-c(10,16,18,19,28,30,32)]) {
 
 # Create seurat object
 datObjs[[x]] <- CreateSeuratObject(counts = rawDat[[x]][['Gene Expression']],
                                    project = x)
 
 # Add antibody to new assay object
 datObjs[[x]][['ADT']] <- CreateAssayObject(counts = rawDat[[x]][['Antibody Capture']])
 
 datObjs[[x]]@meta.data$sample <- x
 
}


