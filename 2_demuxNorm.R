library(Seurat)
library(ggplot2)
library(patchwork)

# Set working directory
setwd("~/Seurat")

# Load data list
mylist <- readRDS("~/Seurat/before_soupx/datObj.rds")

# Split list into LIVER and SVP objects
names(mylist)

# Subset LIVER
grep('LIVER', names(mylist))
LIVERlist <- mylist[grep('LIVER', names(mylist))]

# Subset SVP
grep('SAT-VAT-PBMC', names(mylist))
SVPlist <- mylist[grep('SAT-VAT-PBMC', names(mylist))]

#### Add metadata, normalise, FindVariableFeatures and scale LIVER data ####
# Name LIVER list
names(LIVERlist)
LIVER <- names(LIVERlist)[1]
LIVERlist[[LIVER]]

for (x in names(LIVERlist)) {
 
 # Create metadata columns for Tissue type and Patient_ID 
 LIVERlist[[x]]$Tissue <- 'LIVER'
 id <- gsub('-LIVER', '', x)
 id <- gsub('GC-WL-', '', id)
 LIVERlist[[x]]$Patient_ID <- id
 
 # Normalise ADT data
 DefaultAssay(LIVERlist[[x]]) <- 'ADT'
 LIVERlist[[x]] <- NormalizeData(LIVERlist[[x]], normalization.method = "CLR", margin = 2)
 
 # Normalise and scale RNA data
 DefaultAssay(LIVERlist[[x]]) <- "RNA"
 LIVERlist[[x]] <- NormalizeData(LIVERlist[[x]])
 LIVERlist[[x]] <- FindVariableFeatures(LIVERlist[[x]])
 LIVERlist[[x]] <- ScaleData(LIVERlist[[x]])
}

#### Add metadata, normalise, FindVariableFeatures, scale and demultiplex SVP data ####
# Name SVP list
names(SVPlist)
SVP <- names(SVPlist)[1]
SVPlist[[SVP]]

# Check ADT list
rownames(SVPlist$`GC-WL-10113-1-SAT-VAT-PBMC`@assays$ADT$counts)
rownames(SVPlist$`GC-WL-10113-2-SAT-VAT-PBMC`@assays$ADT$counts)

for (x in names(SVPlist)) {
 # Split ADT and HTO data into separate assays
 SVPlist[[x]][['HTO']] <- CreateAssayObject(counts = SVPlist[[x]]@assays$ADT$counts[c('C0251', 'C0252', 'C0253'), ])
 SVPlist[[x]][['ADTonly']] <- CreateAssayObject(counts = SVPlist[[x]]@assays$ADT$counts[c(1:137), ])

 # Normalise ADT and HTO data
 SVPlist[[x]] <- NormalizeData(SVPlist[[x]], assay = "HTO", normalization.method = "CLR", margin = 2)
 SVPlist[[x]] <- NormalizeData(SVPlist[[x]], assay = "ADTonly", normalization.method = "CLR", margin = 2)

 # Demultiplex cells based on HTO enrichment
 SVPlist[[x]] <- HTODemux(SVPlist[[x]], assay = "HTO", positive.quantile = 0.95)
 
 # Keep only singlets
 SVPlist[[x]] <- subset(SVPlist[[x]], subset = HTO_classification.global == 'Singlet')
 
 # Demux results are saved in HTO classification column under meta data
 # Create tissue vector using HTO classification
 tissue_vector <- c(rep(NA, length(SVPlist[[x]]$HTO_classification)))
 tissue_vector[SVPlist[[x]]$HTO_classification == 'C0251'] <- 'SAT'
 tissue_vector[SVPlist[[x]]$HTO_classification == 'C0252'] <- 'VAT'
 tissue_vector[SVPlist[[x]]$HTO_classification == 'C0253'] <- 'PBMC'

 # Create metadata columns for Tissue type and Patient_ID and add tissue vector in tissue type
 SVPlist[[x]]$Tissue <- tissue_vector
 id <- gsub('-SAT-VAT-PBMC', '', x)
 id <- gsub('GC-WL-', '', id)
 SVPlist[[x]]$Patient_ID <- id
 
 # Normalise ADT data
 DefaultAssay(SVPlist[[x]]) <- 'ADT'
 SVPlist[[x]] <- NormalizeData(SVPlist[[x]], normalization.method = "CLR", margin = 2)
 
 # Normalise RNA data
 DefaultAssay(SVPlist[[x]]) <- "RNA"
 SVPlist[[x]] <- NormalizeData(SVPlist[[x]])
 SVPlist[[x]] <- FindVariableFeatures(SVPlist[[x]])
 SVPlist[[x]] <- ScaleData(SVPlist[[x]])
}

# Merge the liver and SAT, VAT and PBMC list
all.tissue.list <- c(LIVERlist, SVPlist)
SeuObj <- merge(x = all.tissue.list[[1]], 
                y = all.tissue.list[2:30])

saveRDS(SeuObj, 'SeuObj.rds')
