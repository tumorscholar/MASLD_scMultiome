library(Seurat)
library(ggplot2)
library(CellChat)
library(NMF)
library(ggalluvial)
library(ggsci)
library(ComplexHeatmap)
library(patchwork)
options(stringsAsFactors = FALSE)

# Set working directory
setwd("~/cellchat/liver")

#### Data input & processing and initialization of CellChat object ####
# Prepare required input data for CellChat analysis
liverObj <- readRDS("/data/home/hdx044/Seurat/liverObj.rds")

# change cluster name add samples column for cellchat
liverObj$samples <- paste(liverObj$Patient_ID)
liverObj$cluster <- paste0("C", liverObj$seurat_clusters)
Idents(liverObj) <- "cluster"

# Arrange the clusters in increasing order
Idents(liverObj) <- "cluster"
liverObj@active.ident <- factor (liverObj@active.ident,
                                 levels = c('C0','C1','C2','C3','C4','C5','C6','C7','C8','C9','C10','C11','C12','C13','C14', 'C15', 'C16', 'C17', 'C18', 'C19', 'C20', 'C21', 'C22'))

# Keep common clusters of all Stages for comparison
liverObjcellchat <- subset(x = liverObj, idents = c("C19", "C21"), invert = TRUE)

data.input <- liverObjcellchat[["RNA"]]$data # normalized data matrix
labels <- Idents(liverObjcellchat)
# create a dataframe of the cell labels
meta <- data.frame(labels = labels, row.names = names(labels))

# Split obj stage wise for comparison
objList <- SplitObject(liverObjcellchat, split.by = 'Stage')

#### Create a CellChat object for Control ####
cellchat<- objList$H
cellchat <- createCellChat(object = cellchat, group.by = "ident", assay = "RNA")

# Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.human 
showDatabaseCategory(CellChatDB)

# Use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB

# Set the used database in the object
cellchat@DB <- CellChatDB.use

# Preprocessing the expression data for cell-cell communication analysis
# Subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

#Project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
cellchat <- projectData(cellchat, PPI.human)

# Compute the communication probability and infer cellular communication network
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)

# Extract the inferred cellular communication network as a data frame
df.net <- subsetCommunication(cellchat)

# Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

# Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)

# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
netAnalysis_signalingRole_scatter(cellchat)

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")

# Identify global communication patterns
selectK(cellchat, pattern = "outgoing")
nPatterns = 4
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)

# River plot
netAnalysis_river(cellchat, pattern = "outgoing")

# Dot plot
netAnalysis_dot(cellchat, pattern = "outgoing")

# Identify and visualize incoming communication pattern of target cells
selectK(cellchat, pattern = "incoming")
nPatterns = 2
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)

# river plot
netAnalysis_river(cellchat, pattern = "incoming")

# dot plot
netAnalysis_dot(cellchat, pattern = "incoming")

# Identify signaling groups based on their functional similarity
cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")

# Visualization in 2D-space
netVisual_embedding(cellchat, type = "functional", label.size = 3.5)

# Identify signaling groups based on structure similarity
cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")

# Visualization in 2D-space
netVisual_embedding(cellchat, type = "structural", label.size = 3.5)

saveRDS(cellchat, file = "cellchatH.rds")

#### Create a CellChat object for F0 ####
cellchat<- objList$F0
cellchat <- createCellChat(object = cellchat, group.by = "ident", assay = "RNA")

# Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.human 
showDatabaseCategory(CellChatDB)

# Use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB

# Set the used database in the object
cellchat@DB <- CellChatDB.use

# Preprocessing the expression data for cell-cell communication analysis
# Subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

#Project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
cellchat <- projectData(cellchat, PPI.human)

# Compute the communication probability and infer cellular communication network
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)

# Extract the inferred cellular communication network as a data frame
df.net <- subsetCommunication(cellchat)

# Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

# Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)

# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
netAnalysis_signalingRole_scatter(cellchat,vertex.weight = 10)

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")

# Identify global communication patterns
selectK(cellchat, pattern = "outgoing")
nPatterns = 3
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)

# River plot
netAnalysis_river(cellchat, pattern = "outgoing")

# Dot plot
netAnalysis_dot(cellchat, pattern = "outgoing")

# Identify and visualize incoming communication pattern of target cells
selectK(cellchat, pattern = "incoming")
nPatterns = 3
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)

# River plot
netAnalysis_river(cellchat, pattern = "incoming")

# Dot plot
netAnalysis_dot(cellchat, pattern = "incoming")

# Identify signaling groups based on their functional similarity
cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")

# Visualization in 2D-space
netVisual_embedding(cellchat, type = "functional", label.size = 3.5)

# Identify signaling groups based on structure similarity
cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")

# Visualization in 2D-space
netVisual_embedding(cellchat, type = "structural", label.size = 3.5)

saveRDS(cellchat, file = "cellchatF0.rds")

#### Create a CellChat object for F1_2 ####
cellchat<- objList$F1_2
cellchat <- createCellChat(object = cellchat, group.by = "ident", assay = "RNA")

# Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.human 
showDatabaseCategory(CellChatDB)

# Use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB

# Set the used database in the object
cellchat@DB <- CellChatDB.use

# Preprocessing the expression data for cell-cell communication analysis
# Subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

#Project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
cellchat <- projectData(cellchat, PPI.human)

# Compute the communication probability and infer cellular communication network
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)

# Extract the inferred cellular communication network as a data frame
df.net <- subsetCommunication(cellchat)

# Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

# Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)

# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
netAnalysis_signalingRole_scatter(cellchat)

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")

# Identify global communication patterns
selectK(cellchat, pattern = "outgoing")
nPatterns = 4
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)

# River plot
netAnalysis_river(cellchat, pattern = "outgoing")

# Dot plot
netAnalysis_dot(cellchat, pattern = "outgoing")

# Identify and visualize incoming communication pattern of target cells
selectK(cellchat, pattern = "incoming")
nPatterns = 5
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)

# River plot
netAnalysis_river(cellchat, pattern = "incoming")

# Dot plot
netAnalysis_dot(cellchat, pattern = "incoming")

# Identify signaling groups based on their functional similarity
cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")

# Visualization in 2D-space
netVisual_embedding(cellchat, type = "functional", label.size = 3.5)

# Identify signaling groups based on structure similarity
cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")

# Visualization in 2D-space
netVisual_embedding(cellchat, type = "structural", label.size = 3.5)

saveRDS(cellchat, file = "cellchatF1_2.rds")

#### Load CellChat object of H, F0 and F1_2 dataset ####
cellchatH <- readRDS("cellchatH.rds")
cellchatF0 <- readRDS("cellchatF0.rds")
cellchatF1_2 <- readRDS("cellchatF1_2.rds")

# Merge H, F0 and F1_2 datasets for comparison
object.list <- list(H = cellchatH, F0 = cellchatF0, F1_2 = cellchatF1_2)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

# Identify cell populations with significant changes in sending or receiving signals
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
 gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax, label.size = 6, font.size = 14, font.size.title = 14 ) +  scale_y_continuous(limits = c(0,25)) + scale_x_continuous(limits = c(0,25))+
  
  theme(
   axis.text.x = element_text(face = "bold"),  # Bold x-axis text
   axis.text.y = element_text(face = "bold"),  # Bold y-axis text
   axis.title.x = element_text(face = "bold"), # Bold x-axis title
   axis.title.y = element_text(face = "bold")  # Bold y-axis title
  )
}

patchwork::wrap_plots(plots = gg)
ggsave("Differntial_interaction_strength_2D_LIVER.png", dpi = 300, width = 24, height = 8)

# Compare outgoing (or incoming) signaling patterns associated with each cell population
library(ComplexHeatmap)
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i+1]]@netP$pathways, object.list[[i+2]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 10, height = 40, color.heatmap = "Blues")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 10, height = 40, color.heatmap = "Blues")
ht3 = netAnalysis_signalingRole_heatmap(object.list[[i+2]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+2], width = 10, height = 40, color.heatmap = "Blues")

combined_heatmap <- draw(ht1 + ht2 + ht3, ht_gap = unit(0.5, "cm", ))
png(file = "~/plots/cellchat/liver/Incoming_signalling_Heatmap_LIVER.png", width = 5200, height = 5200, res = 300)
draw(combined_heatmap, ht_gap = unit(0.5, "cm"))
dev.off()

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 10, height = 40, color.heatmap = "Oranges")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 10, height = 40, color.heatmap = "Oranges")
ht3 = netAnalysis_signalingRole_heatmap(object.list[[i+2]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 10, height = 40, color.heatmap = "Oranges")

combined_heatmap <- draw(ht1 + ht2 + ht3, ht_gap = unit(0.5, "cm", ))
png(file = "~/plots/cellchat/liver/Outgoing_signalling_Heatmap_LIVER.png", width =5200, height = 5200, res = 300)
draw(combined_heatmap, ht_gap = unit(0.5, "cm"))
dev.off()

# Identify altered signaling with distinct interaction strength
gg <- rankNet(cellchat, mode = "comparison",  comparison = c(2, 3), measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE, color.use = c("blue", "red"))
ggsave("AlteredSignalling_LIVER.png", gg, dpi = 300, width = 6, height = 14)

save(object.list, file = "cellchat_object.list_H_F0_F1_2.RData")
save(cellchat, file = "cellchat_merged_H_F0_F1_2.RData")


#### Load CellChat object of F0 and F1_2 dataset ####
cellchatF0 <- readRDS("cellchatF0.rds")
cellchatF1_2 <- readRDS("cellchatF1_2.rds")

# Merge F0 and F1_2 datasets for comparison
object.list <- list(F0 = cellchatF0, F1_2 = cellchatF1_2)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

# Differential number of interactions among different cell populations across two datasets
par(mfrow = c(1,1), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)

# Identify the up-regulated and down-regulated signaling ligand-receptor pairs
netVisual_bubble(cellchat, sources.use = 12, targets.use = c(1:23),  comparison = c(2, 1), angle.x = 45)

#### Identify dysfunctional signaling in Hepatocytes by comparing the communication probabilities ####

gg1 <- netVisual_bubble(cellchat, sources.use = 12, targets.use = c(1:23),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased hepatocytes signaling in F1_2", angle.x = 45, remove.isolate = T, color.text = c('black', 'blue'), font.size = 16)

plot <- gg1 + theme(
              plot.title = element_text(size=12, margin = margin(b = 5),hjust=0,vjust=0.5, family="Helvetica", face="bold"),
              title = element_text(size = 12, margin = margin(b = 5),hjust=0,vjust=0.5, family="Helvetica", face="bold"),
              legend.text = element_text(size=12,),
              legend.title= element_text(size=12,), 
              legend.key.size=unit(1,"line"),
              plot.margin=unit(c(0.4,0.4,0.4,0.4), "cm"),
              strip.text.x = element_text(size = 12, family="Helvetica", face="bold", vjust=1))
ggsave(
 plot = plot,
 filename = ("~/plots/cellchat/liver/up_hepatocytes_signaling_F1_2_comm_prob.tiff"),
 height = 8,
 width = 10,
 dpi = 300,
 units = 'in'
)

# Write data in table
write.csv(gg1$data, '~/files/cellchat/comm_prob_onlyliver_hep_all_cluster.csv')

gg2 <- netVisual_bubble(cellchat, sources.use = 12, targets.use = c(1:23),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased hepatocytes signaling in F1_2", angle.x = 45, remove.isolate = T, color.text = c('black', 'blue'), font.size = 16)

plot1 <- gg2 + theme(
 plot.title = element_text(size=12, margin = margin(b = 5),hjust=0,vjust=0.5, family="Helvetica", face="bold"),
 title = element_text(size = 12, margin = margin(b = 5),hjust=0,vjust=0.5, family="Helvetica", face="bold"),
 legend.text = element_text(size=12,),
 legend.title= element_text(size=12,), 
 legend.key.size=unit(1,"line"),
 plot.margin=unit(c(0.4,0.4,0.4,0.4), "cm"),
 strip.text.x = element_text(size = 12, family="Helvetica", face="bold", vjust=1))
ggsave(
 plot = plot,
 filename = ("~/plots/cellchat/liver/down_hepatocytes_signaling_F1_2_comm_prob.tiff"),
 height = 8,
 width = 10,
 dpi = 300,
 units = 'in'
)

patchwork::wrap_plots(plots = plot,plot1)

#### Identify dysfunctional signaling in cluster 0, 1, and 2 by comparing the communication probabilities ####

# cluster 0
gg1 <- netVisual_bubble(cellchat, sources.use = 1, targets.use = c(1:23),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased cluster 0 signaling in F1_2", angle.x = 45, remove.isolate = T, color.text = c('black', 'blue'), font.size = 16)

plot <- gg1 + theme(
 plot.title = element_text(size=12, margin = margin(b = 5),hjust=0,vjust=0.5, family="Helvetica", face="bold"),
 title = element_text(size = 12, margin = margin(b = 5),hjust=0,vjust=0.5, family="Helvetica", face="bold"),
 legend.text = element_text(size=12,),
 legend.title= element_text(size=12,), 
 legend.key.size=unit(1,"line"),
 plot.margin=unit(c(0.4,0.4,0.4,0.4), "cm"),
 strip.text.x = element_text(size = 12, family="Helvetica", face="bold", vjust=1))
ggsave(
 plot = plot,
 filename = ("~/plots/cellchat/liver/up_cluster0_signaling_F1_2_comm_prob.tiff"),
 height = 14,
 width = 14,
 dpi = 300,
 units = 'in'
)
# Write data in table
write.csv(gg1$data, '~/files/cellchat/comm_prob_onlyliver_cluster0_all_cluster.csv')

# cluster 1
gg1 <- netVisual_bubble(cellchat, sources.use = 2, targets.use = c(1:23),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased cluster 1 signaling in F1_2", angle.x = 45, remove.isolate = T, color.text = c('black', 'blue'), font.size = 16)

plot <- gg1 + theme(
 plot.title = element_text(size=12, margin = margin(b = 5),hjust=0,vjust=0.5, family="Helvetica", face="bold"),
 title = element_text(size = 12, margin = margin(b = 5),hjust=0,vjust=0.5, family="Helvetica", face="bold"),
 legend.text = element_text(size=12,),
 legend.title= element_text(size=12,), 
 legend.key.size=unit(1,"line"),
 plot.margin=unit(c(0.4,0.4,0.4,0.4), "cm"),
 strip.text.x = element_text(size = 12, family="Helvetica", face="bold", vjust=1))
ggsave(
 plot = plot,
 filename = ("~/plots/cellchat/liver/up_cluster1_signaling_F1_2_comm_prob.tiff"),
 height = 14,
 width = 14,
 dpi = 300,
 units = 'in'
)
# Write data in table
write.csv(gg1$data, '~/files/cellchat/comm_prob_onlyliver_cluster1_all_cluster.csv')

# cluster 2
gg1 <- netVisual_bubble(cellchat, sources.use = 3, targets.use = c(1:23),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased cluster 2 signaling in F1_2", angle.x = 45, remove.isolate = T, color.text = c('black', 'blue'), font.size = 16)

plot <- gg1 + theme(
 plot.title = element_text(size=12, margin = margin(b = 5),hjust=0,vjust=0.5, family="Helvetica", face="bold"),
 title = element_text(size = 12, margin = margin(b = 5),hjust=0,vjust=0.5, family="Helvetica", face="bold"),
 legend.text = element_text(size=12,),
 legend.title= element_text(size=12,), 
 legend.key.size=unit(1,"line"),
 plot.margin=unit(c(0.4,0.4,0.4,0.4), "cm"),
 strip.text.x = element_text(size = 12, family="Helvetica", face="bold", vjust=1))
ggsave(
 plot = plot,
 filename = ("~/plots/cellchat/liver/up_cluster2_signaling_F1_2_comm_prob.tiff"),
 height = 14,
 width = 14,
 dpi = 300,
 units = 'in'
)
# Write data in table
write.csv(gg1$data, '~/files/cellchat/comm_prob_onlyliver_cluster2_all_cluster.csv')

#### Identify dysfunctional signaling by using differential expression analysis ####
# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "F1_2"
# define a char name used for storing the results of differential expression analysis
features.name = paste0(pos.dataset, ".merged")

# perform differential expression analysis 
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.05,thresh.p = 0.05, group.DE.combined = FALSE) 

# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name, variable.all = TRUE)

write.csv(net, '~/files/cellchat/onlyliver_interactions.csv')

# extract the ligand-receptor pairs with upregulated ligands in F1_2
net.up <- subsetCommunication(cellchat, net = net, datasets = "F1_2",ligand.logFC = 0.05, receptor.logFC = NULL)

write.csv(net.up, '~/files/cellchat/onlyliver_up_interactions.csv')

# extract the ligand-receptor pairs with upregulated ligands and upregulated receptors in F0, i.e.,downregulated in F1_2
net.down <- subsetCommunication(cellchat, net = net, datasets = "F0",ligand.logFC = -0.05, receptor.logFC = NULL)

write.csv(net.down, '~/files/cellchat/onlyliver_down_interactions.csv')

gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

# Visualize the identified up-regulated and down-regulated signaling ligand-receptor pairs
pairLR.use.up = net.up[, "interaction_name", drop = F]
gg3 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 12, targets.use = c(1:23), comparison = c(2, 1),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))

# Write data in table
write.csv(gg3$data, '~/files/cellchat/DEG_onlyliver_hep_all_cluster.csv')

pairLR.use.down = net.down[, "interaction_name", drop = F]
gg4 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = 12, targets.use = c(1:23), comparison = c(2, 1),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))

gg3 + gg4

# Visualize the identified up-regulated signaling ligand-receptor pairs using both methods, DEG and communication probabilities
gg1 + gg3

# visualize the enriched ligands in the F0 condition
computeEnrichmentScore(net.down, species = 'human', variable.both = TRUE)

# visualize the enriched ligands in the F1_2 condition
computeEnrichmentScore(net.up, species = 'human', variable.both = TRUE)
