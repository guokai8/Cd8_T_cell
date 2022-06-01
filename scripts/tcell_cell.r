library(CellChat)
###CellChat_1.1.3
library(Seurat)
library(tidyverse)
samf$seurat_clusters<-paste0("C",as.vector(samf$seurat_clusters))
Idents(samf)<-"seurat_clusters"
mockd<-subset(samf,Time=="Mock")
d7d<-subset(samf,Time=="D7")
d14d<-subset(samf,Time=="D14")
####
data.input <- GetAssayData(mockd, assay = "RNA", slot = "data") # normalized data matrix
labels <- Idents(mockd)
meta <- data.frame(labels = labels, row.names = names(labels)) # create a dataframe of the cell labels
mc <- createCellChat(object = data.input, meta = meta, group.by = "labels")
mc <- setIdent(mc, ident.use = "labels") # set "labels" as default cell identity
levels(mc@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(mc@idents)) # number of cells in each cell group
mcDB <- CellChatDB.mouse # use mcDB.human if running on human data
mcDB.use <- CellChatDB.mouse
showDatabaseCategory(mcDB)
#mcDB.use <- subsetDB(mcDB, search = "Secreted Signaling") 
mc@DB<-mcDB.use
mc <- subsetData(mc) # subset the expression data of signaling genes for saving computation cost
mc <- identifyOverExpressedGenes(mc)
mc <- identifyOverExpressedInteractions(mc)
mc <- projectData(mc, PPI.mouse)
mc <- computeCommunProb(mc)
mc <- filterCommunication(mc, min.cells = 10)
mc <- computeCommunProbPathway(mc)
mc <- aggregateNet(mc)
mc <- netAnalysis_computeCentrality(mc, slot.name = "netP") 
pathways.show <- mc@netP$pathways
##
dir.create("mock")
setwd("mock")
vertex.receiver = seq(1,5)
library(igraph)
for (i in 1:length(pathways.show)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(mc, signaling = pathways.show[i], vertex.receiver = vertex.receiver, vertex.weight = groupSize,color.use = pcols,out.format = "pdf" )
  netVisual(mc, signaling = pathways.show[i], vertex.receiver = vertex.receiver, vertex.weight = groupSize,layout ="circle",color.use = pcols,out.format = "pdf" )
  netVisual(mc, signaling = pathways.show[i],vertex.receiver = vertex.receiver, vertex.weight = groupSize, layout = "chord",color.use = pcols,out.format = "pdf" )
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(mc, signaling = pathways.show[i])
  ggsave(filename=paste0(pathways.show[i], "_L-R_contribution_ctrl.pdf"), plot=gg, width = 5, height = 2)
  # Visualize signaling roles of cell groups
  grDevices::pdf(file = paste0("cellSignalingRole_",pathways.show[i],".pdf"))
  netAnalysis_signalingRole_network (mc, pathways.show[i],  width = 10,color.use = pcols)
  dev.off()
}
setwd("../")
####
data.input <- GetAssayData(d7d, assay = "RNA", slot = "data") # normalized data matrix
labels <- Idents(d7d)
meta <- data.frame(labels = labels, row.names = names(labels)) # create a dataframe of the cell labels
d7 <- createCellChat(object = data.input, meta = meta, group.by = "labels")
d7 <- setIdent(d7, ident.use = "labels") # set "labels" as default cell identity
levels(d7@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(d7@idents)) # number of cells in each cell group
d7DB <- CellChatDB.mouse # use d7DB.human if running on human data
d7DB.use <- CellChatDB.mouse
showDatabaseCategory(d7DB)
#d7DB.use <- subsetDB(d7DB, search = "Secreted Signaling") 
d7@DB<-d7DB.use
#CellChatDB.use = CellChatDB.mouse
#d7@DB<-cellchatDB.use
d7 <- subsetData(d7) # subset the expression data of signaling genes for saving computation cost
d7 <- identifyOverExpressedGenes(d7)
d7 <- identifyOverExpressedInteractions(d7)
d7 <- projectData(d7, PPI.mouse)
d7 <- computeCommunProb(d7)
d7 <- filterCommunication(d7, min.cells = 10)
d7 <- computeCommunProbPathway(d7)
d7 <- aggregateNet(d7)
d7 <- netAnalysis_computeCentrality(d7, slot.name = "netP") 
pathways.show <- d7@netP$pathways
##
dir.create("D7")
setwd("D7")
vertex.receiver = seq(1,5)
for (i in 1:length(pathways.show)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(d7, signaling = pathways.show[i], vertex.receiver = vertex.receiver, vertex.weight = groupSize,color.use = pcols,out.format = "pdf" )
  netVisual(d7, signaling = pathways.show[i], vertex.receiver = vertex.receiver, vertex.weight = groupSize,layout ="circle",color.use = pcols,out.format = "pdf" )
  netVisual(d7, signaling = pathways.show[i],vertex.receiver = vertex.receiver, vertex.weight = groupSize, layout = "chord",color.use = pcols,out.format = "pdf" )
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(d7, signaling = pathways.show[i])
  ggsave(filename=paste0(pathways.show[i], "_L-R_contribution_ctrl.pdf"), plot=gg, width = 5, height = 2)
  # Visualize signaling roles of cell groups
  grDevices::pdf(file = paste0("cellSignalingRole_",pathways.show[i],".pdf"))
  netAnalysis_signalingRole_network (d7, pathways.show[i],  width = 10,color.use = pcols)
  dev.off()
}
setwd("../")
###
data.input <- GetAssayData(d14d, assay = "RNA", slot = "data") # normalized data matrix
labels <- Idents(d14d)
meta <- data.frame(labels = labels, row.names = names(labels)) # create a dataframe of the cell labels
d14 <- createCellChat(object = data.input, meta = meta, group.by = "labels")
d14 <- setIdent(d14, ident.use = "labels") # set "labels" as default cell identity
levels(d14@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(d14@idents)) # number of cells in each cell group
d14DB <- CellChatDB.mouse # use d14DB.human if running on human data
d14DB.use <- CellChatDB.mouse
showDatabaseCategory(d14DB)
#d14DB.use <- subsetDB(d14DB, search = "Secreted Signaling") 
d14@DB<-d14DB.use
d14 <- subsetData(d14) # subset the expression data of signaling genes for saving computation cost
d14 <- identifyOverExpressedGenes(d14)
d14 <- identifyOverExpressedInteractions(d14)
d14 <- projectData(d14, PPI.mouse)
d14 <- computeCommunProb(d14)
d14 <- filterCommunication(d14, min.cells = 10)
d14 <- computeCommunProbPathway(d14)
d14 <- aggregateNet(d14)
d14 <- netAnalysis_computeCentrality(d14, slot.name = "netP") 
pathways.show <- d14@netP$pathways
##
dir.create("D14")
setwd("D14")
vertex.receiver = seq(1,5)
for (i in 1:length(pathways.show)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(d14, signaling = pathways.show[i], vertex.receiver = vertex.receiver, vertex.weight = groupSize,color.use = pcols,out.format = "pdf" )
  netVisual(d14, signaling = pathways.show[i], vertex.receiver = vertex.receiver, vertex.weight = groupSize,layout ="circle",color.use = pcols,out.format = "pdf" )
  netVisual(d14, signaling = pathways.show[i],vertex.receiver = vertex.receiver, vertex.weight = groupSize, layout = "chord",color.use = pcols,out.format = "pdf" )
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(d14, signaling = pathways.show[i])
  ggsave(filename=paste0(pathways.show[i], "_L-R_contribution_ctrl.pdf"), plot=gg, width = 5, height = 2)
  # Visualize signaling roles of cell groups
  grDevices::pdf(file = paste0("cellSignalingRole_",pathways.show[i],".pdf"))
  netAnalysis_signalingRole_network (d14, pathways.show[i],  width = 10,color.use = pcols)
  dev.off()
}
setwd("../")
object.list <- list(Mock = mc,D7=d7,D14=d14)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
##
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax,color.use = pcols)+theme_classic()+xlim(0,12)+ylim(0,18)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)
dev.print(pdf,file="all_interaction_v3.pdf")
#####
i = 1
# combining all the identified signaling pathways from different datasets 
#####
library(ComplexHeatmap)
pathway.union <- union(object.list[[1]]@netP$pathways, union(object.list[[2]]@netP$pathways,object.list[[3]]@netP$pathways))
ht1 = netAnalysis_signalingRole_heatmap(object.list[[1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[1], width = 6, height = 12,color.use = pcols)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[2]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[2], width = 6, height = 12,color.use = pcols)
ht3 = netAnalysis_signalingRole_heatmap(object.list[[3]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[3], width = 6, height = 12,color.use = pcols)

draw(ht1 + ht2+ht3, ht_gap = unit(0.2, "cm"))
dev.print(pdf,file="f3B.pdf")
ht1 = netAnalysis_signalingRole_heatmap(object.list[[1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[1],color.use = pcols,width = 6, height = 12,color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[2]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[2],color.use = pcols,width = 6, height = 12,color.heatmap = "GnBu")
ht3 = netAnalysis_signalingRole_heatmap(object.list[[3]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[3],color.use = pcols,width = 6, height = 12,color.heatmap = "GnBu")
draw(ht1 + ht2+ht3, ht_gap = unit(0.2, "cm"))
dev.print(pdf,file="f3C.pdf")
############### 
