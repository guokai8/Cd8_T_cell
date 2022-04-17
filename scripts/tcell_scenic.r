library(Seurat)
library(foreach)
library(SCENIC)
library(Seurat)
library(GENIE3)
library(AUCell)
library(RcisTarget)
#library(loomR) 
exprMat <- samf@assays$RNA@counts
exprMat <- as.matrix(exprMat)
##meta information
cellInfo <- data.frame(samf@meta.data)
colnames(cellInfo)[which(colnames(cellInfo)=="Time")] <- "Group"
cellInfo <- cellInfo[,c("seurat_clusters","Group")]
colnames(cellInfo)[1]<-"CellType"
saveRDS(cellInfo, file="cellInfo.Rds",compress = T)
###download the data mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather and mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather
dbs <- c('mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather', 'mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather')
dbDir='../cisTarget_databases/' # RcisTarget databases location
scenicOptions <- initializeScenic(org="mgi", dbDir="../cisTarget_databases",dbs=dbs, nCores=40)
####
scenicOptions@inputDatasetInfo$cellInfo<-"cellInfo.Rds"
saveRDS(scenicOptions, file="scenicOptions.Rds",compress=T)
###
genesKept <- geneFiltering(exprMat, scenicOptions, minCountsPerGene = 3 * 0.01 * ncol(exprMat), minSamples = ncol(exprMat) * 0.01)
exprMat_filtered <- exprMat[genesKept,]
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1)
runGenie3(exprMat_filtered_log, scenicOptions)
runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions) 
runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered_log) 
runSCENIC_4_aucell_binarize(scenicOptions)
######
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$CellType),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))

###
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellInfo[colnames(regulonAUC), "CellType"])
rss_df<-melt(rss)
colnames(rss_df)<-c("regulon","cell_type","RSS")

rss_df%>%filter(!grepl('extended',regulon))%>%group_by(cell_type)%>%arrange(desc(RSS),.by_group=T)%>%slice_head(n=5)%>%
  ggplot(aes(factor(cell_type),regulon,color=factor(cell_type)))+geom_point(size=3)+scale_color_manual(values=pcols)+
  theme_light(base_size=12)+xlab("")+ylab("")+labs(color="Cluster")
#########
binaryRegulonActivity <- loadInt(scenicOptions, "aucell_binary_nonDupl")
cellInfo<-cellInfo[order(cellInfo$CellType,cellInfo$Group),]
binaryRegulonActivity<-binaryRegulonActivity[,rownames(cellInfo)]
anncol<-samf@meta.data[,c("seurat_clusters","Time","Group")]
ann_colors<-list(seurat_clusters=pcols,Group=ccols, Time=c("Mock" = "#7570B3", "D7" = "#E7298A", "D14" = "#66A61E"))
anncol<-anncol[order(anncol$seurat_clusters,anncol$Group,anncol$Time),]
pheatmap(binaryRegulonActivity[grep('extend',rownames(binaryRegulonActivity),invert = T),rownames(anncol)],cluster_rows = T,cluster_cols = F,color = colorRampPalette(c("white","darkgreen"))(100), breaks=seq(0, 1, length.out = 100),
         treeheight_row=10,  border_color=NA,annotation_col = anncol,show_colnames = F,fontsize_row = 3,annotation_colors = ann_colors)
dev.print(pdf,file="binary_heatmap.pdf")
###########################################################

umap<-Embeddings(samf,reduction = "umap")
umap<-cbind(umap,meta[rownames(umap),c("Group","celltype","Time")])
aucmtx<-t(as.matrix(getAUC(regulonAUC)))
umap<-cbind(umap,aucmtx[rownames(umap),])
thres <- loadInt(scenicOptions, "aucell_thresholds")

thres2<-unlist(lapply(thres, function(x)x[[1]][[1]]))
names(thres2)<-sub('\\.L_k2','',names(thres2))

ggplot(umap,aes(UMAP_1,UMAP_2))+
  geom_point(aes(color=`Runx3 (65g)`),size=0.1,alpha=0.75,shape =19)+
  scale_color_gradient2(high="red",low = "white",mid="lightgrey",midpoint = thres2["Runx3 (65g)"])+facet_wrap(~Time)+
  theme_classic()
#####
ggplot(umap,aes(UMAP_1,UMAP_2))+
geom_point(aes(color=`Runx3 (65g)`),size=0.1,alpha=0.75,shape =19)+
scale_color_gradient2(high="red",low = "ghostwhite",mid="lightgrey",midpoint = thres2["Runx3 (65g)"])+facet_wrap(~Time)+
theme_classic()
dev.print(pdf,file="Runx3_tf.pdf")
ggplot(umap,aes(UMAP_1,UMAP_2))+
geom_point(aes(color=`Runx2 (227g)`),size=0.1,alpha=0.75,shape =19)+
scale_color_gradient2(high="red",low = "ghostwhite",mid="lightgrey",midpoint = thres2["Runx2 (227g)"])+facet_wrap(~Time)+
theme_classic()
dev.print(pdf,file="Runx2_tf.pdf")
ggplot(umap,aes(UMAP_1,UMAP_2))+
geom_point(aes(color=`Irf8 (388g)`),size=0.1,alpha=0.75,shape =19)+
scale_color_gradient2(high="red",low = "ghostwhite",mid="lightgrey",midpoint = thres2["Irf8 (388g)"])+facet_wrap(~Time)+
theme_classic()
dev.print(pdf,file="Irf8_tf.pdf")
ggplot(umap,aes(UMAP_1,UMAP_2))+
geom_point(aes(color=`Irf7 (262g)`),size=0.1,alpha=0.75,shape =19)+
scale_color_gradient2(high="red",low = "ghostwhite",mid="lightgrey",midpoint = thres2["Irf7 (262g)"])+facet_wrap(~Time)+
theme_classic()
dev.print(pdf,file="Irf7_tf.pdf")
ggplot(umap,aes(UMAP_1,UMAP_2))+
geom_point(aes(color=`Stat1 (279g)`),size=0.1,alpha=0.75,shape =19)+
scale_color_gradient2(high="red",low = "ghostwhite",mid="lightgrey",midpoint = thres2["Stat1 (279g)"])+facet_wrap(~Time)+
theme_classic()
dev.print(pdf,file="Stat1_tf.pdf")
ggplot(umap,aes(UMAP_1,UMAP_2))+
geom_point(aes(color=`Mafb (51g)`),size=0.1,alpha=0.75,shape =19)+
scale_color_gradient2(high="red",low = "ghostwhite",mid="lightgrey",midpoint = thres2["Mafb (51g)"])+facet_wrap(~Time)+
theme_classic()
dev.print(pdf,file="Mafb_tf.pdf")
ggplot(umap,aes(UMAP_1,UMAP_2))+
geom_point(aes(color=`Irf5 (214g)`),size=0.1,alpha=0.75,shape =19)+
scale_color_gradient2(high="red",low = "ghostwhite",mid="lightgrey",midpoint = thres2["Irf5 (214g)"])+facet_wrap(~Time)+
theme_classic()
dev.print(pdf,file="Irf5_tf.pdf")
ggplot(umap,aes(UMAP_1,UMAP_2))+
geom_point(aes(color=`E2f3 (484g)`),size=0.1,alpha=0.75,shape =19)+
scale_color_gradient2(high="red",low = "ghostwhite",mid="lightgrey",midpoint = thres2["E2f3 (484g)"])+facet_wrap(~Time)+
theme_classic()
dev.print(pdf,file="E2f3_tf.pdf")
ggplot(umap,aes(UMAP_1,UMAP_2))+
geom_point(aes(color=`E2f2 (367g)`),size=0.1,alpha=0.75,shape =19)+
scale_color_gradient2(high="red",low = "ghostwhite",mid="lightgrey",midpoint = thres2["E2f2 (367g)"])+facet_wrap(~Time)+
theme_classic()
dev.print(pdf,file="E2f2_tf.pdf")
ggplot(umap,aes(UMAP_1,UMAP_2))+
geom_point(aes(color=`Brca1 (380g)`),size=0.1,alpha=0.75,shape =19)+
scale_color_gradient2(high="red",low = "ghostwhite",mid="lightgrey",midpoint = thres2["Brca1 (380g)"])+facet_wrap(~Time)+
theme_classic()
dev.print(pdf,file="Brca1_tf.pdf")
ggplot(umap,aes(UMAP_1,UMAP_2))+
geom_point(aes(color=`E2f1 (348g)`),size=0.1,alpha=0.75,shape =19)+
scale_color_gradient2(high="red",low = "ghostwhite",mid="lightgrey",midpoint = thres2["E2f1 (348g)"])+facet_wrap(~Time)+
theme_classic()
dev.print(pdf,file="E2f1_tf.pdf")
ggplot(umap,aes(UMAP_1,UMAP_2))+
geom_point(aes(color=`Elf1 (6536g)`),size=0.1,alpha=0.75,shape =19)+
scale_color_gradient2(high="red",low = "ghostwhite",mid="lightgrey",midpoint = thres2["Elf1 (6536g)"])+facet_wrap(~Time)+
theme_classic()
dev.print(pdf,file="Elf1_tf.pdf")
ggplot(umap,aes(UMAP_1,UMAP_2))+
geom_point(aes(color=`Mef2a (27g)`),size=0.1,alpha=0.75, shape =19)+
scale_color_gradient2(high="red",low = "ghostwhite",mid="lightgrey",midpoint = thres2["Mef2a (27g)"])+facet_wrap(~Time)+
theme_classic()
dev.print(pdf,file="Mef2a_tf.pdf")
ggplot(umap,aes(UMAP_1,UMAP_2))+
geom_point(aes(color=`Ets1 (5505g)`),size=0.1,alpha=0.75,shape =19)+
scale_color_gradient2(high="red",low = "ghostwhite",mid="lightgrey",midpoint = thres2["Ets1 (5505g)"])+facet_wrap(~Time)+
theme_classic()
dev.print(pdf,file="Ets1_tf.pdf")
ggplot(umap,aes(UMAP_1,UMAP_2))+
geom_point(aes(color=`Spi1 (656g)`),size=0.1,alpha=0.75,shape = 19)+
scale_color_gradient2(high="red",low = "ghostwhite",mid="lightgrey",midpoint = thres2["Spi1 (656g)"])+facet_wrap(~Time)+
theme_classic()
dev.print(pdf,file="Spi1_tf.pdf")
ggplot(umap,aes(UMAP_1,UMAP_2))+
geom_point(aes(color=`Spi1 (656g)`),size=0.1,alpha=0.75, shape =19)+
scale_color_gradient2(high="red",low = "ghostwhite",mid="lightgrey",midpoint = thres2["Spi1 (656g)"])+facet_wrap(~Time)+
theme_classic()




