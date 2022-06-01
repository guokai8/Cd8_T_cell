library(Seurat)
library(tidyverse)
##MAST MAST_1.18.0 Seurat 3.2.2
mock.data<-Read10X("data/mock/filtered_feature_bc_matrix/")
d7.data<-Read10X("data/d7/filtered_feature_bc_matrix/")
d14.data<-Read10X("data/d14/filtered_feature_bc_matrix/")
mock<-CreateSeuratObject(counts = mock.data, project = "Mock",min.cells = 3, min.features = 200)
d7<-CreateSeuratObject(counts = d7.data, project = "D7",min.cells = 3, min.features = 200)
d14<-CreateSeuratObject(counts = d14.data, project = "D14",min.cells = 3, min.features = 200)
### assign group
mock$group<-"Mock"
d7$group<-"D7"
d14$group<-"D14"
##
mock[["percent.mt"]]<-PercentageFeatureSet(mock,pattern = "^mt-")
pp1<-VlnPlot(mock, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
d7[["percent.mt"]]<-PercentageFeatureSet(d7,pattern = "^mt-")
pp2<-VlnPlot(d7, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
d14[["percent.mt"]]<-PercentageFeatureSet(d14,pattern = "^mt-")
pp3<-VlnPlot(d14, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
##
mock <- subset(mock, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 25 )
d7 <- subset(d7, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 25 )
d14 <- subset(d14, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 25 )
#
mock <- NormalizeData(mock)
d7 <- NormalizeData(d7)
d14 <- NormalizeData(d14)
#
mock <- FindVariableFeatures(mock, selection.method = "vst", nfeatures = 2000)
d7 <- FindVariableFeatures(d7, selection.method = "vst", nfeatures = 2000)
d14 <- FindVariableFeatures(d14, selection.method = "vst", nfeatures = 2000)
#
sample.anchors <- FindIntegrationAnchors(object.list = list(mock,d7,d14), dims = 1:50)
sample.combined <- IntegrateData(anchorset = sample.anchors, dims = 1:50)
############
DefaultAssay(sample.combined) <- "integrated"
all.genes <- rownames(sample.combined)
sample.combined <- ScaleData(sample.combined, verbose = TRUE,features = all.genes)
sample.combined <- RunPCA(sample.combined, npcs = 50, features = VariableFeatures(object = sample.combined))
##
pct <- sample.combined[["pca"]]@stdev / sum(sample.combined[["pca"]]@stdev) * 100
# Calculate cumulative percents for each PC
cumu <- cumsum(pct)
# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
co1
##################
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
# last point where change of % of variation is more than 0.1%.
co2
######
pcs <- min(co1, co2)
pcs

###Choose pcs=12
###
sample.combined <- RunUMAP(sample.combined, reduction = "pca", dims = 1:12)
sample.combined <- FindNeighbors(sample.combined, reduction = "pca", dims = 1:12)
sample.combined <- FindClusters(sample.combined, resolution = 0.8)
###
DimPlot(sample.combined, reduction = "umap", split.by = "group",label = T)
library(VennDetail)
mycol<-setcolor(16)
mycol[15]<-"#90C423"
#
mycol<-mycol[c(0,3,2,1,4,5,6,7,11,9,10,8,12,13,14,15)+1]
##
DimPlot(sample.combined, reduction = "umap", label = T,cols=mycol)
#
sample.combined$Time <- sample.combined$group
meta<-sample.combined@meta.data
library(ggrepel)
meta%>%group_by(Time,seurat_clusters)%>%summarise(count=n())%>%mutate(cell=count/sum(count))%>%ggplot(aes(Time,cell,fill=seurat_clusters))+
geom_bar(stat="identity")+
scale_fill_manual(values=mycol)+theme_classic(base_size = 15)+xlab("")+ylab("")
dev.print(pdf,file="proportion_cluster_number.pdf")
meta%>%group_by(Time,seurat_clusters)%>%summarise(count=n())%>%mutate(cell=count/sum(count))%>%ggplot(aes(Time,cell,fill=seurat_clusters))+
geom_bar(stat="identity")+
scale_fill_manual(values=mycol)+theme_classic(base_size = 14)+xlab("")+ylab("")
dev.print(pdf,file="cell_pro.pdf")

#
DefaultAssay(sample.combined) <- "RNA"
sample.markers <- FindAllMarkers(sample.combined, min.pct = 0.25, logfc.threshold = 0.25)
####
####Please make sure the clusters numbers due to the random seed setting by Seurat
## 4,8,9,10,11,13,15

clusters<-c("CD8+ Tn","CD8+ Tn","C2","CD8+ Tem",
            "CD8+ Tem","CD8+ Tn","CD8+ Tem",
            "CD8+ Tcm","CD8+ Tn","C9",'C10',"C11",
            "C12","C13","CD8+ Tcm","C15")
names(clusters) <- levels(sample.combined)
sample.combined <- RenameIdents(sample.combined, clusters)
sample.combined$celltype<-Idents(sample.combined) 
#
DotPlot(sample.combined,features = c("Cd8a","Cd8b1"))+scale_color_viridis_c()
dev.print(pdf,file="Cd8_marker.pdf")
## remove clusters with low or no CD8 expressed
samf<-subset(sample.combined,celltype%in%c("CD8+ Tn","CD8+ Tem","CD8+ Tcm"))
samf$group<-factor(samf$group,levels=c("Mock","D7","D14"))
samf$Time<-samf$group
samf$seurat_clusters<-as.vector(samf$seurat_clusters)
samf$seurat_clusters<-paste0("C",samf$seurat_clusters)
meta<-samf@meta.data
pcols<-mycol[unique(samf$seurat_clusters)]
cd8col<-c("#45B35F","#229885","#9D1E2A")
names(cd8col)<-c("CD8+ Tn","CD8+ Tem","CD8+ Tcm")
Idents(samf)<-"celltype"
meta%>%group_by(Time,celltype)%>%summarise(count=n())%>%mutate(cell=count/sum(count))%>%ggplot(aes(Time,cell,fill=celltype))+
  geom_bar(stat="identity")+geom_text(aes(label = paste0(round(100*cell,1), "%")),position = position_stack())+theme_classic(base_size = 15)+xlab("")+ylab("")+labs(fill="cluster")+
  scale_fill_manual(values=cd8col)
dev.print(pdf,file="f1m.pdf")
#
meta%>%filter(seurat_clusters%in%c("C3","C4","C6"))%>%group_by(Time,seurat_clusters)%>%summarise(count=n())%>%mutate(cell=count/sum(count))%>%ggplot(aes(Time,cell,fill=seurat_clusters))+
geom_bar(stat="identity")+geom_text(aes(label = paste0(round(100*cell,1), "%")),position = position_stack())+
scale_fill_manual(values=pcols[c("C3","C4","C6")])+theme_classic(base_size = 15)+xlab("")+ylab("")+labs(fill="cluster")
dev.print(pdf,file="f1n.pdf")

##
Idents(samf)<-"seurat_clusters"
###
genes<-c("Ccr7","Lef1","Sell","Tcf7","Cd27","Cd28","S1pr1","Il7r","Gzma","Ccl5","Gpr183","Gzmk","Cxcr4","Cxcr3","Cd44")
DotPlot(samf,features = genes,cluster.idents = T)+coord_flip()+theme(axis.text.x=element_text(angle=90,vjust = 0.5,hjust = 1))+scale_color_viridis_c()
Idents(samf)<-"seurat_clusters"
samf$seurat_clusters<-factor(samf$seurat_clusters,levels=c("C0","C1","C5","C8","C7","C14","C3","C4","C6"))
dev.print(pdf,file="f1L.pdf")
###
samf$Group<-ifelse(samf$celltype=="CD8+ Tem","Cxcr3High","Cxcr3Low")
Idents(samf)<-"Group"
DimPlot(samf,cols=c("cyan4","darkred"))
dev.print(pdf,file="f2A.pdf")
Idents(samf)<-"celltype"
VlnPlot(samf,features = "Cxcr3",pt.size = 0,cols=c("grey","darkred","grey"))
dev.print(pdf,file="f2B.pdf")
###
samf$time <- paste0(samf$group,"_",samf$Group)
Idents(samf)<-"time"
### DEGs Cxcr3 high vs Cxcr3 low
cxcrm<-FindMarkers(samf,ident.1 = "Mock_Cxcr3High",ident.2 = "Mock_Cxcr3Low",test.use = "MAST",logfc.threshold = 0)
cxcr7<-FindMarkers(samf,ident.1 = "D7_Cxcr3High",ident.2 = "D7_Cxcr3Low",test.use = "MAST",logfc.threshold = 0)
cxcr14<-FindMarkers(samf,ident.1 = "D14_Cxcr3High",ident.2 = "D14_Cxcr3Low",test.use = "MAST",logfc.threshold = 0)
cxcrmd<-subset(cxcrm,p_val_adj<0.05&abs(avg_log2FC)>0.2)
cxcr7d<-subset(cxcr7,p_val_adj<0.05&abs(avg_log2FC)>0.2)
cxcr14d<-subset(cxcr14,p_val_adj<0.05&abs(avg_log2FC)>0.2)

##devtools::install_github('guokai8/richR')
library(richR)
mmko<-buildAnnot(species="mouse",keytype="SYMBOL",anntype="KEGG",builtin=FALSE)
gse<-function(x,mmko){
  fc<-x$avg_log2FC
  names(fc)<-rownames(x)
  res<-richGSEA(fc,mmko,minSize = 5)
  return(result.GSEAResult(res))
}
### do GSEA
cxcrmg<-gse(cxcrm,mmko)
cxcr7g<-gse(cxcr7,mmko)
cxcr14g<-gse(cxcr14,mmko)
write.csv(cxcrmg,file="cxcr_mock_gsea.csv")
write.csv(cxcr7g,file="cxcr_D7_gsea.csv")
write.csv(cxcr14g,file="cxcr_D14_gsea.csv")
######
samf$cond <- paste0(samf$group,"_",samf$seurat_clusters)
Idents(samf)<-"cond"
deg7dm<-lapply(c("C3","C4","C6"), function(x)FindMarkers(samf,ident.1 = paste("D7",x,sep="_"),ident.2=paste("Mock",x,sep="_"),test.use = "MAST",logfc.threshold = 0))
deg14dm<-lapply(c("C3","C4","C6"), function(x)FindMarkers(samf,ident.1 = paste("D14",x,sep="_"),ident.2=paste("Mock",x,sep="_"),test.use = "MAST",logfc.threshold = 0))
deg14d7<-lapply(c("C3","C4","C6"), function(x)FindMarkers(samf,ident.1 = paste("D14",x,sep="_"),ident.2=paste("D7",x,sep="_"),test.use = "MAST",logfc.threshold = 0))
### write out
deg7dms<-lapply(deg7dm, function(x)subset(x,p_val_adj<0.05&abs(avg_log2FC)>0.2))
deg14dms<-lapply(deg14dm, function(x)subset(x,p_val_adj<0.05&abs(avg_log2FC)>0.2))
deg14d7s<-lapply(deg14d7, function(x)subset(x,p_val_adj<0.05&abs(avg_log2FC)>0.2))
###                 
sapply(names(deg7dms), function(x)write.csv(deg7dms[[x]],file=paste0(x,"_DEG.csv")))
sapply(names(deg14dms), function(x)write.csv(deg14dms[[x]],file=paste0(x,"_DEG.csv")))
sapply(names(deg14d7s), function(x)write.csv(deg14d7s[[x]],file=paste0(x,"_DEG.csv")))
####
deg7dmg<-lapply(deg7dm, function(x)gse(x,mmko))
#############
deg14dmg<-lapply(deg14dm, function(x)gse(x,mmko))
#############
deg14d7g<-lapply(deg14d7, function(x)gse(x,mmko))       
######
library(VennDetail)
d7ven<-venndetail(list(c1=rownames(deg7dms$D7vsMock_1),c4=rownames(deg7dms$D7vsMock_4),c6=rownames(deg7dms$D7vsMock_6)))
d14ven<-venndetail(list(c1=rownames(deg14dms$D14vsMock_1),c4=rownames(deg14dms$D14vsMock_4),c6=rownames(deg14dms$D14vsMock_6)))
d14v7ven<-venndetail(list(c1=rownames(deg14d7s$D14vsD7_1),c4=rownames(deg14d7s$D14vsD7_4),c6=rownames(deg14d7s$D14vsD7_6)))
####
deg7dmgs<-lapply(deg7dmg,function(x)subset(x,pval<0.05))              
deg14dmgs<-lapply(deg14dmg,function(x)subset(x,pval<0.05))
deg14d7gs<-lapply(deg14d7g,function(x)subset(x,pval<0.05))                
####
d7vengsea<-venndetail(lapply(deg7dmgs[c(2,4,6)],function(x)x$pathway))
d14vengsea<-venndetail(lapply(deg14dmgs[c(2,4,6)],function(x)x$pathway))
d14v7vengsea<-venndetail(lapply(deg14d7gs[c(2,4,6)],function(x)x$pathway))
## then plot the venndiagram
dev.off()
plot(d7vengsea)
dev.print(pdf,file="Venndiagram_E1_GSEA.pdf")
dev.off()
plot(d14vengsea)
dev.print(pdf,file="Venndiagram_E4_GSEA.pdf") 
dev.off()
plot(d14v7vengsea)
dev.print(pdf,file="Venndiagram_E6_GSEA.pdf") 
### GSVA
## devtools::install_github('guokai8/scGSVA')
library(scGSVA)
sc <- scgsva(samf,mmko)
##########################
paths<-c("Viral.protein.interaction.with.cytokine.and.cytokine.receptor","Natural.killer.cell.mediated.cytotoxicity","Antigen.processing.and.presentation",
         "Graft.versus.host.disease","Toll.like.receptor.signaling.pathway",
         "Neutrophil.extracellular.trap.formation","Tight.junction","Focal.adhesion","Leukocyte.transendothelial.migration","PD.L1.expression.and.PD.1.checkpoint.pathway.in.cancer",
         "Chemokine.signaling.pathway","Necroptosis","Apoptosis","PI3K.Akt.signaling.pathway",
         "NOD.like.receptor.signaling.pathway","Influenza.A")
Heatmap(sc,features=paths,group_by="seurat_clusters",cluster_rows=F,cluster_cols=F,average = T)
featurePlot(sc,features="Natural.killer.cell.mediated.cytotoxicity",group_by = "group",label="seurat_cluster")

DotPlot(samf,features = c("Ifnar1","Ifnar2","Ifna1","Ifna","Ifna2","Ifna3","Ifna4","Ifnab","Ifnb1","Ifne","Ifnk"))
