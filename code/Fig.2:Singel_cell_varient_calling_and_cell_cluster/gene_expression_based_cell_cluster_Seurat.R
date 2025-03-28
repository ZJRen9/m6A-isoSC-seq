setwd("/home/ZJRen/Single_cell_m6A/result/illumina_10X/Seurat")
#.libPaths(c("/home/ZJRen/R/Single-Cell/Seurat_3.1.0","/usr/lib64/R/library","/usr/share/R/library"))
.libPaths(c("/home/ZJRen/R/Single_Cell_2/Seurat_3.1.0","/usr/bin/R-4.1.3/lib64/R/library"))

library(rlang)
####################
library(Seurat)
library(ggplot2)
library(data.table)
library(dplyr)
library(patchwork)
#set.seed(123)

YTH.obj <- CreateSeuratObject(counts = Read10X("./YTH_fusion_CBE_matrix"),project = "YTH_fusion_CBE",min.cells = 3,min.features = 200)

YTH.obj[["percent.mt"]] <- PercentageFeatureSet(YTH.obj, pattern = "^MT-")

#VlnPlot(YTH.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#plot1 <- FeatureScatter(YTH.obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
#plot2 <- FeatureScatter(YTH.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#plot1 + plot2

median(YTH.obj@meta.data$nFeature_RNA)

YTH.obj <- subset(YTH.obj,subset = nFeature_RNA >= 200 & nFeature_RNA <= 5000 & percent.mt <= 20 & nCount_RNA < 50000)

VlnPlot(YTH.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

YTH.obj <- NormalizeData(YTH.obj,normalization.method = "LogNormalize",scale.factor = 10000)

#YTH.obj <- FindVariableFeatures(YTH.obj,selection.method = "vst",nfeatures = 5000,mean.cutoff = c(0.0125,3),dispersion.cutoff = c(0.5,Inf))
YTH.obj <- FindVariableFeatures(YTH.obj,selection.method = "vst")

top10 <- head(VariableFeatures(YTH.obj),10)
#plot1 <- VariableFeaturePlot(YTH.obj)
#plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
#plot1 + plot2

# APOA2 HepG2
# XIST Hela,BASP1
# BEX2,FAM71E1,HSPA1A HEK293T

YTH.obj <- ScaleData(YTH.obj,vars.to.regress = c("nFeature_RNA","percent.mt"))
YTH.obj <- RunPCA(YTH.obj,features = VariableFeatures(YTH.obj),npcs = 50)
DimPlot(YTH.obj, reduction = "pca")
ElbowPlot(YTH.obj)

YTH.obj <- FindNeighbors(YTH.obj,dims = c(1:10))
#YTH.obj <- FindClusters(object = YTH.obj,resolution=0.3,dims.use = 1:10)
YTH.obj <- FindClusters(object = YTH.obj,resolution=0.05,dims.use = 1:10)
YTH.obj <- RunTSNE(YTH.obj,perplexity=30)
YTH.obj <- RunUMAP(YTH.obj,dims = c(1:10))

###########################
g2m_genes <- cc.genes$g2m.genes
s_genes <- cc.genes$s.genes
g2m_genes <- CaseMatch(search = g2m_genes,match = rownames(YTH.obj))
s_genes <- CaseMatch(search = s_genes,match = rownames(YTH.obj))
YTH.obj <- CellCycleScoring(YTH.obj,g2m.features = g2m_genes,s.features = s_genes)


###########################
new.cluster.ids <- c("Hela","HepG2","HEK293T")
names(new.cluster.ids) <- levels(YTH.obj)
YTH.obj <- RenameIdents(YTH.obj,new.cluster.ids)
print(new.cluster.ids)
###########################
pdf("SeuratPlot1_1.pdf",height = 5,width = 6)
DimPlot(object = YTH.obj, reduction = "umap",label = TRUE)
dev.off()
pdf("SeuratPlot1_2.pdf",height = 5,width = 5)
FeaturePlot(YTH.obj,features = "CGA")
dev.off()
pdf("SeuratPlot1_3.pdf",height = 5,width = 5)
FeaturePlot(YTH.obj,features = "APOA2",min.cutoff = 2)
dev.off()
pdf("SeuratPlot1_4.pdf",height = 5,width = 5)
FeaturePlot(YTH.obj,features = "BEX2",max.cutoff = 3)
dev.off()
pdf("SeuratPlot1_5.pdf",height = 5,width = 5)
FeaturePlot(YTH.obj,features = "HSPA1A",min.cutoff = 0,max.cutoff = 4)
dev.off()
###########################
tiff("SeuratPlot1_1.tiff",height = 530,width = 600,res = 100,pointsize = 5)
DimPlot(object = YTH.obj, reduction = "umap",label = TRUE)
dev.off()
tiff("SeuratPlot1_2.tiff",height = 530,width = 530,res=100,units="px",pointsize = 5)
FeaturePlot(YTH.obj,features = "CGA")
dev.off()
tiff("SeuratPlot1_3.tiff",height = 530,width = 530,res=100,units="px",pointsize = 5)
FeaturePlot(YTH.obj,features = "APOA2",min.cutoff = 2)
dev.off()
tiff("SeuratPlot1_4.tiff",height = 530,width = 530,res=100,units="px",pointsize = 5)
FeaturePlot(YTH.obj,features = "BEX2",max.cutoff = 3)
dev.off()
tiff("SeuratPlot1_5.tiff",height = 530,width = 530,res=100,units="px",pointsize = 5)
FeaturePlot(YTH.obj,features = "HSPA1A",min.cutoff = 0,max.cutoff = 4)
dev.off()
###########################
YTH.obj.markers <- FindAllMarkers(YTH.obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
YTH.obj.markers %>% group_by(cluster) -> YTH.obj.markers.group_by_cluster
top10_markers_genes <- top_n(x=YTH.obj.markers.group_by_cluster,n=20,wt = avg_logFC)
pdf("SeuratPlot2_1.pdf",height = 5,width = 10)
DoHeatmap(YTH.obj, features = top10_markers_genes$gene,slot = "scale.data",label = F) #+ NoLegend()
dev.off()
write.table(top10_markers_genes,sep="\t",quote = FALSE,file = "top20_diff_exp_gene.txt",row.names = TRUE,col.names = TRUE)


cell_group_table <- data.frame(YTH.obj@meta.data)
write.table(cell_group_table,file = "cell_group_table_with_cellcycle.txt",sep = "\t",col.names = TRUE,row.names = TRUE,quote = FALSE)

cell_class_table <- data.frame(Idents(YTH.obj))
write.table(cell_class_table,file = "cell_class_matrix_with_cellcycle.txt",sep = "\t",col.names = TRUE,row.names = TRUE,quote = FALSE)

umap_table <- data.frame(YTH.obj@reductions$umap@cell.embeddings)
write.table(umap_table,file = "umap_table.txt",sep = "\t",col.names = TRUE,row.names = TRUE,quote = FALSE)

saveRDS(YTH.obj,file = "YTH_obj.rds")

