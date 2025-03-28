#############################################
## Nanopore Single Cell
#############################################

### load libraries
library(Seurat)
library(SeuratData)
library(data.table)
library(SeuratWrappers)
library(tidyverse)
library(patchwork)
library(data.table)
library(dplyr)
library(patchwork)
library(readr)
library(ggplot2)
library(RColorBrewer)
######################
### Fig3B
######################
C2T_editing_itst <- read_tsv("nh_c_to_t_snpmolinfos.txt",col_names = T,num_threads = 20) %>% mutate(alt = gsub(".*\\.\\.(\\w)", "\\1", transcriptId))
C2T_editing_itst <- C2T_editing_itst %>% dplyr::rename("barcode" = "cellBC") %>% mutate(read_name = paste0(barcode,"-",UMI,"-",nbReads),chrom=str_split_fixed(geneId,"_",n=4)[,1],snp_start = str_split_fixed(geneId,"_",n=4)[,2], snp_end = str_split_fixed(geneId,"_",n=4)[,3], strand = str_split_fixed(geneId,"_",n=4)[,4]) %>% dplyr::select("chrom","snp_start","snp_end","strand","alt","barcode","read_name")
C2T_editing_itst <- C2T_editing_itst %>% mutate(strand = ifelse(strand == "f","+","-"))
short_barcode <- read.csv("short_barcodes_cell.csv",sep = " ",col.names = c("barcode","cell_line"))

HeLa_long_methy <- C2T_editing_itst %>% filter(barcode %in% (short_barcode %>% filter(cell_line == "HeLa") %>% .$barcode)) %>% group_by(chrom,snp_start,snp_end,strand,alt) %>% summarize(counts=n()) %>% ungroup() %>% pivot_wider(names_from = alt,values_from = counts,values_fill =  0)
HepG2_long_methy <- C2T_editing_itst %>% filter(barcode %in% (short_barcode %>% filter(cell_line == "HepG2") %>% .$barcode)) %>% group_by(chrom,snp_start,snp_end,strand,alt) %>% summarize(counts=n()) %>% ungroup() %>% pivot_wider(names_from = alt,values_from = counts,values_fill =  0)
HEK293T_long_methy <- C2T_editing_itst %>% filter(barcode %in% (short_barcode %>% filter(cell_line == "HEK293T") %>% .$barcode)) %>% group_by(chrom,snp_start,snp_end,strand,alt) %>% summarize(counts=n()) %>% ungroup() %>% pivot_wider(names_from = alt,values_from = counts,values_fill =  0)

HeLa_long_methy_dt <- HeLa_long_methy %>% dplyr::select(-"A",-"G") %>% dplyr::rename("Ref" = "C","Alt" = "T") %>% mutate(ratio = Alt/(Alt + Ref)) %>% dplyr::rename("start" = "snp_start","end" = "snp_end") %>% mutate(start = as.numeric(start),end = as.numeric(end))
HepG2_long_methy_dt <- HepG2_long_methy %>% dplyr::select(-"A",-"G") %>% dplyr::rename("Ref" = "C","Alt" = "T") %>% mutate(ratio = Alt/(Alt + Ref)) %>% dplyr::rename("start" = "snp_start","end" = "snp_end") %>% mutate(start = as.numeric(start),end = as.numeric(end))
HEK293T_long_methy_dt <- HEK293T_long_methy %>% dplyr::select(-"A",-"G") %>% dplyr::rename("Ref" = "C","Alt" = "T") %>% mutate(ratio = Alt/(Alt + Ref))  %>% dplyr::rename("start" = "snp_start","end" = "snp_end") %>% mutate(start = as.numeric(start),end = as.numeric(end))

short_methy <- read.table("total_C2T_editing_with_YTH_V2.txt",header = T)
HeLa_short_methy <- short_methy %>% dplyr::select(chrom,start,end,Hela) %>% mutate(methy = as.numeric(str_split_fixed(Hela,";",n=2)[,2]), unmethy = as.numeric(str_split_fixed(Hela,";",n=2)[,1])) %>% mutate(ratio = methy/(methy+unmethy))
HepG2_short_methy <- short_methy %>% dplyr::select(chrom,start,end,HepG2) %>% mutate(methy = as.numeric(str_split_fixed(HepG2,";",n=2)[,2]), unmethy = as.numeric(str_split_fixed(HepG2,";",n=2)[,1])) %>% mutate(ratio = methy/(methy+unmethy))
HEK293T_short_methy <- short_methy %>% dplyr::select(chrom,start,end,HEK293T) %>% mutate(methy = as.numeric(str_split_fixed(HEK293T,";",n=2)[,2]), unmethy = as.numeric(str_split_fixed(HEK293T,";",n=2)[,1])) %>% mutate(ratio = methy/(methy+unmethy))

# Set color palette for 2D heatmap
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(32)
### HeLa
HeLa_merge_data <- inner_join(HeLa_long_methy_dt,HeLa_short_methy,by=c("chrom","start","end")) %>% filter(Ref + Alt >= 25, methy + unmethy >=25) %>% dplyr::select(chrom,start,end,ratio.x,ratio.y)
c2 <- cor.test(HeLa_merge_data$ratio.x, HeLa_merge_data$ratio.y)
title <- sprintf("N = %d r = %.3f pvalue = %d", nrow(HeLa_merge_data), c2$estimate, c2$p.value)
ggplot(HeLa_merge_data, aes(x= ratio.x,y= ratio.y)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_bin2d(bins=100) +  
  scale_fill_gradientn(colors=r, trans="log10") +
  xlab("Nanopore methylation level") +
  ylab("Bulk-eDART-seq methylation level") +
  theme(text = element_text(size = 15),
        plot.title = element_text(size = 20, vjust = 0.5, hjust = 0.5),
        legend.position = c(0.9,0.25)) +
  xlim(0,1) + ylim(0,1) +
  annotate(geom ="text",x=0.4,y=0.9, 
           label=title,
           vjust=-1,size = 6)+
  ggtitle("HeLa")

######################
### SFig6G
######################
### HEK293T
HEK293T_merge_data <- inner_join(HEK293T_long_methy_dt,HEK293T_short_methy,by=c("chrom","start","end")) %>% filter(Ref + Alt >= 25, methy + unmethy >= 25) %>% dplyr::select(chrom,start,end,ratio.x,ratio.y)
c2 <- cor.test(HEK293T_merge_data$ratio.x, HEK293T_merge_data$ratio.y)
title <- sprintf("N = %d r = %.3f pvalue = %d", nrow(HEK293T_merge_data), c2$estimate, c2$p.value)
ggplot(HEK293T_merge_data, aes(x= ratio.x,y= ratio.y)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_bin2d(bins=100) +
  scale_fill_gradientn(colors=r, trans="log10") +
  xlab("Nanopore methylation level") +
  ylab("Bulk-eDART-seq methylation level") +
  theme(text = element_text(size = 15),
        plot.title = element_text(size = 20, vjust = 0.5, hjust = 0.5),
        legend.position = c(0.9,0.25)) +
  xlim(0,1) + ylim(0,1) +
  annotate(geom ="text",x=0.4,y=0.9, 
           label=title,
           vjust=-1,size = 6)+
  ggtitle("HEK293T")
dev.off()
######################
### SFig6H
######################
### HepG2
HepG2_merge_data <- inner_join(HepG2_long_methy_dt,HepG2_short_methy,by=c("chrom","start","end")) %>% filter(Ref + Alt >= 25, methy + unmethy >=25) %>% dplyr::select(chrom,start,end,ratio.x,ratio.y)
c2 <- cor.test(HepG2_merge_data$ratio.x, HepG2_merge_data$ratio.y)
title <- sprintf("N = %d r = %.3f pvalue = %d", nrow(HepG2_merge_data), c2$estimate, c2$p.value)
ggplot(HepG2_merge_data, aes(x= ratio.x,y= ratio.y)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_bin2d(bins=100) +
  scale_fill_gradientn(colors=r, trans="log10") +
  xlab("Nanopore methylation level") +
  ylab("Bulk-eDART-seq methylation level") +
  theme(text = element_text(size = 15),
        plot.title = element_text(size = 20, vjust = 0.5, hjust = 0.5),
        legend.position = c(0.9,0.25)) +
  xlim(0,1) + ylim(0,1) +
  annotate(geom ="text",x=0.4,y=0.9, 
           label=title,
           vjust=-1,size = 6)+
  ggtitle("HepG2")
dev.off()
######################
### Fig3F;Fig3G
######################
barcode = read_csv("filter_barcodes_cells_750_filter.csv",col_names = F)
m6aLongTrans.markers = read_tsv("m6aLongTrans.markers.tsv",col_names = T)
trans_ID_name <- read_tsv("genes.t2g.mapping",col_names = c("transcriptId","transcriptName","geneID","geneName"))
m6aLongTrans.markers <- m6aLongTrans.markers %>% merge(.,trans_ID_name[,c(1,2)],by.x = "gene",by.y = "transcriptName")

c2t_long_dt = fread("st_combined_ratio.mtx") 
c2t_long_dt <- c2t_long_dt %>% filter(!(transId %in% m6aLongTrans.markers$transcriptId))
rowname = c2t_long_dt$transId
c2t_long_dt = as.matrix(c2t_long_dt[,-1])
rownames(c2t_long_dt) = rowname
c2t_long_dt <- c2t_long_dt[,colnames(c2t_long_dt) %in% barcode$X1]
#########################################
set.seed(12345)
c2t_Long_Ratio = CreateSeuratObject(counts = c2t_long_dt , project = "c2t_Long_Ratio",min.cells = 10)
c2t_Long_Ratio <- subset(c2t_Long_Ratio,subset = nFeature_RNA >= 100)
c2t_Long_Ratio <- NormalizeData(c2t_Long_Ratio)
c2t_Long_Ratio <- RunALRA(c2t_Long_Ratio) 
c2t_Long_Ratio <- FindVariableFeatures(c2t_Long_Ratio, selection.method = "vst", nfeatures = 3000) 
c2t_Long_Ratio <- ScaleData(c2t_Long_Ratio, features = rowname)
c2t_Long_Ratio <- RunPCA(c2t_Long_Ratio)
c2t_Long_Ratio <- FindNeighbors(c2t_Long_Ratio, dims = 1:10) 
c2t_Long_Ratio <- FindClusters(object = c2t_Long_Ratio,resolution=0.08,dims.use = 1:10) # 0.08
c2t_Long_Ratio <- RunUMAP(c2t_Long_Ratio, dims = 1:10, graph.name = "RNA_snn")

DimPlot(object = c2t_Long_Ratio, reduction = "umap",label = TRUE)
c2t.Long.markers <- FindAllMarkers(c2t_Long_Ratio, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
c2t.Long.markers %>% group_by(cluster) -> c2t_Long_Ratio.group_by_cluster
top20_markers_genes <- top_n(x=c2t_Long_Ratio.group_by_cluster,n=20,wt = avg_log2FC)

write.table(top20_markers_genes, "c2t_site_ratio_no_diff_Transcript.Top20markers.tsv", sep = "\t", quote = F, col.names = T, row.names = F)
DoHeatmap(c2t_Long_Ratio, features = top20_markers_genes$gene) + scale_fill_gradientn(colors = c("DarkBlue","Blue","white","Yellow","Yellow3"))

######################
### Fig3H
######################
short_barcode <- read.csv("short_barcodes_cell.csv",sep = " ",col.names = c("barcode","cell_line"))
long_barcode <- read.csv("long_barcodes_cell.csv",sep = " ",col.names = c("barcode","cell_line"))

c2t_cluster <- as.data.frame(c2t_Long_Ratio$seurat_clusters)
colnames(c2t_cluster) <- c("type")
barcode = rownames(c2t_cluster)
c2t_cluster <- data.frame(barcode,c2t_cluster)
c2t_cluster_0 <- c2t_cluster %>% filter(type==0)
cluster0 <- merge(c2t_cluster_0,short_barcode,by = "barcode")
cluster0 <- cluster0 %>% group_by(cell_line) %>% summarize(n=n()) %>% arrange(n) %>% mutate (cluster = "cluster0")
c2t_cluster_1 <- c2t_cluster %>% filter(type==1)
cluster1 <- merge(c2t_cluster_1,short_barcode,by = "barcode")
cluster1 <- cluster1 %>% group_by(cell_line) %>% summarize(n=n()) %>% arrange(n) %>% mutate (cluster = "cluster1")
c2t_cluster_2 <- c2t_cluster %>% filter(type==2)
cluster2 <- merge(c2t_cluster_2,short_barcode,by = "barcode")
cluster2 <- cluster2 %>% group_by(cell_line) %>% summarize(n=n()) %>% arrange(n) %>% mutate (cluster = "cluster2")

cluster_merge <- rbind(cluster0,cluster1,cluster2)
cluster_merge$cluster <-factor(cluster_merge$cluster,
                               levels = c("cluster0","cluster1","cluster2"))
cluster_merge$cell_line <-factor(cluster_merge$cell_line,
                                 levels = c("HeLa","HepG2","HEK293T"))

cluster_merge <- cluster_merge %>% group_by(cluster) %>% mutate(new_col=cumsum(n))
cluster_percent <- cluster_merge %>% group_by(cluster) %>% mutate(percent = n/max(new_col)) %>% filter(n==max(n))

ggplot(cluster_merge,aes(x=cluster,y=n,fill=cell_line))+
  geom_bar(stat="identity",
           position = "stack")+
  scale_fill_manual(values = c("HeLa"="#e99e9c",
                               "HEK293T"="#459943","HepG2"="#88c4e8"))+
  geom_text(data = cluster_percent, aes(x=cluster,y=new_col+100,label=n),
            vjust=1)+
  geom_text(data = cluster_percent, aes(x=cluster,y=new_col+50,label=round(percent,2)),
            vjust=1)+
  labs(x=NULL,y="Barcode Counts")+
  theme_classic()+
  theme(legend.position = c(0.8,0.8),
        plot.title = element_text(size = 20, vjust = 0.5, hjust = 0.5),
        legend.title = element_blank())+
  theme(text = element_text(size = 15)) +
  ggtitle("c2t_site_trans_ratio")

######################
### SFig6A
######################
m6a_long_dt = fread("nh_genematrix.txt")
mt.genes <- m6a_long_dt[grep("^MT-",m6a_long_dt$geneId)]$geneId #remove 13 MT gene
m6a_long_dt <- m6a_long_dt[!(m6a_long_dt$geneId %in% mt.genes)]

rowname = m6a_long_dt$geneId  # name rows by the gene symbol
#convert m6a_long_dt to matrix as seurat accept matrix rather data.frame as input
m6a_long_dt = as.matrix(m6a_long_dt[,-1]) # remove the firt column, which is the gene ID
m6a_long_dt <- m6a_long_dt[,colSums(m6a_long_dt[1:21657,])>750]
rownames(m6a_long_dt) = rowname
table(colSums(m6a_long_dt) >= 100000)
# do seurat analis
set.seed(12345)

m6aLongGene = CreateSeuratObject(counts = m6a_long_dt , project = "m6aLongGene")
m6aLongGene <- NormalizeData(m6aLongGene)
m6aLongGene <- FindVariableFeatures(m6aLongGene, selection.method = "vst", nfeatures = 3000)
m6aLongGene <- ScaleData(m6aLongGene, features = rowname)
m6aLongGene <- RunPCA(m6aLongGene)
m6aLongGene <- FindNeighbors(m6aLongGene, dims = 1:10)
m6aLongGene <- FindClusters(m6aLongGene, resolution = 0.017, graph.name = "RNA_snn")
m6aLongGene <- RunUMAP(m6aLongGene, dims = 1:10, graph.name = "RNA_snn")
DimPlot(object = m6aLongGene, reduction = "umap",label = TRUE)

######################
### SFig6B
######################
m6aLongGene.markers <- FindAllMarkers(m6aLongGene, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,min.diff.pct = 0.2)
m6aLongGene.markers <- FindAllMarkers(m6aLongGene, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
m6aLongGene.markers %>% group_by(cluster) -> m6aLongGene.markers.group_by_cluster

top10_markers_genes <- top_n(x=m6aLongGene.markers.group_by_cluster,n = 20, wt = avg_log2FC)
DoHeatmap(m6aLongGene, features = top10_markers_genes$gene) + scale_fill_gradientn(colors = c("DarkBlue","Blue","white","Yellow","Yellow3"))

######################
### SFig6C
######################
short_barcode <- read.csv("short_barcodes_cell.csv",sep = " ",col.names = c("barcode","cell_line"))

m6a_cluster <- as.data.frame(m6aLongGene$seurat_clusters)
colnames(m6a_cluster) <- c("type")
barcode = rownames(m6a_cluster)
m6a_cluster <- data.frame(barcode,m6a_cluster)
m6a_cluster_0 <- m6a_cluster %>% filter(type==0)
cluster0 <- merge(m6a_cluster_0,short_barcode,by = "barcode")
cluster0 <- cluster0 %>% group_by(cell_line) %>% summarize(n=n()) %>% arrange(n) %>% mutate (cluster = "cluster0")
m6a_cluster_1 <- m6a_cluster %>% filter(type==1)
cluster1 <- merge(m6a_cluster_1,short_barcode,by = "barcode")
cluster1 <- cluster1 %>% group_by(cell_line) %>% summarize(n=n()) %>% arrange(n) %>% mutate (cluster = "cluster1")
m6a_cluster_2 <- m6a_cluster %>% filter(type==2)
cluster2 <- merge(m6a_cluster_2,short_barcode,by = "barcode")
cluster2 <- cluster2 %>% group_by(cell_line) %>% summarize(n=n()) %>% arrange(n) %>% mutate (cluster = "cluster2")
m6a_cluster_3 <- m6a_cluster %>% filter(type==3)
cluster3 <- merge(m6a_cluster_3,short_barcode,by = "barcode")
cluster3 <- cluster3 %>% group_by(cell_line) %>% summarize(n=n()) %>% arrange(n) %>% mutate (cluster = "cluster3")
cluster_merge <- rbind(cluster0,cluster1,cluster2)

cluster_merge$cluster <-factor(cluster_merge$cluster,
                               levels = c("cluster0","cluster1","cluster2"))

cluster_merge$cell_line <-factor(cluster_merge$cell_line,
                                 levels = c("HeLa","HepG2","HEK293T"))

cluster_merge <- cluster_merge %>% 
  group_by(cluster) %>% 
  mutate(new_col=cumsum(n))

cluster_percent <- cluster_merge %>% group_by(cluster) %>% mutate(percent = n/max(new_col)) %>% filter(n==max(n))

ggplot(cluster_merge,aes(x=cluster,y=n,fill=cell_line))+
  geom_bar(stat="identity",
           position = "stack")+
  scale_fill_manual(values = c("HeLa"="#e99e9c",
                               "HEK293T"="#88c4e8","HepG2"="#459943"))+
  geom_text(data = cluster_percent, aes(x=cluster,y=new_col+100,label=n),
            vjust=1)+
  geom_text(data = cluster_percent, aes(x=cluster,y=new_col+50,label=round(percent,2)),
            vjust=1)+
  theme_classic()+
  #scale_y_continuous(expand = expansion(mult = c(0,0)))+
  theme(legend.position = c(0.8,0.8),
        legend.title = element_blank())+
  theme(text = element_text(size = 15)) +
  labs(x=NULL,y="Barcode Counts")+
  ggtitle("nh_genematrix")

######################
### SFig6D
######################
barcode = read_csv("filter_barcodes_cells_750_filter.csv",col_names = F)
m6a_long_dt = fread("/home/hjliang/m6a_sc/nanohunter/nh_isomatrix_transcript.txt") ###filter

rowname = m6a_long_dt$transcriptId  # name rows by the gene symbol
#convert m6a_long_dt to matrix as seurat accept matrix rather data.frame as input
m6a_long_dt = as.matrix(m6a_long_dt[,-1]) # remove the firt column, which is the gene ID
rownames(m6a_long_dt) = rowname
m6a_long_dt <- m6a_long_dt[,colnames(m6a_long_dt) %in% barcode$X1]
# do seurat analis
set.seed(12345)

m6aLongGene = CreateSeuratObject(counts = m6a_long_dt , project = "m6aLongGene")
m6aLongGene <- NormalizeData(m6aLongGene)
#write.table(m6aLongGene [["RNA"]]@data, "nh_isomatrix_normalize_transcript.txt", sep="\t", quote = F, col.names = T, row.names = T)
m6aLongGene <- FindVariableFeatures(m6aLongGene, selection.method = "vst", nfeatures = 3000)
m6aLongGene <- ScaleData(m6aLongGene, features = rowname)
m6aLongGene <- RunPCA(m6aLongGene)
m6aLongGene <- FindNeighbors(m6aLongGene, dims = 1:15)
m6aLongGene <- FindClusters(m6aLongGene, resolution = 0.03, graph.name = "RNA_snn",dims.use = 1:15)

m6aLongGene <- RunUMAP(m6aLongGene, dims = 1:15, graph.name = "RNA_snn")
DimPlot(object = m6aLongGene, reduction = "umap",label = TRUE)

######################
### SFig6E
######################
m6aLongGene.markers <- FindAllMarkers(m6aLongGene, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,min.diff.pct = 0.15)
m6aLongGene.markers <- FindAllMarkers(m6aLongGene, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
m6aLongGene.markers %>% group_by(cluster) -> m6aLongGene.markers.group_by_cluster

top10_markers_genes <- top_n(x=m6aLongGene.markers.group_by_cluster,n = 20, wt = avg_log2FC)
DoHeatmap(m6aLongGene, features = top10_markers_genes$gene) + scale_fill_gradientn(colors = c("DarkBlue","Blue","white","Yellow","Yellow3"))

######################
### SFig6F
######################
m6a_cluster <- as.data.frame(m6aLongGene$seurat_clusters)
colnames(m6a_cluster) <- c("type")
barcode = rownames(m6a_cluster)
m6a_cluster <- data.frame(barcode,m6a_cluster)
m6a_cluster_0 <- m6a_cluster %>% filter(type==0)
cluster0 <- merge(m6a_cluster_0,short_barcode,by = "barcode")
cluster_HeLa <- cluster0 %>% filter(cell_line=="HeLa")
cluster0 <- cluster0 %>% group_by(cell_line) %>% summarize(n=n()) %>% arrange(n) %>% mutate (cluster = "cluster0")
m6a_cluster_1 <- m6a_cluster %>% filter(type==1)
cluster1 <- merge(m6a_cluster_1,short_barcode,by = "barcode")
cluster_HepG2 <- cluster1 %>% filter(cell_line=="HepG2")
cluster1 <- cluster1 %>% group_by(cell_line) %>% summarize(n=n()) %>% arrange(n) %>% mutate (cluster = "cluster1")
m6a_cluster_2 <- m6a_cluster %>% filter(type==2)
cluster2 <- merge(m6a_cluster_2,short_barcode,by = "barcode")
cluster_HEK293T <- cluster2 %>% filter(cell_line=="HEK293T")
cluster2 <- cluster2 %>% group_by(cell_line) %>% summarize(n=n()) %>% arrange(n) %>% mutate (cluster = "cluster2")
m6a_cluster_3 <- m6a_cluster %>% filter(type==3)
cluster3 <- merge(m6a_cluster_3,short_barcode,by = "barcode")
cluster3 <- cluster3 %>% group_by(cell_line) %>% summarize(n=n()) %>% arrange(n) %>% mutate (cluster = "cluster3")

cluster_merge <- rbind(cluster0,cluster1,cluster2)

cluster_merge$cluster <-factor(cluster_merge$cluster,
                               levels = c("cluster0","cluster1","cluster2"))

cluster_merge$cell_line <-factor(cluster_merge$cell_line,
                                 levels = c("HeLa","HepG2","HEK293T"))
cluster_merge <- cluster_merge %>% 
  group_by(cluster) %>% 
  mutate(new_col=cumsum(n))

cluster_percent <- cluster_merge %>% group_by(cluster) %>% mutate(percent = n/max(new_col)) %>% filter(n==max(n))

ggplot(cluster_merge,aes(x=cluster,y=n,fill=cell_line))+
  geom_bar(stat="identity",
           position = "stack")+
  scale_fill_manual(values = c("HeLa"="#e99e9c",
                               "HEK293T"="#88c4e8","HepG2"="#459943"))+
  geom_text(data = cluster_percent, aes(x=cluster,y=new_col+100,label=n),
            vjust=1)+
  geom_text(data = cluster_percent, aes(x=cluster,y=new_col+50,label=round(percent,2)),
            vjust=1)+
  labs(x=NULL,y="Barcode Counts")+
  theme_classic()+
  #scale_y_continuous(expand = expansion(mult = c(0,0)))+
  theme(legend.position = c(0.8,0.8),
        legend.title = element_blank())+
  theme(text = element_text(size = 15)) +
  ggtitle("nh_isomatrix")

######################
### Fig3I
######################
c2t_Long_Ratio@meta.data$barcode <- rownames(c2t_Long_Ratio@meta.data)
HEK293T_barcode <- c2t_Long_Ratio@meta.data %>% filter(seurat_clusters == "1") %>% .$barcode
##################
### m6A
##################
c2t_long_dt = fread("st_combined_ratio.mtx") 
c2t_long_dt <- c2t_long_dt %>% filter(!(transId %in% m6aLongTrans.markers$transcriptId))
rowname = c2t_long_dt$transId
c2t_long_dt = as.matrix(c2t_long_dt[,-1])
rownames(c2t_long_dt) = rowname
c2t_long_dt <- c2t_long_dt[,colnames(c2t_long_dt) %in% barcode$X1]
c2t_long_sub_dt <- c2t_long_dt[,colnames(c2t_long_dt) %in% HEK293T_barcode]
#########################################
set.seed(12345)
c2t_Long_SubCluster_Ratio = CreateSeuratObject(counts = c2t_long_sub_dt , project = "c2t_Long_SubCluster_Ratio",min.cells = 10)
c2t_Long_SubCluster_Ratio <- NormalizeData(c2t_Long_SubCluster_Ratio)
c2t_Long_SubCluster_Ratio <- FindVariableFeatures(c2t_Long_SubCluster_Ratio, selection.method = "vst", nfeatures = 3000) 
c2t_Long_SubCluster_Ratio <- ScaleData(c2t_Long_SubCluster_Ratio, features = rowname)
c2t_Long_SubCluster_Ratio <- RunPCA(c2t_Long_SubCluster_Ratio)
c2t_Long_SubCluster_Ratio <- FindNeighbors(c2t_Long_SubCluster_Ratio, dims = 1:10) 
c2t_Long_SubCluster_Ratio <- FindClusters(object = c2t_Long_SubCluster_Ratio,resolution=0.3,dims.use = 1:10)
c2t_Long_SubCluster_Ratio <- RunUMAP(c2t_Long_SubCluster_Ratio, dims = 1:10, graph.name = "RNA_snn")
SubCluster_cell_class_table <- data.frame(Idents(c2t_Long_SubCluster_Ratio)) %>% dplyr::rename(Cluster = Idents.c2t_Long_SubCluster_Ratio.)
write.table(SubCluster_cell_class_table,file = "HEK293T_subcluster_c2t_site_ratio_cell_class_matrix.txt",sep = "\t",col.names = TRUE,row.names = TRUE,quote = FALSE)

use_cols <-  c("0" = "#FEB24C", #"#fd8d3C",
               "1" =  "#FA9FB5",
               "2" = "#807dba"
)
DimPlot(object = c2t_Long_SubCluster_Ratio, reduction = "umap",label = TRUE, cols = use_cols)

######################
### Fig3J
######################
c2t.Long.SubCluster.markers <- FindAllMarkers(c2t_Long_SubCluster_Ratio, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
c2t.Long.SubCluster.markers %>% group_by(cluster) -> c2t_Long_SubCluster_Ratio.group_by_cluster
c2t_top20_markers_genes <- top_n(x=c2t_Long_SubCluster_Ratio.group_by_cluster,n=20,wt = avg_log2FC)

#### cluster2
c2t_cluster2_gene <- c2t_Long_SubCluster_Ratio.group_by_cluster %>% filter(cluster == 2) %>% inner_join(.,transID2transname,by= c("gene" = "transID")) %>% dplyr::rename(transcript_id = gene,transcript_name = transname, gene_id = geneID, gene_name = genename) %>% ungroup()
c2t_cluster2_genename <- c2t_cluster2_gene$gene_name
library(clusterProfiler)
c2t_cluster2_genename <- bitr(c2t_cluster2_genename, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db"); head(c2t_cluster2_genename)
### GO
library(org.Hs.eg.db)
c2t_cluster2_GO <- enrichGO(
  gene  =  unique(c2t_cluster2_genename$ENTREZID), 
  keyType = "ENTREZID", 
  OrgDb   = org.Hs.eg.db,
  ont     = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE)
write_tsv(c2t_cluster2_GO@result %>% as.data.frame(),"HEK293T_c2t_subcluster2_GO_term.tsv",col_names = T)
library(topGO);library(ggnewscale);library(enrichplot)
ggplot(c2t_cluster2_GO@result %>% head(n=10),aes(x=-log10(p.adjust),y=reorder(str_wrap(Description,50),-log10(p.adjust))))+
  geom_col(aes(fill=p.adjust),width = 0.7) +
  scale_fill_material(name="p.adjust")+
  scale_fill_gradientn(colours = rev(cols),
                       name="FDR") +
  labs(x="-Log10(FDR)",y=NULL,title = 'HEK293T c2t subcluster2 GO BP')+
  theme_bw()+
  scale_x_continuous(expand = c(0.01,0.01))+
  theme(plot.title = element_text(size = 15, vjust = 0.5, hjust = 0.5),
        panel.grid = element_blank(),
        text = element_text(size = 15))

