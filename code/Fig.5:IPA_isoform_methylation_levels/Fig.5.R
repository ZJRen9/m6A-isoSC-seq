######################
###  Fig5J
######################
c2t_long_dt = fread("st_combined_ratio.mtx") 
c2t_long_dt_sub <- c2t_long_dt %>% filter(transId %in% c("ENST00000337432","ENST00000421782")) ## RAD51C
rowname = c2t_long_dt_sub$transId
c2t_long_mtx_sub = as.matrix(c2t_long_dt_sub[,-1])
rownames(c2t_long_mtx_sub) = rowname
c2t_long_mtx_sub <- c2t_long_mtx_sub[,colnames(c2t_long_mtx_sub) %in% barcode$X1]
sums <- colSums(c2t_long_mtx_sub)
c2t_long_mtx_sub <- c2t_long_mtx_sub[, sums != 0]

c2t_Long_Ratio_subset <- c2t_Long_Ratio #### subset ratio
c2t_Long_Ratio_subset@reductions$umap@cell.embeddings <- c2t_Long_Ratio_subset@reductions$umap@cell.embeddings[rownames(c2t_Long_Ratio_subset@reductions$umap@cell.embeddings) %in% colnames(c2t_long_mtx_sub),] 
FeaturePlot(c2t_Long_Ratio_subset,features = c("ENST00000337432","ENST00000421782"),reduction = "umap",keep.scale = "all" , slot = "data",order = T, pt.size = 3) ### RAD51C
FeaturePlot(c2t_Long_Ratio_expr,features = c("ENST00000337432","ENST00000421782"),reduction = "umap",keep.scale = "all" , order = T) ### RAD51C

######################
###  SFig7A
######################
YTH_APOBEC_readid <- read_tsv("all_cons.YTHfusion.readid",col_names = c("read_name")) %>% mutate(barcode = str_split_fixed(read_name,"-",3)[,1]) %>% mutate(gene_name = "APOBEC")
YTH_APOBEC_dt <-  YTH_APOBEC_readid %>% group_by(barcode) %>% summarize(counts = n()) %>% ungroup()
YTH_APOBEC_wider_dt <- YTH_APOBEC_dt %>% filter(barcode %in% colnames(c2t_Long_SubCluster_Ratio@assays$RNA@counts)) %>% 
  pivot_wider(names_from = barcode, values_from = counts, values_fill = list(counts = 0))
YTH_APOBEC_mtx <- YTH_APOBEC_wider_dt %>% as.matrix()
rownames(YTH_APOBEC_mtx) <- "APOBEC"

c2t_Long_SubCluster_APOBEC_expr <- c2t_Long_SubCluster_Ratio ### expression
c2t_Long_SubCluster_APOBEC_expr@assays$RNA@scale.data = YTH_APOBEC_mtx
c2t_Long_SubCluster_APOBEC_expr@assays$RNA@data = YTH_APOBEC_mtx

### Violin plot
APOBEC_count_Seurat_0 <- as.vector(subset(c2t_Long_SubCluster_APOBEC_expr, subset = APOBEC > 0, idents = "0")@assays$RNA@data)
APOBEC_count_Seurat_1 <- as.vector(subset(c2t_Long_SubCluster_APOBEC_expr, subset = APOBEC > 0, idents = "1")@assays$RNA@data)
APOBEC_count_Seurat_2 <- as.vector(subset(c2t_Long_SubCluster_APOBEC_expr, subset = APOBEC > 0, idents = "2")@assays$RNA@data)
APOBEC_count_Seurat_cells <- rbind(tibble(gene_name = "APOBEC", counts = APOBEC_count_Seurat_0,group1 = "0"),
                                   tibble(gene_name = "APOBEC", counts = APOBEC_count_Seurat_1,group1 = "1"),
                                   tibble(gene_name = "APOBEC", counts = APOBEC_count_Seurat_2,group1 = "2"))
APOBEC_count_Seurat_cells_dt <- APOBEC_count_Seurat_cells %>% group_by(group1) %>% mutate(High = quantile(counts, 0.75),Med = median(counts),Low = quantile(counts, 0.25)) %>% ungroup() 

compare_means(ratio  ~ group1,subset_tibble_cells_dt,method = "wilcox.test",paired = F) ###  0.25
ggplot(APOBEC_count_Seurat_cells_dt,aes(x=group1,y= counts,fill = group1))+
  geom_violin(trim = FALSE)+
  geom_linerange(aes(x = group1, ymin = Low, ymax = High),
                 position = position_dodge(width = 0.9))+
  geom_point(aes(x = group1, y = Med,),
             position = position_dodge(width = 0.9), size = 3)+
  scale_fill_manual(values = c("0" = "#FEB24C",#"#fd8d3C",
                               "1" =  "#FA9FB5",
                               "2" = "#807dba"))+
  theme_classic()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = 20, vjust = 0.5, hjust = 0.5),
        legend.position = "none")+  
  geom_signif(comparisons = list(c("0", "1"),c("1","2"),c("0","2")),test = "wilcox.test",
              map_signif_level = function(p) sprintf("P = %.2g",p),
              textsize = 5, 
              y_position = c(6,6.5,7))+
  theme(text = element_text(size = 15)) +
  labs(y=TeX(r"(YTH-APOBEC counts)"),x = "HEK293T Subcluster")+
  scale_y_continuous(limits=c(0,7.5)) +
  ggtitle("APOBEC")

######################
###  SFig7B
######################
c2t_Long_Ratio@meta.data$barcode <- rownames(c2t_Long_Ratio@meta.data)
HEK293T_barcode <- c2t_Long_Ratio@meta.data %>% filter(seurat_clusters == "1") %>% .$barcode ### 1
##################
### expression
##################
m6a_long_dt = fread("nh_isomatrix_transcript.txt")
m6a_long_dt <- m6a_long_dt %>% filter(!(transcriptId %in% m6aLongTrans.markers$transcriptId))
rowname = m6a_long_dt$transcriptId
m6a_long_dt = as.matrix(m6a_long_dt[,-1])
rownames(m6a_long_dt) = rowname
m6a_long_dt <- m6a_long_dt[,colnames(m6a_long_dt) %in% barcode$X1]
m6a_long_dt <- m6a_long_dt[,colnames(m6a_long_dt) %in% HEK293T_barcode]

set.seed(12345)
m6a_long_SubCluster_Expr = CreateSeuratObject(counts = m6a_long_dt , project = "m6a_long_SubCluster_Expr",min.cells = 10)
m6a_long_SubCluster_Expr <- subset(m6a_long_SubCluster_Expr,subset = nFeature_RNA >= 100)
m6a_long_SubCluster_Expr <- NormalizeData(m6a_long_SubCluster_Expr)
m6a_long_SubCluster_Expr <- FindVariableFeatures(m6a_long_SubCluster_Expr, selection.method = "vst", nfeatures = 3000)
m6a_long_SubCluster_Expr <- ScaleData(m6a_long_SubCluster_Expr, features = rowname)
m6a_long_SubCluster_Expr <- RunPCA(m6a_long_SubCluster_Expr)
m6a_long_SubCluster_Expr <- FindNeighbors(m6a_long_SubCluster_Expr, dims = 1:10)
# 0.02 -> 3 clus # 0.04 750 filter
m6a_long_SubCluster_Expr <- FindClusters(m6a_long_SubCluster_Expr, resolution = 0.03, graph.name = "RNA_snn",dims.use = 1:10)
m6a_long_SubCluster_Expr <- RunUMAP(m6a_long_SubCluster_Expr, dims = 1:10, graph.name = "RNA_snn") # 
use_cols <-  c("0" = "#FEB24C",#"#fd8d3C",
               "1" =  "#FA9FB5",
               "2" = "#807dba"
)
DimPlot(object = m6a_long_SubCluster_Expr, reduction = "umap",label = TRUE, cols = use_cols)

######################
###  SFig7C
######################
## barplot
HEK293T_m6A_cluster <- as.data.frame(c2t_Long_SubCluster_Ratio$seurat_clusters)
colnames(HEK293T_m6A_cluster) <- c("type")
barcode = rownames(HEK293T_m6A_cluster)
HEK293T_m6A_cluster <- data.frame(barcode,HEK293T_m6A_cluster)

HEK293T_Expr_cluster <- as.data.frame(m6a_long_SubCluster_Expr$seurat_clusters)
colnames(HEK293T_Expr_cluster) <- c("type")
barcode = rownames(HEK293T_Expr_cluster)
HEK293T_Expr_cluster <- data.frame(barcode,HEK293T_Expr_cluster)

HEK293T_m6A_Expr_cluster <- inner_join(HEK293T_m6A_cluster,HEK293T_Expr_cluster,by=c("barcode")) %>% dplyr::rename(m6A_cluster = type.x, expr_cluster = type.y)

HEK293T_m6A_Expr_cluster_0 <- HEK293T_m6A_Expr_cluster %>% filter(m6A_cluster==0) %>% group_by(expr_cluster) %>% summarize(n=n()) %>% arrange(n) %>% mutate (cluster = "cluster0")
HEK293T_m6A_Expr_cluster_1 <- HEK293T_m6A_Expr_cluster %>% filter(m6A_cluster==1) %>% group_by(expr_cluster) %>% summarize(n=n()) %>% arrange(n) %>% mutate (cluster = "cluster1")
HEK293T_m6A_Expr_cluster_2 <- HEK293T_m6A_Expr_cluster %>% filter(m6A_cluster==2) %>% group_by(expr_cluster) %>% summarize(n=n()) %>% arrange(n) %>% mutate (cluster = "cluster2")

cluster_merge <- rbind(HEK293T_m6A_Expr_cluster_0,HEK293T_m6A_Expr_cluster_1,HEK293T_m6A_Expr_cluster_2)

cluster_merge$cluster <-factor(cluster_merge$cluster,
                               levels = c("cluster0","cluster1","cluster2"))

cluster_merge$expr_cluster <-factor(cluster_merge$expr_cluster,
                                    levels = c("0","1","2"))

cluster_merge <- cluster_merge %>% 
  group_by(cluster) %>% 
  mutate(new_col=cumsum(n))

cluster_percent <- cluster_merge %>% group_by(cluster) %>% mutate(percent = n/max(new_col)) %>% filter(n==max(n))
use_cols <-  c("0" = "#FEB24C",#"#fd8d3C",
               "1" =  "#FA9FB5",
               "2" = "#807dba"
)

ggplot(cluster_merge,aes(x=cluster,y=n,fill=expr_cluster))+
  geom_bar(stat="identity",
           position = "stack")+
  scale_fill_manual(values = c("0"="#FEB24C",
                               "1"="#FA9FB5","2"="#807dba"))+
  geom_text(data = cluster_percent, aes(x=cluster,y=new_col+25,label=n),
            vjust=1, size = 5)+
  geom_text(data = cluster_percent, aes(x=cluster,y=new_col,label=round(percent,2)),
            vjust=1, size = 5)+
  labs(x="m6A cluster",y="Number of cells",fill = "Expression cluster")+
  theme_classic()+
  scale_y_continuous(expand = expansion(mult = c(0.02,0.02)))+
  theme(legend.position = c(0.8,0.8),
        plot.title = element_text(size = 20, vjust = 0.5, hjust = 0.5))+
  theme(text = element_text(size = 15)) +
  ggtitle("HEK293T subcluster")

########################################
###  Fig4K;Fig5I;Fig6G;SFig6P;SFig9h
########################################
### single read plot
source("methylation_R_utils.R")
library(scales)
library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(data.table)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
devtools::install_github("dzhang32/ggtranscript")
library(ggtranscript)
db <- TxDb.Hsapiens.UCSC.hg38.knownGene

gtf <- as.data.frame(rtracklayer::import("genes.gtf"))
gtf <- gtf %>% dplyr::select("seqnames","start","end","strand","type","gene_name","gene_id","gene_type","transcript_id","transcript_name","transcript_type")
Transctipt_Name <- c("ENST00000395119","ENST00000528610")
centers <- gtf %>% filter(transcript_id %in% Transctipt_Name) %>% filter(type == "transcript") %>% dplyr::select(seqnames,start,end,transcript_id,transcript_name) %>% rename(seqnames = "chrom")
write_tsv(centers,"c2t_site_ratio_isoform_Cell_specific.txt",col_names = F)

C2T_editing_itst <- read_tsv("nh_c_to_t_snpmolinfos_trans.txt",col_names = c("chrom","snp_start","snp_end","transcript_id","strand","alt","read_name"),num_threads = 20)

C2T_editing_itst_specific <- C2T_editing_itst %>% filter(transcript_id %in% Transcript_Name) %>% ### Transcript_Name #HeLa_retained_protein_trans$transcript_id
  mutate(snp_end = snp_start + 1) %>%
  mutate(barcode = str_split_fixed(read_name,"-",n=3)[,1]) %>% 
  dplyr::select("chrom","snp_start","snp_end","strand","alt","read_name","barcode")
write_tsv(C2T_editing_itst_specific,"c2t_site_ratio_isoform_Cell_specific_snp.txt",col_names = F )  ### filter snp

```{Linux Bash
  python3 c2t_site_ratio_isoform_Cell_specific_read.py
} 
```

c2t_long_dt = fread("st_combined_ratio.mtx")
barcode = read_csv("filter_barcodes_cells_750_filter.csv",col_names = F)
cluster = read_tsv("c2t_site_ratio_cell_cluster_V2.tsv")
transID2transname <- read_tsv("genes.t2g.mapping",col_names = c("transID","transname","geneID","genename"))
c2t_long_dt <- inner_join(transID2transname,c2t_long_dt,by = c("transID" = "transId")) %>% dplyr::select(-c("transID","geneID","genename"))
rowname = c2t_long_dt$transname
c2t_long_dt = as.matrix(c2t_long_dt[,-1])
rownames(c2t_long_dt) = rowname
c2t_long_dt <- c2t_long_dt[,colnames(c2t_long_dt) %in% barcode$X1]

c2t_long_dt_HeLa <- as.data.frame(c2t_long_dt[,colnames(c2t_long_dt) %in% (cluster %>% filter(cell_line=="HeLa"))$barcode]) %>% mutate(transcript_name = rowname) %>% inner_join(. ,transID2transname %>% dplyr::select(transname,genename),by = c("transcript_name" = "transname")) %>% dplyr::rename("gene_name" = "genename")
c2t_long_dt_HEK293 <- as.data.frame(c2t_long_dt[,colnames(c2t_long_dt) %in% (cluster %>% filter(cell_line=="HEK293T"))$barcode]) %>% mutate(transcript_name = rowname) %>% inner_join(. ,transID2transname %>% dplyr::select(transname,genename),by = c("transcript_name" = "transname")) %>% dplyr::rename("gene_name" = "genename")
c2t_long_dt_HepG2 <- as.data.frame(c2t_long_dt[,colnames(c2t_long_dt) %in% (cluster %>% filter(cell_line=="HepG2"))$barcode]) %>% mutate(transcript_name = rowname) %>% inner_join(. ,transID2transname %>% dplyr::select(transname,genename),by = c("transcript_name" = "transname")) %>% dplyr::rename("gene_name" = "genename")

### HeLa
c2t_long_dt_HeLa_dt <- c2t_long_dt_HeLa %>% pivot_longer(cols = !c(transcript_name,gene_name), names_to = "barcode",values_to = "ratio") 
c2t_long_dt_HeLa_filter_dt <- c2t_long_dt_HeLa_dt %>% filter(ratio > 0)
### HEK293T
c2t_long_dt_HEK293_dt <- c2t_long_dt_HEK293 %>% pivot_longer(cols = !c(transcript_name,gene_name), names_to = "barcode",values_to = "ratio") 
c2t_long_dt_HEK293_filter_dt <- c2t_long_dt_HEK293_dt %>% filter(ratio > 0)
### HepG2
c2t_long_dt_HepG2_dt <- c2t_long_dt_HepG2 %>% pivot_longer(cols = !c(transcript_name,gene_name), names_to = "barcode",values_to = "ratio") 
c2t_long_dt_HepG2_filter_dt <- c2t_long_dt_HepG2_dt %>% filter(ratio > 0)

# rename as regs.sub
# make window 2kb
regs.sub <- centers %>%
  mutate(center = (start + end + 1)/2) %>% merge(.,gtf %>% filter(type=="transcript") %>% dplyr::select(transcript_id,transcript_type),by="transcript_id") %>% arrange(transcript_name) %>% inner_join(.,gtf %>% filter(type == "transcript") %>% dplyr::select("gene_name","transcript_name"))

reads <- read_tsv("c2t_site_ratio_isoform_Cell_specific_reads.tsv",col_names = c("chrom","trans_start","trans_end","strand","rname","transcript_id"),num_threads = 20)
snp <- read_tsv("c2t_site_ratio_isoform_Cell_specific_snp.txt",col_names = c("chrom","snp_start","snp_end","strand","mod","rname","barcode"))
cluster = read_tsv("c2t_site_ratio_cell_cluster_V2.tsv")

reads <- reads %>% mutate(barcode = str_split_fixed(rname,"-",n=3)[,1]) %>%
  mutate(cell_line = case_when(
    barcode %in%  (cluster %>% filter(cell_line=="HeLa"))$barcode ~  "HeLa",
    barcode %in%  (cluster %>% filter(cell_line=="HEK293T"))$barcode ~ "HEK293T",
    barcode %in%  (cluster %>% filter(cell_line=="HepG2"))$barcode ~ "HepG2"))
reads <- na.omit(reads)
reads_HeLa <- reads %>% filter(cell_line == "HeLa")
reads_HEK293T <- reads %>% filter(cell_line == "HEK293T")
reads_HepG2 <- reads %>% filter(cell_line == "HepG2")
reads_snp_HeLa <- inner_join(reads_HeLa,snp,by=c("rname","chrom","strand","barcode")) %>% filter(trans_start <= snp_start, trans_end >= snp_end)
reads_snp_HEK293T <- inner_join(reads_HEK293T,snp,by=c("rname","chrom","strand","barcode")) %>% filter(trans_start <= snp_start, trans_end >= snp_end)
reads_snp_HepG2 <- inner_join(reads_HepG2,snp,by=c("rname","chrom","strand","barcode")) %>% filter(trans_start <= snp_start, trans_end >= snp_end)
reads_snp_merge <- rbind(reads_snp_HeLa,reads_snp_HEK293T,reads_snp_HepG2)

### HeLa HEK293T HepG2
HeLa_reads.list <- list()
HEK293T_reads.list <- list()
HepG2_reads.list <- list()
st <- Sys.time()
for (i in seq(nrow(regs.sub))) { 
  print(i)
  reg <- regs.sub[i,]
  c2t_long_dt_HeLa_sub <- c2t_long_dt_HeLa_filter_dt %>% filter(transcript_name %in% reg$transcript_name)
  c2t_long_dt_HEK293_sub <- c2t_long_dt_HEK293_filter_dt %>% filter(transcript_name %in% reg$transcript_name)
  c2t_long_dt_HepG2_sub <- c2t_long_dt_HepG2_filter_dt %>% filter(transcript_name %in% reg$transcript_name)
  readcalls_HeLa <- reads_snp_HeLa %>% 
    filter(transcript_id %in% reg$transcript_id) %>% 
    filter(reg$chrom == chrom, reg$start <= trans_start, reg$end >= trans_end, trans_start<= snp_start, trans_end >= snp_end) %>%
    filter(barcode %in% c2t_long_dt_HeLa_sub$barcode) %>%  ### Use ratio != 0 barcode
    mutate(mcall = ifelse(mod == "T","1","0")) %>% distinct()
  readcalls_HEK293T <- reads_snp_HEK293T %>%  
    filter(transcript_id %in% reg$transcript_id) %>% 
    filter(reg$chrom == chrom, reg$start <= trans_start, reg$end >= trans_end, trans_start<= snp_start, trans_end >= snp_end) %>%
    filter(barcode %in% c2t_long_dt_HEK293_sub$barcode) %>%  ### Use ratio != 0 barcode 
    mutate(mcall = ifelse(mod == "T","1","0")) %>% distinct()
  readcalls_HepG2 <- reads_snp_HepG2 %>% 
    filter(transcript_id %in% reg$transcript_id) %>% 
    filter(reg$chrom == chrom, reg$start <= trans_start, reg$end >= trans_end, trans_start<= snp_start, trans_end >= snp_end) %>%
    filter(barcode %in% c2t_long_dt_HepG2_sub$barcode) %>%  ### Use ratio != 0 barcode
    mutate(mcall = ifelse(mod == "T","1","0")) %>% distinct()
  # subset reads overlapping this region 
  # modified the distance
  # add info and label, get distance
  #readcalls$distance <- readcalls$trans_start - reg$center
  HeLa_reads.list[[i]] <- readcalls_HeLa %>%
    mutate(lab = paste(rname,chrom,reg$start,sep="_"))
  HEK293T_reads.list[[i]] <- readcalls_HEK293T %>%
    mutate(lab = paste(rname,chrom,reg$start,sep="_"))
  HepG2_reads.list[[i]] <- readcalls_HepG2 %>%
    mutate(lab = paste(rname,chrom,reg$start,sep="_"))
}
Sys.time() - st

#single read plot HEK293T HeLa HepG2
pal <- pal_npg("nrc")(10)
meth_pal <- c(pal[1],pal[2])

keepi <- intersect(which(!sapply(cgdist.list,function(x){is.null(nrow(x))})),
                   which(!sapply(gcdist.list,function(x){is.null(nrow(x))})))

library(ggtranscript)
library(patchwork)
plotpath <- file.path("EEF1D_223_202_singleread.pdf")
#width=5,height=8 Enlarge
pdf(plotpath,width = 5, height = 8, useDingbats = F)
for (i in seq(nrow(regs.sub))){
  print(i)
  reg <- regs.sub[i,]
  HeLa_readdat <- HeLa_reads.list[[i]]
  HEK293T_readdat <- HEK293T_reads.list[[i]]
  HepG2_readdat <- HepG2_reads.list[[i]]
  center <- tibble(center=(reg$start+reg$end + 1)/2, color = "Center")
  # assign protein binding state
  HeLa_readruns <- getRuns_cell_snp(HeLa_readdat,maxGap = 1,pad = 5) %>%
    mutate(m = ifelse(values == 1, "Mutated","Unmutated")) 
  HEK293T_readruns <- getRuns_cell_snp(HEK293T_readdat,maxGap = 1,pad = 5) %>%
    mutate(m = ifelse(values == 1, "Mutated","Unmutated"))
  HepG2_readruns <- getRuns_cell_snp(HepG2_readdat,maxGap = 1,pad = 5) %>%
    mutate(m = ifelse(values == 1, "Mutated","Unmutated"))
  HeLa_order <- order_reads_cell_snp(HeLa_readruns,reads_HeLa %>% filter(transcript_id == reg$transcript_id))
  HEK293T_order <- order_reads_cell_snp(HEK293T_readruns,reads_HEK293T %>% filter(transcript_id == reg$transcript_id))
  HepG2_order <- order_reads_cell_snp(HepG2_readruns,reads_HepG2 %>% filter(transcript_id == reg$transcript_id))

  readruns.both <- bind_rows(HeLa_order$x,HEK293T_order$x,HepG2_order$x)
  readbounds.both <- bind_rows(HeLa_order$bounds,HEK293T_order$bounds,HepG2_order$bounds)
  #bounds.both$start <- bounds.both$start - 50 #looks more beautiful
  #bounds.both$end <- bounds.both$end + 50 #looks more beautiful
  
  # breaks
  if(reg$end - reg$start >= 60000){
    breaks <- seq(floor((reg$start)/8000)*8000,  ## 1000 2000
                  ceiling((reg$end)/8000)*8000,8000) ## 1000 2000
    breaklabs <- breaks/1000
    #if (min(breaklabs) > 1000){  # 1000 2000
    #  morethan1k <- floor(breaklabs/1000) # 100
    #  lessthan1k <- breaklabs - morethan1k*1000 #100
    #  breaklabs <- paste(morethan1k,lessthan1k,sep = ",") 
    #}
  }
  else if(reg$end - reg$start >= 30000  && reg$end - reg$start < 60000){
    breaks <- seq(floor((reg$start)/5000)*5000,  ## 1000 2000
                  ceiling((reg$end)/5000)*5000,5000) ## 1000 2000
    breaklabs <- breaks/1000
    #if (min(breaklabs) > 1000){  # 1000 2000
    #  morethan1k <- floor(breaklabs/1000) # 100
    #  lessthan1k <- breaklabs - morethan1k*1000 #100
    #  breaklabs <- paste(morethan1k,lessthan1k,sep = ",") 
    #}
  }
  else{
    breaks <- seq(floor((reg$start)/2000)*2000,  ## 1000 2000
                  ceiling((reg$end)/2000)*2000,2000) ## 1000 2000
    breaklabs <- breaks/1000
    #if (min(breaklabs) > 1000){  # 1000 2000
    #  morethan1k <- floor(breaklabs/1000) # 100
    #  lessthan1k <- breaklabs - morethan1k*1000 #100
    #  breaklabs <- paste(morethan1k,lessthan1k,sep = ",") 
    #}
  }
  
  p1 <- ggplot(readruns.both,aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax)) + 
    theme_classic() + 
    facet_grid(cell_line ~.,labeller =as_labeller(c(HeLa = "HeLa",HEK293T = "HEK293T",HepG2 = "HepG2")),scales = "free_y",space = "free") +
    geom_rect(data = readbounds.both, fill = "grey90") + 
    geom_rect(aes(fill = m))   +
    scale_fill_manual(name = "State", values = meth_pal) +
    scale_linetype_manual(name = "Binding Site", values = "dashed") +
    labs(x = paste("Coordinate on",reg$chrom[1],"(kb)"), y = "Reads") +
    theme(axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(),
          panel.spacing = unit(0.1, "cm"), ##消除分面间距
          legend.position = "bottom") +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
          panel.grid=element_blank()) +
    theme(text = element_text(size = 15),
          strip.text.y = element_text(size = 12)) + #7
    scale_x_continuous(breaks = breaks, labels = breaklabs) + 
    #xlim(98044500,98045500) ## enlarge
    xlim(143579500,143581000)
  #coord_cartesian(xlim = c(reg$start,reg$end)) ## default
  
  regs.sub_CDS <- gtf %>% filter(transcript_id == reg$transcript_id) %>% filter(type == "CDS") %>% filter(end - start >= 10)
  p2 <- ggplot(data = gtf %>% filter(transcript_id == reg$transcript_id,type == "exon"),aes(xstart = start,xend = end,y = transcript_id)) +
    theme_bw() +
    geom_range(aes(fill = transcript_type),color = NA,height = 0.25) +
    geom_range(data = regs.sub_CDS,aes(fill = transcript_type),color = NA) +
    scale_fill_manual(values = case_when(reg$transcript_type=="protein_coding" ~ "#E0272499",
                                         reg$transcript_type=="nonsense_mediated_decay" ~ "#85B22E",
                                         reg$transcript_type=="retained_intron" ~ "#7985CBFF",
                                         reg$transcript_type=="processed_transcript" ~ "#fd8d3C"), lab = reg$transcript_type) + ### "#7985CBFF" nonsense配色 #E0272499 protein coding
    geom_intron(data = to_intron(gtf %>% filter(transcript_name == reg$transcript_name,type == "exon"),"transcript_name"),aes(strand = strand)) + 
    theme(
          legend.position = c(0.7,0.9),
          legend.text.align = 0,
          legend.title = element_blank(),
          legend.background = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank()) +
    guides(color = FALSE) + # 删除color图例
    labs(x = NULL,y = NULL) +
    theme(text = element_text(size = 13)) + #15
    theme(panel.grid.major =element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(size=12,hjust=0.5)) + #去除背景
    #xlim(98044500,98045500) +
    xlim(143579500,143581000) +
    ggtitle(reg$transcript_name)
  p <- (p2/p1) + plot_layout(heights = c(1.2, 4))
  print(p)
}
dev.off()

