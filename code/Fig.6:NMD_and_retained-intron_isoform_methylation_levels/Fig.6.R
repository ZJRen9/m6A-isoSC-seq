#########################
###  Fig6A
#########################
### gggenomes
library(readxl)
library(gggenomes) #devtools::install_github("thackl/thacklr") devtools::install_github("thackl/gggenomes")
library(ggalluvial)
gtf <- as.data.frame(rtracklayer::import("genes.gtf"))
gtf <- gtf %>% dplyr::select("seqnames","start","end","strand","type","gene_name","transcript_id","transcript_name","transcript_type","gene_type")
HeLa_high_gtf <- gtf %>% 
  filter(transcript_name %in% (HeLa_total_dt %>% filter(group == "high_m6A_transcript"))$transcript_name) %>% 
  filter(type =="transcript") %>% 
  group_by(transcript_type) %>% 
  summarize(n=n()) %>% 
  mutate(group="m6A_high_transcript") %>% 
  mutate(transcript_type=fct_relevel(transcript_type,c("protein_coding","retained_intron","nonsense_mediated_decay","lncRNA"))) %>%
  arrange(transcript_type) %>% 
  mutate(cumsum = cumsum(n)) %>% 
  mutate(start = c(1,dplyr::lag(cumsum)[-1]+1),end = cumsum) %>% 
  mutate(freq = n/sum(n))
HepG2_high_gtf <- gtf %>% filter(transcript_name %in% (HepG2_ratio_expression_dt %>% filter(group == "high_m6A_transcript"))$transcript_name) %>% filter(type =="transcript") %>% group_by(transcript_type) %>% summarize(n=n()) %>% mutate(group="m6A_high_transcript") %>% mutate(transcript_type=fct_relevel(transcript_type,c("protein_coding","lncRNA","retained_intron","nonsense_mediated_decay"))) %>% arrange(transcript_type) %>% mutate(cumsum = cumsum(n)) %>% mutate(start = c(1,dplyr::lag(cumsum)[-1]+1),end = cumsum) %>% mutate(freq = n/sum(n))
HEK293_high_gtf <- gtf %>% filter(transcript_name %in% (HEK293_ratio_expression_dt %>% filter(group == "high_m6A_transcript"))$transcript_name) %>% filter(type =="transcript") %>% group_by(transcript_type) %>% summarize(n=n()) %>% mutate(group="m6A_high_transcript") %>% mutate(transcript_type=fct_relevel(transcript_type,c("protein_coding","lncRNA","retained_intron","nonsense_mediated_decay"))) %>% arrange(transcript_type) %>% mutate(cumsum = cumsum(n)) %>% mutate(start = c(1,dplyr::lag(cumsum)[-1]+1),end = cumsum) %>% mutate(freq = n/sum(n))
HeLa_low_gtf <- gtf %>% filter(transcript_name %in% (HeLa_total_dt %>% filter(group == "low_m6A_transcript"))$transcript_name) %>% 
  filter(type =="transcript") %>% group_by(transcript_type) %>% summarize(n=n()) %>% mutate(group="m6A_low_transcript") %>% 
  mutate(transcript_type=fct_relevel(transcript_type,c("protein_coding","retained_intron","nonsense_mediated_decay","lncRNA"))) %>%
  arrange(transcript_type) %>% mutate(cumsum = cumsum(n)) %>% mutate(start = c(1,dplyr::lag(cumsum)[-1]+1),end = cumsum) %>% mutate(freq = n/sum(n))
HepG2_low_gtf <- gtf %>% filter(transcript_name %in% (HepG2_ratio_expression_dt %>% filter(group == "low_m6A_transcript"))$transcript_name) %>% filter(type =="transcript") %>% group_by(transcript_type) %>% summarize(n=n()) %>% mutate(group="m6A_low_transcript") %>% mutate(transcript_type=fct_relevel(transcript_type,c("protein_coding","lncRNA","retained_intron","nonsense_mediated_decay"))) %>% arrange(transcript_type) %>% mutate(cumsum = cumsum(n)) %>% mutate(start = c(1,dplyr::lag(cumsum)[-1]+1),end = cumsum) %>% mutate(freq = n/sum(n))
HEK293_low_gtf <- gtf %>% filter(transcript_name %in% (HEK293_ratio_expression_dt %>% filter(group == "low_m6A_transcript"))$transcript_name) %>% filter(type =="transcript") %>% group_by(transcript_type) %>% summarize(n=n()) %>% mutate(group="m6A_low_transcript") %>% mutate(transcript_type=fct_relevel(transcript_type,c("protein_coding","lncRNA","retained_intron","nonsense_mediated_decay"))) %>% arrange(transcript_type) %>% mutate(cumsum = cumsum(n)) %>% mutate(start = c(1,dplyr::lag(cumsum)[-1]+1),end = cumsum) %>% mutate(freq = n/sum(n))
HEK293_nonse <- tibble(transcript_type = c("nonsense_mediated_decay"), n = c(0), group = c("m6A_low_transcript"), cumsum = c(82), start = c(82), end = c(82), freq = c(0)) # nonsensed end
HEK293_low_gtf <- rbind(HEK293_low_gtf,HEK293_nonse) %>% mutate(transcript_type=fct_relevel(transcript_type,c("protein_coding","lncRNA","retained_intron","nonsense_mediated_decay"))) %>% arrange(transcript_type)

HeLa_gtf <- rbind(HeLa_low_gtf,HeLa_high_gtf) %>% group_by(transcript_type) %>% mutate(group = factor(group,levels=c("m6A_low_transcript","m6A_high_transcript"))) %>% mutate(transcript_type=fct_relevel(transcript_type,c("nonsense_mediated_decay","retained_intron","protein_coding")))
HepG2_gtf <- rbind(HepG2_low_gtf,HepG2_high_gtf) %>% group_by(transcript_type) %>% mutate(group = factor(group,levels=c("m6A_low_transcript","m6A_high_transcript"))) %>% mutate(transcript_type=fct_relevel(transcript_type,c("nonsense_mediated_decay","retained_intron","lncRNA","protein_coding")))
HEK293_gtf <- rbind(HEK293_low_gtf,HEK293_high_gtf) %>% group_by(transcript_type) %>% mutate(group = factor(group,levels=c("m6A_low_transcript","m6A_high_transcript"))) %>% mutate(transcript_type=fct_relevel(transcript_type,c("nonsense_mediated_decay","retained_intron","lncRNA","protein_coding")))

mycolor = c("#db6968","#4d97cd","#459943","#f8984e")
chisq.test(matrix(c(329,12,7,6,315,52,90,50),nrow=4,ncol=2))$p.value ## 4.6e-23
plotpath <- file.path("HeLa_gggenomes.pdf")
ggplot(HeLa_gtf,aes(x = group, stratum = transcript_type,alluvium = transcript_type,y = freq*100,
                    fill = transcript_type, label = transcript_type)) +
  scale_x_discrete(expand = c(0.2, 0.2)) +
  scale_y_continuous(expand = expansion(mult = c(0,0)))+
  theme_classic() +
  geom_flow(alpha = 0.6) +
  geom_stratum(alpha = 1,width = 0.5) +
  scale_fill_manual(values = c("protein_coding"="#459943",
                               #"lncRNA"= "#db6968",
                               "retained_intron" = "#f8984e",
                               "nonsense_mediated_decay" = "#4d97cd"),
                    labels=c("protein_coding"="Protein-coding",
                             #"lncRNA"= "LncRNA",
                             "retained_intron" = "Retained intron",
                             "nonsense_mediated_decay" = "NMD"))+
  labs(y= "Proportion (%)",x = NULL) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = 20, vjust = 0.5, hjust = 0.5),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        text = element_text(size = 15)) +
  annotate(geom = "text",
           x=2,y=Inf,
           label=TeX(r"(\textit{P} = 4.6e-23)"), # 4.6e-23
           vjust=1,hjust=1,size = 5)+
  ggtitle("HeLa")
dev.off()

#########################
###  SFig9A
#########################
#HEK293T
mycolor = c("#db6968","#4d97cd","#459943","#f8984e")
chisq.test(matrix(c(77,3,2,0,79,14,21,7),nrow=4,ncol=2))$p.value ## 3.7e-05
ggplot(HEK293_gtf,aes(x = group, stratum = transcript_type,alluvium = transcript_type,y = freq*100,
                      fill = transcript_type, label = transcript_type)) +
  scale_x_discrete(expand = c(0.2, 0.2)) +
  scale_y_continuous(expand = expansion(mult = c(0,0)))+
  theme_classic() +
  geom_flow(alpha = 0.6) +
  geom_stratum(alpha = 1,width = 0.5) +
  scale_fill_manual(values = c("protein_coding"="#459943",
                               "lncRNA"= "#db6968",
                               "retained_intron" = "#f8984e",
                               "nonsense_mediated_decay" = "#4d97cd"),
                    labels=c("protein_coding"="Protein-coding",
                             "lncRNA"= "LncRNA",
                             "retained_intron" = "Retained intron",
                             "nonsense_mediated_decay" = "NMD"))+
  labs(y= "Proportion (%)",x = NULL) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = 20, vjust = 0.5, hjust = 0.5),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        text = element_text(size = 15))+
  annotate(geom = "text",
           x=2,y=Inf,
           label=TeX(r"(\textit{P} = 3.7e-05)"),
           vjust=1,hjust=1,size = 5)+
  ggtitle("HEK293T")

#########################
###  SFig9B
#########################
### HepG2
mycolor = c("#db6968","#4d97cd","#459943","#f8984e")
chisq.test(matrix(c(94,3,1,1,93,12,21,11),nrow=4,ncol=2))$p.value ## 7.5e-06
ggplot(HepG2_gtf,aes(x = group, stratum = transcript_type,alluvium = transcript_type,y = freq*100,
                     fill = transcript_type, label = transcript_type)) +
  scale_x_discrete(expand = c(0.2, 0.2)) +
  scale_y_continuous(expand = expansion(mult = c(0,0)))+
  theme_classic() +
  geom_flow(alpha = 0.6) +
  geom_stratum(alpha = 1,width = 0.5) +
  scale_fill_manual(values = c("protein_coding"="#459943",
                               "lncRNA"= "#db6968",
                               "retained_intron" = "#f8984e",
                               "nonsense_mediated_decay" = "#4d97cd"),
                    labels=c("protein_coding"="Protein-coding",
                             "lncRNA"= "LncRNA",
                             "retained_intron" = "Retained intron",
                             "nonsense_mediated_decay" = "NMD"))+
  labs(y= "Proportion (%)",x = NULL) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = 20, vjust = 0.5, hjust = 0.5),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        text = element_text(size = 15)) +
  annotate(geom = "text",
           x=2,y=Inf,
           label=TeX(r"(\textit{P} = 7.5e-06)"),
           vjust=1,hjust=1,size = 5)+
  ggtitle("HepG2")

#########################
###  SFig9C
#########################
barcode = read_csv("filter_barcodes_cells_750_filter.csv",col_names = F)
cluster = read_tsv("c2t_site_ratio_cell_cluster_V2.tsv")

C2T_editing_itst <- read_tsv("/nh_c_to_t_snpmolinfos_trans.txt",col_names = c("chrom","snp_start","snp_end","transcript_id","strand","alt","read_name"),num_threads = 20)
duplicate_gene <- c2t_long_dt_T_HeLa %>% group_by(gene_name) %>% summarise(freq = n()) %>% filter(freq > 1) %>% dplyr::select(gene_name)  

## HeLa
### NMD
nonsense_trans <- gtf %>% filter(gene_name %in% duplicate_gene$gene_name) %>% filter(type == "transcript") %>% filter(transcript_type == "nonsense_mediated_decay")
nonsense_protein_trans <- gtf %>% filter(gene_name %in% nonsense_trans$gene_name) %>% filter(type == "transcript") %>% filter(transcript_type == "protein_coding")

HeLa_nonsense_trans <- c2t_long_dt_merge_HeLa_filter_dt %>% filter(transcript_name %in% nonsense_trans$transcript_name)
HeLa_nonsense_trans_dt <- HeLa_nonsense_trans %>% group_by(transcript_name,gene_name) %>% summarize(mean_ratio = mean(ratio)) %>% mutate(type = "NMD") %>% group_by(gene_name,type) %>% summarize(gene_mean_ratio = mean(mean_ratio)) ### 计算基因内NMD转录本的均值
HeLa_protein_trans <- c2t_long_dt_merge_HeLa_filter_dt %>% filter(transcript_name %in% nonsense_protein_trans$transcript_name) %>% filter(gene_name %in% HeLa_nonsense_trans_dt$gene_name)  ### 筛选配对基因的Protein-coding
HeLa_protein_trans_dt <- HeLa_protein_trans %>% group_by(transcript_name,gene_name) %>% summarize(mean_ratio = mean(ratio)) %>% mutate(type ="Protein-coding") %>% group_by(gene_name,type) %>% summarize(gene_mean_ratio = mean(mean_ratio)) ### 计算基因内PC转录本的均值
HeLa_nonsense_trans_dt <- HeLa_nonsense_trans_dt %>% filter(gene_name %in% HeLa_protein_trans_dt$gene_name)

HeLa_nonsense_protein_merge_trans <- rbind(HeLa_nonsense_trans_dt,HeLa_protein_trans_dt)
HeLa_nonsense_protein_merge_trans %>% group_by(type) %>% summarize(mean(mean_ratio))
HeLa_nonsense_protein_merge_trans_dt <- inner_join(HeLa_nonsense_trans_dt,HeLa_protein_trans_dt,by = "gene_name")

wilcox.test(HeLa_nonsense_protein_merge_trans_dt$mean_ratio.x,HeLa_nonsense_protein_merge_trans_dt$mean_ratio.y)
ggplot(HeLa_nonsense_protein_merge_trans_dt,aes(x = mean_ratio.y, y= mean_ratio.x))+
  geom_point(color = "#459943")+
  theme_classic()+
  theme(panel.border = element_blank(),
        plot.title = element_text(size = 20, vjust = 0.5, hjust = 0.5),
        text = element_text(size = 15),
        legend.justification = c(1,1))+
  scale_x_continuous(limits = c(0,0.5),
                     breaks = seq(0,1,0.1),
                     expand = expansion(mult = c(0,0)))+
  scale_y_continuous(limits = c(0,0.5),
                     breaks = seq(0,1,0.1),
                     expand = expansion(mult = c(0,0)))+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") + 
  labs(y=TeX(r"(Mean m${^6}$A level of NMD transcript)"),x = TeX(r"(Mean m${^6}$A level of Protein-coding transcript)")) +
  ggtitle("NMD vs Protein-coding")

#########################
###  SFig9D
#########################
### RI
retained_trans <- gtf %>% filter(gene_name %in% duplicate_gene$gene_name) %>% filter(type == "transcript") %>% filter(transcript_type == "retained_intron")
retained_protein_trans <- gtf %>% filter(gene_name %in% retained_trans$gene_name) %>% filter(type == "transcript") %>% filter(transcript_type == "protein_coding")

HeLa_retained_trans <- c2t_long_dt_merge_HeLa_filter_dt %>% filter(transcript_name %in% retained_trans$transcript_name)
HeLa_retained_trans_dt <- HeLa_retained_trans %>% group_by(transcript_name,gene_name) %>% summarize(mean_ratio = mean(ratio)) %>% mutate(type = "RI") %>% group_by(gene_name,type) %>% summarize(gene_mean_ratio = mean(mean_ratio))
HeLa_protein_trans <- c2t_long_dt_merge_HeLa_filter_dt %>% filter(transcript_name %in% retained_protein_trans$transcript_name) %>% filter(gene_name %in% HeLa_retained_trans_dt$gene_name)
HeLa_protein_trans_dt <- HeLa_protein_trans %>% group_by(transcript_name,gene_name) %>% summarize(mean_ratio = mean(ratio)) %>% mutate(type ="Protein-coding") %>% group_by(gene_name,type) %>% summarize(gene_mean_ratio = mean(mean_ratio))
HeLa_retained_trans_dt <- HeLa_retained_trans_dt %>% filter(gene_name %in% HeLa_protein_trans_dt$gene_name)

HeLa_retained_protein_merge_trans <- rbind(HeLa_retained_trans_dt,HeLa_protein_trans_dt) 

HeLa_retained_protein_merge_trans %>% group_by(type) %>% summarize(mean(mean_ratio))
HeLa_retained_protein_merge_trans_dt <- inner_join(HeLa_retained_trans_dt,HeLa_protein_trans_dt,by = "gene_name")

wilcox.test(HeLa_retained_protein_merge_trans_dt$mean_ratio.x,HeLa_retained_protein_merge_trans_dt$mean_ratio.y) # 3.163e-11
ggplot(HeLa_retained_protein_merge_trans_dt,aes(x = mean_ratio.y, y= mean_ratio.x))+
  geom_point(color = "#e64b35")+
  theme_classic()+
  theme(panel.border = element_blank(),
        plot.title = element_text(size = 20, vjust = 0.5, hjust = 0.5),
        text = element_text(size = 15),
        legend.justification = c(1,1))+
  scale_x_continuous(limits = c(0,0.5),
                     breaks = seq(0,1,0.1),
                     expand = expansion(mult = c(0,0)))+
  scale_y_continuous(limits = c(0,0.5),
                     breaks = seq(0,1,0.1),
                     expand = expansion(mult = c(0,0)))+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") + 
  labs(y=TeX(r"(Mean m${^6}$A level of RI transcript)"),x = TeX(r"(Mean m${^6}$A level of Protein-coding transcripts)")) +
  ggtitle("RI vs Protein-coding")

#########################
###  SFig9E
#########################
### LncRNA 
lncRNA_trans <- gtf %>% filter(gene_name %in% duplicate_gene$gene_name) %>% filter(type == "transcript") %>% filter(transcript_type == "lncRNA")
lncRNA_protein_trans <- gtf %>% filter(gene_name %in% lncRNA_trans$gene_name) %>% filter(type == "transcript") %>% filter(transcript_type == "protein_coding")

HeLa_lncRNA_trans <- c2t_long_dt_merge_HeLa_filter_dt %>% filter(transcript_name %in% lncRNA_trans$transcript_name)
HeLa_lncRNA_trans_dt <- HeLa_lncRNA_trans %>% group_by(transcript_name,gene_name) %>% summarize(mean_ratio = mean(ratio)) %>% mutate(type = "lncRNA")
HeLa_protein_trans <- c2t_long_dt_merge_HeLa_filter_dt %>% filter(transcript_name %in% lncRNA_protein_trans$transcript_name) %>% filter(gene_name %in% HeLa_lncRNA_trans_dt$gene_name) 
HeLa_protein_trans_dt <- HeLa_protein_trans %>% group_by(transcript_name,gene_name) %>% summarize(mean_ratio = mean(ratio)) %>% mutate(type ="Protein-coding")
HeLa_lncRNA_trans_dt <- HeLa_lncRNA_trans_dt %>% filter(gene_name %in% HeLa_protein_trans_dt$gene_name)

HeLa_lncRNA_protein_merge_trans <- rbind(HeLa_lncRNA_trans_dt,HeLa_protein_trans_dt) #%>
HeLa_lncRNA_protein_merge_trans %>% group_by(type) %>% summarize(mean(mean_ratio))
HeLa_lncRNA_protein_merge_trans_dt <- inner_join(HeLa_lncRNA_trans_dt,HeLa_protein_trans_dt,by = "gene_name")

wilcox.test(HeLa_lncRNA_protein_merge_trans_dt$mean_ratio.x,HeLa_lncRNA_protein_merge_trans_dt$mean_ratio.y)
ggplot(HeLa_lncRNA_protein_merge_trans_dt,aes(x = mean_ratio.y, y= mean_ratio.x))+
  geom_point(color = "#f8984e")+
  theme_classic()+
  theme(panel.border = element_blank(),
        plot.title = element_text(size = 20, vjust = 0.5, hjust = 0.5),
        text = element_text(size = 15),
        legend.justification = c(1,1))+
  scale_x_continuous(limits = c(0,0.5),
                     breaks = seq(0,1,0.1),
                     expand = expansion(mult = c(0,0)))+
  scale_y_continuous(limits = c(0,0.5),
                     breaks = seq(0,1,0.1),
                     expand = expansion(mult = c(0,0)))+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") + 
  labs(y=TeX(r"(Mean m${^6}$A level of lncRNA transcript)"),x = TeX(r"(Mean m${^6}$A level of Protein-coding transcripts)")) +
  ggtitle("LncRNA vs Protein-coding")

#########################
###  Fig6B
#########################
#### Merge boxplot
HeLa_merge_trans_nonsense_retained_lncRNA_dt <- rbind(HeLa_nonsense_protein_merge_trans %>% mutate(group = "NMD vs PC"),
                                                      HeLa_lncRNA_protein_merge_trans %>% mutate(group = "lncRNA vs PC"),
                                                      HeLa_retained_protein_merge_trans %>% mutate(group = "RI vs PC")) %>% mutate(group = factor(group,level = c("NMD vs PC","RI vs PC")),type = factor(type,level = c("NMD","RI","Protein-coding")))
compare_means(gene_mean_ratio ~ type,  data = HeLa_nonsense_protein_merge_trans, method = "wilcox.test") # NMD 0.00072
compare_means(gene_mean_ratio ~ type,  data = HeLa_retained_protein_merge_trans, method = "wilcox.test") # retained 8.9e-05
compare_means(gene_mean_ratio ~ type,  data = HeLa_lncRNA_protein_merge_trans, method = "wilcox.test") # lncRNA 0.00071
####  NMD_RI_LncRNA
ggplot(data=HeLa_merge_trans_nonsense_retained_lncRNA_dt,aes(x=group,y=gene_mean_ratio,fill = type))+
  geom_boxplot(aes(fill=type),outlier.shape = NA)+
  scale_fill_manual(values = c("NMD" = "#459943","Protein-coding" = "#5F80B4", "RI" = "#e64b35", "lncRNA" = "#f8984e"))+ #"#85B22E","#5F80B4" "#8A77A5","#DBE196"
  theme_classic()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = 20, vjust = 0.5, hjust = 0.5),
        legend.position = "none")+  
  annotate(geom = "text",
           x="NMD vs PC",y=0.65,
           label=TeX(r"(\textit{P} = 7.2e-04)"))+ 
  annotate(geom = "text",
           x="RI vs PC",y=0.65,
           label=TeX(r"(\textit{P} = 8.9e-05)"))+ 
  annotate(geom = "text",
           x="lncRNA vs PC",y=0.65,
           label=TeX(r"(\textit{P} = 7.1e-04)"))+ 
  geom_segment(aes(x = 0.75,y = 0.62,xend = 1.25,yend = 0.62), color = "black") +
  geom_segment(aes(x = 1.75,y = 0.62,xend = 2.25,yend = 0.62), color = "black") +
  geom_segment(aes(x = 2.75,y = 0.62,xend = 3.25,yend = 0.62), color = "black") +
  theme(text = element_text(size = 15)) +
  coord_cartesian(ylim = c(0.1,0.7))+
  labs(y=TeX(r"(Mean m${^6}$A level of transcripts)"),x = NULL) +
  ggtitle("HeLa Cell line") 

#########################
###  SFig10B;SFig10E
#########################
gtf <- as.data.frame(rtracklayer::import("gencode.v32.annotation.gtf"))
gtf <- gtf %>% dplyr::select("seqnames","start","end","strand","type","gene_name","gene_id","transcript_id","transcript_name","transcript_type","gene_type") %>% mutate(transcript_id = str_split_fixed(transcript_id,"\\.",n=2)[,1],gene_id = str_split_fixed(gene_id,"\\.",n=2)[,1])

exon_trans <- gtf %>% filter(type == "exon") %>% group_by(gene_name,transcript_id,transcript_name) %>% summarize(n_exons=n())
#### NanoCount
siUPF1_1_counts <- read_tsv("UPF1_task/NanoCount/HEK293T_siUPF1_1_NanoCount.tsv",col_names = T)
siUPF1_2_counts <- read_tsv("UPF1_task/NanoCount/HEK293T_siUPF1_2_NanoCount.tsv",col_names = T)
siCtrl_1_counts <- read_tsv("UPF1_task/NanoCount/HEK293T_siCtrl_1_NanoCount.tsv",col_names = T)
siCtrl_2_counts <- read_tsv("UPF1_task/NanoCount/HEK293T_siCtrl_2_NanoCount.tsv",col_names = T)
DMSO_1_counts <- read_tsv("UPF1_task/NanoCount/HEK293T_DMSO_1_NanoCount.tsv",col_names = T)
DMSO_2_counts <- read_tsv("UPF1_task/NanoCount/HEK293T_DMSO_2_NanoCount.tsv",col_names = T)
SMG1i_1_counts <- read_tsv("UPF1_task/NanoCount/HEK293T_SMG1i_1_NanoCount.tsv",col_names = T)
SMG1i_2_counts <- read_tsv("UPF1_task/NanoCount/HEK293T_SMG1i_2_NanoCount.tsv",col_names = T)

siUPF1_counts_TPM <- full_join(siUPF1_1_counts,siUPF1_2_counts,by = c("transcript_name","transcript_length")) %>% mutate(siUPF1_counts = est_count.x + est_count.y, raw = siUPF1_counts/sum(siUPF1_counts), siUPF1_TPM = 1e6 * raw) %>% mutate(transcript_id = str_split_fixed(transcript_name,"\\|",n=2)[,1],gene_id = str_split_fixed(transcript_name,"\\|",n=2)[,2]) %>% mutate(transcript_id = str_split_fixed(transcript_id,"\\.",n=2)[,1],gene_id = str_split_fixed(gene_id,"\\.",n=2)[,1]) %>% dplyr::select(gene_id,transcript_id,siUPF1_counts,siUPF1_TPM)
siCtrl_counts_TPM <- full_join(siCtrl_1_counts,siCtrl_2_counts,by = c("transcript_name","transcript_length")) %>% mutate(siCtrl_counts = est_count.x + est_count.y, raw = siCtrl_counts/sum(siCtrl_counts), siCtrl_TPM = 1e6 * raw) %>% mutate(transcript_id = str_split_fixed(transcript_name,"\\|",n=2)[,1],gene_id = str_split_fixed(transcript_name,"\\|",n=2)[,2]) %>% mutate(transcript_id = str_split_fixed(transcript_id,"\\.",n=2)[,1],gene_id = str_split_fixed(gene_id,"\\.",n=2)[,1]) %>% dplyr::select(gene_id,transcript_id,siCtrl_counts,siCtrl_TPM)
SMG1i_counts_TPM <- full_join(SMG1i_1_counts,SMG1i_2_counts,by = c("transcript_name","transcript_length")) %>% mutate(SMG1i_counts = est_count.x + est_count.y, raw = SMG1i_counts/sum(SMG1i_counts), SMG1i_TPM = 1e6 * raw) %>% mutate(transcript_id = str_split_fixed(transcript_name,"\\|",n=2)[,1],gene_id = str_split_fixed(transcript_name,"\\|",n=2)[,2]) %>% mutate(transcript_id = str_split_fixed(transcript_id,"\\.",n=2)[,1],gene_id = str_split_fixed(gene_id,"\\.",n=2)[,1]) %>% dplyr::select(gene_id,transcript_id,SMG1i_counts,SMG1i_TPM)
DMSO_counts_TPM <- full_join(DMSO_1_counts,DMSO_2_counts,by = c("transcript_name","transcript_length")) %>% mutate(DMSO_counts = est_count.x + est_count.y, raw = DMSO_counts/sum(DMSO_counts), DMSO_TPM = 1e6 * raw) %>% mutate(transcript_id = str_split_fixed(transcript_name,"\\|",n=2)[,1],gene_id = str_split_fixed(transcript_name,"\\|",n=2)[,2]) %>% mutate(transcript_id = str_split_fixed(transcript_id,"\\.",n=2)[,1],gene_id = str_split_fixed(gene_id,"\\.",n=2)[,1]) %>% dplyr::select(gene_id,transcript_id,DMSO_counts,DMSO_TPM)

siUPF1_siCtrl_counts_TPM <- inner_join(siUPF1_counts_TPM,siCtrl_counts_TPM) %>% filter(siCtrl_TPM >= 10 | siUPF1_TPM >= 10) %>% mutate(TPM_fc = siUPF1_TPM/siCtrl_TPM)
siUPF1_DMSO_counts_TPM <- inner_join(siUPF1_counts_TPM,DMSO_counts_TPM) %>% filter(DMSO_TPM >= 10 | siUPF1_TPM >= 10) %>% mutate(TPM_fc = siUPF1_TPM/DMSO_TPM)
SMG1i_DMSO_counts_TPM <- inner_join(SMG1i_counts_TPM,DMSO_counts_TPM) %>% filter(DMSO_TPM >= 10 | SMG1i_TPM >= 10) %>% mutate(TPM_fc = SMG1i_TPM/DMSO_TPM)

################
### m6Aiso
################
sample_names <- c("HEK293T_DMSO_1","HEK293T_DMSO_2","HEK293T_siUPF1_1","HEK293T_siUPF1_2","HEK293T_siCtrl_1","HEK293T_siCtrl_2","HEK293T_SMG1i_1","HEK293T_SMG1i_2")
merge_transpos_data <- tibble()
for (samplename in sample_names) {
  sub_transpo_data <- read_tsv(paste0("UPF1_task/m6Aiso/",samplename,"_V32.DRACH.tsv"),col_names = T) %>% mutate(group = sub("_[12]$", "", samplename)) %>% mutate(transcript_id = str_split_fixed(transcript_id,"\\.",n=2)[,1])
  merge_transpos_data <- bind_rows(merge_transpos_data,sub_transpo_data)
}
merge_transpos_data <- merge_transpos_data  %>% group_by(chrom, genepos, strand, transcript_id, kmer, group) %>% summarize(depth = sum(depth), modnum = sum(modnum)) %>% ungroup() %>% inner_join(.,gtf %>% filter(type == "transcript") %>% dplyr::select(transcript_id,gene_id,gene_name,transcript_name))
merge_genepos_data <- merge_transpos_data %>% group_by(chrom,genepos,strand,gene_id,gene_name,kmer,group) %>% summarize(depth = sum(depth), modnum = sum(modnum)) %>% ungroup()

siUPF1_genepos_filter_data <- merge_genepos_data %>% filter(group == "HEK293T_siUPF1") %>% filter(modnum >= 5,depth >= 20)
siUPF1_transpos_filter_data <- merge_transpos_data %>% filter(group == "HEK293T_siUPF1") %>% filter(depth >= 10) %>% inner_join(.,siUPF1_genepos_filter_data %>% dplyr::select(chrom,genepos,gene_id)) %>% mutate(site_ratio = modnum/depth)
siUPF1_transpos_filter_data_dt <- siUPF1_transpos_filter_data %>%
  inner_join(.,gtf %>% filter(type == "transcript") %>% dplyr::select("transcript_id","strand","gene_name","transcript_type"))
siUPF1_trans_filter_dt <- siUPF1_transpos_filter_data_dt %>% group_by(chrom,transcript_id,gene_id,gene_name,transcript_type,strand) %>% summarize(depth = sum(depth),modnum = sum(modnum)) %>% ungroup() %>% mutate(ratio = modnum/depth)
siUPF1_trans_sites_number_dt <- siUPF1_transpos_filter_data_dt %>% mutate(group = ifelse(site_ratio>=0.1,"Methy","Unmehty")) %>% group_by(gene_name,gene_id,transcript_id,transcript_name,transcript_type,group) %>% summarize(counts=n()) %>% ungroup() %>% pivot_wider(names_from = "group",values_from = c("counts"),values_fill = 0)

SMG1i_genepos_filter_data <- merge_genepos_data %>% filter(group == "HEK293T_SMG1i") %>% filter(modnum >= 5,depth >= 20)
SMG1i_transpos_filter_data <- merge_transpos_data %>% filter(group == "HEK293T_SMG1i") %>% filter(depth >= 10) %>% inner_join(.,SMG1i_genepos_filter_data %>% dplyr::select(chrom,genepos,gene_id)) %>% mutate(site_ratio = modnum/depth)
SMG1i_transpos_filter_data_dt <- SMG1i_transpos_filter_data %>%
  inner_join(.,gtf %>% filter(type == "transcript") %>% dplyr::select("transcript_id","strand","gene_name","transcript_type"))
SMG1i_trans_filter_dt <- SMG1i_transpos_filter_data_dt %>% group_by(chrom,transcript_id,gene_id,gene_name,transcript_type,strand) %>% summarize(depth = sum(depth),modnum = sum(modnum)) %>% ungroup() %>% mutate(ratio = modnum/depth)
SMG1i_trans_sites_number_dt <- SMG1i_transpos_filter_data_dt %>% mutate(group = ifelse(site_ratio>=0.1,"Methy","Unmehty")) %>% group_by(gene_name,gene_id,transcript_id,transcript_name,transcript_type,group) %>% summarize(counts=n()) %>% ungroup() %>% pivot_wider(names_from = "group",values_from = c("counts"),values_fill = 0)

DMSO_genepos_filter_data <- merge_genepos_data %>% filter(group == "HEK293T_DMSO") %>% filter(modnum >= 5,depth >= 20)
DMSO_transpos_filter_data <- merge_transpos_data %>% filter(group == "HEK293T_DMSO") %>% filter(depth >= 10) %>% inner_join(.,DMSO_genepos_filter_data %>% dplyr::select(chrom,genepos,gene_id)) %>% mutate(site_ratio = modnum/depth)
DMSO_transpos_filter_data_dt <- DMSO_transpos_filter_data %>%
  inner_join(.,gtf %>% filter(type == "transcript") %>% dplyr::select("transcript_id","strand","gene_name","transcript_type"))
DMSO_trans_filter_dt <- DMSO_transpos_filter_data_dt %>% group_by(chrom,transcript_id,gene_id,gene_name,transcript_type,strand) %>% summarize(depth = sum(depth),modnum = sum(modnum)) %>% ungroup() %>% mutate(ratio = modnum/depth)
DMSO_trans_sites_number_dt <- DMSO_transpos_filter_data_dt %>% mutate(group = ifelse(site_ratio>=0.1,"Methy","Unmehty")) %>% group_by(gene_name,gene_id,transcript_id,transcript_name,transcript_type,group) %>% summarize(counts=n()) %>% ungroup() %>% pivot_wider(names_from = "group",values_from = c("counts"),values_fill = 0)

##########################################
siUPF1_trans_m6A_expr_dt <- full_join(siUPF1_trans_data_dt,siUPF1_TPM) %>% mutate(across(.cols = 4:6, ~replace_na(., 0))) %>% mutate(group = case_when(ratio == 0 ~ "0", ratio > 0 & ratio <= 0.2 ~ "(0,0.2]", ratio > 0.2 & ratio <= 0.4 ~ "(0.2,0.4]",ratio > 0.4 & ratio <= 0.6 ~ "(0.4,0.6]",ratio > 0.6 & ratio <= 0.8 ~ "(0.6,0.8]",ratio > 0.8 & ratio <= 1 ~ "(0.8,1]")) ### siUPF1_TPM ### full_join
siUPF1_trans_m6A_expr_dt <- siUPF1_trans_m6A_expr_dt %>% group_by(group) %>%
  mutate(High = quantile(log2(HEK293T_siUPF1_TPM), 0.75), Med = median(log2(HEK293T_siUPF1_TPM)),Low = quantile(log2(HEK293T_siUPF1_TPM), 0.25)) %>% ungroup() %>% mutate(group = factor(group,levels = c("0","(0,0.2]","(0.2,0.4]","(0.4,0.6]","(0.6,0.8]","(0.8,1]"))) 
siUPF1_trans_m6A_fc_dt <- right_join(siUPF1_trans_data_dt,siUPF1_siCtrl_TPM) %>% mutate(across(.cols = 4:6, ~replace_na(., 0))) %>% mutate(group = case_when(ratio == 0 ~ "0", ratio > 0 & ratio <= 0.2 ~ "(0,0.2]", ratio > 0.2 & ratio <= 0.4 ~ "(0.2,0.4]",ratio > 0.4 & ratio <= 0.6 ~ "(0.4,0.6]",ratio > 0.6 & ratio <= 0.8 ~ "(0.6,0.8]",ratio > 0.8 & ratio <= 1 ~ "(0.8,1]")) ### siUPF1_siCtrl_TPM
siUPF1_trans_m6A_fc_dt <- siUPF1_trans_m6A_fc_dt %>% group_by(group) %>%
  mutate(High = quantile(log2(fc), 0.75), Med = median(log2(fc)),Low = quantile(log2(fc), 0.25)) %>% ungroup() %>% mutate(group = factor(group,levels = c("0","(0,0.2]","(0.2,0.4]","(0.4,0.6]","(0.6,0.8]","(0.8,1]"))) 

DMSO_trans_m6A_expr_dt <- full_join(DMSO_trans_data_dt,DMSO_TPM) %>% mutate(across(.cols = 4:6, ~replace_na(., 0))) %>% filter(ratio > 0) %>% mutate(group = case_when(ratio > 0 & ratio <= 0.2 ~ "(0,0.2]", ratio > 0.2 & ratio <= 0.4 ~ "(0.2,0.4]",ratio > 0.4 & ratio <= 0.6 ~ "(0.4,0.6]",ratio > 0.6 & ratio <= 0.8 ~ "(0.6,0.8]",ratio > 0.8 & ratio <= 1 ~ "(0.8,1]"))
DMSO_trans_m6A_expr_dt <- DMSO_trans_m6A_expr_dt %>% group_by(group) %>%
  mutate(High = quantile(log2(HEK293T_DMSO_TPM), 0.75), Med = median(log2(HEK293T_DMSO_TPM)),Low = quantile(log2(HEK293T_DMSO_TPM), 0.25)) %>% ungroup() %>% mutate(group = factor(group,levels = c("(0,0.2]","(0.2,0.4]","(0.4,0.6]","(0.6,0.8]","(0.8,1]"))) 
DMSO_trans_m6A_fc_dt <- right_join(DMSO_trans_data_dt,SMG1i_DMSO_TPM) %>% mutate(across(.cols = 4:6, ~replace_na(., 0))) %>% filter(ratio > 0) %>% mutate(group = case_when(ratio > 0 & ratio <= 0.2 ~ "(0,0.2]", ratio > 0.2 & ratio <= 0.4 ~ "(0.2,0.4]",ratio > 0.4 & ratio <= 0.6 ~ "(0.4,0.6]",ratio > 0.6 & ratio <= 0.8 ~ "(0.6,0.8]",ratio > 0.8 & ratio <= 1 ~ "(0.8,1]"))
DMSO_trans_m6A_fc_dt <- DMSO_trans_m6A_fc_dt %>% group_by(group) %>%
  mutate(High = quantile(log2(fc), 0.75), Med = median(log2(fc)),Low = quantile(log2(fc), 0.25)) %>% mutate(group = factor(group,levels = c("(0,0.2]","(0.2,0.4]","(0.4,0.6]","(0.6,0.8]","(0.8,1]"))) 

SMG1i_trans_m6A_expr_dt <- full_join(SMG1i_trans_data_dt,SMG1i_TPM) %>% mutate(across(.cols = 4:6, ~replace_na(., 0))) %>% filter(ratio > 0) %>% mutate(group = case_when(ratio > 0 & ratio <= 0.2 ~ "(0,0.2]", ratio > 0.2 & ratio <= 0.4 ~ "(0.2,0.4]",ratio > 0.4 & ratio <= 0.6 ~ "(0.4,0.6]",ratio > 0.6 & ratio <= 0.8 ~ "(0.6,0.8]",ratio > 0.8 & ratio <= 1 ~ "(0.8,1]"))
SMG1i_trans_m6A_expr_dt <- SMG1i_trans_m6A_expr_dt %>% group_by(group) %>%
  mutate(High = quantile(log2(HEK293T_SMG1i_TPM), 0.75), Med = median(log2(HEK293T_SMG1i_TPM)),Low = quantile(log2(HEK293T_SMG1i_TPM), 0.25)) %>% ungroup() %>% mutate(group = factor(group,levels = c("(0,0.2]","(0.2,0.4]","(0.4,0.6]","(0.6,0.8]","(0.8,1]"))) 
SMG1i_trans_m6A_fc_dt <- right_join(SMG1i_trans_data_dt,SMG1i_DMSO_TPM) %>%mutate(across(.cols = 4:6, ~replace_na(., 0))) %>% filter(ratio > 0) %>% mutate(group = case_when(ratio > 0 & ratio <= 0.2 ~ "(0,0.2]", ratio > 0.2 & ratio <= 0.4 ~ "(0.2,0.4]",ratio > 0.4 & ratio <= 0.6 ~ "(0.4,0.6]",ratio > 0.6 & ratio <= 0.8 ~ "(0.6,0.8]",ratio > 0.8 & ratio <= 1 ~ "(0.8,1]"))
SMG1i_trans_m6A_fc_dt <- SMG1i_trans_m6A_fc_dt %>% group_by(group) %>%
  mutate(High = quantile(log2(fc), 0.75), Med = median(log2(fc)),Low = quantile(log2(fc), 0.25)) %>% ungroup() %>% mutate(group = factor(group,levels = c("(0,0.2]","(0.2,0.4]","(0.4,0.6]","(0.6,0.8]","(0.8,1]"))) 


### LncRNA m6A
siUPF1_lncRNA_trans_m6A_expr_dt <- siUPF1_trans_m6A_expr_dt %>% filter(transcript_type == "lncRNA") %>% group_by(group) %>%
  mutate(High = quantile(log2(HEK293T_siUPF1_TPM), 0.75), Med = median(log2(HEK293T_siUPF1_TPM)),Low = quantile(log2(HEK293T_siUPF1_TPM), 0.25)) %>% ungroup() %>% mutate(group = factor(group,levels = c("(0,0.2]","(0.2,0.4]","(0.4,0.6]","(0.6,0.8]","(0.8,1]"))) 
siCtrl_lncRNA_trans_m6A_expr_dt <- siCtrl_trans_m6A_expr_dt %>% filter(transcript_type == "lncRNA") %>% group_by(group) %>%
  mutate(High = quantile(log2(HEK293T_siCtrl_TPM), 0.75), Med = median(log2(HEK293T_siCtrl_TPM)),Low = quantile(log2(HEK293T_siCtrl_TPM), 0.25)) %>% ungroup() %>% mutate(group = factor(group,levels = c("(0,0.2]","(0.2,0.4]","(0.4,0.6]","(0.6,0.8]","(0.8,1]"))) 
SMG1i_lncRNA_trans_m6A_expr_dt <- SMG1i_trans_m6A_expr_dt %>% filter(transcript_type == "lncRNA") %>% group_by(group) %>%
  mutate(High = quantile(log2(HEK293T_SMG1i_TPM), 0.75), Med = median(log2(HEK293T_SMG1i_TPM)),Low = quantile(log2(HEK293T_SMG1i_TPM), 0.25)) %>% ungroup() %>% mutate(group = factor(group,levels = c("(0,0.2]","(0.2,0.4]","(0.4,0.6]","(0.6,0.8]","(0.8,1]"))) 
DMSO_lncRNA_trans_m6A_expr_dt <- DMSO_trans_m6A_expr_dt %>% filter(transcript_type == "lncRNA") %>% group_by(group) %>%
  mutate(High = quantile(log2(HEK293T_DMSO_TPM), 0.75), Med = median(log2(HEK293T_DMSO_TPM)),Low = quantile(log2(HEK293T_DMSO_TPM), 0.25)) %>% ungroup() %>% mutate(group = factor(group,levels = c("(0,0.2]","(0.2,0.4]","(0.4,0.6]","(0.6,0.8]","(0.8,1]"))) 

siUPF1_lncRNA_trans_m6A_fc_dt <- siUPF1_trans_m6A_fc_dt %>% filter(transcript_type == "lncRNA") %>% group_by(group) %>%
  mutate(High = quantile(log2(fc), 0.75), Med = median(log2(fc)),Low = quantile(log2(fc), 0.25)) %>% ungroup() %>% mutate(group = factor(group,levels = c("(0,0.2]","(0.2,0.4]","(0.4,0.6]","(0.6,0.8]","(0.8,1]"))) 
siCtrl_lncRNA_trans_m6A_fc_dt <- siCtrl_trans_m6A_fc_dt %>% filter(transcript_type == "lncRNA") %>% group_by(group) %>%
  mutate(High = quantile(log2(fc), 0.75), Med = median(log2(fc)),Low = quantile(log2(fc), 0.25)) %>% ungroup() %>% mutate(group = factor(group,levels = c("(0,0.2]","(0.2,0.4]","(0.4,0.6]","(0.6,0.8]","(0.8,1]"))) 
SMG1i_lncRNA_trans_m6A_fc_dt <- SMG1i_trans_m6A_fc_dt %>% filter(transcript_type == "lncRNA") %>% group_by(group) %>%
  mutate(High = quantile(log2(fc), 0.75), Med = median(log2(fc)),Low = quantile(log2(fc), 0.25)) %>% ungroup() %>% mutate(group = factor(group,levels = c("(0,0.2]","(0.2,0.4]","(0.4,0.6]","(0.6,0.8]","(0.8,1]"))) 
DMSO_lncRNA_trans_m6A_fc_dt <- DMSO_trans_m6A_fc_dt %>% filter(transcript_type == "lncRNA") %>% group_by(group) %>%
  mutate(High = quantile(log2(fc), 0.75), Med = median(log2(fc)),Low = quantile(log2(fc), 0.25)) %>% ungroup() %>% mutate(group = factor(group,levels = c("(0,0.2]","(0.2,0.4]","(0.4,0.6]","(0.6,0.8]","(0.8,1]"))) 

table(DMSO_lncRNA_trans_m6A_expr_dt$group)
###(0,0.2] (0.2,0.4] (0.4,0.6] (0.6,0.8]   (0.8,1] 
### 210       128        97        36        19 
compare_means(HEK293T_DMSO_TPM ~ group,DMSO_lncRNA_trans_m6A_expr_dt,method = "wilcox.test")
compare_means(fc ~ group,DMSO_lncRNA_trans_m6A_fc_dt,method = "wilcox.test")

cum_data <- data.frame(group = c("(0,0.2]","(0.2,0.4]","(0.4,0.6]","(0.6,0.8]","(0.8,1]"), y = rep(-2,5), n_trans = c(210,128,97,36,19))
plot_TPM_violinplot(data = DMSO_lncRNA_trans_m6A_expr_dt %>% dplyr::rename(TPM = "HEK293T_DMSO_TPM"), title = "DMSO_lncRNA",cum_data = cum_data,y_position = c(13,11.5,10,8.5,7),additional_args = coord_cartesian(ylim = c(-3,15)), plot_path = "UPF1_task/Rplot/DMSO_lncRNA_m6Aiso_transcript_m6A_TPM_violinplot.pdf")
cum_data <- data.frame(group = c("(0,0.2]","(0.2,0.4]","(0.4,0.6]","(0.6,0.8]","(0.8,1]"), y = rep(-2,5), n_trans = c(210,127,97,35,18))
plot_fc_violinplot(data = DMSO_lncRNA_trans_m6A_fc_dt, title = "DMSO_lncRNA",cum_data = cum_data,y_position = c(5,4.5,4,3),coord_cartesian(ylim = c(-2.5,5.5)), labs(y=TeX(r"(Log${_2}$ FC(SMG1i/DMSO))"),x = TeX(r"(m${^6}$A level of transcripts)")), plot_path = "UPF1_task/Rplot/DMSO_lncRNA_m6Aiso_transcript_m6A_fc_violinplot.pdf")

#########################
###  SFig10C; SFig10F
#########################
siUPF1_RI_trans_m6A_expr_dt <- siUPF1_trans_m6A_expr_dt %>% filter(transcript_type == "retained_intron") %>% group_by(group) %>%
  mutate(High = quantile(log2(HEK293T_siUPF1_TPM), 0.75), Med = median(log2(HEK293T_siUPF1_TPM)),Low = quantile(log2(HEK293T_siUPF1_TPM), 0.25)) %>% ungroup() %>% mutate(group = factor(group,levels = c("(0,0.2]","(0.2,0.4]","(0.4,0.6]","(0.6,0.8]","(0.8,1]"))) 
siCtrl_RI_trans_m6A_expr_dt <- siCtrl_trans_m6A_expr_dt %>% filter(transcript_type == "retained_intron") %>% group_by(group) %>%
  mutate(High = quantile(log2(HEK293T_siCtrl_TPM), 0.75), Med = median(log2(HEK293T_siCtrl_TPM)),Low = quantile(log2(HEK293T_siCtrl_TPM), 0.25)) %>% ungroup() %>% mutate(group = factor(group,levels = c("(0,0.2]","(0.2,0.4]","(0.4,0.6]","(0.6,0.8]","(0.8,1]"))) 
SMG1i_RI_trans_m6A_expr_dt <- SMG1i_trans_m6A_expr_dt %>% filter(transcript_type == "retained_intron") %>% group_by(group) %>%
  mutate(High = quantile(log2(HEK293T_SMG1i_TPM), 0.75), Med = median(log2(HEK293T_SMG1i_TPM)),Low = quantile(log2(HEK293T_SMG1i_TPM), 0.25)) %>% ungroup() %>% mutate(group = factor(group,levels = c("(0,0.2]","(0.2,0.4]","(0.4,0.6]","(0.6,0.8]","(0.8,1]"))) 
DMSO_RI_trans_m6A_expr_dt <- DMSO_trans_m6A_expr_dt %>% filter(transcript_type == "retained_intron") %>% group_by(group) %>%
  mutate(High = quantile(log2(HEK293T_DMSO_TPM), 0.75), Med = median(log2(HEK293T_DMSO_TPM)),Low = quantile(log2(HEK293T_DMSO_TPM), 0.25)) %>% ungroup() %>% mutate(group = factor(group,levels = c("(0,0.2]","(0.2,0.4]","(0.4,0.6]","(0.6,0.8]","(0.8,1]"))) 

siUPF1_RI_trans_m6A_fc_dt <- siUPF1_trans_m6A_fc_dt %>% filter(transcript_type == "retained_intron") %>% group_by(group) %>%
  mutate(High = quantile(log2(fc), 0.75), Med = median(log2(fc)),Low = quantile(log2(fc), 0.25)) %>% ungroup() %>% mutate(group = factor(group,levels = c("(0,0.2]","(0.2,0.4]","(0.4,0.6]","(0.6,0.8]","(0.8,1]"))) 
siCtrl_RI_trans_m6A_fc_dt <- siCtrl_trans_m6A_fc_dt %>% filter(transcript_type == "retained_intron") %>% group_by(group) %>%
  mutate(High = quantile(log2(fc), 0.75), Med = median(log2(fc)),Low = quantile(log2(fc), 0.25)) %>% ungroup() %>% mutate(group = factor(group,levels = c("(0,0.2]","(0.2,0.4]","(0.4,0.6]","(0.6,0.8]","(0.8,1]"))) 
SMG1i_RI_trans_m6A_fc_dt <- SMG1i_trans_m6A_fc_dt %>% filter(transcript_type == "retained_intron") %>% group_by(group) %>%
  mutate(High = quantile(log2(fc), 0.75), Med = median(log2(fc)),Low = quantile(log2(fc), 0.25)) %>% ungroup() %>% mutate(group = factor(group,levels = c("(0,0.2]","(0.2,0.4]","(0.4,0.6]","(0.6,0.8]","(0.8,1]"))) 
DMSO_RI_trans_m6A_fc_dt <- DMSO_trans_m6A_fc_dt %>% filter(transcript_type == "retained_intron") %>% group_by(group) %>%
  mutate(High = quantile(log2(fc), 0.75), Med = median(log2(fc)),Low = quantile(log2(fc), 0.25)) %>% ungroup() %>% mutate(group = factor(group,levels = c("(0,0.2]","(0.2,0.4]","(0.4,0.6]","(0.6,0.8]","(0.8,1]"))) 

table(DMSO_RI_trans_m6A_expr_dt$group)
### (0,0.2] (0.2,0.4] (0.4,0.6] (0.6,0.8]   (0.8,1] 
## 417       298       124        44        17 
compare_means(HEK293T_DMSO_TPM ~ group,DMSO_RI_trans_m6A_expr_dt,method = "wilcox.test")
compare_means(fc ~ group,DMSO_RI_trans_m6A_fc_dt,method = "wilcox.test")

cum_data <- data.frame(group = c("(0,0.2]","(0.2,0.4]","(0.4,0.6]","(0.6,0.8]","(0.8,1]"), y = rep(-2,5), n_trans = c(417,298,124,44,17))
plot_TPM_violinplot(data = DMSO_RI_trans_m6A_expr_dt %>% dplyr::rename(TPM = "HEK293T_DMSO_TPM"), title = "DMSO_RI",cum_data = cum_data,y_position = c(13,11.5,10,8.5,7),additional_args = coord_cartesian(ylim = c(-3,15)), plot_path = "UPF1_task/Rplot/DMSO_RI_m6Aiso_transcript_m6A_TPM_violinplot.pdf")
cum_data <- data.frame(group = c("(0,0.2]","(0.2,0.4]","(0.4,0.6]","(0.6,0.8]","(0.8,1]"), y = rep(-2,5), n_trans = c(416,298,124,44,17))
plot_fc_violinplot(data = DMSO_RI_trans_m6A_fc_dt, title = "DMSO_RI",cum_data = cum_data,y_position = c(4,3.5,3,2.5,2),coord_cartesian(ylim = c(-3,5)), labs(y=TeX(r"(Log${_2}$ FC(SMG1i/DMSO))"),x = TeX(r"(m${^6}$A level of transcripts)")), plot_path = "UPF1_task/Rplot/DMSO_RI_m6Aiso_transcript_m6A_fc_violinplot.pdf")

#########################
###  SFig10O
#########################
HeLa_GLORI_m6A_rep1 <- read_tsv("GSM6432595_Hela-1_35bp_m2.totalm6A.FDR.csv",col_names = T)
HeLa_GLORI_m6A_rep2 <- read_tsv("GSM6432596_Hela-2_35bp_m2.totalm6A.FDR.csv",col_names = T)
HeLa_GLORI_m6A <- rbind(HeLa_GLORI_m6A_rep1,HeLa_GLORI_m6A_rep2) %>% group_by(Chr,Sites,Strand,Gene,Transcript) %>% summarize(AGcov = sum(AGcov),Acov = sum(Acov)) %>% ungroup() %>% mutate(Ratio = Acov/AGcov) %>% inner_join(.,gtf %>% filter(type == "gene") %>% dplyr::select(gene_name,gene_id), by = c("Gene" = "gene_name"))
HeLa_GLORI_m6A_bed <- HeLa_GLORI_m6A %>% mutate(start = Sites,end = Sites) %>% dplyr::select(Chr,start,end,gene_id,Strand)
write_tsv(HeLa_GLORI_m6A_bed,"HeLa_GLORI_m6A_V43.bed",col_names = F)

HeLa_UPF1_iCLIP <- read_tsv("HeLa_UPF1_iclip_sites_hg38_rep2.bed",col_names = c("chrom","start","end","peak_name","coverage","strand"))
HeLa_UPF1_iCLIP_gr <- GRanges(HeLa_UPF1_iCLIP$chrom,IRanges(start = HeLa_UPF1_iCLIP$start, end = HeLa_UPF1_iCLIP$end),HeLa_UPF1_iCLIP$strand)

HeLa_UPF1_iCLIP_gr_reduce <- reduce(HeLa_UPF1_iCLIP_gr, min.gapwidth = 50)
HeLa_UPF1_iCLIP_ovl <- findOverlaps(query = HeLa_UPF1_iCLIP_gr, subject = HeLa_UPF1_iCLIP_gr_reduce)  
HeLa_UPF1_iCLIP_cluster_counts <- table(subjectHits(HeLa_UPF1_iCLIP_ovl))
mcols(HeLa_UPF1_iCLIP_gr_reduce)$count <- as.integer(HeLa_UPF1_iCLIP_cluster_counts) 
HeLa_UPF1_iCLIP_peaks_dt <- HeLa_UPF1_iCLIP_gr_reduce %>% as.tibble() %>% filter(count >= 3)

HeLa_UPF1_iCLIP_peaks_dt_filter <- inner_join(HeLa_UPF1_iCLIP_peaks_dt,gtf %>% filter(type == "gene"),by = c("seqnames" = "seqnames","strand")) %>% filter(start.x >= start.y, end.x <= end.y)
HeLa_UPF1_iCLIP_peaks_filter_dt <- HeLa_UPF1_iCLIP_peaks_dt_filter %>% dplyr::rename(start = start.x, end = end.x) %>% dplyr::select(seqnames,start,end,gene_id,strand)  ## 38781
write_tsv(HeLa_UPF1_iCLIP_peaks_filter_dt,"/home/hjliang/m6a_sc/RBP_KD/HeLa_UPF1_iCLIP_peaks_V43.bed",col_names = F)
HeLa_UPF1_iCLIP_peaks_filter_dt <- read_tsv("/home/hjliang/m6a_sc/RBP_KD/HeLa_UPF1_iCLIP_peaks_V43.bed",col_names = c("chrom","start","end","gene_id","strand"))

DF2_PAR_CLIP_rep1_hg38_dt <- read_tsv("siYTHDF2_rep1_PARCLIP_hg38.bed",col_names = c("chromosome","cluster_start","cluster_end","strand","gene_name"))
DF2_PAR_CLIP_rep2_hg38_dt <- read_tsv("siYTHDF2_rep2_PARCLIP_hg38.bed",col_names = c("chromosome","cluster_start","cluster_end","strand","gene_name"))
DF2_PAR_CLIP_rep3_hg38_dt <- read_tsv("siYTHDF2_rep3_PARCLIP_hg38.bed",col_names = c("chromosome","cluster_start","cluster_end","strand","gene_name"))
DF2_PAR_CLIP_rep1_hg38_bed <- DF2_PAR_CLIP_rep1_hg38_dt %>% inner_join(.,gtf %>% filter(type == "gene") %>% dplyr::select(gene_name,gene_id)) %>% dplyr::select("chromosome","cluster_start","cluster_end","gene_id","strand")
DF2_PAR_CLIP_rep2_hg38_bed <- DF2_PAR_CLIP_rep2_hg38_dt %>% inner_join(.,gtf %>% filter(type == "gene") %>% dplyr::select(gene_name,gene_id)) %>% dplyr::select("chromosome","cluster_start","cluster_end","gene_id","strand")
DF2_PAR_CLIP_rep3_hg38_bed <- DF2_PAR_CLIP_rep3_hg38_dt %>% inner_join(.,gtf %>% filter(type == "gene") %>% dplyr::select(gene_name,gene_id)) %>% dplyr::select("chromosome","cluster_start","cluster_end","gene_id","strand")

write_tsv(DF2_PAR_CLIP_rep1_hg38_bed,"HeLa_YTHDF2_PAR_CLIP_rep1_V43.bed",col_names = F)
write_tsv(DF2_PAR_CLIP_rep2_hg38_bed,"HeLa_YTHDF2_PAR_CLIP_rep2_V43.bed",col_names = F)
write_tsv(DF2_PAR_CLIP_rep3_hg38_bed,"HeLa_YTHDF2_PAR_CLIP_rep3_V43.bed",col_names = F)

``` {linux}
{
  bash relative_position.sh ### V1.py
  #!/bin/bash ### 30,40,30
  python /home/hjliang/R_code/m6A_peak_quanlity_evaluate2-master/bin/m6A_peak_relative_position_to_transcript.py --peak_format_file /home/hjliang/m6a_sc/RBP_KD/HeLa_UPF1_iCLIP_peaks_V43.bed --max_length_transcript_file /home/hjliang/genomes/hg38/gencode.v43.maxtranscript.gpd --genome_size_file /home/hjliang/R_code/m6A_peak_quanlity_evaluate2-master/data/hg38_genome_size.txt --binSize 100000 --already_know_geneid True --output_file /home/hjliang/m6a_sc/RBP_KD/HeLa_UPF1_iCLIP_peaks_V43.result
  python /home/hjliang/R_code/m6A_peak_quanlity_evaluate2-master/bin/m6A_peak_relative_position_to_transcript.py --peak_format_file /home/hjliang/m6a_sc/RBP_KD/HeLa_YTHDF2_PAR_CLIP_rep3_V43.bed --max_length_transcript_file /home/hjliang/genomes/hg38/gencode.v43.maxtranscript.gpd --genome_size_file /home/hjliang/R_code/m6A_peak_quanlity_evaluate2-master/data/hg38_genome_size.txt --binSize 100000 --already_know_geneid True --output_file /home/hjliang/m6a_sc/RBP_KD/HeLa_YTHDF2_PAR_CLIP_rep3_V43.result
  python /home/hjliang/R_code/m6A_peak_quanlity_evaluate2-master/bin/m6A_peak_relative_position_to_transcript.py --peak_format_file /home/hjliang/m6a_sc/RBP_KD/HeLa_GLORI_m6A_V43.bed --max_length_transcript_file /home/hjliang/genomes/hg38/gencode.v43.maxtranscript.gpd --genome_size_file /home/hjliang/R_code/m6A_peak_quanlity_evaluate2-master/data/hg38_genome_size.txt --binSize 100000 --already_know_geneid True --output_file /home/hjliang/m6a_sc/RBP_KD/HeLa_GLORI_m6A_V43.result
}
```
dataset1 <- read.table("HeLa_UPF1_iCLIP_peaks_V43.result.precentage_for_plot.txt",stringsAsFactors = FALSE,header = TRUE)
dataset2 <- read.table("HeLa_YTHDF2_PAR_CLIP_rep3_V43.result.precentage_for_plot.txt",stringsAsFactors = FALSE,header = TRUE)
dataset3 <- read.table("HeLa_GLORI_m6A_V43.result.precentage_for_plot.txt",stringsAsFactors = FALSE,header = TRUE)

dataset1 %>% summarize(sum=sum(count)) #31768
dataset2 %>% summarize(sum=sum(count)) #32374
dataset3 %>% summarize(sum=sum(count)) #112006

dataset1 <- dataset1 %>% mutate(density = count/31768) %>% mutate(group = "HeLa_UPF1_iCLIP")
dataset2 <- dataset2 %>% mutate(density = count/32374) %>% mutate(group = "HeLa_YTHDF2_PAR_CLIP")
dataset3 <- dataset3 %>% mutate(density = count/112006) %>% mutate(group = "HeLa_GLORI")

dataset_merge <- rbind(dataset1,dataset2,dataset3)

ggplot() +
  geom_line(data = dataset_merge, aes(x = position, y = density, color = group), size = 1.5) +
  theme_classic()+
  ylim(0, 0.06) +
  labs(y="Density",x = NULL)+
  scale_color_manual(values = c("HeLa_UPF1_iCLIP" = "#e64b35","HeLa_YTHDF2_PAR_CLIP" = "#4d97cd","HeLa_GLORI" = "#459943"), 
                     labels = c("HeLa_UPF1_iCLIP" = "HeLa UPF1 iCLIP peaks(31768)","HeLa_YTHDF2_PAR_CLIP" = "HeLa YTHDF2 PAR CLIP(32374)","HeLa_GLORI" = "HeLa GLORI(112006)")) +
  theme(legend.position = c(0.4,0.9),
        plot.title = element_text(size = 20, vjust = 0.5, hjust = 0.5),
        legend.title = element_blank(),
        legend.background = element_blank(),
        text = element_text(size = 15))+
  scale_x_continuous(limits = c(0,100),
                     breaks = c(15, 50, 85),
                     labels = c("5'UTR", "CDS", "3'UTR"))+
  geom_vline(xintercept = c(30, 70), color = "gray", 
             linetype = "dashed", size = 1)+
  ggtitle("HeLa")


