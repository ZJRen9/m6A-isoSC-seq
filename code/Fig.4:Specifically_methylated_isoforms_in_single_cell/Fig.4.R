


######################
### Fig4A
######################




c2t_long_T_dt = fread("st_combined_sites_T.mtx")
c2t_long_C_dt = fread("st_combined_sites_C.mtx")
m6a_long_dt = fread("nh_isomatrix_normalize_transcript.txt") %>% dplyr::rename(transcript_name = V1)
cluster = read_tsv("c2t_site_ratio_cell_cluster_V2.tsv")
barcode = read_csv("filter_barcodes_cells_750_filter.csv",col_names = F)
transID2transname <- read_tsv("genes.t2g.mapping",col_names = c("transID","transname","geneID","genename"))
### alt
c2t_long_T_dt <- merge(transID2transname,c2t_long_T_dt,by.x = "transID",by.y = "transId") %>% dplyr::select(-c("transID","geneID","genename"))
rowname = c2t_long_T_dt$transname
c2t_long_T_dt = as.matrix(c2t_long_T_dt[,-1])
rownames(c2t_long_T_dt) = rowname
c2t_long_T_dt <- c2t_long_T_dt[,colnames(c2t_long_T_dt) %in% barcode$X1]
### ref
c2t_long_C_dt <- merge(transID2transname,c2t_long_C_dt,by.x = "transID",by.y = "transId") %>% dplyr::select(-c("transID","geneID","genename"))
rowname = c2t_long_C_dt$transname
c2t_long_C_dt = as.matrix(c2t_long_C_dt[,-1])
rownames(c2t_long_C_dt) = rowname
c2t_long_C_dt <- c2t_long_C_dt[,colnames(c2t_long_C_dt) %in% barcode$X1]

c2t_long_T_dt <-  as.data.frame(c2t_long_T_dt) %>% mutate(transcript_name = rowname) %>% inner_join(. ,transID2transname %>% dplyr::select(transname,genename),by = c("transcript_name" = "transname")) %>% dplyr::rename("gene_name" = "genename") %>% pivot_longer(cols = !c(transcript_name,gene_name), names_to = "barcode",values_to = "alt_counts") %>% filter(alt_counts > 0)
c2t_long_C_dt <-  as.data.frame(c2t_long_C_dt) %>% mutate(transcript_name = rowname) %>% inner_join(. ,transID2transname %>% dplyr::select(transname,genename),by = c("transcript_name" = "transname")) %>% dplyr::rename("gene_name" = "genename") %>% pivot_longer(cols = !c(transcript_name,gene_name), names_to = "barcode",values_to = "ref_counts") %>% filter(ref_counts > 0)

### expression
####  All Gene
m6a_long_dt = fread("nh_isomatrix_normalize_transcript.txt") %>% dplyr::rename(transcript_name = V1)
rowname = m6a_long_dt$transcript_name
m6a_long_dt = as.matrix(m6a_long_dt[,-1])
rownames(m6a_long_dt) = rowname
m6a_long_dt <- m6a_long_dt[,colnames(m6a_long_dt) %in% barcode$X1]
m6a_long_dt <-  as.data.frame(m6a_long_dt) %>% mutate(transcript_name = rowname) %>% inner_join(. ,transID2transname %>% dplyr::select(transname,genename),by = c("transcript_name" = "transname")) %>% dplyr::rename("gene_name" = "genename") %>% pivot_longer(cols = !c(transcript_name,gene_name), names_to = "barcode",values_to = "expression") %>% filter(expression > 0)

### HeLa 
c2t_long_merge_expr_dt <- full_join(c2t_long_T_dt,c2t_long_C_dt) %>% mutate(across(.cols = 4:5, ~replace_na(., 0))) %>% filter(alt_counts > 0)
c2t_long_merge_expr_HeLa_dt <- c2t_long_merge_expr_dt %>% inner_join(.,cluster) %>% filter(cell_line == "HeLa") %>% inner_join(.,m6a_long_dt)
c2t_long_filter_trans_HeLa_dt <- c2t_long_T_dt %>% inner_join(.,cluster) %>% filter(cell_line == "HeLa") %>% group_by(transcript_name,gene_name) %>% summarize(methy_cells =n()) %>% filter(methy_cells>=10) %>% ungroup()
c2t_long_merge_expr_filter_HeLa_dt <- c2t_long_merge_expr_HeLa_dt %>% filter(transcript_name %in% c2t_long_filter_trans_HeLa_dt$transcript_name) %>% mutate(ratio = alt_counts/(alt_counts + ref_counts))
c2t_long_merge_expr_filter_trans_HeLa_dt <- c2t_long_merge_expr_filter_HeLa_dt %>% group_by(transcript_name,gene_name) %>% summarize(mean_ratio = mean(ratio)) %>% inner_join(.,m6a_long_trans_dt) %>% inner_join(.,c2t_long_filter_trans_HeLa_dt) %>% mutate(methy_ratio = methy_cells/total_cells) %>% ungroup()
c2t_long_merge_expr_filter_trans_HeLa_group_dt <- c2t_long_merge_expr_filter_trans_HeLa_dt %>% mutate(group = ifelse(mean_ratio >= 0.75, ">=0.75","<0.75"))

c2t_long_merge_expr_iso_HeLa_dt <- tibble()

for (i in c2t_long_filter_trans_HeLa_dt$transcript_name){
  print(i)
  c2t_long_filter_sub <- c2t_long_merge_expr_filter_HeLa_dt %>% filter(transcript_name == i)
  c2t_long_filter_sub_stat <- cor.test(c2t_long_filter_sub$ratio,c2t_long_filter_sub$expression,method = "spearman")
  c2t_long_filter_sub_dt <- tibble(gene_name = unique(c2t_long_filter_sub$gene_name),transcript_name = unique(c2t_long_filter_sub$transcript_name),
                                   n = nrow(c2t_long_filter_sub), R = c2t_long_filter_sub_stat$estimate %>% as.numeric(), P.value = c2t_long_filter_sub_stat$p.value)
  c2t_long_merge_expr_iso_HeLa_dt <- rbind(c2t_long_merge_expr_iso_HeLa_dt,c2t_long_filter_sub_dt)
}
c2t_long_merge_expr_iso_HeLa_dt <- c2t_long_merge_expr_iso_HeLa_dt %>% mutate(padj = p.adjust(P.value,method = "fdr")) 

c2t_long_merge_expr_iso_HeLa_dt <- c2t_long_merge_expr_iso_HeLa_dt %>% mutate(group = case_when(R > 0 & padj < 0.05 ~ "Up transcripts",
                                                                                                R < 0 & padj < 0.05 ~ "Down transcripts",
                                                                                                TRUE ~ "Not Sig"))
table(c2t_long_merge_expr_iso_HeLa_dt$group)
### Down transcripts    Not Sig   Up transcripts 
###     516             3115               25 
ggplot(data=c2t_long_merge_expr_iso_HeLa_dt,aes(x=R,y=-log(padj,10)))+
  geom_point(size = 0.5,aes(color = group))+
  geom_hline(yintercept = -log(0.05,10),linetype = "dashed") +
  scale_color_manual(values = c(#"Sig" = "#e64b35",
    "Up transcripts"="#e64b35", #467897
    "Down transcripts"= "#2873B3",
    "Not Sig"="#aaaaaa"),
    labels=c(#"Sig" = "Significant: 861", # FDR < 0.05
      "Up transcripts"= "Positive (25)",
      "Down transcripts"= "Negative (516)",
      "Not Sig"="Non-significant (3115)"))+
  theme_classic()+
  theme(legend.position = c(0.7,0.8),
        plot.title = element_text(size = 20, vjust = 0.5, hjust = 0.5),
        legend.title = element_blank(),
        legend.background = element_blank())+
  theme(text = element_text(size = 15)) +
  scale_x_continuous(limits=c(-1,1)) +
  labs(x="R",
       y=TeX(r"(-Log${_1}{_0}$ {(}\textit{FDR}{)})")) + ### P value
  ggtitle("HeLa Single cell")

######################
### SFig7D
######################
#### HEK293T
c2t_long_merge_expr_dt <- full_join(c2t_long_T_dt,c2t_long_C_dt) %>% mutate(across(.cols = 4:5, ~replace_na(., 0))) %>% filter(alt_counts > 0)
c2t_long_merge_expr_HEK293T_dt <- c2t_long_merge_expr_dt %>% inner_join(.,cluster) %>% filter(cell_line == "HEK293T") %>% inner_join(.,m6a_long_dt)
c2t_long_filter_trans_HEK293T_dt <- c2t_long_T_dt %>% inner_join(.,cluster) %>% filter(cell_line == "HEK293T") %>% group_by(transcript_name,gene_name) %>% summarize(methy_cells =n()) %>% filter(methy_cells>=10) %>% ungroup()
c2t_long_merge_expr_filter_HEK293T_dt <- c2t_long_merge_expr_HEK293T_dt %>% filter(transcript_name %in% c2t_long_filter_trans_HEK293T_dt$transcript_name) %>% mutate(ratio = alt_counts/(alt_counts + ref_counts))
c2t_long_merge_expr_filter_trans_HEK293T_dt <- c2t_long_merge_expr_filter_HEK293T_dt %>% group_by(transcript_name,gene_name) %>% summarize(mean_ratio = mean(ratio)) %>% inner_join(.,m6a_long_trans_dt) %>% inner_join(.,c2t_long_filter_trans_HEK293T_dt) %>% mutate(methy_ratio = methy_cells/total_cells) %>% ungroup()
c2t_long_merge_expr_filter_trans_HEK293T_group_dt <- c2t_long_merge_expr_filter_trans_HEK293T_dt %>% mutate(group = ifelse(mean_ratio >= 0.75, ">=0.75","<0.75"))

c2t_long_merge_expr_iso_HEK293T_dt <- tibble()

for (i in c2t_long_filter_trans_HEK293T_dt$transcript_name){
  print(i)
  c2t_long_filter_sub <- c2t_long_merge_expr_filter_HEK293T_dt %>% filter(transcript_name == i)
  c2t_long_filter_sub_stat <- cor.test(c2t_long_filter_sub$ratio,c2t_long_filter_sub$expression,method = "spearman")
  c2t_long_filter_sub_dt <- tibble(gene_name = unique(c2t_long_filter_sub$gene_name),transcript_name = unique(c2t_long_filter_sub$transcript_name),
                                   n = nrow(c2t_long_filter_sub), R = c2t_long_filter_sub_stat$estimate %>% as.numeric(), P.value = c2t_long_filter_sub_stat$p.value)
  c2t_long_merge_expr_iso_HEK293T_dt <- rbind(c2t_long_merge_expr_iso_HEK293T_dt,c2t_long_filter_sub_dt)
}
c2t_long_merge_expr_iso_HEK293T_dt <- c2t_long_merge_expr_iso_HEK293T_dt %>% mutate(padj = p.adjust(P.value,method = "fdr")) 

c2t_long_merge_expr_iso_HEK293T_dt <- c2t_long_merge_expr_iso_HEK293T_dt %>% mutate(group = case_when(R > 0 & padj < 0.05 ~ "Up transcripts",
                                                                                                      R < 0 & padj < 0.05 ~ "Down transcripts",
                                                                                                      TRUE ~ "Not Sig"))
table(c2t_long_merge_expr_iso_HEK293T_dt$group)
## Down transcripts          Not Sig   Up transcripts 
##         126             1367               22
ggplot(data=c2t_long_merge_expr_iso_HEK293T_dt,aes(x=R,y=-log(padj,10)))+
  geom_point(size = 0.5,aes(color = group))+
  geom_hline(yintercept = -log(0.05,10),linetype = "dashed") +
  scale_color_manual(values = c(#"Sig" = "#e64b35",
    "Up transcripts"="#e64b35", #467897
    "Down transcripts"= "#2873B3",
    "Not Sig"="#aaaaaa"),
    labels=c(#"Sig" = "Significant: 861", # FDR < 0.05
      "Up transcripts"= "Positive (22)",
      "Down transcripts"= "Negative (126)",
      "Not Sig"="Non-significant (1367)"))+
  theme_classic()+
  theme(legend.position = c(0.7,0.8),
        plot.title = element_text(size = 20, vjust = 0.5, hjust = 0.5),
        legend.title = element_blank(),
        legend.background = element_blank())+
  theme(text = element_text(size = 15)) +
  scale_x_continuous(limits=c(-1,1)) +
  labs(x="R",
       y=TeX(r"(-Log${_1}{_0}$ {(}\textit{FDR}{)})")) + ### P value
  ggtitle("HEK293T Single cell")

######################
### SFig7E
######################
### HepG2
c2t_long_merge_expr_dt <- full_join(c2t_long_T_dt,c2t_long_C_dt) %>% mutate(across(.cols = 4:5, ~replace_na(., 0))) %>% filter(alt_counts > 0)
c2t_long_merge_expr_HepG2_dt <- c2t_long_merge_expr_dt %>% inner_join(.,cluster) %>% filter(cell_line == "HepG2") %>% inner_join(.,m6a_long_dt)
c2t_long_filter_trans_HepG2_dt <- c2t_long_T_dt %>% inner_join(.,cluster) %>% filter(cell_line == "HepG2") %>% group_by(transcript_name,gene_name) %>% summarize(methy_cells =n()) %>% filter(methy_cells>=10) %>% ungroup()
c2t_long_merge_expr_filter_HepG2_dt <- c2t_long_merge_expr_HepG2_dt %>% filter(transcript_name %in% c2t_long_filter_trans_HepG2_dt$transcript_name) %>% mutate(ratio = alt_counts/(alt_counts + ref_counts))
c2t_long_merge_expr_filter_trans_HepG2_dt <- c2t_long_merge_expr_filter_HepG2_dt %>% group_by(transcript_name,gene_name) %>% summarize(mean_ratio = mean(ratio)) %>% inner_join(.,m6a_long_trans_dt) %>% inner_join(.,c2t_long_filter_trans_HepG2_dt) %>% mutate(methy_ratio = methy_cells/total_cells) %>% ungroup()
c2t_long_merge_expr_filter_trans_HepG2_group_dt <- c2t_long_merge_expr_filter_trans_HepG2_dt %>% mutate(group = ifelse(mean_ratio >= 0.75, ">=0.75","<0.75"))

c2t_long_merge_expr_iso_HepG2_dt <- tibble()

for (i in c2t_long_filter_trans_HepG2_dt$transcript_name){
  print(i)
  c2t_long_filter_sub <- c2t_long_merge_expr_filter_HepG2_dt %>% filter(transcript_name == i)
  c2t_long_filter_sub_stat <- cor.test(c2t_long_filter_sub$ratio,c2t_long_filter_sub$expression,method = "spearman")
  c2t_long_filter_sub_dt <- tibble(gene_name = unique(c2t_long_filter_sub$gene_name),transcript_name = unique(c2t_long_filter_sub$transcript_name),
                                   n = nrow(c2t_long_filter_sub), R = c2t_long_filter_sub_stat$estimate %>% as.numeric(), P.value = c2t_long_filter_sub_stat$p.value)
  c2t_long_merge_expr_iso_HepG2_dt <- rbind(c2t_long_merge_expr_iso_HepG2_dt,c2t_long_filter_sub_dt)
}
c2t_long_merge_expr_iso_HepG2_dt <- c2t_long_merge_expr_iso_HepG2_dt %>% mutate(padj = p.adjust(P.value,method = "fdr")) 

c2t_long_merge_expr_iso_HepG2_dt <- c2t_long_merge_expr_iso_HepG2_dt %>% mutate(group = case_when(R > 0 & padj < 0.05 ~ "Up transcripts",
                                                                                                  R < 0 & padj < 0.05 ~ "Down transcripts",
                                                                                                  TRUE ~ "Not Sig"))
table(c2t_long_merge_expr_iso_HepG2_dt$group)
## Down transcripts          Not Sig   Up transcripts 
##           101             1807               27
ggplot(data=c2t_long_merge_expr_iso_HepG2_dt,aes(x=R,y=-log(padj,10)))+
  geom_point(size = 0.5,aes(color = group))+
  geom_hline(yintercept = -log(0.05,10),linetype = "dashed") +
  scale_color_manual(values = c(#"Sig" = "#e64b35",
    "Up transcripts"="#e64b35", #467897
    "Down transcripts"= "#2873B3",
    "Not Sig"="#aaaaaa"),
    labels=c(#"Sig" = "Significant: 861", # FDR < 0.05
      "Up transcripts"= "Positive (27)",
      "Down transcripts"= "Negative (101)",
      "Not Sig"="Non-significant (1807)"))+
  theme_classic()+
  theme(legend.position = c(0.7,0.8),
        plot.title = element_text(size = 20, vjust = 0.5, hjust = 0.5),
        legend.title = element_blank(),
        legend.background = element_blank())+
  theme(text = element_text(size = 15)) +
  scale_x_continuous(limits=c(-1,1)) +
  labs(x="R",
       y=TeX(r"(-Log${_1}{_0}$ {(}\textit{FDR}{)})")) + ### P value
  ggtitle("HepG2 Single cell")

######################
### Fig4B
######################
### m6A
c2t_long_merge_expr_filter_HeLa_trans_dt <- c2t_long_merge_expr_filter_HeLa_dt %>% group_by(transcript_name,gene_name,cell_line) %>% summarize(expression = mean(expression),ratio = mean(ratio)) %>% ungroup()
c2t_long_merge_expr_iso_HeLa_trans_dt <- inner_join(c2t_long_merge_expr_iso_HeLa_dt,c2t_long_merge_expr_filter_HeLa_trans_dt) %>% filter(group != "Up transcripts") %>% mutate(group = ifelse(group == "Not Sig","Non-significant","Negative correlation"))

c2t_long_merge_expr_iso_HeLa_trans_gene_name <- rbind(c2t_long_merge_expr_iso_HeLa_trans_dt %>% filter(group == "Non-significant") %>% distinct(gene_name), c2t_long_merge_expr_iso_HeLa_trans_dt %>% filter(group == "Negative correlation") %>% distinct(gene_name)) %>% group_by(gene_name) %>% summarize(n=n()) %>% filter(n>1)
c2t_long_merge_expr_iso_HeLa_trans_filter_dt <- c2t_long_merge_expr_iso_HeLa_trans_dt %>% filter(gene_name %in% c2t_long_merge_expr_iso_HeLa_trans_gene_name$gene_name) %>% na.omit()

compare_means(ratio ~ group,c2t_long_merge_expr_iso_HeLa_trans_dt,method = "wilcox.test") # 3.2e-15
ggplot(c2t_long_merge_expr_iso_HeLa_trans_dt,aes(x=group,y=ratio,fill = group))+
  stat_boxplot(geom = "errorbar",width=0.2)+
  geom_boxplot(aes(fill=group),outlier.shape = NA)+ 
  scale_fill_manual(values = c("Negative correlation" = "#3288bd", "Non-significant" = "#f46d43"),
                    labels = c("Negative correlation" = "Negative correlation","Non-significant" = "Non-significant"))+
  theme_classic()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = 20, vjust = 0.5, hjust = 0.5),
        legend.position = "none")+  
  geom_signif(comparisons = list(c("Negative correlation", "Non-significant")),test = "wilcox.test",
              map_signif_level = function(p) sprintf("P = %.2g",p),
              textsize = 5,
              y_position = c(0.7))+
  theme(text = element_text(size = 15)) +
  labs(y=TeX(r"(Isoform m${^6}$A level)"),x = NULL)+
  coord_cartesian(ylim = c(0,1))+
  ggtitle("HeLa Neagative correlation transcripts") 

######################
### Fig4C
######################
c2t_long_merge_expr_filter_HeLa_dt <- c2t_long_merge_expr_filter_HeLa_dt %>% mutate(group = case_when(ratio > 0 & ratio <= 0.2 ~ "(0,0.2]",
                                                                                                      ratio > 0.2 & ratio <= 0.4 ~ "(0.2,0.4]",
                                                                                                      ratio > 0.4 & ratio <= 0.6 ~ "(0.4,0.6]",
                                                                                                      ratio > 0.6 & ratio <= 1 ~ "(0.6,1]"))
library(RColorBrewer)
blue_range <- colorRampPalette(c("#a7bfe8","#6190e8"))
compare_means(expression ~ group,  data = c2t_long_merge_expr_filter_HeLa_dt %>% filter(transcript_name == "ISG15-203"), method = "wilcox.test")
ggplot(data=c2t_long_merge_expr_filter_HeLa_dt %>% filter(transcript_name == "ISG15-203"),aes(x=group,y=expression))+  #sub_total_wilcox
  stat_boxplot(geom = "errorbar",width=0.2)+
  geom_boxplot(aes(fill=group),outlier.shape = NA)+
  scale_fill_manual(values = blue_range(4))+
  theme_classic()+
  theme(legend.position = "none",
        plot.title = element_text(size = 20, vjust = 0.5, hjust = 0.5),
        legend.title = element_blank(),
        legend.background = element_blank())+
  theme(text = element_text(size = 15)) +
  scale_y_continuous(limits=c(1,7.5)) +
  labs(x="m6A level of isoforms",
       y = "Expression of isoforms") + ### P value
  geom_signif(comparisons = list(c("(0,0.2]", "(0.2,0.4]"),
                                 c("(0.2,0.4]","(0.4,0.6]"),
                                 c("(0.4,0.6]","(0.6,1]")),test = "wilcox.test",
              map_signif_level = function(p) sprintf("P = %.2g",p),
              textsize = 5,
              y_position = c(7,5,4))+
  ggtitle("HeLa_ISG15-203") 

######################
### Fig4D
######################
### exon number
c2t_long_merge_expr_iso_HeLa_exon_number_dt <- c2t_long_merge_expr_iso_HeLa_dt %>% inner_join(.,exon_trans) %>% filter(group != "Up transcripts") %>% mutate(group = ifelse(group == "Not Sig","Non-significant","Negative correlation"))
c2t_long_merge_expr_iso_HeLa_exon_number_gene_name <- rbind(c2t_long_merge_expr_iso_HeLa_exon_number_dt %>% filter(group == "Non-significant") %>% distinct(gene_name), c2t_long_merge_expr_iso_HeLa_exon_number_dt %>% filter(group == "Negative correlation") %>% distinct(gene_name)) %>% group_by(gene_name) %>% summarize(n=n()) %>% filter(n>1)
c2t_long_merge_expr_iso_HeLa_exon_number_filter_dt <- c2t_long_merge_expr_iso_HeLa_exon_number_dt %>% filter(gene_name %in% c2t_long_merge_expr_iso_HeLa_exon_number_gene_name$gene_name) %>% na.omit()

compare_means(n_exons ~ group,c2t_long_merge_expr_iso_HeLa_exon_number_dt,method = "wilcox.test") #  9.2e-07
ggplot(c2t_long_merge_expr_iso_HeLa_exon_number_dt,aes(x=group,y=n_exons,fill = group))+
  stat_boxplot(geom = "errorbar",width=0.2)+
  geom_boxplot(aes(fill=group),outlier.shape = NA)+
  scale_fill_manual(values = c("Negative correlation" = "#4d97cd", "Non-significant" = "#459943"),
                    labels = c("Negative correlation" = "Negative correlation","Non-significant" = "Non-significant"))+
  theme_classic()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = 20, vjust = 0.5, hjust = 0.5),
        legend.position = "none")+  
  geom_signif(comparisons = list(c("Negative correlation", "Non-significant")),test = "wilcox.test",
              map_signif_level = function(p) sprintf("P = %.2g",p),
              textsize = 5,
              y_position = c(12))+
  theme(text = element_text(size = 15)) +
  labs(y="Exon number",x = NULL)+
  coord_cartesian(ylim = c(1,15))+
  ggtitle("HeLa Neagative correlation transcripts") 

######################
### Fig4E
######################
transcript_HeLa_dt <- read_tsv("HeLa_diff_ratio_based_iso_specific_all_m6a.txt",col_names = T)
transcript_HeLa_dt_pvalue_df <- transcript_HeLa_dt %>%
  mutate(change=case_when(
    padj < 0.05 & diff_ratio > 0 ~ "Hypermethylated transcripts",
    padj < 0.05 & diff_ratio < 0 ~ "Hypomethylated transcripts",
    padj > 0.05 ~ "Not Sig"
  ))

ggplot(data=transcript_HeLa_dt_pvalue_df,aes(x=diff_ratio,y=-log(padj,10)))+
  geom_point(aes(color=change),size = 1)+
  geom_hline(yintercept = -log(0.05,10),linetype = "dashed") +
  scale_color_manual(values = c(#"Sig" = "#e64b35"
    "Hypermethylated transcripts"="#e64b35", #467897
    "Hypomethylated transcripts"= "#2873B3",
    "Not Sig"="#aaaaaa"),
    labels=c(#"Sig" = "Significant: 861", # FDR < 0.05
      "Hypermethylated transcripts"= "Highly methylated transcripts (507)",
      "Hypomethylated transcripts"= "Lowly methylated  transcripts (354)",
      "Not Sig"="Non-significant (1112)"))+
  theme_classic()+
  theme(legend.position = c(0.9,0.9),
        plot.title = element_text(size = 20, vjust = 0.5, hjust = 0.5),
        legend.title = element_blank(),
        legend.background = element_blank())+
  theme(text = element_text(size = 15)) +
  labs(x=TeX(r"(Mean m${^6}$A ratio)"),
       y=TeX(r"(-Log${_1}{_0}$ {(}\textit{FDR}{)})")) + ### P value FDR
  ggtitle("HeLa")

######################
### Fig4F
######################
m6a_long_dt = fread("nh_isomatrix_normalize_transcript.txt") %>% dplyr::rename(transcript_name = V1)
rowname = m6a_long_dt$transcript_name   # name rows by the gene symbol
m6a_long_dt = as.matrix(m6a_long_dt[,-1]) # remove the first column, which is the gene ID
rownames(m6a_long_dt) = rowname
m6a_long_dt <- m6a_long_dt[,colnames(m6a_long_dt) %in% barcode$X1]

m6a_long_dt_HeLa <- as.data.frame(m6a_long_dt[,colnames(m6a_long_dt) %in% (cluster %>% filter(cell_line=="HeLa"))$barcode]) %>% mutate(transcript_name = rowname) %>% inner_join(. ,transID2transname %>% dplyr::select(transname,genename),by = c("transcript_name" = "transname")) %>% dplyr::rename("gene_name" = "genename")
m6a_long_dt_HEK293 <- as.data.frame(m6a_long_dt[,colnames(m6a_long_dt) %in% (cluster %>% filter(cell_line=="HEK293T"))$barcode]) %>% mutate(transcript_name = rowname) %>% inner_join(. ,transID2transname %>% dplyr::select(transname,genename),by = c("transcript_name" = "transname")) %>% dplyr::rename("gene_name" = "genename")
m6a_long_dt_HepG2 <- as.data.frame(m6a_long_dt[,colnames(m6a_long_dt) %in% (cluster %>% filter(cell_line=="HepG2"))$barcode]) %>% mutate(transcript_name = rowname) %>% inner_join(. ,transID2transname %>% dplyr::select(transname,genename),by = c("transcript_name" = "transname")) %>% dplyr::rename("gene_name" = "genename")

### HeLa
m6a_long_dt_HeLa_dt <- m6a_long_dt_HeLa %>% pivot_longer(cols = !c(transcript_name,gene_name), names_to = "barcode",values_to = "expression") 
m6a_long_dt_HeLa_filter_dt <- m6a_long_dt_HeLa_dt %>% filter(expression > 0)
### HEK293T
m6a_long_dt_HEK293_dt <- m6a_long_dt_HEK293 %>% pivot_longer(cols = !c(transcript_name,gene_name), names_to = "barcode",values_to = "expression") 
m6a_long_dt_HEK293_filter_dt <- m6a_long_dt_HEK293_dt %>% filter(expression > 0)
### HepG2
m6a_long_dt_HepG2_dt <- m6a_long_dt_HepG2 %>% pivot_longer(cols = !c(transcript_name,gene_name), names_to = "barcode",values_to = "expression") 
m6a_long_dt_HepG2_filter_dt <- m6a_long_dt_HepG2_dt %>% filter(expression > 0)

HeLa_ratio_expression_dt <- m6a_long_dt_HeLa_filter_dt %>% group_by(transcript_name,gene_name) %>% summarize(mean_expression = mean(expression)) %>% ungroup()
HeLa_ratio_expression_dt <- inner_join(HeLa_ratio_expression_dt,HeLa_total_dt)

library(paletteer)
library(latex2exp)
compare_means(mean_expression ~ group,  data = HeLa_ratio_expression_dt, method = "wilcox.test") # 4.2e-32
ggplot(data=HeLa_ratio_expression_dt %>% mutate(group = ifelse(group == "high_m6A_transcript","high","low")),aes(x=group,y=mean_expression))+
  geom_boxplot(aes(fill=group))+
  scale_fill_manual(values = c("#85B22E","#5F80B4"))+
  geom_jitter(width = 0.2)+
  theme_classic()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = 20, vjust = 0.5, hjust = 0.5),
        legend.position = "none")+  
  geom_signif(comparisons = list(c("high", "low")),test = "t.test",
              map_signif_level = F,
              textsize = 5, 
              annotations = c("P = 4.2e-32"))+
  theme(text = element_text(size = 15)) +
  labs(y="Mean of expression",x = TeX(r"(Isoform m${^6}$A level)"))+
  ggtitle("HeLa") 

######################
### SFig7F
######################
### HepG2
transcript_HepG2_dt <- read_tsv("HepG2_diff_ratio_based_iso_specific_all_m6a.txt",col_names = T)
transcript_HepG2_dt_pvalue_df <- transcript_HepG2_dt %>%
  mutate(change=case_when(
    padj < 0.05 & diff_ratio > 0 ~ "Hypermethylated transcripts",
    padj < 0.05 & diff_ratio < 0 ~ "Hypomethylated transcripts",
    padj > 0.05 ~ "Not Sig"
  ))

ggplot(data=transcript_HepG2_dt_pvalue_df,aes(x=diff_ratio,y=-log(padj,10)))+
  geom_point(aes(color=change))+
  geom_hline(yintercept = -log(0.05,10),linetype = "dashed") +
  #geom_vline(xintercept = -0.1,linetype = "dashed") +
  scale_color_manual(values = c(#"Sig" = "#e64b35"
    "Hypermethylated transcripts"="#e64b35", #467897
    "Hypomethylated transcripts"= "#2873B3",
    "Not Sig"="#aaaaaa"),
    labels=c(#"Sig" = "Significant: 861", # FDR < 0.05
      "Hypermethylated transcripts"= "Highly methylated transcripts (137)",
      "Hypomethylated transcripts"= "Lowly methylated  transcripts (99)",
      "Not Sig"="Non-significant (518)"))+
  theme_classic()+
  theme(legend.position = c(0.9,0.9),
        plot.title = element_text(size = 20, vjust = 0.5, hjust = 0.5),
        #legend.text.align = 0,
        legend.title = element_blank(),
        legend.background = element_blank())+
  theme(text = element_text(size = 15)) +
  labs(x=TeX(r"(m${^6}$A Mean ratio)"),
       y=TeX(r"(-Log${_1}{_0}$ {(}\textit{FDR}{)})")) + ### P value
  ggtitle("HepG2")

######################
### SFig7G
######################
HEK293_ratio_expression_dt <- m6a_long_dt_HEK293_filter_dt %>% group_by(transcript_name,gene_name) %>% summarize(mean_expression = mean(expression)) %>% ungroup()
HEK293_ratio_expression_dt <- inner_join(HEK293_ratio_expression_dt,HEK293_total_dt)

compare_means(mean_expression ~ group,  data = HEK293_ratio_expression_dt, method = "wilcox.test") # 5.6e-14
ggplot(data=HEK293_ratio_expression_dt %>% mutate(group = ifelse(group == "high_m6A_transcript","high","low")),aes(x=group,y=mean_expression))+
  geom_boxplot(aes(fill=group))+
  scale_fill_manual(values = c("#85B22E","#5F80B4"))+
  geom_jitter(width = 0.2)+
  theme_classic()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = 20, vjust = 0.5, hjust = 0.5),

        legend.position = "none")+  
  geom_signif(comparisons = list(c("high", "low")),test = "t.test",
              map_signif_level = F,
              textsize = 5, 
              annotations = c("P = 5.6e-14"))+
  theme(text = element_text(size = 15)) +
  labs(y="Mean of expression",x = TeX(r"(Isoform m${^6}$A level)"))+
  ggtitle("HEK293T") 

######################
### SFig7H
######################
### HEK293T
transcript_HEK293_dt <- read_tsv("HEK293T_diff_ratio_based_iso_specific_all_m6a.txt",col_names = T)
transcript_HEK293_dt_pvalue_df <- transcript_HEK293_dt %>%
  mutate(change=case_when(
    padj < 0.05 & diff_ratio > 0 ~ "Hypermethylated transcripts",
    padj < 0.05 & diff_ratio < 0 ~ "Hypomethylated transcripts",
    padj > 0.05 ~ "Not Sig"
  ))

ggplot(data=transcript_HEK293_dt_pvalue_df,aes(x=diff_ratio,y=-log(padj,10)))+
  geom_point(aes(color=change))+
  geom_hline(yintercept = -log(0.05,10),linetype = "dashed") +
  #geom_vline(xintercept = -0.1,linetype = "dashed") +
  scale_color_manual(values = c(#"Sig" = "#e64b35",
    "Hypermethylated transcripts"="#e64b35", #467897
    "Hypomethylated transcripts"= "#2873B3",
    "Not Sig"="#aaaaaa"),
    labels=c(#"Sig" = "Significant: 861", # FDR < 0.05
      "Hypermethylated transcripts"= "Highly methylated transcripts (121)",
      "Hypomethylated transcripts"= "Lowly methylated  transcripts (82)",
      "Not Sig"="Non-significant (423)"))+
  theme_classic()+
  theme(legend.position = c(0.9,0.9),
        #legend.text.align = 0,
        plot.title = element_text(size = 20, vjust = 0.5, hjust = 0.5),
        legend.title = element_blank(),
        legend.background = element_blank())+
  theme(text = element_text(size = 15)) +
  labs(x=TeX(r"(m${^6}$A Mean ratio)"),
       y=TeX(r"(-Log${_1}{_0}$ {(}\textit{FDR}{)})")) + ### P value
  ggtitle("HEK293T")

######################
### SFig7I
######################
HepG2_ratio_expression_dt <- m6a_long_dt_HepG2_filter_dt %>% group_by(transcript_name,gene_name) %>% summarize(mean_expression = mean(expression)) %>% ungroup()
HepG2_ratio_expression_dt <- inner_join(HepG2_ratio_expression_dt,HepG2_total_dt) 

compare_means(mean_expression ~ group,  data = HepG2_ratio_expression_dt, method = "wilcox.test") # 5.30e-16
ggplot(data=HepG2_ratio_expression_dt %>% mutate(group = ifelse(group == "high_m6A_transcript","high","low")),aes(x=group,y=mean_expression))+
  geom_boxplot(aes(fill=group))+
  scale_fill_manual(values = c("#85B22E","#5F80B4"))+
  geom_jitter(width = 0.2)+
  theme_classic()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = 20, vjust = 0.5, hjust = 0.5),
        legend.position = "none")+  
  geom_signif(comparisons = list(c("high", "low")),test = "t.test",
              map_signif_level = F,
              textsize = 5, 
              annotations = c("P = 5.3e-16"))+
  theme(text = element_text(size = 15)) +
  labs(y="Mean of expression",x = TeX(r"(Isoform m${^6}$A level)"))+
  ggtitle("HepG2") 

######################
###  Fig4G
######################
### eulerr
##  phyper  Use FDR

HeLa_HepG2_high_m6A_trans <- unique(intersect(HeLa_high_m6A_trans,HepG2_high_m6A_trans))
length(HeLa_HepG2_high_m6A_trans) # 120
HeLa_HEK293_high_m6A_trans <- unique(intersect(HeLa_high_m6A_trans,HEK293_high_m6A_trans))
length(HeLa_HEK293_high_m6A_trans) # 111
HEK293_HepG2_high_m6A_trans <- unique(intersect(HEK293_high_m6A_trans,HepG2_high_m6A_trans))
length(HEK293_HepG2_high_m6A_trans) # 70
c2t_long_dt_T_HeLa %>% nrow() ## 24225 background transcript

HeLa_high_m6A_trans <- unique((HeLa_ratio_expression_dt %>% filter(group == "high_m6A_transcript"))$transcript_name)
length(HeLa_high_m6A_trans) # 507
HepG2_high_m6A_trans <- unique((HepG2_ratio_expression_dt %>% filter(group == "high_m6A_transcript"))$transcript_name)
length(HepG2_high_m6A_trans) # 137
HEK293_high_m6A_trans <- unique((HEK293_ratio_expression_dt %>% filter(group == "high_m6A_transcript"))$transcript_name)
length(HEK293_high_m6A_trans) # 121

HeLa_HepG2_phyper <- phyper(120-1,507, 24225-507, 137, lower.tail = F)
HeLa_HepG2_phyper # 1.444869e-187
HeLa_HEK293T_phyper <- phyper(111-1,507, 24225-507, 121, lower.tail = F)
HeLa_HEK293T_phyper  # 1.218773e-178
HEK293_HepG2_phyper <- phyper(70-1,137, 24225-137, 121, lower.tail = F)
HEK293_HepG2_phyper # 7.268166e-133
library(eulerr)
Venn <- function(lst.A, lst.B, lst.C){
  A = lst.A %>% length()
  B = lst.B %>% length()
  C = lst.C %>% length()
  AB = intersect(lst.A, lst.B) %>% length()
  AC = intersect(lst.A, lst.C) %>% length() 
  BC = intersect(lst.B, lst.C) %>% length()
  ABC = intersect(lst.A, lst.B) %>% intersect(., lst.C) %>% length()
  euler(c(
    "A" = A - AB - AC + ABC,
    "B" = B - AB - BC + ABC,
    "C" = C - AC - BC + ABC,
    "A&B" = AB - ABC,
    "A&C" = AC - ABC,
    "B&C" = BC - ABC,
    "A&B&C" = ABC
  )) %>% return()
} 

vd <- Venn(HeLa_high_m6A_trans,HepG2_high_m6A_trans,HEK293_high_m6A_trans)
plot(vd,
     labels = list(
       labels = c('HeLa', 'HepG2', 'HEK293T'),
       col = "black", font = 2
     ), 
     edges = list(col="white",alpha=0),
     fills = list(fill = c("dodgerblue", "seagreen3", "orchid3"), alpha = 0.4),
     quantities = list(cex=1, col='black'),
     main = list(label=c("Highly methylated transcripts in each celltype"),cex=2),
     legend = list(labels=c("A:HeLa","B:HepG2","C:HEK293T"),cex=1.2)
)
grid.text("P(AvsB)=1.4e-187", x=0.8, y=0.25, rot=1,gp=gpar(fontsize=13, col="black"))
grid.text("P(AvsC)=1.2e-178", x=0.8, y=0.2, rot=1,gp=gpar(fontsize=13, col="black"))
grid.text("P(BvsC)=7.3e-133", x=0.8, y=0.15, rot=1,gp=gpar(fontsize=13, col="black"))

######################
###  Fig4H
######################
### eulerr
##  phyper  Use FDR

HeLa_low_m6A_trans <- unique((HeLa_ratio_expression_dt %>% filter(group == "low_m6A_transcript"))$transcript_name)
length(HeLa_low_m6A_trans) #354
HepG2_low_m6A_trans <- unique((HepG2_ratio_expression_dt %>% filter(group == "low_m6A_transcript"))$transcript_name)
length(HepG2_low_m6A_trans) #99
HEK293_low_m6A_trans <- unique((HEK293_ratio_expression_dt %>% filter(group == "low_m6A_transcript"))$transcript_name)
length(HEK293_low_m6A_trans) #82

HeLa_HepG2_low_m6A_trans <- unique(intersect(HeLa_low_m6A_trans,HepG2_low_m6A_trans))
length(HeLa_HepG2_low_m6A_trans) # 89
HeLa_HEK293_low_m6A_trans <- unique(intersect(HeLa_low_m6A_trans,HEK293_low_m6A_trans))
length(HeLa_HEK293_low_m6A_trans) # 74
HEK293_HepG2_low_m6A_trans <- unique(intersect(HEK293_low_m6A_trans,HepG2_low_m6A_trans))
length(HEK293_HepG2_low_m6A_trans) # 51
c2t_long_dt_T_HeLa %>% nrow() ## 24225 background transcript

HeLa_HepG2_phyper <- phyper(89-1,354, 24225-354, 99, lower.tail = F)
HeLa_HepG2_phyper # 4.107468e-156
HeLa_HEK293T_phyper <- phyper(74-1,354, 24225-354, 82, lower.tail = F)
HeLa_HEK293T_phyper  # 1.516836e-129
HEK293_HepG2_phyper <- phyper(51-1,99, 24225-99, 82, lower.tail = F)
HEK293_HepG2_phyper # 7.0245e-107
library(eulerr)
Venn <- function(lst.A, lst.B, lst.C){
  A = lst.A %>% length()
  B = lst.B %>% length()
  C = lst.C %>% length()
  AB = intersect(lst.A, lst.B) %>% length()
  AC = intersect(lst.A, lst.C) %>% length() 
  BC = intersect(lst.B, lst.C) %>% length()
  ABC = intersect(lst.A, lst.B) %>% intersect(., lst.C) %>% length()
  euler(c(
    "A" = A - AB - AC + ABC,
    "B" = B - AB - BC + ABC,
    "C" = C - AC - BC + ABC,
    "A&B" = AB - ABC,
    "A&C" = AC - ABC,
    "B&C" = BC - ABC,
    "A&B&C" = ABC
  )) %>% return()
} 

vd <- Venn(HeLa_low_m6A_trans,HepG2_low_m6A_trans,HEK293_low_m6A_trans)
plot(vd,
     labels = list(
       labels = c('HeLa', 'HepG2', 'HEK293T'),
       col = "black", font = 2
     ), 
     edges = list(col="white",alpha=0),
     fills = list(fill = c("dodgerblue", "seagreen3", "orchid3"), alpha = 0.4),
     quantities = list(cex=1, col='black'),
     main = list(label=c("Lowly methylated transcripts in each celltype"),cex=2),
     legend = list(labels=c("A:HeLa","B:HepG2","C:HEK293T"),cex=1.2)
)
grid.text("P(AvsB)=4.1e-156", x=0.8, y=0.25, rot=1,gp=gpar(fontsize=13, col="black"))
grid.text("P(AvsC)=1.5e-129", x=0.8, y=0.2, rot=1,gp=gpar(fontsize=13, col="black"))
grid.text("P(BvsC)=7.0e-107", x=0.8, y=0.15, rot=1,gp=gpar(fontsize=13, col="black"))

######################
###  Fig4I
######################
HeLa_gene <- unique(HeLa_ratio_expression_dt$gene_name)
HepG2_gene <- unique(HepG2_ratio_expression_dt$gene_name)
HEK293_gene <- unique(HEK293_ratio_expression_dt$gene_name)
HeLa_setdiff_gene <- setdiff(HeLa_gene,c(HepG2_gene,HEK293_gene)) ### only HeLa gene
HeLa_setdiff_high_trans <- setdiff(HeLa_high_m6A_trans,c(HepG2_high_m6A_trans,HEK293_high_m6A_trans))
HeLa_setdiff_low_trans <- setdiff(HeLa_low_m6A_trans,c(HepG2_low_m6A_trans,HEK293_low_m6A_trans))

####  All Gene
c2t_long_dt_HeLa_low_sign <- c2t_long_dt_HeLa %>% filter(transcript_name %in% HeLa_setdiff_low_trans) ### filter transcript
c2t_long_dt_HEK293_low_sign <- c2t_long_dt_HEK293 %>% filter(transcript_name %in% HeLa_setdiff_low_trans)
c2t_long_dt_HepG2_low_sign <- c2t_long_dt_HepG2 %>% filter(transcript_name %in% HeLa_setdiff_low_trans)

c2t_long_dt_HeLa_high_sign <- c2t_long_dt_HeLa %>% filter(transcript_name %in% HeLa_setdiff_high_trans) ### filter transcript
c2t_long_dt_HEK293_high_sign <- c2t_long_dt_HEK293 %>% filter(transcript_name %in% HeLa_setdiff_high_trans)
c2t_long_dt_HepG2_high_sign <- c2t_long_dt_HepG2 %>% filter(transcript_name %in% HeLa_setdiff_high_trans)

m6a_long_dt_HeLa_low_sign <- m6a_long_dt_HeLa %>% filter(transcript_name %in% HeLa_setdiff_low_trans)
m6a_long_dt_HEK293_low_sign <- m6a_long_dt_HEK293 %>% filter(transcript_name %in% HeLa_setdiff_low_trans)
m6a_long_dt_HepG2_low_sign <- m6a_long_dt_HepG2 %>% filter(transcript_name %in% HeLa_setdiff_low_trans)

m6a_long_dt_HeLa_high_sign <- m6a_long_dt_HeLa %>% filter(transcript_name %in% HeLa_setdiff_high_trans)
m6a_long_dt_HEK293_high_sign <- m6a_long_dt_HEK293 %>% filter(transcript_name %in% HeLa_setdiff_high_trans)
m6a_long_dt_HepG2_high_sign <- m6a_long_dt_HepG2 %>% filter(transcript_name %in% HeLa_setdiff_high_trans)

##### low m6A modification transcript
transcript_low_trans_ratio_expression_dt <- bind_rows(mclapply(mc.cores=10,seq(HeLa_setdiff_low_trans),function(i){
  print(i)
  HeLa_ratio <- c2t_long_dt_HeLa_low_sign %>% filter(transcript_name == HeLa_setdiff_low_trans[i]) %>% dplyr::select(-transcript_name) %>% unlist() %>% as.numeric()
  HeLa_ratio_Mean <- mean(HeLa_ratio[HeLa_ratio!=0])
  HepG2_ratio <- c2t_long_dt_HepG2_low_sign %>% filter(transcript_name == HeLa_setdiff_low_trans[i]) %>% dplyr::select(-transcript_name) %>% unlist() %>% as.numeric()
  HepG2_ratio_Mean <- mean(HepG2_ratio[HepG2_ratio!=0])
  HEK293_ratio <- c2t_long_dt_HEK293_low_sign %>% filter(transcript_name == HeLa_setdiff_low_trans[i]) %>% dplyr::select(-transcript_name) %>% unlist() %>% as.numeric()
  HEK293_ratio_Mean <- mean(HEK293_ratio[HEK293_ratio!=0])
  HeLa_rowname <- c2t_long_dt_HeLa_low_sign %>% filter(transcript_name == HeLa_setdiff_low_trans[i]) %>% dplyr::select(-transcript_name) %>% select_if(function(x) any(x > 0)) %>% colnames()
  HeLa_expression <-  m6a_long_dt_HeLa_low_sign %>% filter(transcript_name == HeLa_setdiff_low_trans[i]) %>% dplyr::select(all_of(HeLa_rowname)) %>% as.numeric() ### Use ratio != 0 barcode
  HeLa_expression_Mean <- mean(HeLa_expression) 
  HepG2_rowname <- c2t_long_dt_HepG2_low_sign %>% filter(transcript_name == HeLa_setdiff_low_trans[i]) %>% dplyr::select(-transcript_name) %>% select_if(function(x) any(x > 0)) %>% colnames()
  HepG2_expression <-  m6a_long_dt_HepG2_low_sign %>% filter(transcript_name == HeLa_setdiff_low_trans[i]) %>% dplyr::select(all_of(HepG2_rowname)) %>% as.numeric() ### Use ratio != 0 barcode
  HepG2_expression_Mean <- mean(HepG2_expression) 
  HEK293_rowname <- c2t_long_dt_HEK293_low_sign %>% filter(transcript_name == HeLa_setdiff_low_trans[i]) %>% dplyr::select(-transcript_name) %>% select_if(function(x) any(x > 0)) %>% colnames()
  HEK293_expression <-  m6a_long_dt_HEK293_low_sign %>% filter(transcript_name == HeLa_setdiff_low_trans[i]) %>% dplyr::select(all_of(HEK293_rowname)) %>% as.numeric() ### Use ratio != 0 barcode
  HEK293_expression_Mean <- mean(HEK293_expression)
  HeLa_sub_dt <- tibble(transcript_name = HeLa_setdiff_low_trans[i],ratio = HeLa_ratio_Mean,expression = HeLa_expression_Mean,type = "HeLa", length = length(HeLa_rowname), group = "low_m6A_transcript")
  HepG2_sub_dt <- tibble(transcript_name = HeLa_setdiff_low_trans[i],ratio = HepG2_ratio_Mean,expression = HepG2_expression_Mean,type = "HepG2", length = length(HepG2_rowname), group = "low_m6A_transcript")
  HEK293_sub_dt <- tibble(transcript_name = HeLa_setdiff_low_trans[i],ratio = HEK293_ratio_Mean,expression = HEK293_expression_Mean,type = "HEK293T", length = length(HEK293_rowname), group = "low_m6A_transcript")
  merge_dt <- rbind(HeLa_sub_dt,HepG2_sub_dt,HEK293_sub_dt)
}
))

transcript_low_trans_ratio_expression_filter_dt <- transcript_low_trans_ratio_expression_dt %>% na.omit() %>% filter(length >= 10)

##### high m6A modification transcript
transcript_high_trans_ratio_expression_dt <- bind_rows(mclapply(mc.cores=10,seq(HeLa_setdiff_high_trans),function(i){
  print(i)
  HeLa_ratio <- c2t_long_dt_HeLa_high_sign %>% filter(transcript_name == HeLa_setdiff_high_trans[i]) %>% dplyr::select(-transcript_name) %>% unlist() %>% as.numeric()
  HeLa_ratio_Mean <- mean(HeLa_ratio[HeLa_ratio!=0])
  HepG2_ratio <- c2t_long_dt_HepG2_high_sign %>% filter(transcript_name == HeLa_setdiff_high_trans[i]) %>% dplyr::select(-transcript_name) %>% unlist() %>% as.numeric()
  HepG2_ratio_Mean <- mean(HepG2_ratio[HepG2_ratio!=0])
  HEK293_ratio <- c2t_long_dt_HEK293_high_sign %>% filter(transcript_name == HeLa_setdiff_high_trans[i]) %>% dplyr::select(-transcript_name) %>% unlist() %>% as.numeric()
  HEK293_ratio_Mean <- mean(HEK293_ratio[HEK293_ratio!=0])
  HeLa_rowname <- c2t_long_dt_HeLa_high_sign %>% filter(transcript_name == HeLa_setdiff_high_trans[i]) %>% dplyr::select(-transcript_name) %>% select_if(function(x) any(x > 0)) %>% colnames()
  HeLa_expression <-  m6a_long_dt_HeLa_high_sign %>% filter(transcript_name == HeLa_setdiff_high_trans[i]) %>% dplyr::select(all_of(HeLa_rowname)) %>% as.numeric() ### Use ratio != 0 barcode
  HeLa_expression_Mean <- mean(HeLa_expression) 
  HepG2_rowname <- c2t_long_dt_HepG2_high_sign %>% filter(transcript_name == HeLa_setdiff_high_trans[i]) %>% dplyr::select(-transcript_name) %>% select_if(function(x) any(x > 0)) %>% colnames()
  HepG2_expression <-  m6a_long_dt_HepG2_high_sign %>% filter(transcript_name == HeLa_setdiff_high_trans[i]) %>% dplyr::select(all_of(HepG2_rowname)) %>% as.numeric() ### Use ratio != 0 barcode
  HepG2_expression_Mean <- mean(HepG2_expression) 
  HEK293_rowname <- c2t_long_dt_HEK293_high_sign %>% filter(transcript_name == HeLa_setdiff_high_trans[i]) %>% dplyr::select(-transcript_name) %>% select_if(function(x) any(x > 0)) %>% colnames()
  HEK293_expression <-  m6a_long_dt_HEK293_high_sign %>% filter(transcript_name == HeLa_setdiff_high_trans[i]) %>% dplyr::select(all_of(HEK293_rowname)) %>% as.numeric() ### Use ratio != 0 barcode
  HEK293_expression_Mean <- mean(HEK293_expression)
  HeLa_sub_dt <- tibble(transcript_name = HeLa_setdiff_high_trans[i],ratio = HeLa_ratio_Mean,expression = HeLa_expression_Mean,type = "HeLa", length = length(HeLa_rowname), group = "high_m6A_transcript")
  HepG2_sub_dt <- tibble(transcript_name = HeLa_setdiff_high_trans[i],ratio = HepG2_ratio_Mean,expression = HepG2_expression_Mean,type = "HepG2", length = length(HepG2_rowname), group = "high_m6A_transcript")
  HEK293_sub_dt <- tibble(transcript_name = HeLa_setdiff_high_trans[i],ratio = HEK293_ratio_Mean,expression = HEK293_expression_Mean,type = "HEK293T", length = length(HEK293_rowname), group = "high_m6A_transcript")
  merge_dt <- rbind(HeLa_sub_dt,HepG2_sub_dt,HEK293_sub_dt)
}
))

transcript_high_trans_ratio_expression_filter_dt <- transcript_high_trans_ratio_expression_dt %>% na.omit() %>% filter(length >= 10)

#### merge 
transcript_merge_trans_ratio_expression_dt <- rbind(transcript_low_trans_ratio_expression_filter_dt,transcript_high_trans_ratio_expression_filter_dt) %>% 
  mutate(group = ifelse(group == "high_m6A_transcript","high","low"))
transcript_merge_trans_ratio_expression_dt <- transcript_merge_trans_ratio_expression_dt%>% inner_join(.,gtf %>% filter(type == "transcript") %>% dplyr::select(gene_name,transcript_name,transcript_id))
transcript_merge_trans_ratio_expression_group_dt <- transcript_merge_trans_ratio_expression_dt %>% group_by(transcript_id) %>% filter(n() ==3) %>% group_by(gene_name) %>% filter(n() >= 6) %>% arrange(gene_name)
filter_genename <- transcript_merge_trans_ratio_expression_group_dt %>% group_by(gene_name,group) %>% summarize(n=n()) %>% summarize(n=n()) %>% filter(n==2)
transcript_merge_trans_ratio_expression_group_dt <- transcript_merge_trans_ratio_expression_group_dt %>% filter(gene_name%in% filter_genename$gene_name)

### ratio
ggplot(data=transcript_merge_trans_ratio_expression_dt,aes(x=group,y=ratio))+
  stat_boxplot(geom = "errorbar",width=0.2)+
  geom_boxplot(aes(fill=group),width = 0.7, outlier.shape = NA)+
  facet_wrap(~type)+
  scale_fill_manual(values = c("#3288bd","#f46d43"))+ #"#85B22E","#5F80B4" "#8A77A5","#DBE196"
  theme_classic()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = 20, vjust = 0.5, hjust = 0.5),
        legend.position = "none")+  
  geom_signif(comparisons = list(c("high","low")),test = "wilcox.test",
              map_signif_level = F,
              textsize = 5,
              y_position = c(0.7))+
  theme(text = element_text(size = 15)) +
  coord_cartesian(ylim = c(0,0.8)) +
  labs(y=TeX(r"(Mean m${^6}$A ratio of transcripts)"),x = TeX(r"(Isoform m${^6}$A level)"))

######################
###  Fig4J
######################
### expression
ggplot(data=transcript_merge_trans_ratio_expression_dt,aes(x=group,y=expression))+
  geom_boxplot(aes(fill=group))+
  facet_wrap(~type)+
  scale_fill_manual(values = c("#85B22E","#5F80B4"))+ #"#85B22E","#5F80B4" "#8A77A5","#DBE196"
  geom_jitter(width = 0.2)+
  theme_classic()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = 20, vjust = 0.5, hjust = 0.5),

        legend.position = "none")+  
  geom_signif(comparisons = list(c("high","low")),test = "wilcox.test",
              map_signif_level = F,
              textsize = 5)+
  theme(text = element_text(size = 15)) +
  labs(y=TeX(r"(Mean of expression)"),x = TeX(r"(Isoform m${^6}$A level)"))

######################
###  Fig4M
######################
c2t_long_dt = fread("st_combined_ratio.mtx") 
c2t_long_dt_sub <- c2t_long_dt %>% filter(transId %in% c("ENST00000395119","ENST00000528610")) ## EEF1D
rowname = c2t_long_dt_sub$transId
c2t_long_mtx_sub = as.matrix(c2t_long_dt_sub[,-1])
rownames(c2t_long_mtx_sub) = rowname
c2t_long_mtx_sub <- c2t_long_mtx_sub[,colnames(c2t_long_mtx_sub) %in% barcode$X1]
sums <- colSums(c2t_long_mtx_sub)
c2t_long_mtx_sub <- c2t_long_mtx_sub[, sums != 0]

c2t_Long_Ratio_subset <- c2t_Long_Ratio #### subset ratio
c2t_Long_Ratio_subset@reductions$umap@cell.embeddings <- c2t_Long_Ratio_subset@reductions$umap@cell.embeddings[rownames(c2t_Long_Ratio_subset@reductions$umap@cell.embeddings) %in% colnames(c2t_long_mtx_sub),] 
FeaturePlot(c2t_Long_Ratio_subset,features = c("ENST00000395119","ENST00000528610"),reduction = "umap") ### EEF1D
FeaturePlot(c2t_Long_Ratio_expr,features = c("ENST00000395119","ENST00000528610"),reduction = "umap",keep.scale = "all" , order = T) ### EEF1D

######################
###  Fig4L
######################
### Violinplot 
subset_ratio_Seurat_1 <- subset(c2t_Long_Ratio_subset, subset = ENST00000395119 > 0)
subset_ratio_1 <- as.vector(subset_ratio_Seurat_1@assays$RNA@data[rownames(subset_ratio_Seurat_1@assays$RNA@data) %in% "ENST00000395119",])
subset_ratio_Seurat_2 <- subset(c2t_Long_Ratio_subset, subset = ENST00000528610 > 0)
subset_ratio_2 <- as.vector(subset_ratio_Seurat_2@assays$RNA@data[rownames(subset_ratio_Seurat_2@assays$RNA@data) %in% "ENST00000528610",])
subset_tibble <- rbind(tibble(transcript_name = "ENST00000395119", ratio = subset_ratio_1,group = "low"),tibble(transcript_name = "ENST00000528610", ratio = subset_ratio_2, group = "high"))
subset_tibble_dt <- subset_tibble %>% group_by(group) %>% mutate(High = quantile(ratio, 0.75),Med = median(ratio),Low = quantile(ratio, 0.25)) %>% ungroup() %>% mutate(transcript_name = factor(transcript_name,levels = c("ENST00000528610","ENST00000395119")))

ggplot(subset_tibble_dt,aes(x=transcript_name,y= ratio,fill = group))+
  geom_violin(trim = FALSE)+
  geom_linerange(aes(x = transcript_name, ymin = Low, ymax = High),
                 position = position_dodge(width = 0.9))+
  geom_point(aes(x = transcript_name, y = Med,),
             position = position_dodge(width = 0.9), size = 3)+
  scale_fill_manual(values = c("#3288bd","#f46d43"))+
  theme_classic()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = 20, vjust = 0.5, hjust = 0.5),
        legend.position = "none")+  
  geom_signif(comparisons = list(c("ENST00000528610", "ENST00000395119")),test = "wilcox.test",
              map_signif_level = function(p) sprintf("P = %.2g",p),
              textsize = 5, 
              y_position = c(1.15))+
  theme(text = element_text(size = 15)) +
  labs(y=TeX(r"(Mean m${^6}$A level of transcripts)"),x = NULL)+
  scale_y_continuous(limits=c(0,1.2),breaks = c(0.25,0.5,0.75,1)) +
  ggtitle("EEF1D")

######################
###  SFig7J
######################
high_m6A_trans <- "ENST00000528610"
low_m6A_trans <- "ENST00000395119"
subset_ratio_Seurat_1_HeLa <- subset(c2t_Long_Ratio_subset, subset = ENST00000528610 > 0, idents = "HeLa")
subset_ratio_Seurat_1_HEK293T <- subset(c2t_Long_Ratio_subset, subset = ENST00000528610 > 0, idents = "HEK293T")
subset_ratio_Seurat_1_HepG2 <- subset(c2t_Long_Ratio_subset, subset = ENST00000528610 > 0, idents = "HepG2")
subset_ratio_Seurat_2_HeLa <- subset(c2t_Long_Ratio_subset, subset = ENST00000395119 > 0, idents = "HeLa")
subset_ratio_Seurat_2_HEK293T <- subset(c2t_Long_Ratio_subset, subset = ENST00000395119 > 0, idents = "HEK293T")
subset_ratio_Seurat_2_HepG2 <- subset(c2t_Long_Ratio_subset, subset = ENST00000395119 > 0, idents = "HepG2")

subset_expr_1_HeLa <- as.vector(c2t_Long_Ratio_expr@assays$RNA@data[high_m6A_trans, colnames(c2t_Long_Ratio_expr@assays$RNA@data) %in% colnames(subset_ratio_Seurat_1_HeLa@assays$RNA@data)])
subset_expr_1_HEK293T <- as.vector(c2t_Long_Ratio_expr@assays$RNA@data[high_m6A_trans,colnames(c2t_Long_Ratio_expr@assays$RNA@data) %in% colnames(subset_ratio_Seurat_1_HEK293T@assays$RNA@data)])
subset_expr_1_HepG2 <- as.vector(c2t_Long_Ratio_expr@assays$RNA@data[high_m6A_trans,colnames(c2t_Long_Ratio_expr@assays$RNA@data) %in% colnames(subset_ratio_Seurat_1_HepG2@assays$RNA@data)])
subset_expr_2_HeLa <- as.vector(c2t_Long_Ratio_expr@assays$RNA@data[low_m6A_trans,colnames(c2t_Long_Ratio_expr@assays$RNA@data) %in% colnames(subset_ratio_Seurat_2_HeLa@assays$RNA@data)])
subset_expr_2_HEK293T <- as.vector(c2t_Long_Ratio_expr@assays$RNA@data[low_m6A_trans,colnames(c2t_Long_Ratio_expr@assays$RNA@data) %in% colnames(subset_ratio_Seurat_2_HEK293T@assays$RNA@data)])
subset_expr_2_HepG2 <- as.vector(c2t_Long_Ratio_expr@assays$RNA@data[low_m6A_trans,colnames(c2t_Long_Ratio_expr@assays$RNA@data) %in% colnames(subset_ratio_Seurat_2_HepG2@assays$RNA@data)])

subset_expr_tibble_cells <- rbind(tibble(transcript_name = high_m6A_trans, expr = subset_expr_1_HeLa,group1 = "high",group2 = "HeLa"),
                                  tibble(transcript_name = high_m6A_trans, expr = subset_expr_1_HEK293T,group1 = "high",group2 = "HEK293T"),
                                  tibble(transcript_name = high_m6A_trans, expr = subset_expr_1_HepG2,group1 = "high",group2 = "HepG2"),
                                  tibble(transcript_name = low_m6A_trans, expr = subset_expr_2_HeLa,group1 = "low",group2 = "HeLa"),
                                  tibble(transcript_name = low_m6A_trans, expr = subset_expr_2_HEK293T,group1 = "low",group2 = "HEK293T"),
                                  tibble(transcript_name = low_m6A_trans, expr = subset_expr_2_HepG2,group1 = "low",group2 = "HepG2"))
subset_expr_tibble_cells_dt <- subset_expr_tibble_cells %>% group_by(group1,group2) %>% mutate(High = quantile(expr, 0.75),Med = median(expr),Low = quantile(expr, 0.25)) %>% ungroup() %>% mutate(group2 = factor(group2,levels = c("HeLa","HEK293T","HepG2"))) %>% 
  mutate(transcript_name = factor(transcript_name,levels = c(high_m6A_trans,low_m6A_trans)))

compare_means(expr  ~ group1,subset_expr_tibble_cells_dt,method = "wilcox.test",paired = F) ###  0.25
ggplot(subset_expr_tibble_cells_dt,aes(x=transcript_name,y= expr,fill = group1))+
  geom_violin(trim = FALSE)+
  facet_wrap(~group2)+
  geom_linerange(aes(x = transcript_name, ymin = Low, ymax = High, group = group2),
                 position = position_dodge(width = 0.9))+
  geom_point(aes(x = transcript_name, y = Med, group = group2),
             position = position_dodge(width = 0.9), size = 3)+
  scale_fill_manual(values = c("#85B22E","#5F80B4"))+
  theme_classic()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = 20, vjust = 0.5, hjust = 0.5),
        legend.position = "none")+  
  geom_signif(comparisons = list(c(high_m6A_trans, low_m6A_trans)),test = "wilcox.test",
              map_signif_level = function(p) sprintf("P = %.2g",p),
              textsize = 5, 
              y_position = c(4.5))+
  theme(text = element_text(size = 15)) +
  labs(y=TeX(r"(Mean normalization expression of transcripts)"),x = NULL)+
  scale_y_continuous(limits=c(-0.5,5)) +
  ggtitle("EEF1D")


