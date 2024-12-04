rm(list = ls())
library(pheatmap)
library(data.table)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(xlsx)
library(cols4all)
library(ggpubr)
library(smplot2)
library(magrittr)
library(ggalluvial)

set.seed(406206)
setwd("submit_source_data/source data")

all_genes_tmm <- read.table("salmon_abundance_estimates_to_matrix.isoform.TMM.EXPR.matrix",header = T,row.names = 1,sep = "\t")
colnames(all_genes_tmm) <- gsub("[.]","-",colnames(all_genes_tmm))
con_mean <- rowMeans(all_genes_tmm[,1:5])
st_mean <- rowMeans(all_genes_tmm[,6:10])
lt_mean <- rowMeans(all_genes_tmm[,11:15])
all_tmm_mean <- data.frame(gene = rownames(all_genes_tmm),con_mean = con_mean,st_mean = st_mean, lt_mean = lt_mean)

######################## Fig.2B ########################

degs.total <- read.table("diffExpr.P5e-2_C1.matrix",header = T,row.names = 1)
colnames(degs.total) <- gsub("[.]","-",colnames(degs.total))

degs.log <- log(degs.total+1)

rowMeans(degs.log %>% dplyr::select(starts_with("Con"))) ->conmean
rowMeans(degs.log %>% dplyr::select(starts_with("ST"))) ->stmean
rowMeans(degs.log %>% dplyr::select(starts_with("LT"))) ->ltmean
data.frame(row.names = row.names(degs.log),Con = conmean,ST = stmean,LT = ltmean) -> degs.mean


degs.phm <- pheatmap::pheatmap(degs.mean,scale = "row",show_rownames = F,treeheight_row = 0,
                   treeheight_col = 10,cluster_cols = F,cutree_rows = 4,gaps_col = 1:2,
                   color =colorRampPalette(rev(brewer.pal(5,"RdBu")))(500),show_colnames = T,
                   border_color = "black",cellwidth = 40)

######################## Fig.2B ########################

go.enrichment.res <- read.xlsx("Fig.2C.xlsx",sheetIndex = 1)
go.enrichment.res$Group %<>% factor(., levels = c("ST vs Con","LT vs Con","LT vs ST"))
go.enrichment.res$Category <- factor(go.enrichment.res$Category,levels = c("Homeostasis","Immune","Mitochondrial metabolic processes"))
go.enrichment.res$Displayed_Rank <- go.enrichment.res$Displayed_Rank %>% as.numeric()

## for better display, change the up-regulated terms' log2p * -1

go.enrichment.res <- go.enrichment.res %>%
  mutate(log2p = ifelse(Trend == "Up", log2p * -1, log2p))


go.plot <- ggplot(data = go.enrichment.res,aes(x = reorder(Description,Displayed_Rank,decreasing = T), y = Group,color = log2p))+ #
  geom_point(data = go.enrichment.res,aes(size = GeneRatio))+
  scale_color_continuous_c4a_div("red_blue_diverging",reverse = T)+
  theme_bw()+
  labs(x="GO terms",y = " Group")+
  theme(axis.line = element_line(colour = "black"), 
        axis.text = element_text(color = "black",size = 10), #,angle = 90, hjust = 1
        legend.text = element_text(size = 9),
        legend.title=element_text(size=9),
        panel.grid.major = element_blank(),
        axis.title.x = element_text(size = 11),
        axis.text.x = element_text(angle = 90,hjust = 0),
        axis.title.y = element_text(size = 11),
        strip.text.x = element_text(size=8, color = "white"),
        strip.background = element_rect(fill="black"))+ 
  coord_flip()


######################## Fig.2D ########################
## GO enrichment results from all comparisons

go.all.res <- readRDS("C:/Users/amin/OneDrive/mussel/submit/add_metagenome/for cb/revision3/Supplementary_files/submit_source_data/source data/go.all.res.rds")

go.enrichment.res$dot_group <- paste0(go.enrichment.res$Group,"_",go.enrichment.res$Trend) %>% gsub(" ","_",.)
care_genes_df_list <- list()

### extract gene in care Go terms (Fig.2C) and their expressions

for (i in 1:length(go.enrichment.res$ID)){
  
  care_genes <- go.all.res[[go.enrichment.res[i,"dot_group"]]] %>%
    dplyr::filter(ID == go.enrichment.res[i,1]) %>%
    .$geneID %>% strsplit(.,"/") %>% unlist() 
  care_genes_df <- data.frame(ID = rep(go.enrichment.res[i,1],length(care_genes)),
                              category = rep(go.enrichment.res[i,"Category"],length(care_genes)),
                              gene = care_genes)
  care_genes_df_list[[i]] <- care_genes_df
  
}

care_df <- do.call(rbind,care_genes_df_list) %>% unique
care_tmm_df <- left_join(care_df,all_tmm_mean)
colnames(care_tmm_df)[4:6] <- c("Con","ST","LT")

### 3 GO categories' genes boxplot

# Homeostasis

care_pathway <- care_tmm_df %>% dplyr::filter(category == "Homeostasis") %>% .[,3:6] %>% unique
df_long <- reshape2::melt(care_pathway, id.vars = "gene", variable.name = "Condition", value.name = "Value")

Homeostasis.p <- ggplot(df_long, aes(x = Condition, y = log2(Value+1), fill = Condition,color = Condition)) +
  sm_raincloud(boxplot.params = list(outlier.shape = NA),sep_level = 1)+
  stat_compare_means(comparisons = list(c("Con", "ST"), c("ST", "LT"), c("Con", "LT")), method = "wilcox.test") +   #,paired = T
  scale_color_discrete_c4a_cat("blue") +
  scale_fill_discrete_c4a_cat("blue")+
  theme_classic2()+
  theme(axis.title.x = element_blank()) +
  labs(y = "log2(TPM)")

# Immune

care_pathway <- care_tmm_df %>% dplyr::filter(category == "Immune") %>% .[,3:6] %>% unique
df_long <- reshape2::melt(care_pathway, id.vars = "gene", variable.name = "Condition", value.name = "Value")

Immune.p <- ggplot(df_long, aes(x = Condition, y = log2(Value+1), fill = Condition,color = Condition)) +
  sm_raincloud(boxplot.params = list(outlier.shape = NA),sep_level = 1)+
  stat_compare_means(comparisons = list(c("Con", "ST"), c("ST", "LT"), c("Con", "LT")), method = "wilcox.test") +   #
  scale_color_discrete_c4a_cat("green") +
  scale_fill_discrete_c4a_cat("green")+
  theme_classic2()+
  theme(axis.title.x = element_blank()) +
  labs(y = "log2(TPM)")

# Mitochondrial metabolic processes


care_pathway <- care_tmm_df %>% dplyr::filter(category == "Mitochondrial metabolic processes") %>% .[,3:6] %>% unique
df_long <- reshape2::melt(care_pathway, id.vars = "gene", variable.name = "Condition", value.name = "Value")

mmp.p <- ggplot(df_long, aes(x = Condition, y = log2(Value+1), fill = Condition,color = Condition)) +
  sm_raincloud(boxplot.params = list(outlier.shape = NA),sep_level = 1)+
  stat_compare_means(comparisons = list(c("Con", "ST"), c("ST", "LT"), c("Con", "LT")), method = "wilcox.test") +   #,paired = T
  scale_color_discrete_c4a_cat("orange") +
  scale_fill_discrete_c4a_cat("orange")+
  theme_classic2()+
  theme(axis.title.x = element_blank()) +
  labs(y = "log2(TPM)")


######################## Fig.2E ########################

degs_phm_clusters<- fread("Fig.2E.txt") %>% .[,-2:-4]
colnames(degs_phm_clusters)[1] <- "gene"

care_tmm_cluster <- left_join(care_tmm_df,degs_phm_clusters)
care_tmm_cluster[,c(2,7)] %>% table() %>% data.frame() -> care_freq
care_freq$Cluster <- paste0("C",care_freq$Cluster)

## Alluvium plot

care_freq$category <- factor(care_freq$category,levels = c("Homeostasis","Immune","Mitochondrial metabolic processes"))

p <- ggplot(care_freq,aes(y = Freq, axis1 = category, axis2 = Cluster)) +
  scale_x_discrete(limits = c("GO terms", "Cluster"), expand = c(.2, .05)) +
  geom_alluvium(aes(fill = Cluster)) +
  geom_stratum() + 
  scale_fill_manual(values = c("C1" = "#A30545",
                               "C2" = "#7ECBA4", 
                               "C3" = "#F26E42",
                               "C4" = "#4A65AE")) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), angle = 90) +
  theme_classic()+
  theme(
    legend.position = 'none',
    axis.title.x = element_blank(), 
    axis.title.y = element_blank(), 
    axis.text.x = element_blank(), 
    axis.text.y = element_blank(), 
    axis.ticks = element_blank(), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank() 
  )

