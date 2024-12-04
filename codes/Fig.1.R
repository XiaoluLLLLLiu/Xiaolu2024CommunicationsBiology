rm(list = ls())
library(factoextra)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(xlsx)

set.seed(406206)
setwd("submit_source_data/source data")


all_genes_tmm <- read.table("salmon_abundance_estimates_to_matrix.isoform.TMM.EXPR.matrix",header = T,row.names = 1,sep = "\t")
colnames(all_genes_tmm) <- gsub("[.]","-",colnames(all_genes_tmm))
con_mean <- rowMeans(all_genes_tmm[,1:5])
st_mean <- rowMeans(all_genes_tmm[,6:10])
lt_mean <- rowMeans(all_genes_tmm[,11:15])
all_tmm_mean <- data.frame(gene = rownames(all_genes_tmm),con_mean = con_mean,st_mean = st_mean, lt_mean = lt_mean)

all_genes_tmm_t <- t(all_genes_tmm)
all_genes_tmm_use <- all_genes_tmm_t[,colSums(all_genes_tmm_t) > 10] %>% data.frame()

######################## Fig.1B ########################
rownames(all_genes_tmm_use) <- rownames(all_genes_tmm_use) %>% 
  gsub("SCK","Con",.) %>% gsub("SP","ST",.) %>% gsub("LP","LT",.)
pca_dat <- prcomp(all_genes_tmm_use, scale = TRUE,retx=T) 
group = factor(c(rep('Con',5),rep('ST',5),rep('LT',5)), level = c("Con", "ST", "LT"))

summary(pca_dat)
pca_Variance <- round((pca_dat$sdev^2/sum(pca_dat$sdev^2)) * 100,2)

p <- fviz_pca_ind(pca_dat, repel = T,geom.ind = c('point','text'),
                  col.ind=group, mean.point=F, addEllipses = T, legend.title="Groups", 
                  ellipse.type="confidence", ellipse.level=0.9, 
                  palette = c("#8FAADC", "#A9D18E", "#203864"))+ 
  theme(panel.border = element_rect(fill=NA,color="black", linewidth = 1, linetype="solid"))+#加个边框
  theme_classic()+
  labs(x = paste0("PC1 (",pca_Variance[1],"% explained var.)"), y = paste0("PC2 (",pca_Variance[2],"% explained var.)" ))+
  theme_test()+
  ggtitle('PCA plot')


######################## Fig.1C ########################
degs.total <- read.table("diffExpr.P5e-2_C1.matrix",header = T,row.names = 1)
colnames(degs.total) <- gsub("[.]","-",colnames(degs.total))
colnames(degs.total)[1:5] <- paste0("Con-",1:5)
tmm.cor <- stats::cor(degs.total,method = "spearman")

pheatmap(as.matrix(tmm.cor),col = colorRampPalette(brewer.pal(9, "Blues")) (255),
         border_color = NA,treeheight_row = 0,treeheight_col = 10,cutree_rows = 3,cutree_cols = 3) -> cor.p

######################## Fig.1D ########################
degs.num <- read.xlsx("Fig.1D.xlsx",sheetIndex = 1,header = T)
degs.num <- reshape2::melt(degs.num,variable.name = "Group",value.name = "Num")
degs.num$NA. <- factor(degs.num$NA.,levels = c("ST vs Con", "LT vs Con","LT vs ST"))

degs.stats.p <- ggplot(degs.num,aes(x = NA., y = Num,fill = Group))+
  geom_bar(stat="identity",position = position_dodge(width = 0.8), width = 0.6)+
  geom_text(aes(label = Num),position=position_dodge(width = 0.5),size = 3,vjust = -0.2)+
  theme_classic()+
  scale_fill_manual(values = c("#E3170D","#1E90FF"))+
  labs(x = "Comparsion",y = "", title = "DEGs stats")

