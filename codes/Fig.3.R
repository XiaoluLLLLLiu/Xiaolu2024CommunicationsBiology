rm(list = ls())
library(ggtree)
library(ggplot2)
library(dplyr)
library(sva)
library(DESeq2)
library(pheatmap)
library(VennDiagram)
library(xlsx)
library(data.table)
library(cols4all)
library(factoextra)
library(RColorBrewer)

set.seed(406206)
setwd("submit_source_data/source data")

######################## Fig.3B ########################

tr <- read.tree("Fig.3B.txt")
ggtree(tr,branch.length = "none",layout = 'circular') -> tr.p
tr.p + geom_text2(aes(label=label)) +
  geom_hilight(node = 9, fill = "#08306B", alpha = 0.2) +
  geom_hilight(node = 11, fill = "#8FAADC", alpha = 0.2)+
  geom_tippoint(color = "#1F77B4FF", alpha = 0.5, size = 1.5) -> tr.p.f


### Comparative phylotranscriptomic analysis data processing

single.list <- fread("Orthogroups_SingleCopyOrthologues.txt",header = F) %>% .$V1
ortho.df <- fread("Orthogroups.tsv",header = T)

single.df <- ortho.df[ortho.df$Orthogroup %in% single.list ,] %>% as.data.frame()
single.df <- apply(single.df, c(1, 2), function(x) {
  gsub("\\.[pP]\\d+","", x)
}) %>% data.frame()
colnames(single.df) <- gsub("_last_orfs","",colnames(single.df))


quant.list <- readRDS("quant.list.rds")


###extract each species expressions
quant.list[grepl("Ba",names(quant.list))] %>% do.call(cbind,.) %>% dplyr::select(c(1,ends_with("NumReads")))-> ba.orig.counts
quant.list[grepl("Bm",names(quant.list))] %>% do.call(cbind,.) %>% dplyr::select(c(1,ends_with("NumReads"))) -> bm.orig.counts
quant.list[grepl("Gp",names(quant.list))] %>% do.call(cbind,.) %>% dplyr::select(c(1,ends_with("NumReads"))) -> gp.orig.counts
quant.list[grepl("Mc",names(quant.list))] %>% do.call(cbind,.) %>% dplyr::select(c(1,ends_with("NumReads"))) -> mc.orig.counts
quant.list[grepl("Mg",names(quant.list))] %>% do.call(cbind,.) %>% dplyr::select(c(1,ends_with("NumReads"))) -> mg.orig.counts
quant.list[grepl("Pv",names(quant.list))] %>% do.call(cbind,.) %>% dplyr::select(c(1,ends_with("NumReads"))) -> pv.orig.counts
quant.list[c(10:14,23:27)]%>% do.call(cbind,.) %>% dplyr::select(c(1,ends_with("NumReads"))) -> mussel.orig.counts

###change each species transcripts prefix
ba.orig.counts$`Ba-1.Name` %>% gsub("TRINITY","BA",.) -> ba.orig.counts$`Ba-1.Name` 
bm.orig.counts$`Bm-1.Name` %>% gsub("TRINITY","BM",.) -> bm.orig.counts$`Bm-1.Name`
gp.orig.counts$`Gp-1.Name` %>% gsub("TRINITY","GP",.) -> gp.orig.counts$`Gp-1.Name` 
mc.orig.counts$`Mc-1.Name` %>% gsub("TRINITY","MC",.) -> mc.orig.counts$`Mc-1.Name`
mg.orig.counts$`Mg-1.Name` %>% gsub("TRINITY","MG",.) -> mg.orig.counts$`Mg-1.Name` 
pv.orig.counts$`Pv-1.Name` %>% gsub("TRINITY","PV",.) -> pv.orig.counts$`Pv-1.Name` 

####change colnames for left_join species tpm and orthogroup
colnames(ba.orig.counts)[1] <- "Ba"
colnames(bm.orig.counts)[1] <- "Bm"
colnames(gp.orig.counts)[1] <- "Gp"
colnames(mc.orig.counts)[1] <- "Mc"
colnames(mg.orig.counts)[1] <- "Mg"
colnames(pv.orig.counts)[1] <- "Pv"
colnames(mussel.orig.counts)[1] <- "Mussel"

###alternate each transcripts ID to orthologous ID
left_join(ba.orig.counts,data.frame(Ba = single.df$Ba,Orthogroup = single.df$Orthogroup),by = "Ba") %>% na.omit() -> ba.orthogroup.rc
left_join(bm.orig.counts,data.frame(Bm = single.df$Bm,Orthogroup = single.df$Orthogroup),by = "Bm") %>% na.omit() -> bm.orthogroup.rc
left_join(gp.orig.counts,data.frame(Gp = single.df$Gp,Orthogroup = single.df$Orthogroup),by = "Gp") %>% na.omit() -> gp.orthogroup.rc
left_join(mc.orig.counts,data.frame(Mc = single.df$Mc,Orthogroup = single.df$Orthogroup),by = "Mc") %>% na.omit() -> mc.orthogroup.rc
left_join(mg.orig.counts,data.frame(Mg = single.df$Mg,Orthogroup = single.df$Orthogroup),by = "Mg") %>% na.omit() -> mg.orthogroup.rc
left_join(pv.orig.counts,data.frame(Pv = single.df$Pv,Orthogroup = single.df$Orthogroup),by = "Pv") %>% na.omit() -> pv.orthogroup.rc
left_join(mussel.orig.counts,data.frame(Mussel = single.df$Mussel,Orthogroup = single.df$Orthogroup),by = "Mussel") %>% na.omit() -> mussel.orthogroup.rc

### merge orthogroup rc matrix
left_join(mussel.orthogroup.rc, mc.orthogroup.rc, by = "Orthogroup") %>% 
  left_join(., mg.orthogroup.rc, by = "Orthogroup") %>% 
  left_join(., pv.orthogroup.rc, by = "Orthogroup") %>% 
  left_join(., ba.orthogroup.rc, by = "Orthogroup") %>% 
  left_join(., bm.orthogroup.rc, by = "Orthogroup") %>% 
  left_join(., gp.orthogroup.rc, by = "Orthogroup") -> merged_data


merged_data %>% dplyr::select(c("Orthogroup",ends_with("NumReads"))) -> merged_data
merged_data <- data.frame(row.names = merged_data$Orthogroup,merged_data[,-1])
colnames(merged_data) <- colnames(merged_data) %>% gsub("[.]NumReads","",.) %>% gsub("[.]","-",.)

#####prepared the length matrix
quant.list[grepl("Ba",names(quant.list))] %>% do.call(cbind,.) %>% dplyr::select(c(1,ends_with("EffectiveLength")))-> ba.orig.el
quant.list[grepl("Bm",names(quant.list))] %>% do.call(cbind,.) %>% dplyr::select(c(1,ends_with("EffectiveLength"))) -> bm.orig.el
quant.list[grepl("Gp",names(quant.list))] %>% do.call(cbind,.) %>% dplyr::select(c(1,ends_with("EffectiveLength"))) -> gp.orig.el
quant.list[grepl("Mc",names(quant.list))] %>% do.call(cbind,.) %>% dplyr::select(c(1,ends_with("EffectiveLength"))) -> mc.orig.el
quant.list[grepl("Mg",names(quant.list))] %>% do.call(cbind,.) %>% dplyr::select(c(1,ends_with("EffectiveLength"))) -> mg.orig.el
quant.list[grepl("Pv",names(quant.list))] %>% do.call(cbind,.) %>% dplyr::select(c(1,ends_with("EffectiveLength"))) -> pv.orig.el
quant.list[c(10:14,23:27)]%>% do.call(cbind,.) %>% dplyr::select(c(1,ends_with("EffectiveLength"))) -> mussel.orig.el

ba.orig.el$`Ba-1.Name` %>% gsub("TRINITY","BA",.) -> ba.orig.el$`Ba-1.Name` 
bm.orig.el$`Bm-1.Name` %>% gsub("TRINITY","BM",.) -> bm.orig.el$`Bm-1.Name`
gp.orig.el$`Gp-1.Name` %>% gsub("TRINITY","GP",.) -> gp.orig.el$`Gp-1.Name` 
mc.orig.el$`Mc-1.Name` %>% gsub("TRINITY","MC",.) -> mc.orig.el$`Mc-1.Name`
mg.orig.el$`Mg-1.Name` %>% gsub("TRINITY","MG",.) -> mg.orig.el$`Mg-1.Name` 
pv.orig.el$`Pv-1.Name` %>% gsub("TRINITY","PV",.) -> pv.orig.el$`Pv-1.Name` 

colnames(ba.orig.el)[1] <- "Ba"
colnames(bm.orig.el)[1] <- "Bm"
colnames(gp.orig.el)[1] <- "Gp"
colnames(mc.orig.el)[1] <- "Mc"
colnames(mg.orig.el)[1] <- "Mg"
colnames(pv.orig.el)[1] <- "Pv"
colnames(mussel.orig.el)[1] <- "Mussel"

left_join(ba.orig.el,data.frame(Ba = single.df$Ba,Orthogroup = single.df$Orthogroup),by = "Ba") %>% na.omit() -> ba.orthogroup.el
left_join(bm.orig.el,data.frame(Bm = single.df$Bm,Orthogroup = single.df$Orthogroup),by = "Bm") %>% na.omit() -> bm.orthogroup.el
left_join(gp.orig.el,data.frame(Gp = single.df$Gp,Orthogroup = single.df$Orthogroup),by = "Gp") %>% na.omit() -> gp.orthogroup.el
left_join(mc.orig.el,data.frame(Mc = single.df$Mc,Orthogroup = single.df$Orthogroup),by = "Mc") %>% na.omit() -> mc.orthogroup.el
left_join(mg.orig.el,data.frame(Mg = single.df$Mg,Orthogroup = single.df$Orthogroup),by = "Mg") %>% na.omit() -> mg.orthogroup.el
left_join(pv.orig.el,data.frame(Pv = single.df$Pv,Orthogroup = single.df$Orthogroup),by = "Pv") %>% na.omit() -> pv.orthogroup.el
left_join(mussel.orig.el,data.frame(Mussel = single.df$Mussel,Orthogroup = single.df$Orthogroup),by = "Mussel") %>% na.omit() -> mussel.orthogroup.el

left_join(mussel.orthogroup.el, mc.orthogroup.el, by = "Orthogroup") %>% 
  left_join(., mg.orthogroup.el, by = "Orthogroup") %>% 
  left_join(., pv.orthogroup.el, by = "Orthogroup") %>% 
  left_join(., ba.orthogroup.el, by = "Orthogroup") %>% 
  left_join(., bm.orthogroup.el, by = "Orthogroup") %>% 
  left_join(., gp.orthogroup.el, by = "Orthogroup") -> merged_length

merged_length %>% dplyr::select(c("Orthogroup",ends_with("EffectiveLength"))) -> merged_length
merged_length <- data.frame(row.names = merged_length$Orthogroup,merged_length[,-1])
colnames(merged_length) <- colnames(merged_length) %>%  gsub("[.]","-",.)

###raw merged data pca   ---- with strong effect in each species
condition <- factor(c(rep('LP',5),rep('Con',5),(rep('Mc',3)),rep('Mg',3),rep('Pv',2),rep('Ba',3),rep('Bm',3),rep('Gp',3)),levels = c("LP","Con","Mc","Mg","Pv","Ba",'Bm',"Gp"))
env <- factor(c(rep("Deep",5),rep("Shallow",5),rep("Shallow",8),rep("Deep",9)),levels = c("Deep","Shallow"))
platform <- factor(c(rep("NovaSeq",13),rep("Hiseq",3),rep("Genome_Analyzer_IIx",2),rep("NovaSeq",3),rep("Hiseq",6)))
project <- factor(c(rep("Self",10),rep("PRJNA770460",3),rep("PRJNA824625",3),rep("PRJNA254094",2),
                    rep("PRJEB34925",3),rep("PRJNA360359",6)))
species <- factor(c(rep('Mussel',10),(rep('Mc',3)),rep('Mg',3),rep('Pv',2),rep('Ba',3),rep('Bm',3),rep('Gp',3)),levels = c("Mussel","Mc","Mg","Pv","Ba",'Bm',"Gp"))
resource <- factor(c(rep('Self',10),rep('Shallow',8),rep('Deep',9)),levels = c('Self','Shallow','Deep'))
colData <- data.frame(row.names=colnames(merged_data), condition,platform,project,env,species,resource)


### raw readcount input to combat-seq
count_m <- round(merged_data) %>% as.matrix()
combat_Expr <- ComBat_seq(count_m,batch = colData$project,group = colData$env) 



####explore public DESeq2, re-build the DESeq2 subject
combat_expr_pub <- combat_Expr[,-1:-10]
colData_pub <- colData[-1:-10,]

dds <- DESeqDataSetFromMatrix(combat_expr_pub,
                              colData = colData_pub,
                              design = ~ resource)
dds$resource <- relevel(dds$resource,ref = "Shallow")
dds <- dds[rowSums(counts(dds)) > 1, ] 
colnames(merged_length) <- colnames(merged_length) %>% gsub("-EffectiveLength","",.)
public_merged_length <- merged_length[,-1:-10]
public_merged_length <- as.matrix(public_merged_length)
normFactors <- public_merged_length / exp(rowMeans(log(public_merged_length)))
normalizationFactors(dds) <- normFactors


dds <- DESeq(dds)
res <- results(dds, contrast=c("resource","Deep","Shallow"))
baseMean_Deep <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$resource == "Deep"])
baseMean_Shallow <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$resource == "Shallow"])
res <- cbind(baseMean_Deep, baseMean_Shallow, as.data.frame(res))
res$padj[is.na(res$padj)]  <- 1
res <- as.data.frame(res[order(res$pvalue),])

res_sig <- res[abs(res$log2FoldChange) > 1 & res$padj < 0.05,]
res_sig <- res_sig %>% mutate(sig_trend = case_when(log2FoldChange > 1 ~ "Up",
                                                    log2FoldChange < -1 ~ "Down",
                                                    TRUE ~ "No change"))


######################## Fig.3C ########################
public_res <- read.xlsx("Fig.3C.xlsx",sheetIndex = 1)
public_res$Trend <- factor(public_res$Trend,levels = c("Up","Down","None"))

ma_p <- ggplot() +
  geom_point(data = public_res[public_res$Trend == "None", ],
             aes(x = log10(baseMean), y = log2FoldChange, color = Trend),
             size = 1, alpha = 0.5) +  
  geom_point(data = public_res[public_res$Trend %in% c("Up", "Down"), ],
             aes(x = log10(baseMean), y = log2FoldChange, color = Trend),
             size = 1, alpha = 1) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("#F2F2F2","#C5243D", "#6191BB" )) +
  scale_alpha(range = c(0.1, 0.5))

######################## Fig.3D ########################

ogs_go <- read.xlsx("Fig.3D.xlsx",sheetIndex = 1)
ogs_go$log10padjust <- ifelse(ogs_go$Trend == 'Up', -log10(ogs_go$p.adjust), -(-log10(ogs_go$p.adjust))) ## for better illustration

level <- c(ogs_go$Description[ogs_go$Trend=='Up'],
           ogs_go$Description[ogs_go$Trend=='Down'])
ogs_go$Description <- factor(ogs_go$Description, levels = rev(level))

ogs_go.p <- ggplot(ogs_go,aes(x = log10padjust, y = Description)) +
  scale_y_discrete( labels=function(x) str_wrap(x, width=20))+
  theme_classic() +
  geom_point(aes( colour = log10padjust,size = Count)) + 
  geom_col(aes(fill = log10padjust), width = 0.1)+
  scale_size(range = c(2, 7)) +
  scale_color_continuous_c4a_div("red_blue_diverging",mid = 0,reverse = T) +
  scale_fill_continuous_c4a_div("red_blue_diverging",mid = 0,reverse = T)+
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 11),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        plot.margin = margin(t = 5.5, r = 10, l = 5.5, b = 5.5))+
  
  ylab('')


######################## Fig.3E ########################
pca_dat <- prcomp(t(combat_Expr), scale = TRUE,retx=T) 
pca_Variance <- round((pca_dat$sdev^2/sum(pca_dat$sdev^2)) * 100,2)
adjust.pca <- fviz_pca_ind(pca_dat, repel = T,geom.ind = c('point','text'),
                           col.ind=colData$env, mean.point=F,  legend.title="Groups", 
                           ellipse.type="confidence", ellipse.level=0.9,addEllipses = T,
                           palette = c("#F8766D", "#1f497d","orange"))+
  #palette = c("8eb4e3", "#00BA38", "#619CFF","purple","orange","black","royalblue","firebrick"))+ 
  theme(panel.border = element_rect(fill=NA,color="black", linewidth = 1, linetype="solid"))+#加个边框
  theme_classic()+
  labs(x = paste0("PC1 (",pca_Variance[1],"% explained var.)"), y = paste0("PC2 (",pca_Variance[2],"% explained var.)" ))+
  theme_test()+
  geom_point(aes(shape = colData$species,color = colData$env),size = 3)+
  scale_shape_manual(values = c("Ba" = 15,"Bm" = 16, "Gp" = 17, "Mc" = 18,"Mg" = 19, "Mussel" = 20, "Pv" = 4))+
  ggtitle('PCA plot')


#### explore the samples' correlation, using full data
ddsFull <- DESeqDataSetFromMatrix(combat_Expr,colData = colData,design = ~ resource)
rld<- rlog(ddsFull, blind = FALSE)

######################## Fig.3F ########################
tmm.cor <- stats::cor(assay(rld),method = "spearman")

pheatmap(as.matrix(tmm.cor),col = colorRampPalette(brewer.pal(9, "Blues")) (255),
         border_color = NA,treeheight_row = 0,treeheight_col = 10,cutree_rows = 2,cutree_cols = 2) -> cor.p


######################## self data, all SCOGs, re-build the DESeq2 subject

combat_expr_mussel <- combat_Expr[,1:10]
colData_mussel <- colData[1:10,]

dds2 <- DESeqDataSetFromMatrix(combat_expr_mussel,
                               colData = colData_mussel,
                               design = ~ env)
dds2$env <- relevel(dds2$env,ref = "Shallow")
dds <- dds[rowSums(counts(dds)) > 1, ] 
colnames(merged_length) <- colnames(merged_length) %>% gsub("-EffectiveLength","",.)
self_merged_length <- merged_length[,1:10]
self_merged_length <- as.matrix(self_merged_length)
normFactors <- self_merged_length / exp(rowMeans(log(self_merged_length)))
normalizationFactors(dds2) <- normFactors


dds2 <- DESeq(dds2)
res2 <- results(dds2, contrast=c("env","Deep","Shallow"))
baseMean_LP <- rowMeans(counts(dds2, normalized=TRUE)[,colData(dds2)$env == "Deep"])
baseMean_Con <- rowMeans(counts(dds2, normalized=TRUE)[,colData(dds2)$env == "Shallow"])
res2 <- cbind(baseMean_LP, baseMean_Con, as.data.frame(res2))
res2$padj[is.na(res2$padj)]  <- 1
res2 <- as.data.frame(res2[order(res2$pvalue),])

res2_sig <- res2[abs(res2$log2FoldChange) > 1 & res2$padj < 0.05,]
res2_sig <- res2_sig %>% mutate(sig_trend = case_when(log2FoldChange > 1 ~ "Up",
                                                      log2FoldChange < -1 ~ "Down",
                                                      TRUE ~ "No change"))

######################## Fig.3G ########################
res$scogs <- row.names(res)
res2$scogs <- row.names(res2)
merged_double <- left_join(res,res2,by = "scogs")
print(colSums(is.na(merged_double)))
merged_double$log2FoldChange.y[is.na(merged_double$log2FoldChange.y)] <- 0
#write.csv(merged_double,paste0(dir,"comparative/OrthoFinder_Res_new/x-y_plot_final.csv"),row.names = F)

merged_double <- merged_double %>% dplyr::select("scogs","log2FoldChange.x","log2FoldChange.y")
merged_double <- merged_double %>%
  mutate(Trend = case_when(
    log2FoldChange.x > 0 & log2FoldChange.y > 0 ~ "Same Up",
    log2FoldChange.x < 0 & log2FoldChange.y < 0 ~ "Same Down",
    log2FoldChange.x * log2FoldChange.y <= 0 ~ "Reverse",
    TRUE ~ "Other"
  )) %>% na.omit()


merged_double$Trend <- factor(merged_double$Trend,levels = c("Same Up","Same Down","Reverse"))
cor.test(merged_double$log2FoldChange.x,merged_double$log2FoldChange.y) -> x.y.cor
print(x.y.cor$estimate)

xy.cor.p <- ggplot(merged_double,aes(x = log2FoldChange.y,y = log2FoldChange.x,color = Trend))+
  geom_point(size = 0.5)+
  ylim(-6,6)+
  xlim(-6,6)+
  theme_classic()+
  scale_color_manual(values = c("#C5243D","#6191BB","#d8d8d8"))+ 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_smooth(method = "lm", color = "Black",size=0.2)+
  annotate("text", x = -3.5, y = 6, label = expression(italic("R = 0.821")))+
  xlab("SCOGs log2FoldChange of LT vs Con")+
  ylab("SCOGs log2FoldChange of Public Species (Deep vs Shallow)")


######################## Fig.S6 ########################
descogs_xy <- merged_double[merged_double$scogs %in% row.names(res_sig),]
descogs_xy <- descogs_xy[descogs_xy$scogs %in% row.names(res2_sig),]

descogs_xy$Trend <- factor(descogs_xy$Trend,levels = c("Same Up","Same Down","Reverse"))
cor.test(descogs_xy$log2FoldChange.x,descogs_xy$log2FoldChange.y)
xy.cor.p <- ggplot(descogs_xy,aes(x = log2FoldChange.y,y = log2FoldChange.x,color = Trend))+
  geom_point(size = 0.5)+
  ylim(-6,6)+
  xlim(-10,10)+
  theme_classic()+
  scale_color_manual(values = c("#C5243D","#6191BB","#d8d8d8"))+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_smooth(method = "lm", color = "Black",size=0.2)+
  annotate("text", x = -3.5, y = 6, label = expression(italic("R = 0.884145")))+
  xlab("Shared DESCOGs log2FC of LT vs Con")+
  ylab("Shared DESCOGs log2FC of Public Species (Deep vs Shallow)")


######################## Fig.3H ########################

venn.up.p <- venn.diagram(
  x = list(rownames(res_sig[res_sig$sig_trend == "Up",]),
           rownames(res2_sig[res2_sig$sig_trend == "Up",])),
  category.names = c("Deep vs Shallow" , "LT vs Con"),
  filename =NULL ,output=F,
  resolution = 400,fill = c("#C00000","#DC8EAA"),
  lwd = 2);grid.draw(venn.up.p)

mat.up <- matrix(c(as.numeric(venn.up.p[[7]]$label),
                   as.numeric(venn.up.p[[5]]$label),
                   as.numeric(venn.up.p[[6]]$label),
                   4878 - sum(as.numeric(venn.up.p[[5]]$label),as.numeric(venn.up.p[[6]]$label)) + as.numeric(venn.up.p[[7]]$label)), nrow=2)
up.p_value <- fisher.test(mat.up)$p.value
print(up.p_value)


venn.dn.p <- venn.diagram(
  x = list(rownames(res_sig[res_sig$sig_trend == "Down",]),
           rownames(res2_sig[res2_sig$sig_trend == "Down",])),
  category.names = c("Deep vs Shallow" , "LT vs Con"),
  filename =NULL ,output=F,
  resolution = 400,fill = c("#B3D5ED","#2F5597"),
  lwd = 2);grid.draw(venn.dn.p)


mat.dn <- matrix(c(as.numeric(venn.dn.p[[7]]$label),
                   as.numeric(venn.dn.p[[5]]$label),
                   as.numeric(venn.dn.p[[6]]$label),
                   4878 - sum(as.numeric(venn.dn.p[[5]]$label),as.numeric(venn.dn.p[[6]]$label)) + as.numeric(venn.dn.p[[7]]$label)), nrow=2)
dn.p_value <- fisher.test(mat.dn)$p.value
print(dn.p_value)

