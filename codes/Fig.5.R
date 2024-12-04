rm(list=ls())
library(xlsx)
library(dplyr)
library(metaphlanToPhyloseq)
library(phyloseq)
library(ape)
library(reshape2)
library(stringr)
library(ggplot2)
library(ggpubr)
library(amplicon)
library(ggtree)
library(ggstance)
library(tibble)
library(preprocessCore)
library(data.table)
library(WGCNA)
library(pheatmap)
library(cols4all)
library(colortools)

set.seed(406206)
setwd("submit_source_data/source data")


#### edited fuction
plot_bar2 <- function (physeq, x = "Sample", y = "Abundance", fill = NULL, 
                       title = NULL, facet_grid = NULL,factor_levels = NULL) 
{
  mdf = psmelt(physeq)
  mdf$Sample <- factor(mdf$Sample,levels = factor_levels)
  p = ggplot(mdf, aes_string(x = x, y = y, fill = fill))
  p = p + geom_bar(stat = "identity", position = "stack", 
                   color = "black")
  p = p + theme(axis.text.x = element_text(angle = -90, hjust = 0))
  if (!is.null(facet_grid)) {
    p <- p + facet_grid(facet_grid)
  }
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  return(p)
}

tax_stack_clust_2 <- function (otu = NULL, map = NULL, tax = NULL, dist = "bray", 
                               Group = "Group", j = "Phylum", rep = 6, Top = 10, tran = TRUE, 
                               hcluter_method = "complete", cuttree = 3,Genus_levels = NULL) 
{
  vegan_otu = function(physeq) {
    OTU = otu_table(physeq)
    if (taxa_are_rows(OTU)) {
      OTU = t(OTU)
    }
    return(as(OTU, "matrix"))
  }
  vegan_tax <- function(physeq) {
    tax <- tax_table(physeq)
    return(as(tax, "matrix"))
  }
  idx = rownames(otu) %in% rownames(tax)
  otu = otu[idx, ]
  tax = tax[rownames(otu), ]
  map <- map[Group]
  colnames(map) = "Group"
  map$ID = row.names(map)
  ps = phyloseq(sample_data(map), otu_table(as.matrix(otu), 
                                            taxa_are_rows = TRUE), tax_table(as.matrix(tax)))
  random_tree_raw = rtree(ntaxa(ps), rooted=TRUE, tip.label=taxa_names(ps))
  ps <- merge_phyloseq(ps,random_tree_raw)
  
  ps1_rela = phyloseq::transform_sample_counts(ps, function(x) x/sum(x))
  otu = as.data.frame(t(vegan_otu(ps1_rela)))
  unif = phyloseq::distance(ps1_rela, method = dist)
  hc <- stats::hclust(unif, method = hcluter_method)
  clus <- cutree(hc, cuttree)
  d = data.frame(label = names(clus), member = factor(clus))
  map = as.data.frame(sample_data(ps))
  dd = merge(d, map, by = "row.names", all = F)
  row.names(dd) = dd$Row.names
  dd$Row.names = NULL
  p = ggtree(hc) %<+% dd + geom_tippoint(size = 5, shape = 21, 
                                         aes(fill = Group, x = x)) + geom_tiplab(aes(color = Group, 
                                                                                     x = x * 1.2), hjust = 1)
  psdata = ps1_rela %>% tax_glom(taxrank = j)
  if (tran == TRUE) {
    psdata = psdata %>% transform_sample_counts(function(x) {
      x/sum(x)
    })
  }
  otu = otu_table(psdata)
  tax = tax_table(psdata)
  for (i in 1:dim(tax)[1]) {
    if (row.names(tax)[i] %in% names(sort(rowSums(otu), 
                                          decreasing = TRUE)[1:Top])) {
      tax[i, j] = tax[i, j]
    }
    else {
      tax[i, j] = "Other"
    }
  }
  tax_table(psdata) = tax
  Taxonomies <- psdata %>% psmelt()
  Taxonomies$Abundance = Taxonomies$Abundance * 100
  Taxonomies$OTU = NULL
  colnames(Taxonomies)[1] = "id"
  p <- p + ggnewscale::new_scale_fill()
  p1 <- facet_plot(p, panel = "Stacked Barplot", data = Taxonomies, 
                   geom = geom_barh, mapping = aes(x = Abundance, fill = Species), 
                   color = "black", stat = "identity")
  grotax <- Taxonomies %>% group_by(Group, Species) %>% summarise(Abundance = mean(Abundance))
  
  combined_data <- data.frame()
  for (i in unique(grotax$Group)) {
    tmp <- grotax %>% 
      filter(Group == i) %>%
      mutate(Abundance = if_else(Species == "Other", 
                                 100 - sum(Abundance[Species != "Other"]), 
                                 Abundance)
      )
    combined_data <- bind_rows(combined_data, tmp)
  }
  
  grotax <- combined_data
  ps1_rela = phyloseq::transform_sample_counts(ps, function(x) x/sum(x))
  hc = phyloseq::merge_samples(ps1_rela, "Group", fun = mean) %>% 
    phyloseq::distance(method = dist) %>% stats::hclust(method = hcluter_method)
  clus <- cutree(hc, cuttree)
  d = data.frame(label = names(clus), member = factor(clus))
  map = as.data.frame(sample_data(phyloseq::merge_samples(ps1_rela, 
                                                          "Group", fun = mean)))
  dd = merge(d, map, by = "row.names", all = F)
  row.names(dd) = dd$Row.names
  dd$Row.names = NULL
  p3 = ggtree(hc) %<+% dd + geom_tippoint(size = 5, shape = 21, show.legend = FALSE, 
                                          aes(fill = member, x = x)) + geom_tiplab(aes(color = member, 
                                                                                       x = x * 1.2), hjust = 1)
  p3 <- p3 + ggnewscale::new_scale_fill()
  grotax$Species <- factor(grotax$Species,levels = Genus_levels)
  p4 <- facet_plot(p3, panel = "Stacked Barplot", data = grotax, 
                   geom = geom_barh, mapping = aes(x = Abundance, fill = Species), 
                   color = "black", stat = "identity")+
  scale_fill_manual(values = c("#0070C0","#FCB777","#FFFFDF","#AFDEAC","#FFD966","grey"))
  return(list(p, p1, p3, p4, Taxonomies))
}


#### analysis 

bac_abd <- read.xlsx("Fig.5A、B、C、D_S8A.xlsx",sheetIndex = 1)
colnames(bac_abd) <- colnames(bac_abd) %>% gsub("[.]","-",.)

merged_profiles <- filter_taxa_lvl(df = bac_abd, taxa_lvl = "s")
physeq_merged <- metaphlan_to_phyloseq(
  merged_profiles,
  taxa_lvl = "Species"
)
sample_data(physeq_merged) <- data.frame(row.names =  colnames(bac_abd)[-1],
                                         treatment = c(rep("Con",5),rep("LT",5),rep("DMSG",11),rep("Seawater",4),rep("Sediment",3)),
                                         raw_group = c(rep("Con",5),rep("LP",5),rep("DMSG",11),rep("Seawater",4),rep("Sediment",3)))

random_tree_raw = rtree(ntaxa(physeq_merged), rooted=TRUE, tip.label=taxa_names(physeq_merged))
phy_tree(physeq_merged) <- random_tree_raw

### analysis

#### abd starck plot
class_profiles <- filter_taxa_lvl(df = bac_abd, taxa_lvl = "c")
class_profiles2 <- data.frame(row.names = unlist(lapply(strsplit(class_profiles$clade_name,"c__"),function(x)x[2])),
                              class_profiles[,-1])
colnames(class_profiles2) <- colnames(class_profiles2) %>% gsub("[.]","-",.)


class.ave <- apply(class_profiles2, 1, FUN=mean)
class.2 <- cbind(class_profiles2, class.ave)[order(-class.ave),] #排个序
class.2 <- subset(class.2[1:10,], select=-class.ave)
class.2 <- rbind(class.2, Others=apply(class.2, 2, function(x){100-sum(x)}))
class.2 <- cbind(classID=row.names(class.2), class.2)

class.2$classID <- factor(class.2$classID, levels = rev(class.2$classID))
class.gg <- melt(class.2, id.vars="classID", variable.name="SampleID", value.name="Abundance")
class.gg$group <- unlist(lapply(str_split(class.gg$SampleID,"[-]"),function(x)x[1]))
class.gg <- class.gg %>% 
  mutate(Abundance = ifelse(SampleID %in% c("Con_1", "Con_2") & classID == "Others", 0, Abundance))

class.gg$group <- factor(class.gg$group,levels = c("Con","LT","DMSG","Seawater","Sediment"))

color_set <- rev(wheel("skyblue3")[1:10])
color_set[9:10] <- c("#FFD966","#0070C0")

######################## Fig.5A ########################
all_abd_stark <- ggbarplot(class.gg, x = "SampleID", y="Abundance", color="black", fill="classID",
                           legend="right",legend.title="Class")+
  theme_bw() +
  rotate_x_text() + 
  scale_fill_manual(values=c("gray",color_set)) + # 颜色设置
  facet_grid(~ group, scales = "free_x", space='free') + 
  labs(x = NULL, y = "Relative abundance (100%)") + 
  theme(axis.ticks.x=element_blank()) 


### Manually erase the barplot for Con-1 and Con-2 in the output


######################## Fig.5B ########################

physeq_merged <- subset_samples(physeq_merged,!sample_names(physeq_merged) %in% c("Con-1","Con-2"))

p_r <- plot_richness(physeq_merged, "treatment", measures="Shannon")
p_r$data$treatment <- factor(p_r$data$treatment,levels = c("Con","LT","DMSG","Seawater","Sediment"))
r_plot <- ggplot(data = p_r$data,aes(x = treatment, y = value, fill = treatment)) +
  geom_jitter(aes(color = treatment))+
  geom_boxplot(outlier.shape = NA)+
  theme_bw()+
  labs(x = NULL, y = "Shannon diversity index") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


######################## Fig.5C ########################
pcoa_result <- ordinate(physeq_merged, method = "PCoA",distance = "unifrac")
pcoa.plot <- plot_ordination(physeq_merged, pcoa_result, color="treatment")+
  theme_bw()

pcoa.plot$data$treatment <- factor(pcoa.plot$data$treatment ,levels = c("Con","LT","DMSG","Seawater","Sediment"))
pcoa.plot <- pcoa.plot + geom_point(size=3,aes(color = treatment)) +
  labs(title="Unifrac Distance")+
  theme_bw()+
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("#F8766D","#00BF7D","#A3A500","#E76BF3","#00B0F6"))+
  ggforce::geom_mark_ellipse(aes(fill = treatment), fill = NA, linetype = "dashed", alpha = 0.1, expand = 0.01)+
  ggtitle("PCoA based on Unifrac dissimilarity")


######################## Fig.5D ########################

p_habitat_details = tax_stack_clust_2(otu = as.data.frame(physeq_merged@otu_table), 
                                      map = data.frame(physeq_merged@sam_data), 
                                      tax = data.frame(physeq_merged@tax_table),
                                      Group ="treatment",tran=TRUE,Top = 5,j = "Species",
                                      Genus_levels = c("Bathymodiolus_platifrons_methanotrophic_gill_symbiont",
                                                       "Bradyrhizobium_diazoefficiens","GGB60957_SGB82997",
                                                       "Marine_Group_I_thaumarchaeote","Sphingomonas_sp_3F27F9",
                                                       "Other"),
                                      dist="unifrac", hcluter_method="average",cuttree = 2)

print(p_habitat_details[[4]])


######################## Fig.5F ########################
# use own microbiome data
self_data <- merged_profiles[,c(1:11)]
self_data <- tibble::column_to_rownames(self_data,"clade_name")
self_data <- self_data[rowSums(self_data)!=0,]
self_data[,1:10] <- normalize.quantiles(as.matrix(self_data))

#for model try
#write.csv(self_data_add,"D:/works/io-rna/hanjia_analysis/metagenome/data/quantiles.csv")


core_changed_bac <- self_data %>% filter(grepl("s__Bathymodiolus_platifrons_methanotrophic_gill_symbiont", rownames(self_data))) %>% as.matrix()
rownames(core_changed_bac) <- "Methanotrophic_gill_symbiont"

gene_tmm <- fread("salmon_abundance_estimates_to_matrix.isoform.TMM.EXPR.matrix")
gene_tmm <- data.frame(row.names = gene_tmm$V1,gene_tmm[,c(-1,-7:-11)]) %>% as.matrix()
colnames(gene_tmm) <- colnames(core_changed_bac)

WGCNA::cor(t(gene_tmm),t(core_changed_bac),use = "everything") -> core_cor 
WGCNA::corPvalueStudent(core_cor, nrow(t(gene_tmm))) %>% as.data.frame() -> cor.p
methanotrophic_cor <- cbind(core_cor,cor.p)
colnames(methanotrophic_cor) <- c("Correlation","Pvalue")  ## all correlation results, as same as Supplementary Data 2 


cor_data <- fread("Fig.5F.txt")

### gene expression heatmap
rna_seq_exp_phm <- pheatmap(data.frame(row.names = cor_data$Class[-25],cor_data[-25,6:15]),
                            cluster_rows = F,cluster_cols = F,scale = "row",color = cols4all::c4a("blue_red3",255),
                            gaps_col = 5,border_color = "black",cellheight = 10,cellwidth = 10)


### methanotrphic data for single row heatmap
single_methanotrophic <- cor_data[25,6:15]
single_methanotrophic_phm <- pheatmap(single_methanotrophic,cluster_rows = F,cluster_cols = F,gaps_col = 5,border_color = "black",show_rownames = F,
                                      cellwidth = 10,cellheight = 10,color = colorRampPalette(c("white", "#0070C0"))(100))


### correlation barplot
cor_bar_plot <- cor_data[-25,] %>% 
  ggplot(aes(x = Cor, y = reorder(Class,Rank,decreasing = T),fill = -log2(P))) +
  geom_bar(stat = "identity")+
  xlim(c(0,1)) +
  scale_fill_continuous_c4a_seq("or_rd")+
  theme_classic()




