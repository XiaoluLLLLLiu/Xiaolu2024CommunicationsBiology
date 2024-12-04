rm(list = ls())
library(WGCNA)
library(tibble)
library(xlsx)
library(clusterProfiler)
library(pheatmap)

set.seed(406206)
setwd("submit_source_data/source data")


#data preprocess
rlog.m <- assay(rld)  ### rlg data from Fig.3
rlog.m <- rlog.m[,c(6:18,1:5,19:27)]
mussel.ogs.df <- rownames_to_column(data.frame(rlog.m), "OG_names")
data.frame(OG_names = mussel.orthogroup.rc$Orthogroup,Mussel_id = mussel.orthogroup.rc$Mussel) -> og2mussel

left_join(mussel.ogs.df,og2mussel,by = "OG_names") -> tmp1
data.frame(row.names = tmp1$Mussel_id,tmp1[2:28]) -> rlog.df
colnames(rlog.df) <- colnames(rlog.df) %>% gsub("[.]","-",.)


######################## WGCNA ########################
enableWGCNAThreads(nThreads = 5)

m.mad <- apply(rlog.df, 1, mad)
dataExprVar <- rlog.df[which(m.mad > max(quantile(m.mad, probs=seq(0, 1, 0.25))[1],0)),];dim(dataExprVar)

dataExpr <- as.data.frame(t(dataExprVar))
gsg = goodSamplesGenes(dataExpr, verbose = 3)

if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:",
                     paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:",
                     paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")));
  # Remove the offending genes and samples from the data:
  dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
}

sampleTree = hclust(dist(dataExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="");abline(h = 40000, col = "red") 

clust = cutreeStatic(sampleTree, cutHeight = 40000, minSize = 1)
table(clust)

keepSamples = (clust==1)
dataExpr = dataExpr[keepSamples, ]


#phenotypic data
phenoData <- data.frame(Sample = colnames(rlog.m),
                        Con = c(rep(1,5),rep(0,22)),
                        Shallow = c(rep(0,5),rep(1,8),rep(0,14)),
                        LT = c(rep(0,13),rep(1,5),rep(0,9)),
                        Deep = c(rep(0,18),rep(1,9)))

Samples = rownames(dataExpr) 
traitRows = match(Samples, phenoData$Sample)
datTraits = as.data.frame(phenoData[traitRows, -1])
rownames(datTraits) = phenoData[traitRows, 1]

row.names(datTraits) <- row.names(datTraits)


sampleTree2 = hclust(dist(dataExpr,method = "manhattan"), method = "average")

traitColors = numbers2colors(datTraits, signed = FALSE); 
plotDendroAndColors(sampleTree2, 
                    traitColors,
                    groupLabels = names(phenoData[,-1]), 
                    main = "Sample dendrogram and trait heatmap")


## soft threshold select
powers = c(c(1:15))
sft = pickSoftThreshold(dataExpr, powerVector=powers,
                        networkType="unsigned", verbose=5)

par(mfrow = c(1,2)) 
cex1 = 0.6 

##### Fig.S7A
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2", type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers, cex=cex1, col="red")
abline(h=0.85, col="black", lty=4, lwd=1) 


plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")

plot1 <- recordPlot()


power = sft$powerEstimate
adjacency = adjacency(dataExpr, power = power)


TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM


geneTree = hclust(as.dist(dissTOM), method = "average")
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)

minModuleSize = 30
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)


dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")


MEList = moduleEigengenes(dataExpr, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average")
plot(METree, main = "Clustering of module eigengenes",xlab = "", sub = "") 


MEDissThres = 0.3 
abline(h=MEDissThres, col = "red")


merge = mergeCloseModules(dataExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs


plotDendroAndColors(geneTree, mergedColors,"Module Colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Module Colors")


moduleColors = mergedColors
table(moduleColors)


colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs

## calculate the correlation within the modules

nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(300),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

### calculate similarity between modules and genes

modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(dataExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")


### calculate similarity between datTraits and genes
traitNames=colnames(datTraits)
geneTraitSignificance = as.data.frame(cor(dataExpr, datTraits, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", traitNames, sep="")
names(GSPvalue) = paste("p.GS.", traitNames, sep="")



### output the members in each module
InterestedModule = datTraits$LT
InterestedModule = as.data.frame(InterestedModule);
names(InterestedModule) = "InterestedModule"

modNames = substring(names(MEs), 3) 

geneModuleMembership = as.data.frame(cor(dataExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(dataExpr, InterestedModule, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(InterestedModule), sep="");
names(GSPvalue) = paste("p.GS.", names(InterestedModule), sep="");


anno <- fread("nr.blast-fmt6.tsv")
accession_ids <- gsub("^(\\S+).*", "\\1",  anno$V7)
gene_names <- gsub(".*? (\\S[^\\[]+) \\[.*\\]$", "\\1", anno$V7)
species_names <- gsub(".*\\[([^\\]]+)\\]$", "\\1", anno$V7)

anno_full <- data.frame(anno,accession = accession_ids,gene = gene_names,sp_name = species_names)


probes = names(dataExpr)
probes2anno = match(probes,anno_full$V1)


geneInfo0 = data.frame(gene_id = colnames(dataExpr),
                       geneanno = anno_full$V7[probes2anno],
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)

######################## Fig.4B ########################
all_genes_tmm <- read.table("salmon_abundance_estimates_to_matrix.isoform.TMM.EXPR.matrix",header = T,row.names = 1,sep = "\t")
colnames(all_genes_tmm) <- gsub("[.]","-",colnames(all_genes_tmm))


### perform the GSEA

#go annotations file
go.annotations.files <- fread('de_lck_go_annotation.txt.Parsed.Go2Gene.xls', header=T, stringsAsFactors=F)

goterm2gene <- go.annotations.files[,c(1,4)]
goterm2name <- go.annotations.files[,c(1,3)]

blue_list <- geneInfo0 %>% dplyr::filter(moduleColor == "blue") %>% .$gene_id
blue_expr <- all_genes_tmm[row.names(all_genes_tmm) %in% blue_list,]
blue_expr$log2FoldChange <- apply(blue_expr,1,function(x) log2(sum(x[11:15])/sum(x[1:5])))
blue_expr[is.infinite(blue_expr$log2FoldChange), "log2FoldChange"] <- 10

blue.gsea.list <- blue_expr %>% dplyr::select(.,"log2FoldChange") 
blue.gsea.list$geneid <- row.names(blue.gsea.list)
blue.gsea.list <- blue.gsea.list[order(blue.gsea.list$log2FoldChange,decreasing = T),]
blue.gene_list <- blue.gsea.list$log2FoldChange
names(blue.gene_list) <- blue.gsea.list$geneid

blue.gesa.go.res <- GSEA(gene = blue.gene_list,
                         TERM2GENE = goterm2gene,
                         TERM2NAME = goterm2name,
                         pvalueCutoff = 1,
                         pAdjustMethod = 'BH',
                         nPermSimple = 10000)
green_gsea_res_filter <- blue.gesa.go.res@result %>% dplyr::filter(pvalue < 0.05)


######################## Fig.4C ########################


fig4c_raw <- read.xlsx("Fig.4B、C.xlsx",sheetIndex = 1)
fig4c_hub_info <- fig4c_raw[2:12,11:24]

single.df # Used in Fig.3


single.df[single.df$Mussel %in% fig4c_hub_info$Hub.genes.identified.by.cytoHubba,]$Orthogroup  -> hub_ogs
hub_phm_df <- rlog.m[rownames(rlog.m) %in% hub_ogs,]
rownames(hub_phm_df) <- c("SOX1S","COX11","TPMT","KLF2","NDUFB11","RAD51L3","NDUFAF1","NDUFB5","EEF1B","RP-S11e","RP-L27Ae")

pheatmap(t(hub_phm_df), scale = "column",treeheight_col = 30,
         border_color = "Black",cluster_row = F, gaps_row = 13,cutree_rows = 2)  -> hub.phm


