#!/usr/bin/env Rscript
#   Script to do the WGCNA analysis
# Written by Li Lei, 01-27-2020, Berkley
#getwd();
# If necessary, change the path below to the directory where the data files are stored.
# "." means current directory. On Windows use a forward slash / instead of the usual \.
#workingDir = ".";
#setwd(workingDir);
# Load the WGCNA package
#install.packages("impute")

#install.packages("BiocInstaller")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("impute")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("preprocessCore")

#BiocManager::install("org.Mm.eg.db")
#BiocManager::install("data.table")

#install.packages("WGCNA")
source("https://bioconductor.org/biocLite.R")
biocLite(c("AnnotationDbi", "impute","GO.db", "preprocessCore"))
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
install.packages(c("WGCNA", "stringr", "reshape2"), repos=site)
#dependencies 'robustbase', 'mvtnorm', 'pcaPP' are not available for package 'rrcov'
biocLite("textshape")
biocLite("ggfortify")

install.packages("dynamicTreeCut", "cluster", "flashClust", "Hmisc", "reshape", "foreach", "doParallel" ) 
install.packages("igraph","ggfortify")
install.packages("ggfortify")
library(igraph)
library(textshape)
library(WGCNA)
library(reshape2)
library(stringr)
library(tidyr)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

type = "unsigned"
#type = "signed"
corType = "pearson"
#corFnc = ifelse(corType=="pearson", cor, bicor)
maxPOutliers = ifelse(corType=="pearson",1,0.05)
robustY = ifelse(corType=="pearson",T,F)

tpm <- read.delim("/global/projectb/scratch/llei2019/CSP_Kranthi/sus_beni/normalized_VST_final_sus_beni_count_formal.txt",sep="\t",header=T)
#metaData <- read.delim("/global/projectb/scratch/llei2019/CSP_Kranthi/sus_beni/Meta_data_sus_beni.txt",header = T)
#Since I tested different types of the data, and found that only focusing on the leaves and the pathogen treatments.
metaData <- read.delim("/global/projectb/scratch/llei2019/CSP_Kranthi/sus_beni/pathogen_meta.txt",header = T)
#metaData <- read.delim("/global/projectb/scratch/llei2019/CSP_Kranthi/sus_beni/meta_pathogen_plus_root.txt",header = T)

#extract the column I needed!!!
col.num <- which( colnames(tpm) %in% metaData$code)
row.num <- which( metaData$code %in% colnames(tpm) )
tpm <- tpm[, col.num]
metaData <- metaData[row.num,]
ncol(tpm)
head(tpm)
head(metaData)
#code are th333same with the code in the data
#header <- paste(metaData$accessions,metaData$experiment,metaData$tissue,metaData$treatment, metaData$timepoint, metaData$replicate,sep='_')
#Since only focus on the leaves, we do not need to add tissues
header <- paste(metaData$accessions,metaData$experiment,metaData$treatment, metaData$timepoint, metaData$replicate,sep='_')

metaData$group <- paste(metaData$accessions,metaData$experiment,metaData$treatment, sep='_')
#metaData$group <- paste(metaData$accessions,metaData$experiment,metaData$tissue,metaData$treatment,metaData$timepoint, sep='_')
#metaData$group <- paste(metaData$accessions,metaData$experiment,metaData$tissue,metaData$treatment, sep='_')

datTraits <- data.frame(header,metaData$group)
str(datTraits)
colnames(datTraits) <- c("ID","group")
head(datTraits)
#tpm <- read.delim("/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/VST_normalized_count_abiotic.txt",sep="\t",header=T)
#"/global/cscratch1/sd/llei2019/SeanGordon_Distachyon_gene_counts/counts.txt"
#use log(tpm+1) to do a little bit conversion
head(tpm)
str(tpm)
#I compare the order of the code in the tmp file and the metadata file and found that they are teh same
colnames(tpm) <- header
#tpm[1,]
#hist(tpm$Shoot_CK_1h.s1)
#femData <- sapply(tpm,function(x) log(x+1))
femData <- tpm
#column_to_rownames(tpm,"X")
#row.names(femData) <- tpm$GID
femData <-data.frame(femData)

### select genes with MAD >85%
dim(femData) #[1] 19982    78 #pathogen[1] 19982    42 this got picked#[1] 19982    54
#25189/3 = 8396.333
head(femData)
m.mad <- apply(femData,1,mad)#median absolute deviation (MAD)
head(m.mad)
#femData <- femData[which(m.mad > 0),]
#This is to take the top75% genes
#femData <- femData[which(m.mad > 
#colors()[1:17]   
quantile(m.mad, probs=seq(0, 1, 0.05))
#0%        5%       10%       15%       20%       25%       30%       35%       40% 
#0.0000000 0.1947251 0.2402565 0.2774659 0.3134411 0.3482913 0.3846300 0.4248659 0.4686480 
#45%       50%       55%       60%       65%       70%       75%       80%       85% 
#0.5150579 0.5667942 0.6299805 0.6932431 0.7738405 0.8664585 0.9766338 1.1200309 1.2992371 
#90%       95%      100% 
#1.5503801 1.9688731 5.6023635 
#42 samples:
#0%        5%       10%       15%       20%       25%       30%       35%       40% 
#0.0000000 0.1315377 0.1791380 0.2097212 0.2365621 0.2631993 0.2886189 0.3143882 0.3447190 
#45%       50%       55%       60%       65%       70%       75%       80%       85% 
#0.3744002 0.4066705 0.4430818 0.4833315 0.5290972 0.5870723 0.6513397 0.7363367 0.8407120 
#90%       95%      100% 
#0.9902778 1.2206653 5.5611807 
#above 75%
femData <- femData[which(m.mad > 
                           max(quantile(m.mad, probs=seq(0, 1, 0.05))[6],0.01)),]
#femData <- femData[which(m.mad > 
#                           max(quantile(m.mad, probs=seq(0, 1, 0.05))[16],0.01)),]

###
datExpr0 <- as.data.frame(t(femData))
head(datExpr0[,1:8])
dim(datExpr0) #[1]    78 14986 #[1]    42 14986 #[1]    54 14986
rownames(datExpr0)
colnames(datExpr0)

#This function iteratively identifies samples and genes with too many missing entries and genes with zero variance. 
#Iterations may be required since excluding samples effectively changes criteria on genes and vice versa. 
#The process is repeated until the lists of good samples and genes are stable. 
gsg = goodSamplesGenes(datExpr0, verbose = 3)

gsg$allOK
### remove the offending genes and samples from the data
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
nGenes = ncol(datExpr0)#75%[1] 14986 #[1] 14986
nSamples = nrow(datExpr0)#78 #42#54#pick 42 set 
head(datExpr0[,1:8])
sampleTree = hclust(dist(datExpr0), method = "average")
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small. sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
#par(cex = 0.6);
par(mar = c(0,5,2,0))
pdf(file = "/global/projectb/scratch/llei2019/CSP_Kranthi/sus_beni/VST_sus_42_patho_tree_mad75_trait.pdf",width = 25,height = 20)
plot(sampleTree,main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 2,cex.axis = 1.5, cex.main = 2)
abline(h = 150, col = "red");
myCol = c("pink1", "violet", "mediumpurple1", "slateblue1", "purple", "purple3",
          "turquoise2", "skyblue", "steelblue", "blue2", "navyblue",
          "orange", "tomato", "coral2", "palevioletred", "violetred", "red2")

sample_colors <- numbers2colors(as.numeric(factor(datTraits$group)), 
                                colors = myCol, signed = FALSE)
par(mar = c(1,4,3,1),cex=0.8)
plotDendroAndColors(sampleTree, sample_colors,
                    groupLabels = colnames(sample),
                    cex.dendroLabels = 0.8,
                    marAll = c(1, 4, 3, 1),
                    cex.rowText = 0.01,
                    main = "Sample dendrogram and trait heatmap")

dev.off()
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 150, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust!=0)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)#14986
nSamples = nrow(datExpr)#78

#all of the samples has passed!
#do a PCA analysis for the samples to see what will happen:
ir.pca <- prcomp(datExpr,
                 center = TRUE,
                 scale. = TRUE) 
print(ir.pca)
#meta_data <- read.delim("/global/u2/l/llei2019/Brachypodium/Sylvaticum/scripts/Meta_data_syl.txt",header = T)
#meta_data <- read.delim("/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/distachyon_metadata_minor1.txt",header = T)
#head(meta_data)
#nrow(meta_data)
#84
#abiotic_shoot <- meta_data[(meta_data$tissue == "shoot" & (meta_data$experiment == "timecourse" | meta_data$experiment == "abiotic")),]
#abiotic_shoot$group <- paste(abiotic_shoot$treatment, abiotic_shoot$timepoint, sep='_')

library(ggfortify)
#pdf(file = "/global/projectb/scratch/llei2019/CSP_Kranthi/VST_sus_all_sample_pca_mad75.pdf",width = 10,height = 10)
autoplot(ir.pca,data = datTraits, colour = 'group',shape = FALSE, label.size = 3)
ggsave("/global/projectb/scratch/llei2019/CSP_Kranthi/sus_beni/VST_sus_42_patho_mad75.pdf", width = 30, height = 25)
dev.off()
#####
#sampleTree = hclust(dist(datExpr), method = "average")
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small. sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
#par(cex = 0.6);
#par(mar = c(0,4,2,0))
#pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/VST_combined_abiotic_distachyon_minor1.pdf",width = 25,height = 20)
#plot(sampleTree,main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)
#abline(h = 160, col = "red");
#dev.off()


###Automatic block-wise network
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
head(sft)
power = sft$powerEstimate
power
#6
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9
pdf(file = "/global/projectb/scratch/llei2019/CSP_Kranthi/sus_beni/power_soft_thresh_VST_sus_42_patho_mad75.pdf",width = 10,height = 10)
# Scale-free topology fit index as a function of the soft-thresholding power plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
plot(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n", main = "Scale independence")

text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red")
abline(h=0.8,col="red")
#abline(h=0.9,col="red")

# this line corresponds to using an R^2 cut-off of h abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power 
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
#abline(h=0.85,col="red")
dev.off()

#power = 16 is the best according to the plot, I set the mergeCutHeight = 0.3
#try 6 to see how it look like
#test to see if signed is better since this suggested https://peterlangfelder.com/2018/11/25/signed-or-unsigned-which-network-type-is-preferable/
#power = 12 is the best for pathogen (only leaves) datasets
#power = 8 is the best for pathogen datasets

bwnet = blockwiseModules(datExpr,maxBlockSize = nGenes,
                         power = 12, TOMType = type, minModuleSize = 30,
                         reassignThreshold = 0, mergeCutHeight = 0.25, numericLabels = TRUE,pamRespectsDendro = FALSE,
                         saveTOMs = TRUE,corType = corType,
                         maxPOutliers=maxPOutliers, loadTOMs=TRUE,
                         saveTOMFileBase = "BrachySus_TOM-blockwise", verbose = 3)
head(bwnet)
bwLabels = bwnet$colors
table(bwnet$colors)
#0    1    2    3    4    5    6    7    8    9   10   11   12   13 
#405 4941 4756 1716  954  650  507  315  221  173  140  124   48   36 
#power=6
# 0    1    2    3    4    5    6    7    8    9   10 
#53 6453 5552  962  636  577  385  134   88   88   58 
#pathogen power=12 only for leaves alll#####
#0    1    2    3    4    5    6    7    8    9   10   11 
#62 5321 4380 2352 1364  511  487  168  120  105   83   33 
#########################################
#pathogen+root power=8
#0    1    2    3    4    5    6    7    8    9   10   11 
#36 6723 4642 1427  747  646  244  182   91   85   82   81
bwModuleColors = labels2colors(bwLabels)
table(bwModuleColors)
#black        blue       brown       green greenyellow        grey     magenta 
#315        4756        1716         650         124         405         173 
#pink      purple         red      salmon         tan   turquoise      yellow 
#221         140         507          36          48        4941         954 

#Power=6
#black      blue     brown     green      grey   magenta      pink    purple       red 
#134      5552       962       577        53        88        88        58       385 
#turquoise    yellow 
#6453       636 
#pathogen power=12:#######
#bwModuleColors
#black        blue       brown       green greenyellow        grey     magenta 
#168        4380        2352         511          33          62         105 
#pink      purple         red   turquoise      yellow 
#120          83         487        5321        1364 
#############################
#plut root +pathogen power=8
#black        blue       brown       green greenyellow        grey     magenta 
#182        4642        1427         646          81          36          85 
#pink      purple         red   turquoise      yellow 
#91          82         244        6723         747 

#try to merge:
#merge = mergeCloseModules(datExpr, bwModuleColors[bwnet$blockGenes[[1]]], cutHeight = 0.3, verbose = 3)
#mergedColors = merge$colors
#moduleLabels = merge$colors
#cbind(moduleLabels, bwnet$colors)
#table(merge$colors)
#black        blue       green greenyellow        grey     magenta        pink 
#315        6472         650         124         405         173         221 
#purple         red      salmon      yellow 
#188        5448          36         954

# open a graphics window
sizeGrWindow(6,6)
pdf(file = "/global/projectb/scratch/llei2019/CSP_Kranthi/sus_beni/Dendro_VST_sus_patho_mad75_cut0.25_p12.pdf",width = 25,height = 25)
# Plot the dendrogram and the module colors underneath for block 1
plotDendroAndColors(bwnet$dendrograms[[1]], bwModuleColors[bwnet$blockGenes[[1]]],
                    "Module colors", main = "Gene dendrogram and module colors in block 1",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwModuleColors[bwnet$blockGenes[[1]]], mergedColors),
#                    c("Dynamic Tree Cut", "Merged dynamic"), main = "Gene dendrogram and module colors in block 1",
#                    dendroLabels = FALSE, hang = 0.03,
#                    addGuide = TRUE, guideHang = 0.05)
dev.off()
#Power=6 seems to bebetter
#head(bwnet) Since there are not huge different between the merged and unmerged,so still keep the unmerged
#mergedMEs = merge$newMEs
MEs = bwnet$MEs
MEs_col = MEs
#MEs_col = mergedMEs
colnames(MEs_col) = paste0("ME", labels2colors(as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)
pdf(file = "/global/projectb/scratch/llei2019/CSP_Kranthi/sus_beni/eigennetwork_VST_sus_patho_mads75_0.25_power12.pdf",width = 15,height = 13)
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                      xLabelsAngle = 90)
dev.off()

head(MEs_col)
graph<-wgcna2igraph(net = bwnet, datExpr = datExpr,
                    modules2plot = c("blue","green","turquoise","brown"),
                    colors2plot = c("orange","darkred","cyan","cornflowerblue"),
                    kME.threshold = 0.5, adjacency.threshold = 0.1,
                    adj.power = 12, verbose = T,
                    node.size = 0, frame.color = NA, node.color = NA,
                    edge.alpha = .5, edge.width =1)
plot(graph)

#hubs    = chooseTopHubInEachModule(datExpr, 
#                                   bwModuleColors,omitColors = "grey", 
#                                   power = 12, 
#                                 type = "unsigned" )
#This function can only find the single hub genes
#hubs
#trait module association
#In R:
#Plot the relationships between modules and traits
str(datTraits)
design = model.matrix(~0+ as.factor(datTraits$group))
design <- design[,c(5,6,1,2,4,3)]
#colnames(design) = c("Bd21_fungi_leaves_Mock","Bd21_fungi_leaves_PCA",
#                     "Bd21-3_bacteria_leaves_Mock","Bd21-3_bacteria_leaves_XT",
#                     "Bd21-3_virus_leaves_Mock","Bd21-3_virus_leaves_BSMV",
#                     "Bd21-3_virus_root_Mock", "Bd21-3_virus_root_BSMV",
#                     "Bd21-3_Endophytes_shoot_Mock","Bd21-3_Endophytes_shoot_SV",
#                     "Bd21-3_Endophytes_root_Mock", "Bd21-3_Endophytes_root_SV")
#colnames(design) = levels(as.factor(datTraits$group))
colnames(design) = c("Bd21_fungi_Mock","Bd21_fungi_PCA",
                     "Bd21-3_bacteria_Mock","Bd21-3_bacteria_XT",
                     "Bd21-3_virus_Mock", "Bd21-3_virus_BSMV")
#View(design)
#export the design matrix and see what they look like
#write.table(x = design,
#            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/sus_beni/design.txt",
#            quote = FALSE,
#            sep = "\t",
#            eol = "\n",
#            col.names = TRUE,
#            row.names = TRUE)

moduleColors = labels2colors(bwnet$colors)
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, design, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

#Display the correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
pdf("/global/projectb/scratch/llei2019/CSP_Kranthi/sus_beni/pheno-module_time_combined_VST_suspatho_mads75_encoding1_0.25power12.pdf",width = 15,height = 15)
#Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(design),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.6,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

#trait <- "/global/projectb/scratch/llei2019/CSP_Kranthi/pheno_sus_encoding1.txt"

#if(trait != "") {
# traitData <- read.table(file=trait, sep='\t', header=T, row.names=1,
#                         check.names=FALSE, comment='',quote="")
# head(traitData)
# sampleName = rownames(datExpr)
# traitData = traitData[match(sampleName, rownames(traitData)), ]
#}
#nrow(MEs_col)
#nrow(traitData)
#robustY = ifelse(corType=="pearson",T,F)
#if (corType=="pearson") {
#  modTraitCor = cor(MEs_col, traitData, use = "p")
#  modTraitP = corPvalueStudent(modTraitCor, nSamples)
#} else {
#  modTraitCorP = bicorAndPvalue(MEs_col, traitData, robustY=robustY)
#  modTraitCor = modTraitCorP$bicor
#  modTraitP   = modTraitCorP$p
#}

#textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
#dim(textMatrix) = dim(modTraitCor)

#sizeGrWindow(13, 9)
#par(mar = c(5, 30, 5, 5),mfrow = c(1,1))
#cex1 = 1
#pdf("/global/projectb/scratch/llei2019/CSP_Kranthi/pheno-module_combined_abiotic_VST_sus_all_sample_mads75_encoding1_0.25p6.pdf",width = 15,height = 15)
#labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(traitData), 
#               yLabels = colnames(MEs_col), 
#               cex.lab = 1, 
#               ySymbols = colnames(MEs_col), colorLabels = FALSE, 
#               colors = blueWhiteRed(50), 
#               textMatrix = textMatrix, setStdMargins = FALSE, 
#               cex.text = 0.8, zlim = c(-1,1),
#               main = paste("Module-trait relationships"))
#dev.off()

if (corType=="pearson") {
  geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
  MMPvalue = as.data.frame(corPvalueStudent(
    as.matrix(geneModuleMembership), nSamples))
} else {
  geneModuleMembershipA = bicorAndPvalue(datExpr, MEs, robustY=robustY)
  geneModuleMembership = geneModuleMembershipA$bicor
  MMPvalue   = geneModuleMembershipA$p
}
head(geneModuleMembership)
head(datExpr[,1:8])
if (corType=="pearson") {
  geneTraitCor = as.data.frame(cor(datExpr, design, use = "p"))
  geneTraitP = as.data.frame(corPvalueStudent(
    as.matrix(geneTraitCor), nSamples))
} else {
  geneTraitCorA = bicorAndPvalue(datExpr, design, robustY=robustY)
  geneTraitCor = as.data.frame(geneTraitCorA$bicor)
  geneTraitP   = as.data.frame(geneTraitCorA$p)
}

###Read correlation cyan for heat positive regulation
module = "green"
pheno = "Bd21_fungi_PCA"
modNames = substring(colnames(MEs_col), 3)
module_column = match(module, modNames)
pheno_column = match(pheno,colnames(design))
moduleGenes = bwModuleColors == module
sizeGrWindow(7, 7)
par(mfrow = c(1,1))
#??geneTraitSignificance
pdf("/global/projectb/scratch/llei2019/CSP_Kranthi/sus_beni/MEgreen_Bd21_fungi_PCA.pdf",width = 15,height = 15)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                   abs(geneTraitCor[moduleGenes, pheno_column]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", pheno),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()
###green
module = "green"
pheno = "Bd21-3_bacteria_XT"
modNames = substring(colnames(MEs_col), 3)
module_column = match(module, modNames)
pheno_column = match(pheno,colnames(design))
moduleGenes = bwModuleColors == module
sizeGrWindow(7, 7)
par(mfrow = c(1,1))
#??geneTraitSignificance
pdf("/global/projectb/scratch/llei2019/CSP_Kranthi/sus_beni/MEgreen_Bd21_3_bacteria_XT.pdf",width = 15,height = 15)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                   abs(geneTraitCor[moduleGenes, pheno_column]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", pheno),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()
###green
module = "green"
pheno = "Bd21-3_bacteria_Mock"
modNames = substring(colnames(MEs_col), 3)
module_column = match(module, modNames)
pheno_column = match(pheno,colnames(design))
moduleGenes = bwModuleColors == module
sizeGrWindow(7, 7)
par(mfrow = c(1,1))
#??geneTraitSignificance
pdf("/global/projectb/scratch/llei2019/CSP_Kranthi/sus_beni/MEgreen_Bd21_3_bacteria_Mock.pdf",width = 15,height = 15)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                   abs(geneTraitCor[moduleGenes, pheno_column]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", pheno),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()
#####
#####
###Read correlation turquoise for heat positive regulation
#module = "turquoise"
#pheno = "drought"
#modNames = substring(colnames(MEs_col), 3)
#module_column = match(module, modNames)
#pheno_column = match(pheno,colnames(traitData))
#moduleGenes = bwModuleColors == module
#sizeGrWindow(7, 7)
#par(mfrow = c(1,1))
#??geneTraitSignificance
#pdf("/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/MEblue_quoise_drought_combined_abiotic_VST_distachyon.pdf",width = 15,height = 15)
#verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
#                   abs(geneTraitCor[moduleGenes, pheno_column]),
#                   xlab = paste("Module Membership in", module, "module"),
#                   ylab = paste("Gene significance for", pheno),
#                   main = paste("Module membership vs. gene significance\n"),
#                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
#dev.off()

###
#In R:
#Intramodular connectivity, module membership, and screening for intramodular hub genes
#Intramodular connectivity
connet=abs(cor(datExpr,use="p"))^6
Alldegrees1=intramodularConnectivity(connet, moduleColors)
head(Alldegrees1)
#Relationship between gene significance and intramodular connectivity
#Here is doing the module and trait association
which.module="green"
fungi= as.data.frame(design[,2])
names(fungi) = "Bd21_fungi_PCA"
GS1 = as.numeric(cor(datExpr, fungi, use = "p"))
GeneSignificance=abs(GS1)
quantile(GS1,probs = seq(0, 1, 0.01))
#0%           1%           2%           3%           4%           5% 
#  -0.860621312 -0.680473427 -0.632524721 -0.599909001 -0.576758801 -0.553280183 
#6%           7%           8%           9%          10%          11% 
#  -0.532407545 -0.515274425 -0.498998121 -0.484476876 -0.471777964 -0.458154755 
#12%          13%          14%          15%          16%          17% 
#  -0.447399986 -0.435097063 -0.425101957 -0.414745522 -0.404469041 -0.393667085 
#18%          19%          20%          21%          22%          23% 
#  -0.384112454 -0.375080571 -0.365985227 -0.358345646 -0.350855860 -0.343757757 
#24%          25%          26%          27%          28%          29% 
#  -0.335326619 -0.327417369 -0.319564153 -0.310266186 -0.301980841 -0.292801432 
#30%          31%          32%          33%          34%          35% 
#  -0.282707652 -0.272281856 -0.263613259 -0.251584785 -0.239527249 -0.228824444 
#36%          37%          38%          39%          40%          41% 
#  -0.215679891 -0.205012394 -0.194503833 -0.181083914 -0.168739846 -0.152572212 
#42%          43%          44%          45%          46%          47% 
#  -0.138596297 -0.124663802 -0.109933016 -0.094193193 -0.079921754 -0.063916805 
#48%          49%          50%          51%          52%          53% 
#  -0.047510572 -0.033035777 -0.016934462 -0.001916816  0.013781407  0.030353664 
#54%          55%          56%          57%          58%          59% 
#  0.047223886  0.061720862  0.074843595  0.091533373  0.105647387  0.120853517 
#60%          61%          62%          63%          64%          65% 
#  0.135181480  0.150818169  0.167183783  0.182492514  0.196454954  0.209812093 
#66%          67%          68%          69%          70%          71% 
#  0.221511849  0.233468408  0.246052950  0.258843558  0.271944233  0.283692324 
#72%          73%          74%          75%          76%          77% 
#  0.295382015  0.306414866  0.317715223  0.329277747  0.341022148  0.352744655 
#78%          79%          80%          81%          82%          83% 
#  0.362281608  0.373539763  0.383691800  0.395522581  0.406976344  0.417271094 
#84%          85%          86%          87%          88%          89% 
#  0.428111633  0.438319772  0.448350652  0.458980732  0.469836131  0.484371444 
#90%          91%          92%          93%          94%          95% 
#  0.498566359  0.513265843  0.527360519  0.542749454  0.559086116  0.577608472 
#96%          97%          98%          99%         100% 
#0.595565851  0.619017562  0.647854695  0.698435090  0.857100006 

#Generalizing intramodular connectivity for all genes on the array
datKME=signedKME(datExpr, MEs, outputColumnName="MM.")
#Display the first few rows of the data frame
head(datKME)
#Finding genes with high gene significance and high intramodular connectivity in specific modules
#abs(GS1)>.8 #adjust parameter based on actual situations
#abs(datKME$MM.black)>.8 #at least larger than 0.8
quantile(abs(datKME$MM.green),probs=seq(0,1,0.01))
#quantile(abs(datKME$MM.green),probs=seq(0,1,0.01))
#0%           1%           2%           3%           4%           5% 
#  2.031462e-05 6.330906e-03 1.227046e-02 1.801202e-02 2.382275e-02 3.097812e-02 
#6%           7%           8%           9%          10%          11% 
#  3.789011e-02 4.398869e-02 4.981650e-02 5.590288e-02 6.305847e-02 6.921211e-02 
#12%          13%          14%          15%          16%          17% 
#  7.511029e-02 8.148059e-02 8.804477e-02 9.432843e-02 1.000621e-01 1.062745e-01 
#18%          19%          20%          21%          22%          23% 
#  1.113425e-01 1.180336e-01 1.238950e-01 1.294604e-01 1.357057e-01 1.419054e-01 
#24%          25%          26%          27%          28%          29% 
#  1.498401e-01 1.562455e-01 1.619815e-01 1.686345e-01 1.748645e-01 1.811216e-01 
#30%          31%          32%          33%          34%          35% 
#  1.885857e-01 1.954869e-01 2.032078e-01 2.104219e-01 2.173745e-01 2.250208e-01 
#36%          37%          38%          39%          40%          41% 
#  2.317862e-01 2.389813e-01 2.453869e-01 2.517888e-01 2.578363e-01 2.634387e-01 
#42%          43%          44%          45%          46%          47% 
#  2.707447e-01 2.776263e-01 2.842557e-01 2.910413e-01 2.982143e-01 3.039836e-01 
#48%          49%          50%          51%          52%          53% 
#  3.111153e-01 3.174785e-01 3.236871e-01 3.310809e-01 3.370617e-01 3.437081e-01 
#54%          55%          56%          57%          58%          59% 
#  3.501841e-01 3.579712e-01 3.657251e-01 3.726287e-01 3.803208e-01 3.870735e-01 
#60%          61%          62%          63%          64%          65% 
#  3.945300e-01 4.018860e-01 4.085529e-01 4.155821e-01 4.224601e-01 4.298464e-01 
#66%          67%          68%          69%          70%          71% 
#  4.360888e-01 4.419362e-01 4.487321e-01 4.549778e-01 4.624692e-01 4.697983e-01 
#72%          73%          74%          75%          76%          77% 
#  4.768282e-01 4.831869e-01 4.898871e-01 4.976992e-01 5.047800e-01 5.107301e-01 
#78%          79%          80%          81%          82%          83% 
#  5.169809e-01 5.238996e-01 5.302545e-01 5.367206e-01 5.436206e-01 5.494383e-01 
#84%          85%          86%          87%          88%          89% 
#  5.562694e-01 5.648972e-01 5.724731e-01 5.802064e-01 5.891093e-01 5.981515e-01 
#90%          91%          92%          93%          94%          95% 
#  6.092596e-01 6.203902e-01 6.327628e-01 6.458615e-01 6.662434e-01 6.865891e-01 
#96%          97%          98%          99%         100% 
#7.117339e-01 7.501331e-01 7.966119e-01 8.715988e-01 9.793678e-01 

#> FilterGenes= abs(GS1)>0.31 & abs(datKME$MM.green)>0.66
#> table(FilterGenes)
#FilterGenes
#FALSE  TRUE 
#18825    66 
#GS1 98% quantile && datKME$MM.green 98%
FilterGenes= abs(GS1)>0.65 & abs(datKME$MM.green)>0.8
#GS1 90% quantile && datKME$MM.green 90%
#FilterGenes= abs(GS1)>0.45 & abs(datKME$MM.green)>0.68
table(FilterGenes)
##FilterGenes
#FALSE  TRUE 
#14961    25 
#FilterGenes
#FALSE  TRUE 
#14916    70 
hubgenes <- rownames(datKME)[FilterGenes]
hubgenes
#[1] "Bradi1g13080.v3.2" "Bradi1g19230.v3.2" "Bradi1g23710.v3.2"
#[4] "Bradi1g43880.v3.2" "Bradi1g48780.v3.2" "Bradi1g53947.v3.2"
#[7] "Bradi1g54760.v3.2" "Bradi1g57580.v3.2" "Bradi1g60110.v3.2"
#[10] "Bradi1g70730.v3.2" "Bradi2g01030.v3.2" "Bradi2g01130.v3.2"
#[13] "Bradi2g12580.v3.2" "Bradi2g26752.v3.2" "Bradi2g41110.v3.2"
#[16] "Bradi2g44320.v3.2" "Bradi2g51570.v3.2" "Bradi2g53700.v3.2"
#[19] "Bradi3g01779.v3.2" "Bradi3g27500.v3.2" "Bradi3g39980.v3.2"
#[22] "Bradi3g59830.v3.2" "Bradi4g05040.v3.2" "Bradi4g24590.v3.2"
#[25] "Bradi5g16770.v3.2"
#?chooseTopHubInEachModule
#check the green and fungi control
which.module="green"
fungimock= as.data.frame(design[,1])
names(fungimock) = "Bd21_fungi_Mock"
GS1 = as.numeric(cor(datExpr, fungimock, use = "p"))
GeneSignificance=abs(GS1)
quantile(GS1,probs = seq(0, 1, 0.01))
#0%           1%           2%           3%           4%           5%           6% 
#-0.784972454 -0.647182899 -0.604516903 -0.578255938 -0.553595513 -0.531284203 -0.511217152 
#7%           8%           9%          10%          11%          12%          13% 
#-0.495355450 -0.475825936 -0.461859529 -0.448657768 -0.433995969 -0.418949943 -0.405360654 
#14%          15%          16%          17%          18%          19%          20% 
#-0.393505368 -0.382445169 -0.371044723 -0.358139448 -0.345621213 -0.336342290 -0.325862518 
#21%          22%          23%          24%          25%          26%          27% 
#-0.315656855 -0.305916260 -0.295358695 -0.284547540 -0.274635196 -0.265361774 -0.253983754 
#28%          29%          30%          31%          32%          33%          34% 
#-0.244189683 -0.235340615 -0.225906002 -0.216729659 -0.207443603 -0.197442370 -0.186445196 
#35%          36%          37%          38%          39%          40%          41% 
#-0.175462687 -0.164277519 -0.154828940 -0.144153956 -0.133509371 -0.123890144 -0.113070326 
#42%          43%          44%          45%          46%          47%          48% 
#-0.099454266 -0.087705411 -0.075794989 -0.063226987 -0.048370732 -0.036808439 -0.022952491 
#49%          50%          51%          52%          53%          54%          55% 
#-0.009210474  0.005590815  0.022066130  0.036836436  0.051460976  0.066516668  0.082317671 
#56%          57%          58%          59%          60%          61%          62% 
#0.097907195  0.112762707  0.130745328  0.147332983  0.162425368  0.179659470  0.196173139 
#63%          64%          65%          66%          67%          68%          69% 
#0.211954404  0.227904677  0.244251685  0.260787514  0.277320933  0.293872949  0.311889543 
#70%          71%          72%          73%          74%          75%          76% 
#0.329024671  0.345630795  0.360448191  0.373357712  0.386869793  0.399922695  0.413092481 
#77%          78%          79%          80%          81%          82%          83% 
#0.427674385  0.439366074  0.452235173  0.461934803  0.473607130  0.484934145  0.496960309 
#84%          85%          86%          87%          88%          89%          90% 
#0.509589546  0.520579575  0.530499976  0.542269749  0.553271737  0.565025926  0.576497572 
#91%          92%          93%          94%          95%          96%          97% 
#0.589237638  0.601269956  0.613197820  0.625801079  0.639068240  0.653714552  0.672503682 
#98%          99%         100% 
#0.691196102  0.725115372  0.842729286 
#GS1 98% quantile && datKME$MM.green 98%
FilterGenes= abs(GS1)>0.69 & abs(datKME$MM.green)>0.8
table(FilterGenes)
#FilterGenes
#FALSE 
#14986 
#hubgenes <- rownames(datKME)[FilterGenes]
#####conlucion: there is no hub genes which need to excluded from the file

#Meganta:
which.module="magenta"
bacteria = as.data.frame(design[,4])
names(bacteria) = "bacteria"
GS1 = as.numeric(cor(datExpr, bacteria, use = "p"))
GeneSignificance=abs(GS1)
quantile(GS1,probs = seq(0, 1, 0.01))
#0%           1%           2%           3%           4%           5%           6% 
#-0.815852986 -0.662823786 -0.627545927 -0.603020680 -0.582123987 -0.562601417 -0.544843726 
#7%           8%           9%          10%          11%          12%          13% 
#-0.530985504 -0.516979837 -0.505339789 -0.491321517 -0.478674158 -0.467593341 -0.455646566 
#14%          15%          16%          17%          18%          19%          20% 
#-0.446259189 -0.434339683 -0.425552961 -0.415371255 -0.406889591 -0.396429249 -0.387131135 
#21%          22%          23%          24%          25%          26%          27% 
#-0.378087529 -0.369234430 -0.361169744 -0.353106474 -0.345011180 -0.335851530 -0.327606448 
#28%          29%          30%          31%          32%          33%          34% 
#-0.319784314 -0.310978546 -0.303012004 -0.296431292 -0.287438635 -0.280066817 -0.272467374 
#35%          36%          37%          38%          39%          40%          41% 
#-0.264322135 -0.256916453 -0.250035117 -0.242560984 -0.233916819 -0.224143722 -0.215374575 
#42%          43%          44%          45%          46%          47%          48% 
#-0.206228773 -0.196344515 -0.186045293 -0.177414609 -0.168280468 -0.157809552 -0.147268167 
#49%          50%          51%          52%          53%          54%          55% 
#-0.136483849 -0.123775716 -0.110426350 -0.098754281 -0.086600044 -0.073433345 -0.061453385 
#56%          57%          58%          59%          60%          61%          62% 
#-0.047662780 -0.033885328 -0.019918921 -0.005056398  0.009032613  0.026442699  0.041705601 
#63%          64%          65%          66%          67%          68%          69% 
#0.056258421  0.071924315  0.087929257  0.104025887  0.120188442  0.139276024  0.158425132 
#70%          71%          72%          73%          74%          75%          76% 
#0.177457557  0.196858042  0.217584395  0.232878484  0.249937380  0.268255322  0.287130444 
#77%          78%          79%          80%          81%          82%          83% 
#0.305062722  0.321533439  0.338584651  0.354385731  0.369500644  0.387764944  0.400011716 
#84%          85%          86%          87%          88%          89%          90% 
#0.414616540  0.429842445  0.444382206  0.460734707  0.476664505  0.493533443  0.507730847 
#91%          92%          93%          94%          95%          96%          97% 
#0.522684183  0.539032432  0.553843572  0.569668960  0.587493495  0.608048669  0.630428387 
#98%          99%         100% 
#0.655556884  0.690507387  0.983568181  

#Generalizing intramodular connectivity for all genes on the array
datKME=signedKME(datExpr, MEs, outputColumnName="MM.")
#Display the first few rows of the data frame
head(datKME)
#Finding genes with high gene significance and high intramodular connectivity in specific modules
#abs(GS1)>.8 #adjust parameter based on actual situations
#abs(datKME$MM.black)>.8 #at least larger than 0.8
quantile(abs(datKME$MM.magenta),probs=seq(0,1,0.01))
#0%           1%           2%           3%           4%           5%           6% 
#5.604077e-05 6.976416e-03 1.432359e-02 2.115791e-02 2.845793e-02 3.507661e-02 4.230527e-02 
#7%           8%           9%          10%          11%          12%          13% 
#4.906058e-02 5.573578e-02 6.246657e-02 6.910184e-02 7.556644e-02 8.135781e-02 8.808256e-02 
#14%          15%          16%          17%          18%          19%          20% 
#9.387990e-02 1.003347e-01 1.066769e-01 1.125770e-01 1.188312e-01 1.258133e-01 1.323055e-01 
#21%          22%          23%          24%          25%          26%          27% 
#1.387783e-01 1.458320e-01 1.522584e-01 1.578819e-01 1.639333e-01 1.696789e-01 1.758055e-01 
#28%          29%          30%          31%          32%          33%          34% 
#1.817226e-01 1.878153e-01 1.942802e-01 1.992811e-01 2.048051e-01 2.112953e-01 2.175698e-01 
#35%          36%          37%          38%          39%          40%          41% 
#2.230384e-01 2.279338e-01 2.334062e-01 2.394399e-01 2.451134e-01 2.515994e-01 2.571578e-01 
#42%          43%          44%          45%          46%          47%          48% 
#2.635748e-01 2.697900e-01 2.757226e-01 2.823502e-01 2.878089e-01 2.936573e-01 2.991444e-01 
#49%          50%          51%          52%          53%          54%          55% 
#3.044644e-01 3.103963e-01 3.171185e-01 3.229876e-01 3.290688e-01 3.365143e-01 3.430791e-01 
#56%          57%          58%          59%          60%          61%          62% 
#3.496742e-01 3.556594e-01 3.623643e-01 3.688611e-01 3.766011e-01 3.827116e-01 3.886293e-01 
#63%          64%          65%          66%          67%          68%          69% 
#3.948348e-01 4.011689e-01 4.072192e-01 4.142200e-01 4.204574e-01 4.262851e-01 4.320967e-01 
#70%          71%          72%          73%          74%          75%          76% 
#4.388984e-01 4.448820e-01 4.517318e-01 4.587285e-01 4.647725e-01 4.709801e-01 4.773668e-01 
#77%          78%          79%          80%          81%          82%          83% 
#4.838050e-01 4.915321e-01 4.983472e-01 5.054055e-01 5.134566e-01 5.204896e-01 5.278000e-01 
#84%          85%          86%          87%          88%          89%          90% 
#5.359800e-01 5.455407e-01 5.552490e-01 5.645265e-01 5.747009e-01 5.837893e-01 5.947253e-01 
#91%          92%          93%          94%          95%          96%          97% 
#6.050728e-01 6.170533e-01 6.312128e-01 6.430819e-01 6.600168e-01 6.814066e-01 7.023100e-01 
#98%          99%         100% 
#7.345881e-01 7.840931e-01 9.453976e-01 

#> FilterGenes= abs(GS1)>0.31 & abs(datKME$MM.green)>0.66
#> table(FilterGenes)
#FilterGenes
#FALSE  TRUE 
#18825    66 
#GS1 95% quantile && datKME$MM.green 95%
FilterGenes= abs(GS1)>0.69 & abs(datKME$MM.green)>0.78
#GS1 90% quantile && datKME$MM.green 90%
#FilterGenes= abs(GS1)>0.45 & abs(datKME$MM.green)>0.68
table(FilterGenes)
#FilterGenes
#FALSE  TRUE 
#14966    20 
hubgenes <- rownames(datKME)[FilterGenes]
hubgenes
#[1] "Bradi1g25737.v3.2" "Bradi1g30860.v3.2" "Bradi1g32785.v3.2" "Bradi1g45150.v3.2"
#[5] "Bradi2g00486.v3.2" "Bradi2g20127.v3.2" "Bradi2g43910.v3.2" "Bradi2g45220.v3.2"
#[9] "Bradi2g47014.v3.2" "Bradi3g22240.v3.2" "Bradi3g34490.v3.2" "Bradi3g59835.v3.2"
#[13] "Bradi4g06280.v3.2" "Bradi4g08340.v3.2" "Bradi4g32830.v3.2" "Bradi4g39611.v3.2"
#[17] "Bradi5g01167.v3.2"

#check the magenta and bateria control
which.module="magenta"
bactariamock= as.data.frame(design[,3])
names(bactariamock) = "Bd21-3_bacteria_Mock"
GS1 = as.numeric(cor(datExpr, bactariamock, use = "p"))
GeneSignificance=abs(GS1)
quantile(GS1,probs = seq(0, 1, 0.01))
#0%           1%           2%           3%           4%           5%           6% 
#-0.700710811 -0.577966686 -0.544828492 -0.521422556 -0.500179058 -0.483897080 -0.468782323 
#7%           8%           9%          10%          11%          12%          13% 
#-0.456587872 -0.445742272 -0.433518817 -0.423503987 -0.413941789 -0.403631342 -0.393460781 
#14%          15%          16%          17%          18%          19%          20% 
#-0.384216144 -0.376120435 -0.368000149 -0.359679642 -0.351636563 -0.344684006 -0.338152352 
#21%          22%          23%          24%          25%          26%          27% 
#-0.331172476 -0.324406165 -0.317567687 -0.311031271 -0.305108093 -0.298013561 -0.291237659 
#28%          29%          30%          31%          32%          33%          34% 
#-0.285516936 -0.279402694 -0.272147047 -0.266148596 -0.259674775 -0.252345326 -0.245124454 
#35%          36%          37%          38%          39%          40%          41% 
#-0.238796148 -0.231478213 -0.224867256 -0.217100727 -0.208897687 -0.201970944 -0.192850146 
#42%          43%          44%          45%          46%          47%          48% 
#-0.185114665 -0.176061633 -0.167450580 -0.159174864 -0.148988281 -0.139137791 -0.129716313 
#49%          50%          51%          52%          53%          54%          55% 
#-0.118742801 -0.109019093 -0.099780690 -0.090422414 -0.078805763 -0.068265591 -0.055267052 
#56%          57%          58%          59%          60%          61%          62% 
#-0.044529034 -0.033638007 -0.018794918 -0.005495950  0.006627746  0.018866770  0.030681251 
#63%          64%          65%          66%          67%          68%          69% 
#0.043777044  0.057692931  0.072598650  0.087069835  0.098909467  0.113774829  0.128626579 
#70%          71%          72%          73%          74%          75%          76% 
#0.142845215  0.156179963  0.169984118  0.186212403  0.202282864  0.219476385  0.234415992 
#77%          78%          79%          80%          81%          82%          83% 
#0.251882893  0.265263101  0.283302163  0.298620238  0.313323820  0.326362389  0.341812371 
#84%          85%          86%          87%          88%          89%          90% 
#0.357251786  0.372784478  0.387798516  0.402939132  0.418752009  0.436360239  0.453420998 
#91%          92%          93%          94%          95%          96%          97% 
#0.468077075  0.484890382  0.505249072  0.521751786  0.539679001  0.560493640  0.584150711 
#98%          99%         100% 
#0.604371923  0.638943624  0.833825187 
#GS1 98% quantile && datKME$MM.green 98%
FilterGenes= abs(GS1)>0.54 & abs(datKME$MM.green)>0.78
table(FilterGenes)
#FilterGenes
#FALSE 
#14986 
####so all of the genes listed are the hub genes

#green yellow
which.module="greenyellow"
virus = as.data.frame(design[,6])
names(virus) = "virus"
GS1 = as.numeric(cor(datExpr, virus, use = "p"))
GeneSignificance=abs(GS1)
quantile(GS1,probs = seq(0, 1, 0.01))
#0%           1%           2%           3%           4%           5%           6% 
#-0.784171863 -0.615713985 -0.588820091 -0.567655953 -0.551550832 -0.537999747 -0.525060121 
#7%           8%           9%          10%          11%          12%          13% 
#-0.513420212 -0.503634979 -0.491852939 -0.480071054 -0.470163445 -0.459881639 -0.451573017 
#14%          15%          16%          17%          18%          19%          20% 
#-0.442952448 -0.432803465 -0.423551590 -0.413452212 -0.404059415 -0.394486866 -0.385556476 
#21%          22%          23%          24%          25%          26%          27% 
#-0.377038706 -0.368393976 -0.358912810 -0.350281418 -0.341663239 -0.332132123 -0.322264074 
#28%          29%          30%          31%          32%          33%          34% 
#-0.312323189 -0.301790000 -0.291269127 -0.279666365 -0.269387960 -0.258425155 -0.247037779 
#35%          36%          37%          38%          39%          40%          41% 
#-0.236812901 -0.224531598 -0.211683723 -0.198125940 -0.186034946 -0.172158979 -0.158667809 
#42%          43%          44%          45%          46%          47%          48% 
#-0.146164490 -0.134494802 -0.120416017 -0.107394028 -0.092413710 -0.078946302 -0.067005633 
#49%          50%          51%          52%          53%          54%          55% 
#-0.049147670 -0.034124935 -0.017262506 -0.001956911  0.012067859  0.029350594  0.047177680 
#56%          57%          58%          59%          60%          61%          62% 
#0.064364010  0.081463820  0.098667177  0.116796901  0.135482299  0.152146794  0.170877958 
#63%          64%          65%          66%          67%          68%          69% 
#0.188745642  0.205438786  0.221913882  0.239044789  0.255884352  0.272604368  0.289558952 
#70%          71%          72%          73%          74%          75%          76% 
#0.304758110  0.317764096  0.331613629  0.346161286  0.359514503  0.371564009  0.382858156 
#77%          78%          79%          80%          81%          82%          83% 
#0.395753440  0.406354084  0.416550851  0.428589765  0.440575098  0.451758246  0.463263700 
#84%          85%          86%          87%          88%          89%          90% 
#0.473394024  0.483991242  0.495892460  0.505852314  0.514427164  0.524760344  0.536510371 
#91%          92%          93%          94%          95%          96%          97% 
#0.546867039  0.558067673  0.569350706  0.580553542  0.593943821  0.606478753  0.620411717 
#98%          99%         100% 
#0.637649555  0.662293839  0.813758659 
#Generalizing intramodular connectivity for all genes on the array
datKME=signedKME(datExpr, MEs, outputColumnName="MM.")
#Display the first few rows of the data frame
head(datKME)
#Finding genes with high gene significance and high intramodular connectivity in specific modules
#abs(GS1)>.8 #adjust parameter based on actual situations
#abs(datKME$MM.black)>.8 #at least larger than 0.8
quantile(abs(datKME$MM.greenyellow),probs=seq(0,1,0.01))
#0%           1%           2%           3%           4%           5%           6% 
#5.217572e-05 5.466091e-03 1.099943e-02 1.674087e-02 2.289324e-02 2.845480e-02 3.392066e-02 
#7%           8%           9%          10%          11%          12%          13% 
#3.945838e-02 4.555590e-02 5.205939e-02 5.818062e-02 6.330760e-02 6.949220e-02 7.505088e-02 
#14%          15%          16%          17%          18%          19%          20% 
#8.057767e-02 8.646487e-02 9.204471e-02 9.799614e-02 1.037870e-01 1.092472e-01 1.145801e-01 
#21%          22%          23%          24%          25%          26%          27% 
#1.194886e-01 1.240432e-01 1.296598e-01 1.347004e-01 1.399830e-01 1.451295e-01 1.496390e-01 
#28%          29%          30%          31%          32%          33%          34% 
#1.550940e-01 1.599843e-01 1.651947e-01 1.705593e-01 1.751354e-01 1.803788e-01 1.852455e-01 
#35%          36%          37%          38%          39%          40%          41% 
#1.906641e-01 1.952971e-01 2.003363e-01 2.050140e-01 2.092253e-01 2.137390e-01 2.186647e-01 
#42%          43%          44%          45%          46%          47%          48% 
#2.231873e-01 2.270518e-01 2.317932e-01 2.359197e-01 2.401834e-01 2.444336e-01 2.483332e-01 
#49%          50%          51%          52%          53%          54%          55% 
#2.524759e-01 2.570352e-01 2.612638e-01 2.657247e-01 2.702995e-01 2.744006e-01 2.790406e-01 
#56%          57%          58%          59%          60%          61%          62% 
#2.838520e-01 2.877381e-01 2.921087e-01 2.969974e-01 3.010980e-01 3.055218e-01 3.098052e-01 
#63%          64%          65%          66%          67%          68%          69% 
#3.147684e-01 3.190463e-01 3.231902e-01 3.278176e-01 3.322722e-01 3.365117e-01 3.416466e-01 
#70%          71%          72%          73%          74%          75%          76% 
#3.461207e-01 3.509406e-01 3.554997e-01 3.602005e-01 3.649549e-01 3.692906e-01 3.740202e-01 
#77%          78%          79%          80%          81%          82%          83% 
#3.795264e-01 3.845216e-01 3.901514e-01 3.954662e-01 4.012002e-01 4.074418e-01 4.135949e-01 
#84%          85%          86%          87%          88%          89%          90% 
#4.206298e-01 4.267697e-01 4.346869e-01 4.428522e-01 4.512602e-01 4.586362e-01 4.683348e-01 
#91%          92%          93%          94%          95%          96%          97% 
#4.784536e-01 4.893027e-01 5.014823e-01 5.125620e-01 5.269054e-01 5.451025e-01 5.721172e-01 
#98%          99%         100% 
#6.052210e-01 6.588426e-01 9.548486e-01 
#top1% threshold
FilterGenes= abs(GS1)>0.66 & abs(datKME$MM.green)>0.66
#GS1 90% quantile && datKME$MM.green 90%
#FilterGenes= abs(GS1)>0.45 & abs(datKME$MM.green)>0.68
table(FilterGenes)
#FilterGenes
#FALSE  TRUE 
#14964    22 
hubgenes <- rownames(datKME)[FilterGenes]
hubgenes
#top1% #[1] "Bradi1g74310.v3.2" "Bradi2g22180.v3.2" "Bradi3g36357.v3.2" "Bradi4g23695.v3.2"
#top2%
#[1] "Bradi1g04080.v3.2" "Bradi1g20500.v3.2" "Bradi1g29267.v3.2" "Bradi1g42556.v3.2"
#[5] "Bradi1g58290.v3.2" "Bradi1g63520.v3.2" "Bradi1g74310.v3.2" "Bradi1g75580.v3.2"
#[9] "Bradi2g22180.v3.2" "Bradi2g34470.v3.2" "Bradi2g40410.v3.2" "Bradi2g58750.v3.2"
#[13] "Bradi2g60710.v3.2" "Bradi3g03520.v3.2" "Bradi3g36357.v3.2" "Bradi3g38670.v3.2"
#[17] "Bradi4g04270.v3.2" "Bradi4g17950.v3.2" "Bradi4g23695.v3.2" "Bradi4g38895.v3.2"
#[21] "Bradi4g42980.v3.2" "Bradi5g00640.v3.2"
#check the magenta and bateria control
which.module="greenyellow"
virusmock= as.data.frame(design[,5])
names(virusmock) = "Bd21-3_virus_Mock"
GS1 = as.numeric(cor(datExpr, bactariamock, use = "p"))
GeneSignificance=abs(GS1)
quantile(GS1,probs = seq(0, 1, 0.01))
#0%           1%           2%           3%           4%           5%           6% 
#-0.700710811 -0.577966686 -0.544828492 -0.521422556 -0.500179058 -0.483897080 -0.468782323 
#7%           8%           9%          10%          11%          12%          13% 
#-0.456587872 -0.445742272 -0.433518817 -0.423503987 -0.413941789 -0.403631342 -0.393460781 
#14%          15%          16%          17%          18%          19%          20% 
#-0.384216144 -0.376120435 -0.368000149 -0.359679642 -0.351636563 -0.344684006 -0.338152352 
#21%          22%          23%          24%          25%          26%          27% 
#-0.331172476 -0.324406165 -0.317567687 -0.311031271 -0.305108093 -0.298013561 -0.291237659 
#28%          29%          30%          31%          32%          33%          34% 
#-0.285516936 -0.279402694 -0.272147047 -0.266148596 -0.259674775 -0.252345326 -0.245124454 
#35%          36%          37%          38%          39%          40%          41% 
#-0.238796148 -0.231478213 -0.224867256 -0.217100727 -0.208897687 -0.201970944 -0.192850146 
#42%          43%          44%          45%          46%          47%          48% 
#-0.185114665 -0.176061633 -0.167450580 -0.159174864 -0.148988281 -0.139137791 -0.129716313 
#49%          50%          51%          52%          53%          54%          55% 
#-0.118742801 -0.109019093 -0.099780690 -0.090422414 -0.078805763 -0.068265591 -0.055267052 
#56%          57%          58%          59%          60%          61%          62% 
#-0.044529034 -0.033638007 -0.018794918 -0.005495950  0.006627746  0.018866770  0.030681251 
#63%          64%          65%          66%          67%          68%          69% 
#0.043777044  0.057692931  0.072598650  0.087069835  0.098909467  0.113774829  0.128626579 
#70%          71%          72%          73%          74%          75%          76% 
#0.142845215  0.156179963  0.169984118  0.186212403  0.202282864  0.219476385  0.234415992 
#77%          78%          79%          80%          81%          82%          83% 
#0.251882893  0.265263101  0.283302163  0.298620238  0.313323820  0.326362389  0.341812371 
#84%          85%          86%          87%          88%          89%          90% 
#0.357251786  0.372784478  0.387798516  0.402939132  0.418752009  0.436360239  0.453420998 
#91%          92%          93%          94%          95%          96%          97% 
#0.468077075  0.484890382  0.505249072  0.521751786  0.539679001  0.560493640  0.584150711 
#98%          99%         100% 
#0.604371923  0.638943624  0.833825187  
#GS1 99% quantile && datKME$MM.greenyellow 99%
FilterGenes= abs(GS1)>0.64 & abs(datKME$MM.green)>0.66
table(FilterGenes)
#FilterGenes
#FALSE 
#14986 



#export the green module
#In R:
#Export the network
#Recalculate topological overlap
TOM = TOMsimilarityFromExpr(datExpr, power = 12) 
# Select module
module = "grey" #change to green
# Select module probes
probes = colnames(datExpr)
inModule = (moduleColors==module)
modProbes = probes[inModule] 
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)

#Export the network into edge and node list files Cytoscape can read
#default threshold = 0.5, we could adjust parameter based on actual situations or in Cytoscape
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("/global/projectb/scratch/llei2019/CSP_Kranthi/sus_beni/CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("/global/projectb/scratch/llei2019/CSP_Kranthi/sus_beni/CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes, 
                               nodeAttr = moduleColors[inModule])

#Screen the top genes
nTop = 20
IMConn = softConnectivity(datExpr[, modProbes])
top = (rank(-IMConn) <= nTop)
filter <- modTOM[top, top]

cyt = exportNetworkToCytoscape(filter,
                               edgeFile = paste("/global/projectb/scratch/llei2019/CSP_Kranthi/sus_beni/CytoscapeInput-edges-filter-top20-", paste(module, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("/global/projectb/scratch/llei2019/CSP_Kranthi/sus_beni/CytoscapeInput-nodes-filter-top20-", paste(module, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = rownames(filter), 
                               nodeAttr = moduleColors[inModule][1:nTop])

#In R:
#Extract gene IDs in specific module
#Select module
module = "green"
#Select module probes (gene ID)
probes = colnames(datExpr)
inModule = (moduleColors == module)
modProbes = probes[inModule]
write.table(modProbes,file="/global/projectb/scratch/llei2019/CSP_Kranthi/green_0.25p6",sep="\t",quote=F,row.names=F,col.names=F)


#####extract top 50 hub genes for each module!!!
MEs <- moduleEigengenes(datExpr, bwModuleColors)$eigengenes
head(MEs)
kMEs <- signedKME(datExpr, MEs)
head(kMEs)
# rank the genes for each module on kMEs
rankGenes <- function(x){
  kMErank <- rank(-kMEs[ ,x])
  genes <- rownames(kMEs)
  genes <- genes[order(kMErank)]
  #genes[1:50]#top 50 hub genes!!!
}

topGenes <- lapply(1:ncol(kMEs), rankGenes)
nrow(topGenes)
# Get the top results in a data.frame
topGenes <- do.call(cbind, topGenes)
colnames(topGenes) <- substr(colnames(kMEs), start=4, stop=30)
nrow(topGenes)
head(topGenes)
hubgenes <- data.frame(topGenes)
write.csv(hubgenes, file = "/global/projectb/scratch/llei2019/CSP_Kranthi/sus_beni/allhubgenes_list_VST_suspatho_mads75_encoding1_0.25power12.csv")

#
#adj_mat<-adjacency(datExpr,power=12)
adj <- TOM
adj[adj > 0.1] = 1
adj[adj != 1] = 0
network <- graph.adjacency(adj)
network <- simplify(network)  # removes self-loops
#results <- blockwiseModules(data, power=6, TOMType="unsigned", networkType="unsigned")
V(network)$color <- bwnet$colors
par(mar=c(0,0,0,0))
# remove unconnected nodes
network <- delete.vertices(network, degree(network)==0)
plot(network, layout=layout.fruchterman.reingold(network), edge.arrow.size = 0.2)

###Plot the network for all of the modules
top.n.edges = 25000#3000#5000#15000#25000
min.edge = 2
adj_mat<-adjacency(datExpr,power=12)
head(adj_mat)
nrow(adj_mat)
message("number of genes present = ", nrow(adj_mat))#number of genes present = 18891 
if(!is.na(top.n.edges) & nrow(adj_mat) >= 100){
  adjacency.threshold = sort(as.numeric(adj_mat), decreasing=T)[(top.n.edges*2)+nrow(adj_mat)]
}
if(nrow(adj_mat) < 100)
  adjacency.threshold = 0.01

message("adjacency.threshold = ", adjacency.threshold)
#adjacency.threshold = 0.893484769264947#if edge=5000,0.879482033215977；edge=10000，adjacency.threshold = 0.856292297792383,edge15000,adjacency.threshold = 0.839325416948576
#25000，adjacency.threshold = 0.814426058718555,20000,adjacency.threshold = 0.826031854913583
adj_mat[adj_mat > adjacency.threshold] <- 1
adj_mat[adj_mat < adjacency.threshold] <- 0
diag(adj_mat) <- 0
rs<-rowSums(adj_mat)
if(verbose) cat("removing unconnected nodes\n")
adj_mat<-adj_mat[rs>min.edge,rs>min.edge]
message("number of genes present = ", nrow(adj_mat))#356 #461#661#839#1064#962
test <- data.frame(adj_mat)
nrow(test)#356#461#661#839#1064
#row.names(test)[modProbes,]

#turquoise
probes = colnames(datExpr)
turquoise_modules = c("turquoise")
turquoise_inModule = is.finite(match(bwModuleColors, turquoise_modules))
turquoise_modProbes = probes[turquoise_inModule]
turquoise_list <- subset(adj_mat, rownames(adj_mat) %in% turquoise_modProbes)

#pink
probes = colnames(datExpr)
pink_modules = c("pink")
pink_inModule = is.finite(match(bwModuleColors, pink_modules))
pink_modProbes = probes[pink_inModule]
pink_list <- subset(adj_mat, rownames(adj_mat) %in% pink_modProbes)#0

#purple
probes = colnames(datExpr)
purple_modules = c("purple")
purple_inModule = is.finite(match(bwModuleColors, purple_modules))
purple_modProbes = probes[purple_inModule]
purple_list <- subset(adj_mat, rownames(adj_mat) %in% purple_modProbes)

#red
probes = colnames(datExpr)
red_modules = c("red")
red_inModule = is.finite(match(bwModuleColors, red_modules))
red_modProbes = probes[red_inModule]
red_list <- subset(adj_mat, rownames(adj_mat) %in% red_modProbes)

#magenta
probes = colnames(datExpr)
magenta_modules = c("magenta")
magenta_inModule = is.finite(match(bwModuleColors, magenta_modules))
magenta_modProbes = probes[magenta_inModule]
magenta_list <- subset(adj_mat, rownames(adj_mat) %in% magenta_modProbes)#0

#blue
probes = colnames(datExpr)
blue_modules = c("blue")
blue_inModule = is.finite(match(bwModuleColors, blue_modules))
blue_modProbes = probes[blue_inModule]
blue_list <- subset(adj_mat, rownames(adj_mat) %in% blue_modProbes)
#green
probes = colnames(datExpr)
green_modules = c("green")
green_inModule = is.finite(match(bwModuleColors, green_modules))
green_modProbes = probes[green_inModule]
green_list <- subset(adj_mat, rownames(adj_mat) %in% green_modProbes)

#brown
probes = colnames(datExpr)
brown_modules = c("brown")
brown_inModule = is.finite(match(bwModuleColors, brown_modules))
brown_modProbes = probes[brown_inModule]
brown_list <- subset(adj_mat, rownames(adj_mat) %in% brown_modProbes)

#yellow
probes = colnames(datExpr)
yellow_modules = c("yellow")
yellow_inModule = is.finite(match(bwModuleColors, yellow_modules))
yellow_modProbes = probes[yellow_inModule]
yellow_list <- subset(adj_mat, rownames(adj_mat) %in% yellow_modProbes)

#black
probes = colnames(datExpr)
black_modules = c("black")
black_inModule = is.finite(match(bwModuleColors, black_modules))
black_modProbes = probes[black_inModule]
black_list <- subset(adj_mat, rownames(adj_mat) %in% black_modProbes)#0

#greenyellow
probes = colnames(datExpr)
greeny_modules = c("greenyellow")
greeny_inModule = is.finite(match(bwModuleColors, greeny_modules))
greeny_modProbes = probes[greeny_inModule]
greeny_list <- subset(adj_mat, rownames(adj_mat) %in% greeny_modProbes)#0

#grey
probes = colnames(datExpr)
grey_modules = c("grey")
grey_inModule = is.finite(match(bwModuleColors, grey_modules))
grey_modProbes = probes[grey_inModule]
grey_list <- subset(adj_mat, rownames(adj_mat) %in% grey_modProbes)#0

#lightcyan
#probes = colnames(datExpr)
#lightcyan_modules = c("lightcyan")
#lightcyan_inModule = is.finite(match(bwModuleColors, lightcyan_modules))
#lightcyan_modProbes = probes[lightcyan_inModule]
#lightcyan_list <- subset(adj_mat, rownames(adj_mat) %in% lightcyan_modProbes)#0


#midnightblue
#probes = colnames(datExpr)
#midnightblue_modules = c("midnightblue")
#midnightblue_inModule = is.finite(match(bwModuleColors, midnightblue_modules))
#midnightblue_modProbes = probes[midnightblue_inModule]
#midnightblue_list <- subset(adj_mat, rownames(adj_mat) %in% midnightblue_modProbes)#0



#salmon
#probes = colnames(datExpr)
#salmon_modules = c("salmon")
#salmon_inModule = is.finite(match(bwModuleColors, salmon_modules))
#salmon_modProbes = probes[salmon_inModule]
#salmon_list <- subset(adj_mat, rownames(adj_mat) %in% salmon_modProbes)#0

#tan
#probes = colnames(datExpr)
#tan_modules = c("tan")
#tan_inModule = is.finite(match(bwModuleColors, tan_modules))
#tan_modProbes = probes[tan_inModule]
#tan_list <- subset(adj_mat, rownames(adj_mat) %in% tan_modProbes)#0



#darkgreen
#probes = colnames(datExpr)
#darkgreen_modules = c("darkgreen")
#darkgreen_inModule = is.finite(match(bwModuleColors, darkgreen_modules))
#darkgreen_modProbes = probes[darkgreen_inModule]
#darkgreen_list <- subset(adj_mat, rownames(adj_mat) %in% darkgreen_modProbes)
#darkgrey
#probes = colnames(datExpr)
#darkgrey_modules = c("darkgrey")
#darkgrey_inModule = is.finite(match(bwModuleColors, darkgrey_modules))
#darkgrey_modProbes = probes[darkgrey_inModule]
#darkgrey_list <- subset(adj_mat, rownames(adj_mat) %in% darkgrey_modProbes)
#darkred
#probes = colnames(datExpr)
#darkred_modules = c("darkred")
#darkred_inModule = is.finite(match(bwModuleColors, darkred_modules))
#darkred_modProbes = probes[darkred_inModule]
#darkred_list <- subset(adj_mat, rownames(adj_mat) %in% darkred_modProbes)
#darkturquoise
#probes = colnames(datExpr)
#darkturquoise_modules = c("darkturquoise")
#darkturquoise_inModule = is.finite(match(bwModuleColors, darkturquoise_modules))
#darkturquoise_modProbes = probes[darkturquoise_inModule]
#darkturquoise_list <- subset(adj_mat, rownames(adj_mat) %in% darkturquoise_modProbes)
#grey60
#probes = colnames(datExpr)
#grey60_modules = c("grey60")
#grey60_inModule = is.finite(match(bwModuleColors, grey60_modules))
#grey60_modProbes = probes[grey60_inModule]
#grey60_list <- subset(adj_mat, rownames(adj_mat) %in% grey60_modProbes)
#lightgreen
#probes = colnames(datExpr)
#lightgreen_modules = c("lightgreen")
#lightgreen_inModule = is.finite(match(bwModuleColors, lightgreen_modules))
#lightgreen_modProbes = probes[lightgreen_inModule]
#lightgreen_list <- subset(adj_mat, rownames(adj_mat) %in% lightgreen_modProbes)
#lightyellow
#probes = colnames(datExpr)
#lightyellow_modules = c("lightyellow")
#lightyellow_inModule = is.finite(match(bwModuleColors, lightyellow_modules))
#lightyellow_modProbes = probes[lightyellow_inModule]
#lightyellow_list <- subset(adj_mat, rownames(adj_mat) %in% lightyellow_modProbes)

#orange
#probes = colnames(datExpr)
#orange_modules = c("orange")
#orange_inModule = is.finite(match(bwModuleColors, orange_modules))
#orange_modProbes = probes[orange_inModule]
#orange_list <- subset(adj_mat, rownames(adj_mat) %in% orange_modProbes)
#royalblue
#probes = colnames(datExpr)
#royalblue_modules = c("royalblue")
#royalblue_inModule = is.finite(match(bwModuleColors, royalblue_modules))
#royalblue_modProbes = probes[royalblue_inModule]
#royalblue_list <- subset(adj_mat, rownames(adj_mat) %in% royalblue_modProbes)

#rownames(yellow_list) <- subset(yellow_list, colnames(yellow_list) %in% modProbes)
#nrow(yellow_list)
#yellow_list$color <- rep("yellow", nrow(yellow_list))
#yellow_sub <- data.frame(row.names(yellow_list),yellow_list$color)
#network <- graph.adjacency(yellow_list)
#head(network)
#network <- simplify(network)  # removes self-loops
#results <- blockwiseModules(data, power=6, TOMType="unsigned", networkType="unsigned")
#V(network)$color <- results$colors
#par(mar=c(0,0,0,0))
# remove unconnected nodes
#network <- delete.vertices(network, degree(network)==0)
#pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/network_all_abiotic_VST#.pdf",width = 25,height = 25)
#plot(network, vertex.color="yellow", vertex.label="",vertex.size=2,layout=layout_in_circle(network), edge.arrow.size = 0.2)
#plot(network, vertex.color="yellow",vertex.label="",vertex.size=2, layout=layout.fruchterman.reingold(network), edge.arrow.size = 0.2)

#dev.off()
####
network <- graph.adjacency(adj_mat,mode = "undirected",weighted = TRUE,)
head(network)
E(network)
V(network)[rownames(turquoise_list)]$color <- "turquoise"
V(network)[rownames(pink_list)]$color <- "pink"
V(network)[rownames(purple_list)]$color <- "purple"
V(network)[rownames(red_list)]$color <- "red"
V(network)[rownames(magenta_list)]$color <- "magenta"
V(network)[rownames(blue_list)]$color <- "blue"
V(network)[rownames(green_list)]$color <- "green"
V(network)[rownames(brown_list)]$color <- "brown"
V(network)[rownames(yellow_list)]$color <- "yellow"
V(network)[rownames(black_list)]$color <- "black"
V(network)[rownames(greeny_list)]$color <- "greenyellow"
V(network)[rownames(grey)]$color <- "grey"
#V(network)[rownames(cyan_list)]$color <- "cyan"
#V(network)[rownames(grey60_list)]$color <- "grey60"
#V(network)[rownames(tan_list)]$color <- "tan"
#V(network)[rownames(midnightblue_list)]$color <- "midnightblue"
#V(network)[rownames(lightcyan_list)]$color <- "lightcyan"
#V(network)[rownames(salmon_list)]$color <- "salmon"
#V(network)[rownames(darkgreen_list)]$color <- "darkgreen"
#V(network)[rownames(darkgrey_list)]$color <- "darkgrey"
#V(network)[rownames(darkred_list)]$color <- "darkred"
#V(network)[rownames(darkturquoise_list)]$color <- "darkturquoise"
#V(network)[rownames(lightgreen_list)]$color <- "lightgreen"
#V(network)[rownames(lightyellow_list)]$color <- "lightyellow"
#V(network)[rownames(orange_list)]$color <- "orange"
#V(network)[rownames(royalblue_list)]$color <- "royalblue"
#all[rownames(yellow_list)]
#all[rownames(yellow_list)]
network <- simplify(network)  # removes self-loops
#results <- blockwiseModules(data, power=6, TOMType="unsigned", networkType="unsigned")
#V(network)$color <- results$colors
sizeGrWindow(10,10)
par(mar=c(0,0,0,0))
# remove unconnected nodes
network <- delete.vertices(network, degree(network)==0)
#minC <- rep(-Inf, vcount(network))
minC <- rep(-Inf, vcount(network))
maxC <- rep(Inf, vcount(network))
minC[1] <- maxC[1] <- 0
#co <- layout_with_fr(network, minx=minC, maxx=maxC,
#                     miny=minC, maxy=maxC)
co <- layout.fruchterman.reingold(network, minx=minC, maxx=maxC,
                                  miny=minC, maxy=maxC)
#layout.fruchterman.reingold
co[1,]
pdf(file = "/global/projectb/scratch/llei2019/CSP_Kranthi/network_sus_brachy_dges25000.pdf",width = 25,height = 25)
#plot(network, vertex.color="yellow", vertex.label="",vertex.size=2,layout=layout_in_circle(network), edge.arrow.size = 0.2)
#plot(network,layout=co,pch=19,arrow.mode=0,mark.border=NA, edge.color = "#d1d1e0", vertex.label="",vertex.size=2.5, layout=layout.fruchterman.reingold(network), edge.arrow.size = 0.2)
plot(network,layout=co,pch=19,arrow.mode=0,mark.border=NA, edge.color = "#d1d1e0", vertex.label="",vertex.size=2.5,edge.arrow.size = 0.2)

dev.off()

#Since the turquoise (Drought) and cyan (heat) associated with the certain abiotic treatment,
#so I have to do the GO analysis and the plot the network for those two modules

###Find all of the genes in the module:
gs<-colnames(datExpr)
#  cols<-net[[1]]
cols <- bwModuleColors
names(cols)<-gs #assign each gene into different module
head(gs)

kme <- signedKME(datExpr = datExpr, datME = MEs, outputColumnName = NULL)
head(kme)

#row.names(kme)<-gs
kmes<-sapply(1:length(gs), function(x) abs(kme[gs[x],colnames(kme) == cols[[x]]]))
kmes<-data.frame(genes = gs, cols = cols, kme = kmes)
head(kmes)
nrow(kmes)
#check the gene list for each module
#turquoise
turquoise_kmes <- kmes[kmes$cols=="turquoise",]
head(turquoise_kmes)
nrow(turquoise_kmes)
#6244
turquoise_kmes <-turquoise_kmes[order(-turquoise_kmes$kme),]
write.csv(turquoise_kmes, file = "/global/projectb/scratch/llei2019/CSP_Kranthi/sus_brachy_turquoise_kmes.csv")
###pink
pink_kmes <- kmes[kmes$cols=="pink",]
head(pink_kmes)
nrow(pink_kmes)
#347
pink_kmes <-pink_kmes[order(-pink_kmes$kme),]
write.csv(pink_kmes, file = "/global/projectb/scratch/llei2019/CSP_Kranthi/sus_brachy_pink_kmes.csv")

###purple
purple_kmes <- kmes[kmes$cols=="purple",]
head(purple_kmes)
nrow(purple_kmes)
#178
purple_kmes <-purple_kmes[order(-purple_kmes$kme),]
write.csv(purple_kmes, file = "/global/projectb/scratch/llei2019/CSP_Kranthi/sus_brachy_purple_kmes.csv")

###red
red_kmes <- kmes[kmes$cols=="red",]
head(red_kmes)
nrow(red_kmes)
#444
red_kmes <-red_kmes[order(-red_kmes$kme),]
write.csv(red_kmes, file = "/global/projectb/scratch/llei2019/CSP_Kranthi/sus_brachy_red_kmes.csv")

###Magenta
magenta_kmes <- kmes[kmes$cols=="magenta",]
head(magenta_kmes)
nrow(magenta_kmes)
#265
magenta_kmes <-magenta_kmes[order(-magenta_kmes$kme),]
write.csv(magenta_kmes, file = "/global/projectb/scratch/llei2019/CSP_Kranthi/sus_brachy_magenta_kmes.csv")

###Blue
blue_kmes <- kmes[kmes$cols=="blue",]
head(blue_kmes)
nrow(blue_kmes)
#4203
blue_kmes <-blue_kmes[order(-blue_kmes$kme),]
write.csv(blue_kmes, file = "/global/projectb/scratch/llei2019/CSP_Kranthi/sus_brachy_blue_kmes.csv")

###green
green_kmes <- kmes[kmes$cols=="green",]
head(green_kmes)
nrow(green_kmes)
#587
green_kmes <-green_kmes[order(-green_kmes$kme),]
write.csv(green_kmes, file = "/global/projectb/scratch/llei2019/CSP_Kranthi/sus_brachy_green_kmes.csv")

###brown
brown_kmes <- kmes[kmes$cols=="brown",]
head(brown_kmes)
nrow(brown_kmes)
#3801
brown_kmes <-brown_kmes[order(-brown_kmes$kme),]
write.csv(brown_kmes, file = "/global/projectb/scratch/llei2019/CSP_Kranthi/sus_brachy_brown_kmes.csv")

###yellow
yellow_kmes <- kmes[kmes$cols=="yellow",]
head(yellow_kmes)
nrow(yellow_kmes)
#2316
yellow_kmes <-yellow_kmes[order(-yellow_kmes$kme),]
write.csv(yellow_kmes, file = "/global/projectb/scratch/llei2019/CSP_Kranthi/sus_brachy_yellow_kmes.csv")

###black
black_kmes <- kmes[kmes$cols=="black",]
head(black_kmes)
nrow(black_kmes)
#354
black_kmes <-black_kmes[order(-black_kmes$kme),]
write.csv(black_kmes, file = "/global/projectb/scratch/llei2019/CSP_Kranthi/sus_brachy_black_kmes.csv")

###greenyellow
greenyellow_kmes <- kmes[kmes$cols=="greenyellow",]
head(greenyellow_kmes)
nrow(greenyellow_kmes)
#52
greenyellow_kmes <-greenyellow_kmes[order(-greenyellow_kmes$kme),]
write.csv(greenyellow_kmes, file = "/global/projectb/scratch/llei2019/CSP_Kranthi/sus_brachy_greenyellow_kmes.csv")

###grey
grey_kmes <- kmes[kmes$cols=="grey",]
head(grey_kmes)
nrow(grey_kmes)
#100
grey_kmes <-grey_kmes[order(-grey_kmes$kme),]
write.csv(grey_kmes, file = "/global/projectb/scratch/llei2019/CSP_Kranthi/sus_brachy_grey_kmes.csv")

#################
turquoise_gene_list <- rownames(turquoise_kmes)

turquoise_kmes <-turquoise_kmes[order(-turquoise_kmes$kme),]
turquoise_top20_gene_list <- turquoise_kmes[1:20,]$genes

gs<-kmes$genes[kmes$cols=="turquoise"]
#cols<-cols[kmes$kme>=kME.threshold]
datExpr_turquoise = datExpr[,gs]
head(datExpr_turquoise[,1:8])
ncol(datExpr_turquoise)
turquoise_gene_list <- data.frame(colnames(datExpr_turquoise))
#}#yellow module we have 3587 genes!!! 3951 genes for turquoise
###write into the file:
write.csv(turquoise_gene_list, file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/abiotic_turtoise_list_VST_distachyon.csv")

adj_mat<-adjacency(datExpr_turquoise,power=9)
###
head(adj_mat[,1:8])
######### Get the top edges of the blue module:
top.n.edges = 3000
min.edge = 2
#head(adj_mat)
message("number of genes present = ", nrow(adj_mat))#number of genes present = 3587 Distachyon: 3951
if(!is.na(top.n.edges) & nrow(adj_mat) >= 100){
  adjacency.threshold = sort(as.numeric(adj_mat), decreasing=T)[(top.n.edges*2)+nrow(adj_mat)]
}
if(nrow(adj_mat) < 100)
  adjacency.threshold = 0.01

message("adjacency.threshold = ", adjacency.threshold)#adjacency.threshold = 0.690966673259787#distachyon:0.511680160681778
adj_mat[adj_mat > adjacency.threshold] <- 1
adj_mat[adj_mat < adjacency.threshold] <- 0
diag(adj_mat) <- 0
rs<-rowSums(adj_mat)
if(verbose) cat("removing unconnected nodes\n")
adj_mat<-adj_mat[rs>min.edge,rs>min.edge]
message("number of genes present = ", nrow(adj_mat))#number of genes present = 343
#distachyon: 351
#test <- data.frame(adj_mat)
#row.names(test)[modProbes,])
#yellow
#yellow_tophub <- subset(adj_mat, rownames(adj_mat) %in% yellow_top20_gene_list)
#`%notin%` <- Negate(`%in%`)#build a function not in
#yellow_NOTtophub <- subset(adj_mat, rownames(adj_mat) %notin% yellow_top20_gene_list)

#row.names(yellow_tophub)
#nrow(yellow_tophub)
#yellow_tophub_size <- cbind(row.names(yellow_tophub),rep(7,nrow(yellow_tophub)))
#yellow_NOTtophub_size <- cbind(row.names(yellow_NOTtophub),rep(2.5,nrow(yellow_NOTtophub)))
#yellow_size <-rbind(yellow_tophub_size,yellow_NOTtophub_size)
#yellow_tophub_size<-setNames(row.names(yellow_tophub),rep(7,nrow(yellow_tophub)))

#adj_mat
network <- graph.adjacency(adj_mat,mode = "undirected",weighted = TRUE)
head(network)
E(network)
V(network)
deg <- degree(network, mode="all")
V(network)$vertex_degree <-  deg

network <- simplify(network)  # removes self-loops
#results <- blockwiseModules(data, power=6, TOMType="unsigned", networkType="unsigned")
#V(network)$color <- results$colors
V(network)$color <- "turquoise"
scale_factor <- 0.1
par(mar=c(0,0,0,0))
# remove unconnected nodes
network <- delete.vertices(network, degree(network)==0)
pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/network_all_abiotic_VST#_turtoise.pdf",width = 25,height = 25)
#plot(network, vertex.color="yellow", vertex.label="",vertex.size=2,layout=layout_in_circle(network), edge.arrow.size = 0.2)
#plot(network,pch=19,arrow.mode=0,mark.border=NA, edge.color = "#d1d1e0", vertex.label="",vertex.size=2.5, layout=layout.fruchterman.reingold(network), edge.arrow.size = 0.2)
plot(network,pch=19,arrow.mode=0,mark.border=NA, edge.color = "#d1d1e0", vertex.label="",vertex.size=V(network)$vertex_degree*scale_factor, layout=layout.fruchterman.reingold(network), edge.arrow.size = 0.2)
dev.off()

#cyan
gs<-colnames(datExpr)
cols <- bwModuleColors
names(cols)<-gs #assign each gene into different module

kme <- signedKME(datExpr = datExpr, datME = MEs, outputColumnName = NULL)
head(kme)
kmes<-sapply(1:length(gs), function(x) abs(kme[gs[x],colnames(kme) == cols[[x]]]))
kmes<-data.frame(genes = gs, cols = cols, kme = kmes)
head(kmes)
cyan_kmes <- kmes[kmes$cols=="cyan",]
head(cyan_kmes)
cyan_kmes <-cyan_kmes[order(-cyan_kmes$kme),]
cyan_top20_gene_list <- cyan_kmes[1:20,]$genes

gs<-kmes$genes[kmes$cols=="cyan"]
#cols<-cols[kmes$kme>=kME.threshold]
datExpr_cyan = datExpr[,gs]
head(datExpr_cyan[,1:8])
ncol(datExpr_cyan)
cyan_gene_list <- data.frame(colnames(datExpr_cyan))
#}#yellow module we have 3587 genes!!! 3951 genes for turquoise
###write into the file:
write.csv(cyan_gene_list, file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/abiotic_cyan_list_VST_distachyon.csv")

adj_mat<-adjacency(datExpr_cyan,power=9)
###
head(adj_mat[,1:8])
######### Get the top edges of the blue module:
top.n.edges = 3000
min.edge = 2
#head(adj_mat)
message("number of genes present = ", nrow(adj_mat))#number of genes present = 3587 Distachyon: 248
if(!is.na(top.n.edges) & nrow(adj_mat) >= 100){
  adjacency.threshold = sort(as.numeric(adj_mat), decreasing=T)[(top.n.edges*2)+nrow(adj_mat)]
}
if(nrow(adj_mat) < 100)
  adjacency.threshold = 0.01

message("adjacency.threshold = ", adjacency.threshold)#adjacency.threshold = 0.690966673259787#distachyon:0.0952776666923627
adj_mat[adj_mat > adjacency.threshold] <- 1
adj_mat[adj_mat < adjacency.threshold] <- 0
diag(adj_mat) <- 0
rs<-rowSums(adj_mat)
if(verbose) cat("removing unconnected nodes\n")
adj_mat<-adj_mat[rs>min.edge,rs>min.edge]
message("number of genes present = ", nrow(adj_mat))#number of genes present = 158
#distachyon: 351
#adj_mat
network <- graph.adjacency(adj_mat,mode = "undirected",weighted = TRUE)
head(network)
E(network)
V(network)
deg <- degree(network, mode="all")
V(network)$vertex_degree <-  deg

network <- simplify(network)  # removes self-loops
#results <- blockwiseModules(data, power=6, TOMType="unsigned", networkType="unsigned")
#V(network)$color <- results$colors
V(network)$color <- "cyan"
scale_factor <- 0.1
par(mar=c(0,0,0,0))
# remove unconnected nodes
network <- delete.vertices(network, degree(network)==0)
pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/network_all_abiotic_VST#_cyan.pdf",width = 25,height = 25)
#plot(network, vertex.color="yellow", vertex.label="",vertex.size=2,layout=layout_in_circle(network), edge.arrow.size = 0.2)
#plot(network,pch=19,arrow.mode=0,mark.border=NA, edge.color = "#d1d1e0", vertex.label="",vertex.size=2.5, layout=layout.fruchterman.reingold(network), edge.arrow.size = 0.2)
plot(network,pch=19,arrow.mode=0,mark.border=NA, edge.color = "#d1d1e0", vertex.label="",vertex.size=V(network)$vertex_degree*scale_factor, layout=layout.fruchterman.reingold(network), edge.arrow.size = 0.2)
dev.off()
