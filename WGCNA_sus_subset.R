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
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

type = "unsigned"
corType = "pearson"
#corFnc = ifelse(corType=="pearson", cor, bicor)
maxPOutliers = ifelse(corType=="pearson",1,0.05)
robustY = ifelse(corType=="pearson",T,F)

tpm <- read.delim("/global/projectb/scratch/llei2019/CSP_Kranthi/normalized_VST_final_susceptible_count_formal.txt",sep="\t",header=T)
metaData <- read.delim("/global/projectb/scratch/llei2019/CSP_Kranthi/Meta_data_susceptible.txt",header = T)
head(metaData)
sub_meta <- metaData[(metaData$experiment == "BSMV_MMMV_WSMV"),]
head(sub_meta)
nrow(sub_meta)
col.num <- which(sub_meta$code  %in%  colnames(tpm))
tpm <- tpm[, col.num]
ncol(tpm)
header <- paste(sub_meta$tissue,sub_meta$treatment, sub_meta$timepoint, sub_meta$replicate,sep='_')
colnames(tpm) <- header
head(tpm)
femData <- tpm
femData <-data.frame(femData)
### select genes with MAD >85%
dim(femData) #[1] 25189   74
head(femData)
m.mad <- apply(femData,1,mad)#median absolute deviation (MAD)
head(m.mad)
#femData <- femData[which(m.mad > 0),]
#This is to take the 75% genes
#femData <- femData[which(m.mad > 
#above 85%
femData <-  femData[which(m.mad > 
                           max(quantile(m.mad, probs=seq(0, 1, 0.05))[6],0.01)),]
dim(femData)#[1] 18891    74
###
datExpr0 <- as.data.frame(t(femData))
head(datExpr0[,1:8])
dim(datExpr0) #[1]   74 18891
rownames(datExpr0)
colnames(datExpr0)
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
nGenes = ncol(datExpr0)#18891
nSamples = nrow(datExpr0)#74
head(datExpr0[,1:8])
sampleTree = hclust(dist(datExpr0), method = "average")
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small. sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
pdf(file = "/global/projectb/scratch/llei2019/CSP_Kranthi/VST_sus_Virus_sample_tree_mad75.pdf",width = 25,height = 20)
plot(sampleTree,main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)
abline(h = 250, col = "red");
dev.off()
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 250, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust!=0)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)#18891
nSamples = nrow(datExpr)#74

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
#abiotic_shoot <- meta_data
##head(abiotic_shoot)
#head(ir.pca)
head(sub_meta)
sub_meta$group <- paste(sub_meta$tissue,sub_meta$treatment,sub_meta$timepoint, sep='_')

#abiotic_shoot$group <- paste(abiotic_shoot$treatment, abiotic_shoot$timepoint, sep='_')

library(ggfortify)
pdf(file = "/global/projectb/scratch/llei2019/CSP_Kranthi/VST_sus_virus_pca_mad75.pdf",width = 10,height = 10)
autoplot(ir.pca,data = sub_meta, colour = 'group')
autoplot(ir.pca)
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
pdf(file = "/global/projectb/scratch/llei2019/CSP_Kranthi/power_soft_thresh_VST_sus_virus_sample_mad75.pdf",width = 10,height = 10)
# Scale-free topology fit index as a function of the soft-thresholding power plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
plot(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n", main = "Scale independence")

text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red")
abline(h=0.9,col="red")
# this line corresponds to using an R^2 cut-off of h abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power 
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
#abline(h=0.85,col="red")
dev.off()

#power = 12 is the best according to the plot, I set the mergeCutHeight = 0.3
bwnet = blockwiseModules(datExpr,maxBlockSize = nGenes,
                         power = 12, TOMType = type, minModuleSize = 30,
                         reassignThreshold = 0, mergeCutHeight = 0.35, numericLabels = TRUE,pamRespectsDendro = FALSE,
                         saveTOMs = TRUE,corType = corType,
                         maxPOutliers=maxPOutliers, loadTOMs=TRUE,
                         saveTOMFileBase = "allAin_TOM-blockwise", verbose = 3)
head(bwnet)
bwLabels = bwnet$colors
table(bwnet$colors)

###mergeCutHeight = 0.35
#0    1    2    3    4    5    6    7 
#129 5704 5274 4306 1721  841  802  114 
bwModuleColors = labels2colors(bwLabels)
#try to merge:

#merge = mergeCloseModules(datExpr, bwLabels, cutHeight = 0.3, verbose = 3)
#?matchLabels
#moduleLabels = merge$colors
#cbind(moduleLabels, bwnet$colors)
#table(merge$colors)

# open a graphics window
sizeGrWindow(6,6)
pdf(file = "/global/projectb/scratch/llei2019/CSP_Kranthi/Dendro_VST_sus_virus_sample_mad75_merge0.35.pdf",width = 25,height = 25)
# Plot the dendrogram and the module colors underneath for block 1
plotDendroAndColors(bwnet$dendrograms[[1]], bwModuleColors[bwnet$blockGenes[[1]]],
                    "Module colors", main = "Gene dendrogram and module colors in block 1",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

#head(bwnet)
MEs = bwnet$MEs
MEs_col = MEs
colnames(MEs_col) = paste0("ME", labels2colors(as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)
pdf(file = "/global/projectb/scratch/llei2019/CSP_Kranthi/eigennetwork_VST_sus_virus_sample_mads75#_0.35.pdf",width = 15,height = 13)
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                      xLabelsAngle = 90)
dev.off()

head(datExpr)
write.table(x = rownames(datExpr),
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/sample_id_virus.txt",
            quote = FALSE,
            sep = "\t",
            eol = "\n",
            col.names = TRUE,
            row.names = TRUE)

#associate trait
#trait <- "/global/projectb/scratch/llei2019/CSP_Kranthi/pheno_sus_final.txt"
trait <- "/global/projectb/scratch/llei2019/CSP_Kranthi/virus_pheno.txt"

if(trait != "") {
  traitData <- read.table(file=trait, sep='\t', header=T, row.names=1,
                          check.names=FALSE, comment='',quote="")
  head(traitData)
  sampleName = rownames(datExpr)
  traitData = traitData[match(sampleName, rownames(traitData)), ]
}
nrow(MEs_col)
nrow(traitData)
robustY = ifelse(corType=="pearson",T,F)
if (corType=="pearson") {
  modTraitCor = cor(MEs_col, traitData, use = "p")
  modTraitP = corPvalueStudent(modTraitCor, nSamples)
} else {
  modTraitCorP = bicorAndPvalue(MEs_col, traitData, robustY=robustY)
  modTraitCor = modTraitCorP$bicor
  modTraitP   = modTraitCorP$p
}

textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)

sizeGrWindow(13, 9)
par(mar = c(5, 30, 5, 5),mfrow = c(1,1))
cex1 = 1
pdf("/global/projectb/scratch/llei2019/CSP_Kranthi/pheno-module_virus_mads75_encoding1_0.35.pdf",width = 15,height = 15)
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(traitData), 
               yLabels = colnames(MEs_col), 
               cex.lab = 1, 
               ySymbols = colnames(MEs_col), colorLabels = FALSE, 
               colors = blueWhiteRed(50), 
               textMatrix = textMatrix, setStdMargins = FALSE, 
               cex.text = 0.8, zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

if (corType=="pearson") {
  geneModuleMembership = as.data.frame(cor(datExpr, MEs_col, use = "p"))
  MMPvalue = as.data.frame(corPvalueStudent(
    as.matrix(geneModuleMembership), nSamples))
} else {
  geneModuleMembershipA = bicorAndPvalue(datExpr, MEs_col, robustY=robustY)
  geneModuleMembership = geneModuleMembershipA$bicor
  MMPvalue   = geneModuleMembershipA$p
}
head(geneModuleMembership)

if (corType=="pearson") {
  geneTraitCor = as.data.frame(cor(datExpr, traitData, use = "p"))
  geneTraitP = as.data.frame(corPvalueStudent(
    as.matrix(geneTraitCor), nSamples))
} else {
  geneTraitCorA = bicorAndPvalue(datExpr, traitData, robustY=robustY)
  geneTraitCor = as.data.frame(geneTraitCorA$bicor)
  geneTraitP   = as.data.frame(geneTraitCorA$p)
}

###Read correlation cyan for heat positive regulation
module = "brown"
pheno = "leaves_BSMV"
modNames = substring(colnames(MEs_col), 3)
module_column = match(module, modNames)
pheno_column = match(pheno,colnames(traitData))
moduleGenes = bwModuleColors == module
sizeGrWindow(7, 7)
par(mfrow = c(1,1))
#??geneTraitSignificance
pdf("/global/projectb/scratch/llei2019/CSP_Kranthi/MEbrown_leave_BSMV.pdf",width = 15,height = 15)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                   abs(geneTraitCor[moduleGenes, pheno_column]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", pheno),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()
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

#extract top 50 hub genes for each module!!!
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
write.csv(hubgenes, file = "/global/projectb/scratch/llei2019/CSP_Kranthi/allhubgenes_list_brachy_sus_virus_0.35.csv")

###Plot the network for all of the modules
top.n.edges = 20000#3000#5000#15000#25000
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
#25000，adjacency.threshold = 0.80996000130966
#20000,adjacency.threshold = 0.821448186597572
adj_mat[adj_mat > adjacency.threshold] <- 1
adj_mat[adj_mat < adjacency.threshold] <- 0
diag(adj_mat) <- 0
rs<-rowSums(adj_mat)
if(verbose) cat("removing unconnected nodes\n")
adj_mat<-adj_mat[rs>min.edge,rs>min.edge]
message("number of genes present = ", nrow(adj_mat))#1405 #1177
test <- data.frame(adj_mat)
nrow(test)#1177
#row.names(test)[modProbes,]
#brown
probes = colnames(datExpr)
brown_modules = c("brown")
brown_inModule = is.finite(match(bwModuleColors, brown_modules))
brown_modProbes = probes[brown_inModule]
brown_list <- subset(adj_mat, rownames(adj_mat) %in% brown_modProbes)
#green
probes = colnames(datExpr)
green_modules = c("green")
green_inModule = is.finite(match(bwModuleColors, green_modules))
green_modProbes = probes[green_inModule]
green_list <- subset(adj_mat, rownames(adj_mat) %in% green_modProbes)
#black
probes = colnames(datExpr)
black_modules = c("black")
black_inModule = is.finite(match(bwModuleColors, black_modules))
black_modProbes = probes[black_inModule]
black_list <- subset(adj_mat, rownames(adj_mat) %in% black_modProbes)#0
#turquoise
probes = colnames(datExpr)
turquoise_modules = c("turquoise")
turquoise_inModule = is.finite(match(bwModuleColors, turquoise_modules))
turquoise_modProbes = probes[turquoise_inModule]
turquoise_list <- subset(adj_mat, rownames(adj_mat) %in% turquoise_modProbes)
#yellow
probes = colnames(datExpr)
yellow_modules = c("yellow")
yellow_inModule = is.finite(match(bwModuleColors, yellow_modules))
yellow_modProbes = probes[yellow_inModule]
yellow_list <- subset(adj_mat, rownames(adj_mat) %in% yellow_modProbes)
#blue
probes = colnames(datExpr)
blue_modules = c("blue")
blue_inModule = is.finite(match(bwModuleColors, blue_modules))
blue_modProbes = probes[blue_inModule]
blue_list <- subset(adj_mat, rownames(adj_mat) %in% blue_modProbes)
#red
probes = colnames(datExpr)
red_modules = c("red")
red_inModule = is.finite(match(bwModuleColors, red_modules))
red_modProbes = probes[red_inModule]
red_list <- subset(adj_mat, rownames(adj_mat) %in% red_modProbes)

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
V(network)[rownames(red_list)]$color <- "red"
V(network)[rownames(blue_list)]$color <- "blue"
V(network)[rownames(green_list)]$color <- "green"
V(network)[rownames(brown_list)]$color <- "brown"
V(network)[rownames(yellow_list)]$color <- "yellow"
V(network)[rownames(black_list)]$color <- "black"
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
minC <- rep(-Inf, vcount(network))
minC <- rep(-Inf, vcount(network))
maxC <- rep(Inf, vcount(network))
minC[1] <- maxC[1] <- 0
#co <- layout_with_fr(network, minx=minC, maxx=maxC,
#                     miny=minC, maxy=maxC)
co <- layout.fruchterman.reingold(network, minx=minC, maxx=maxC,
                                  miny=minC, maxy=maxC)
#layout.fruchterman.reingold
#co[1,]
pdf(file = "/global/projectb/scratch/llei2019/CSP_Kranthi/network_sus_virus_brachy_dges20000.pdf",width = 25,height = 25)
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
#5704
turquoise_kmes <-turquoise_kmes[order(-turquoise_kmes$kme),]
write.csv(turquoise_kmes, file = "/global/projectb/scratch/llei2019/CSP_Kranthi/sus_brachy_virus_turquoise_kmes.csv")

###red
red_kmes <- kmes[kmes$cols=="red",]
head(red_kmes)
nrow(red_kmes)
#802
red_kmes <-red_kmes[order(-red_kmes$kme),]
write.csv(red_kmes, file = "/global/projectb/scratch/llei2019/CSP_Kranthi/sus_brachy_virus_red_kmes.csv")

###Blue
blue_kmes <- kmes[kmes$cols=="blue",]
head(blue_kmes)
nrow(blue_kmes)
#5274
blue_kmes <-blue_kmes[order(-blue_kmes$kme),]
write.csv(blue_kmes, file = "/global/projectb/scratch/llei2019/CSP_Kranthi/sus_brachy_virus_blue_kmes.csv")

###green
green_kmes <- kmes[kmes$cols=="green",]
head(green_kmes)
nrow(green_kmes)
#841
green_kmes <-green_kmes[order(-green_kmes$kme),]
write.csv(green_kmes, file = "/global/projectb/scratch/llei2019/CSP_Kranthi/sus_brachy_virus_green_kmes.csv")

###brown
brown_kmes <- kmes[kmes$cols=="brown",]
head(brown_kmes)
nrow(brown_kmes)
#4306
brown_kmes <-brown_kmes[order(-brown_kmes$kme),]
write.csv(brown_kmes, file = "/global/projectb/scratch/llei2019/CSP_Kranthi/sus_brachy_virus_brown_kmes.csv")

###yellow
yellow_kmes <- kmes[kmes$cols=="yellow",]
head(yellow_kmes)
nrow(yellow_kmes)
#1721
yellow_kmes <-yellow_kmes[order(-yellow_kmes$kme),]
write.csv(yellow_kmes, file = "/global/projectb/scratch/llei2019/CSP_Kranthi/sus_brachy_virus_yellow_kmes.csv")

###black
black_kmes <- kmes[kmes$cols=="black",]
head(black_kmes)
nrow(black_kmes)
#114
black_kmes <-black_kmes[order(-black_kmes$kme),]
write.csv(black_kmes, file = "/global/projectb/scratch/llei2019/CSP_Kranthi/sus_brachy_virus_black_kmes.csv")

###grey
grey_kmes <- kmes[kmes$cols=="grey",]
head(grey_kmes)
nrow(grey_kmes)
#129
grey_kmes <-grey_kmes[order(-grey_kmes$kme),]
write.csv(grey_kmes, file = "/global/projectb/scratch/llei2019/CSP_Kranthi/sus_brachy_virus_grey_kmes.csv")

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
