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
header <- paste(metaData$accessions,metaData$experiment,metaData$tissue,metaData$treatment, metaData$timepoint, metaData$replicate,sep='_')
metaData$group <- paste(metaData$accessions,metaData$experiment,metaData$tissue,metaData$treatment, sep='_')
datTraits <- data.frame(header,metaData$group)
str(datTraits)
colnames(datTraits) <- c("ID","group")

#tpm <- read.delim("/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/VST_normalized_count_abiotic.txt",sep="\t",header=T)
#"/global/cscratch1/sd/llei2019/SeanGordon_Distachyon_gene_counts/counts.txt"
#use log(tpm+1) to do a little bit conversion
head(tpm)
str(tpm)
colnames(tpm) <- header
#tpm[1,]
#hist(tpm$Shoot_CK_1h.s1)
#femData <- sapply(tpm,function(x) log(x+1))
femData <- tpm
#column_to_rownames(tpm,"X")
#row.names(femData) <- tpm$GID
femData <-data.frame(femData)

### select genes with MAD >85%
dim(femData) #[1] 25189   125
#25189/3 = 8396.333
head(femData)
m.mad <- apply(femData,1,mad)#median absolute deviation (MAD)
head(m.mad)
#femData <- femData[which(m.mad > 0),]
#This is to take the 75% genes
#femData <- femData[which(m.mad > 
#colors()[1:17]   
quantile(m.mad, probs=seq(0, 1, 0.05))
#0%        5%       10%       15%       20%       25%       30%       35%       40% 
#0.0000000 0.2512806 0.3235558 0.3746929 0.4209740 0.4653237 0.5119409 0.5624571 0.6160490 
#45%       50%       55%       60%       65%       70%       75%       80%       85% 
#0.6743627 0.7392820 0.8127742 0.8983165 0.9943977 1.0983623 1.2126194 1.3699733 1.5670278 
#90%       95%      100% 
#1.8540406 2.3308159 6.6080134

#above 75%
#femData <- femData[which(m.mad > 
#                           max(quantile(m.mad, probs=seq(0, 1, 0.05))[6],0.01)),]
femData <- femData[which(m.mad > 
                           max(quantile(m.mad, probs=seq(0, 1, 0.05))[16],0.01)),]

###
datExpr0 <- as.data.frame(t(femData))
head(datExpr0[,1:8])
dim(datExpr0) #[1]   125 18891 #[1]  125 6297
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
nGenes = ncol(datExpr0)#75%[1] 18891 #[1] top 25% 6297
nSamples = nrow(datExpr0)#125
head(datExpr0[,1:8])
sampleTree = hclust(dist(datExpr0), method = "average")
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small. sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
#par(cex = 0.6);
par(mar = c(0,5,2,0))
pdf(file = "/global/projectb/scratch/llei2019/CSP_Kranthi/VST_sus_all_sample_tree_mad25_trait.pdf",width = 25,height = 20)
plot(sampleTree,main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 2,cex.axis = 1.5, cex.main = 2)
abline(h = 250, col = "red");
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
clust = cutreeStatic(sampleTree, cutHeight = 250, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust!=0)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)#18891
nSamples = nrow(datExpr)#125

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
ggsave("/global/projectb/scratch/llei2019/CSP_Kranthi/VST_sus_all_sample_pca_mad25.pdf", width = 30, height = 25)
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
pdf(file = "/global/projectb/scratch/llei2019/CSP_Kranthi/power_soft_thresh_VST_sus_all_sample_mad25.pdf",width = 10,height = 10)
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

#power = 8 is the best according to the plot, I set the mergeCutHeight = 0.3
bwnet = blockwiseModules(datExpr,maxBlockSize = nGenes,
                         power = 8, TOMType = type, minModuleSize = 30,
                         reassignThreshold = 0, mergeCutHeight = 0.25, numericLabels = TRUE,pamRespectsDendro = FALSE,
                         saveTOMs = TRUE,corType = corType,
                         maxPOutliers=maxPOutliers, loadTOMs=TRUE,
                         saveTOMFileBase = "BrachySus_TOM-blockwise", verbose = 3)
head(bwnet)
bwLabels = bwnet$colors
table(bwnet$colors)
#mergeCutHeight = 0.3, softpower=6
#0    1    2    3    4    5    6    7    8    9   10   11 
#149 8115 3480 3143 2000  565  437  381  375  143   65   38 

#mergeCutHeight = 0.3
#0    1    2    3    4    5    6    7    8    9   10   11   12   13 
#100 6244 4079 3635 2316  587  444  354  347  265  178  166  124   52

###mergeCutHeight = 0.35
#0    1    2    3    4    5    6    7    8    9   10   11 
#100 6244 4203 3801 2316  587  444  354  347  265  178   52 
bwModuleColors = labels2colors(bwLabels)
table(bwModuleColors)
#MDA top25
#bwModuleColors
#blue     brown     green      grey       red turquoise 
#1222       845       134        24        62      3175 
#yellow 
#835 
#bwModuleColors
#black        blue       brown       green greenyellow        grey 
#381        3480        3143         565          38         149 
#magenta        pink      purple         red   turquoise      yellow 
#143         375          65         437        8115        2000 

#try to merge:

#merge = mergeCloseModules(datExpr, bwLabels, cutHeight = 0.3, verbose = 3)
#?matchLabels
#moduleLabels = merge$colors
#cbind(moduleLabels, bwnet$colors)
#table(merge$colors)

# open a graphics window
sizeGrWindow(6,6)
pdf(file = "/global/projectb/scratch/llei2019/CSP_Kranthi/Dendro_VST_sus_all_sample_mad25_merge0.25_p6.pdf",width = 25,height = 25)
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
pdf(file = "/global/projectb/scratch/llei2019/CSP_Kranthi/eigennetwork_VST_sus_all_sample_mads25_0.25_p6.pdf",width = 15,height = 13)
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                      xLabelsAngle = 90)
dev.off()

head(MEs_col)


#trait module association
#In R:
#Plot the relationships between modules and traits
str(datTraits)
design = model.matrix(~0+ as.factor(datTraits$group))
colnames(design) = levels(as.factor(datTraits$group))
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
pdf("/global/projectb/scratch/llei2019/CSP_Kranthi/pheno-module_combined_abiotic_VST_sus_all_sample_mads25_encoding1_0.25p6.pdf",width = 15,height = 15)
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

#if (corType=="pearson") {
#  geneModuleMembership = as.data.frame(cor(datExpr, MEs_col, use = "p"))
#  MMPvalue = as.data.frame(corPvalueStudent(
#    as.matrix(geneModuleMembership), nSamples))
#} else {
#  geneModuleMembershipA = bicorAndPvalue(datExpr, MEs_col, robustY=robustY)
#  geneModuleMembership = geneModuleMembershipA$bicor
#  MMPvalue   = geneModuleMembershipA$p
#}
#head(geneModuleMembership)

#if (corType=="pearson") {
#  geneTraitCor = as.data.frame(cor(datExpr, traitData, use = "p"))
#  geneTraitP = as.data.frame(corPvalueStudent(
#    as.matrix(geneTraitCor), nSamples))
#} else {
#  geneTraitCorA = bicorAndPvalue(datExpr, traitData, robustY=robustY)
#  geneTraitCor = as.data.frame(geneTraitCorA$bicor)
#  geneTraitP   = as.data.frame(geneTraitCorA$p)
#}

###Read correlation cyan for heat positive regulation
module = "green"
pheno = "Bd21.3_BSMV_MMMV_WSMV_root_MMMV"
modNames = substring(colnames(MEs_col), 3)
module_column = match(module, modNames)
pheno_column = match(pheno,colnames(traitData))
moduleGenes = bwModuleColors == module
sizeGrWindow(7, 7)
par(mfrow = c(1,1))
#??geneTraitSignificance
pdf("/global/projectb/scratch/llei2019/CSP_Kranthi/MEgreen_Bd21.3_BSMV_MMMV_WSMV_root_MMMV.pdf",width = 15,height = 15)
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
#In R:
#Intramodular connectivity, module membership, and screening for intramodular hub genes
#Intramodular connectivity
connet=abs(cor(datExpr,use="p"))^6
Alldegrees1=intramodularConnectivity(connet, moduleColors)

#Relationship between gene significance and intramodular connectivity
which.module="green"
Luminal= as.data.frame(design[,3])
names(Luminal) = "Luminal"
GS1 = as.numeric(cor(datExpr, Luminal, use = "p"))
GeneSignificance=abs(GS1)
quantile(GS1,probs = seq(0, 1, 0.01))
#0%           1%           2%           3%           4%           5% 
#  -0.494433915 -0.337367511 -0.308488246 -0.290013329 -0.276168227 -0.262586210 
#6%           7%           8%           9%          10%          11% 
#  -0.251044324 -0.240797497 -0.230927874 -0.223031186 -0.214806144 -0.206891319 
#12%          13%          14%          15%          16%          17% 
#  -0.198788490 -0.191688931 -0.185279267 -0.178717028 -0.172118183 -0.166403681 
#18%          19%          20%          21%          22%          23% 
#  -0.161202796 -0.155189415 -0.149564319 -0.144665963 -0.138614236 -0.133531168 
#24%          25%          26%          27%          28%          29% 
#  -0.128732736 -0.123110621 -0.117809100 -0.111991138 -0.106055351 -0.100718833 
#30%          31%          32%          33%          34%          35% 
#  -0.095359333 -0.090528798 -0.084877997 -0.079862472 -0.074056942 -0.068386608 
#36%          37%          38%          39%          40%          41% 
#  -0.061841808 -0.056839969 -0.051181800 -0.046086231 -0.040808162 -0.035686829 
#42%          43%          44%          45%          46%          47% 
#  -0.031378358 -0.026900004 -0.021678292 -0.017022610 -0.012426036 -0.007276159 
#48%          49%          50%          51%          52%          53% 
#  -0.002650443  0.002146609  0.006868507  0.011638022  0.016942317  0.021389867 
#54%          55%          56%          57%          58%          59% 
#  0.026610182  0.031450458  0.036467715  0.042017375  0.047345644  0.052770887 
#60%          61%          62%          63%          64%          65% 
#  0.057524335  0.063255878  0.068434295  0.073599269  0.079021836  0.085073493 
#66%          67%          68%          69%          70%          71% 
#  0.091056623  0.097452508  0.102782985  0.108363166  0.114738246  0.120386871 
#72%          73%          74%          75%          76%          77% 
#  0.126031048  0.132365920  0.138467413  0.144612782  0.151038455  0.157085827 
#78%          79%          80%          81%          82%          83% 
#  0.164652354  0.172499456  0.179490196  0.186130296  0.192694300  0.199054671 
#84%          85%          86%          87%          88%          89% 
#  0.206625938  0.213828723  0.220637970  0.228990236  0.236947385  0.245188160 
#90%          91%          92%          93%          94%          95% 
#  0.253711858  0.262694770  0.274624575  0.284039131  0.296350606  0.309379690 
#96%          97%          98%          99%         100% 
#0.325594746  0.344155905  0.366851639  0.398978903  0.564329968 
#Generalizing intramodular connectivity for all genes on the array
datKME=signedKME(datExpr, MEs, outputColumnName="MM.")
#Display the first few rows of the data frame
head(datKME)
#Finding genes with high gene significance and high intramodular connectivity in specific modules
#abs(GS1)>.8 #adjust parameter based on actual situations
#abs(datKME$MM.black)>.8 #at least larger than 0.8
quantile(abs(datKME$MM.green),probs=seq(0,1,0.01))
quantile(abs(datKME$MM.green),probs=seq(0,1,0.01))
#0%           1%           2%           3%           4%           5% 
#  7.912615e-05 6.231307e-03 1.292665e-02 1.946377e-02 2.685373e-02 3.389659e-02 
#6%           7%           8%           9%          10%          11% 
#  4.116373e-02 4.845548e-02 5.505601e-02 6.187055e-02 6.797368e-02 7.542844e-02 
#12%          13%          14%          15%          16%          17% 
#  8.285638e-02 8.962828e-02 9.616547e-02 1.016459e-01 1.081674e-01 1.145858e-01 
#18%          19%          20%          21%          22%          23% 
#  1.212710e-01 1.275527e-01 1.342355e-01 1.399556e-01 1.464549e-01 1.527050e-01 
#24%          25%          26%          27%          28%          29% 
#  1.596801e-01 1.659851e-01 1.729109e-01 1.791395e-01 1.849279e-01 1.906919e-01 
#30%          31%          32%          33%          34%          35% 
#  1.956302e-01 2.018347e-01 2.083480e-01 2.142807e-01 2.203696e-01 2.265614e-01 
#36%          37%          38%          39%          40%          41% 
#  2.322735e-01 2.381245e-01 2.436093e-01 2.498846e-01 2.551255e-01 2.611000e-01 
#42%          43%          44%          45%          46%          47% 
#  2.667859e-01 2.719285e-01 2.775364e-01 2.828574e-01 2.877896e-01 2.927875e-01 
#48%          49%          50%          51%          52%          53% 
#  2.982537e-01 3.033007e-01 3.091663e-01 3.138837e-01 3.192176e-01 3.244997e-01 
#54%          55%          56%          57%          58%          59% 
#  3.298423e-01 3.354146e-01 3.407841e-01 3.460525e-01 3.521590e-01 3.585840e-01 
#60%          61%          62%          63%          64%          65% 
#  3.643866e-01 3.698403e-01 3.755638e-01 3.814261e-01 3.871539e-01 3.920749e-01 
#66%          67%          68%          69%          70%          71% 
#  3.978397e-01 4.032612e-01 4.088734e-01 4.145457e-01 4.203921e-01 4.258846e-01 
#72%          73%          74%          75%          76%          77% 
#  4.311343e-01 4.373582e-01 4.432215e-01 4.493601e-01 4.553296e-01 4.613267e-01 
#78%          79%          80%          81%          82%          83% 
#  4.676641e-01 4.738324e-01 4.807605e-01 4.878547e-01 4.945590e-01 5.018162e-01 
#84%          85%          86%          87%          88%          89% 
#  5.101032e-01 5.186887e-01 5.278821e-01 5.366233e-01 5.467722e-01 5.584461e-01 
#90%          91%          92%          93%          94%          95% 
#  5.720441e-01 5.872551e-01 6.023320e-01 6.199040e-01 6.395227e-01 6.620948e-01 
#96%          97%          98%          99%         100% 
#6.878042e-01 7.218301e-01 7.585260e-01 8.049063e-01 9.395972e-01 
#> FilterGenes= abs(GS1)>0.31 & abs(datKME$MM.green)>0.66
#> table(FilterGenes)
#FilterGenes
#FALSE  TRUE 
#18825    66 

FilterGenes= abs(GS1)>0.35 & abs(datKME$MM.green)>0.69
table(FilterGenes)
#FilterGenes
#FALSE  TRUE 
#18888     3 
hubgenes <- rownames(datKME)[FilterGenes]
hubgenes
#[1] "Bradi2g32220.v3.2" "Bradi4g30850.v3.2" "Bradi5g18020.v3.2"
#"TRICHOME BIREFRINGENCE-LIKE 25" "glutamate receptor","cis-zeatin O-glucosyltransferase"
#[1] "Bradi1g34385.v3.2 (unknown)" "Bradi1g61840.v3.2 (Mo25 family protein)" "Bradi2g32220.v3.2 (TRICHOME BIREFRINGENCE-LIKE 25)"
#[4] "Bradi3g02770.v3.2 (seven in absentia protein family domain containing protein)" "Bradi3g22707.v3.2 (unknown)" "Bradi3g45667.v3.2 (unknown)"
#[7] "Bradi4g08830.v3.2 (indole-3-glycerol phosphate synthase, chloroplast precursor)" "Bradi4g15425.v3.2 (unknown)" "Bradi4g30850.v3.2(glutamate receptor)"
#[10] "Bradi4g30950.v3.2 (hypro1)" "Bradi5g18020.v3.2 (cis-zeatin O-glucosyltransferase"
#[1] "Bradi1g34385.v3.2 (unknown)" "Bradi1g61840.v3.2 (Mo25 family protein)" "Bradi2g32220.v3.2 (TRICHOME BIREFRINGENCE-)"
#extract top 50 hub genes for each module!!!
#export the green module
#In R:
#Export the network
#Recalculate topological overlap
TOM = TOMsimilarityFromExpr(datExpr, power = 6) 
# Select module
module = "green"
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
                               edgeFile = paste("/global/projectb/scratch/llei2019/CSP_Kranthi/CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("/global/projectb/scratch/llei2019/CSP_Kranthi/CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes, 
                               nodeAttr = moduleColors[inModule])

#Screen the top genes
nTop = 15
IMConn = softConnectivity(datExpr[, modProbes])
top = (rank(-IMConn) <= nTop)
filter <- modTOM[top, top]

cyt = exportNetworkToCytoscape(filter,
                               edgeFile = paste("/global/projectb/scratch/llei2019/CSP_Kranthi/CytoscapeInput-edges-filter-top15-", paste(module, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("/global/projectb/scratch/llei2019/CSP_Kranthi/CytoscapeInput-nodes-filter-top15-", paste(module, collapse="-"), ".txt", sep=""),
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


####
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
write.csv(hubgenes, file = "/global/projectb/scratch/llei2019/CSP_Kranthi/allhubgenes_list_brachy_sus_0.35.csv")

###Plot the network for all of the modules
top.n.edges = 25000#3000#5000#15000#25000
min.edge = 2
adj_mat<-adjacency(datExpr,power=8)
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
