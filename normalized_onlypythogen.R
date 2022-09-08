library(ggplot2)
library(edgeR)
library(DESeq2)
library(tidyr)
library(ggrepel)
library(EnhancedVolcano)

counts <- read.delim("/global/projectb/scratch/llei2019/CSP_Kranthi/sus_beni/final_susceptible_benificial_count_formal.txt", row.names = 1)
head(counts)
nrow(counts)
#[1] 32439
metaData <- read.delim("/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/meta_data_pathogene_only_del2_reps.txt",header = T)
head(metaData)
nrow(metaData)
#102
metaData$header <- paste0(metaData$accessions, ".", metaData$experiment, ".", metaData$tissue, ".", metaData$treatment, ".", metaData$timepoint, ".", metaData$replicate)
head(metaData)
#extract the column I needed!!!
col.num <- which( colnames(counts) %in% metaData$code)
row.num <- which( metaData$code %in% colnames(counts) )
counts <- counts[, col.num]
mata <- metaData[row.num,]
ncol(counts)
#102
head(counts)
colnames(counts) <- metaData$header
head(counts)
metaData$group <- paste0(metaData$accessions, ".", metaData$experiment, ".", metaData$tissue, ".", metaData$treatment, ".", metaData$timepoint)
metaData$subgroup <- paste0( metaData$accessions, ".", metaData$experiment, ".",  metaData$treatment)
head(metaData)
metaData$accessions <- NULL
metaData$experiment <- NULL
metaData$treatment <- metaData$subgroup
metaData$subgroup <- NULL
head(metaData)
#Find the best cutoff for the downstream analysis
#This is to convert a dataframe into single column.
counts_long <- gather(counts,code,counts)
head(counts_long)
#nrow(counts_long)
#this is to find the reasonable cutoff for the threshold setting
quantile(counts_long$counts,probs=seq(0,1,0.01))
#0%     1%     2%     3%     4%     5%     6%     7%     8%     9%    10%    11%    12% 
#0      0      0      0      0      0      0      0      0      0      0      0      0 
#13%    14%    15%    16%    17%    18%    19%    20%    21%    22%    23%    24%    25% 
#0      0      0      0      0      0      0      1      1      1      1      2      2 
#26%    27%    28%    29%    30%    31%    32%    33%    34%    35%    36%    37%    38% 
#3      4      4      5      6      8      9     11     12     15     17     20     23 
#39%    40%    41%    42%    43%    44%    45%    46%    47%    48%    49%    50%    51% 
#26     30     35     40     45     51     58     65     73     81     91    101    112 
#52%    53%    54%    55%    56%    57%    58%    59%    60%    61%    62%    63%    64% 
#124    136    150    164    180    197    216    235    256    279    303    329    356 
#65%    66%    67%    68%    69%    70%    71%    72%    73%    74%    75%    76%    77% 
#385    416    449    483    521    560    602    648    696    748    803    863    927 
#78%    79%    80%    81%    82%    83%    84%    85%    86%    87%    88%    89%    90% 
#995   1069   1149   1237   1333   1436   1550   1676   1816   1972   2148   2347   2579 
#91%    92%    93%    94%    95%    96%    97%    98%    99%   100% 
#  2848   3170   3569   4068   4721   5621   6956   9272  14803 808320 

count.cutoff = 5 # 71% of the data were above 3
bioreplicates.cutoff = 3# 73% of the data were above 5
## RAW COUNTS ##
keep <- rowSums(counts >= count.cutoff) >= bioreplicates.cutoff
counts <- counts[keep, ]
#?rowSums
nrow(counts)
#[1] 28566
#CPM filtering
normalized.counts <- cpm(counts)
head(normalized.counts)
###below is for find the propriate threshold to do filtering based on the cpm normalization
normalized.counts.re <- data.frame(normalized.counts)
head(normalized.counts.re)
normalized.counts_long <- gather(normalized.counts.re,code,counts)
head(normalized.counts_long)
#nrow(counts_long)
#this is to find the reasonable cutoff for the threshold setting
quantile(normalized.counts_long$counts,probs=seq(0,1,0.01))
#0%           1%           2%           3%           4%           5%           6% 
#0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 
#7%           8%           9%          10%          11%          12%          13% 
#0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 2.383607e-02 2.933522e-02 4.216069e-02 
#14%          15%          16%          17%          18%          19%          20% 
#5.583918e-02 6.942965e-02 8.483927e-02 1.051748e-01 1.258272e-01 1.488777e-01 1.760252e-01 
#21%          22%          23%          24%          25%          26%          27% 
#2.085622e-01 2.394949e-01 2.835406e-01 3.254970e-01 3.792673e-01 4.365193e-01 5.024622e-01 
#28%          29%          30%          31%          32%          33%          34% 
#5.736175e-01 6.547789e-01 7.425686e-01 8.505492e-01 9.535743e-01 1.079825e+00 1.209631e+00 
#35%          36%          37%          38%          39%          40%          41% 
#1.356567e+00 1.516916e+00 1.689923e+00 1.876417e+00 2.079631e+00 2.295550e+00 2.534884e+00 
#42%          43%          44%          45%          46%          47%          48% 
#2.790258e+00 3.059986e+00 3.353844e+00 3.662883e+00 3.992502e+00 4.341844e+00 4.722612e+00 
#49%          50%          51%          52%          53%          54%          55% 
#5.117514e+00 5.544131e+00 5.994964e+00 6.479831e+00 6.986218e+00 7.517663e+00 8.081014e+00 
#56%          57%          58%          59%          60%          61%          62% 
#8.672763e+00 9.300201e+00 9.955766e+00 1.064120e+01 1.136231e+01 1.213361e+01 1.292504e+01 
#63%          64%          65%          66%          67%          68%          69% 
#1.376682e+01 1.465465e+01 1.559473e+01 1.657881e+01 1.761713e+01 1.871546e+01 1.987519e+01 
#70%          71%          72%          73%          74%          75%          76% 
#2.110242e+01 2.240178e+01 2.380513e+01 2.526587e+01 2.681957e+01 2.846919e+01 3.021651e+01 
#77%          78%          79%          80%          81%          82%          83% 
#3.208623e+01 3.410481e+01 3.623865e+01 3.855539e+01 4.104341e+01 4.373465e+01 4.667600e+01 
#84%          85%          86%          87%          88%          89%          90% 
#4.988456e+01 5.340227e+01 5.729920e+01 6.163825e+01 6.648493e+01 7.211650e+01 7.853809e+01 
#91%          92%          93%          94%          95%          96%          97% 
#8.611091e+01 9.517912e+01 1.062578e+02 1.201775e+02 1.382994e+02 1.632033e+02 2.000912e+02 
#98%          99%         100% 
#2.638431e+02 4.142158e+02 2.800492e+04 
count.cutoff = 1 #83% of the data is above 1
bioreplicates.cutoff = 3
keep <- rowSums(normalized.counts >= count.cutoff) >= bioreplicates.cutoff
#abiotic_shoot_count(keep)
normalized.counts <- normalized.counts[keep, ]
nrow(normalized.counts)
#[1] 24567
nrow(counts)
#[1] 28566
#keep the datapoint pass the threshold
counts <- counts[row.names(normalized.counts),]
nrow(counts)
#24567
grep("Bradi4g39317", rownames(counts))#check if this NBSLRR gene members got filtered from this link: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3368180/pdf/CFG2012-418208.pdf
#[1] 21631
#it suggested that the filtering seems good!
#24567/36927=0.6652856 only keep around 67% genes

## VST instead of voom
#creat a group 
design <- model.matrix(~0 + group,metaData)
dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, colData = metaData, design = design)
head(counts(dds))
#colSums(counts(dds)) %>% barplot
expr.vst <- assay(DESeq2::vst(dds))
head(expr.vst)
#str(expr.vst)
expr.vst <- as.data.frame(expr.vst)
#This file is for WGCNA
write.table(x = expr.vst,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/normalized_VST_pathogene_noUB_del2reps.txt",
            sep = "\t",
            eol = "\n",
            col.names = TRUE,
            row.names = TRUE
)

ddsTC <- DESeq(dds)
#heatmap:
library("pheatmap")
rld <- vst( ddsTC )
pcaData <- plotPCA(rld, intgroup = c(  "tissue", "treatment", "timepoint"),returnData = TRUE) # vsd and plotPCA are part of DESeq2 package, nothing with my example below. 
head(pcaData)
percentVar <- round(100 * attr(pcaData, "percentVar")) 
pdf(file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/PCA_pathogene_noUB_del2reps.pdf",width = 20,height = 15)
par(mar=c(1,1,1,1))
ggplot(pcaData, aes(x = PC1, y = PC2, color = factor(treatment), shape = factor(tissue))) + 
  geom_point(size =6, aes(fill=factor(treatment))) + 
  geom_point(size =6) + 
  scale_shape_manual(values=c(21,22,23,24)) + 
  #scale_alpha_manual(values=c("Bd21"=0, "Bd21-3"=1)) + 
  xlab(paste0("PC1: ", percentVar[1], "% variance")) + 
  ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
  ggtitle("PCA of all genes, no covariate adjusted")+
  theme(text = element_text(size=25))
dev.off()
#Then we can plot the other PCs 
# The function is the basically the same as https://github.com/mikelove/DESeq2/blob/master/R/plots.R.
# Inspired by https://www.biostars.org/p/243695/. 
# I added two variables, pp1 and pp2 to let user chose which principle to plot. pp1 and pp2 only take integer values.
#Test PC1 vs PC3
plotPCA.DESeqTransform = function(pp1=1, pp2=3, object, intgroup="condition", ntop=500, returnData=FALSE)
{
  # calculate the variance for each gene
  rv <- rowVars(assay(object))
  
  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  
  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(object)[select,]))
  
  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop=FALSE])
  
  # add the intgroup factors together to create a new grouping factor
  group <- if (length(intgroup) > 1) {
    factor(apply( intgroup.df, 1, paste, collapse=":"))
  } else {
    colData(object)[[intgroup]]
  }
  
  # assembly the data for the plot
  d <- data.frame(PC1=pca$x[,pp1], PC2=pca$x[,pp2], group=group, intgroup.df, name=colnames(object))
  
  if (returnData) {
    attr(d, "percentVar") <- percentVar[pp1:pp2]
    return(d)
  }
  
  ggplot(data=d, aes_string(x="PC1", y="PC3", color="group")) + geom_point(size=3) + 
    xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
    ylab(paste0("PC3: ",round(percentVar[2] * 100),"% variance"))
}

pcaData <- plotPCA.DESeqTransform(object=rld, intgroup = c(  "tissue", "treatment", "timepoint"),returnData = TRUE) # vsd and plotPCA are part of DESeq2 package, nothing with my example below. 
head(pcaData)
percentVar <- round(100 * attr(pcaData, "percentVar")) 
pdf(file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/PCA_pathogene_noUB_del2reps_pc1_pc3.pdf",width = 20,height = 15)
par(mar=c(1,1,1,1))
ggplot(pcaData, aes(x = PC1, y = PC2, color = factor(treatment), shape = factor(tissue))) + 
  geom_point(size =6, aes(fill=factor(treatment))) + 
  geom_point(size =6) + 
  scale_shape_manual(values=c(21,22,23,24)) + 
  #scale_alpha_manual(values=c("Bd21"=0, "Bd21-3"=1)) + 
  xlab(paste0("PC1: ", percentVar[1], "% variance")) + 
  ylab(paste0("PC3: ", percentVar[3], "% variance")) + 
  ggtitle("PCA of all genes, no covariate adjusted")+
  theme(text = element_text(size=25))
dev.off()

head(assay(rld))
library( "genefilter" )
#install.packages("gplots")
library(gplots)
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 1000 )
head(topVarGenes)
#dev.off()
pdf(file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/heatmap_1000_DGE_pathogene_noUB_del2reps.pdf",width = 15,height = 15)

par(mar=c(1,1,1,1))
#heatmap(assay(rld)[ topVarGenes, ],
#        labRow = "",
#        scale = "row")
palette <- colorRampPalette(c("red","white","blue"))(256)
reorder_rld <- assay(rld)[,c(4,5,6,1,2,3,10,11,12,7,8,9,16,17,18,13,14,15,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,67,100,101,102)]
head(reorder_rld)
heatmap.2( reorder_rld[ topVarGenes, ], scale="row", Colv=FALSE, 
           trace="none", dendrogram='row',
           col = palette,
           lhei = c(0.8,7),
           ColSideColors = c( "Bd21.fungi.Mock"="#cc561e","Bd21.fungi.PCA" ="#aa2b1d", "Bd21-3.bacteria.Mock" = "#ef8d32","Bd21-3.bacteria.XT" = "#beca5c", "Bd21-3.virus.Mock"="#75cfb8","Bd21-3.virus.BSMV"="#bbdfc8","Bd21-3.virus.MMMV"="#f0e5d8","Bd21-3.virus.WSMV"="#ffc478")[
             colData(rld)$treatment] )
dev.off()
##All seems good!
resultsNames(ddsTC)
Bd21.fungi.leaves.PCA.2.dpi <- results(ddsTC, contrast=list("groupBd21.fungi.leaves.PCA.2.dpi", "groupBd21.fungi.leaves.Mock.2.dpi"))
EnhancedVolcano(Bd21.fungi.leaves.PCA.2.dpi,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Bd21.fungi.leaves.PCA.2.dpi',
                xlim = c(-25, 25),
                ylim = c(0,50),
                FCcutoff=1,
                pCutoff = 0.05,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/Bd21.fungi.leaves.PCA.2.dpi.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)
#strict:
EnhancedVolcano(Bd21.fungi.leaves.PCA.2.dpi,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Bd21.fungi.leaves.PCA.2.dpi',
                xlim = c(-25, 25),
                ylim = c(0,50),
                FCcutoff=5,
                pCutoff = 0.00001,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/Bd21.fungi.leaves.PCA.2.dpi.strict.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)
# pca 4 dpi
Bd21.fungi.leaves.PCA.4.dpi <- results(ddsTC, contrast=list("groupBd21.fungi.leaves.PCA.4.dpi", "groupBd21.fungi.leaves.Mock.4.dpi"))
EnhancedVolcano(Bd21.fungi.leaves.PCA.4.dpi,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Bd21.fungi.leaves.PCA.4.dpi',
                xlim = c(-25, 25),
                ylim = c(0,50),
                FCcutoff=1,
                pCutoff = 0.05,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/Bd21.fungi.leaves.PCA.4.dpi.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)
#strict:
EnhancedVolcano(Bd21.fungi.leaves.PCA.4.dpi,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Bd21.fungi.leaves.PCA.4.dpi',
                xlim = c(-25, 25),
                ylim = c(0,50),
                FCcutoff=5,
                pCutoff = 0.00001,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/Bd21.fungi.leaves.PCA.4.dpi.strict.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)
#Bd21.fungi.leaves.PCA.6.dpi
Bd21.fungi.leaves.PCA.6.dpi <- results(ddsTC, contrast=list("groupBd21.fungi.leaves.PCA.6.dpi", "groupBd21.fungi.leaves.Mock.6.dpi"))
EnhancedVolcano(Bd21.fungi.leaves.PCA.6.dpi,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Bd21.fungi.leaves.PCA.6.dpi',
                xlim = c(-25, 25),
                ylim = c(0,50),
                FCcutoff=1,
                pCutoff = 0.05,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/Bd21.fungi.leaves.PCA.6.dpi.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)
#strict:
EnhancedVolcano(Bd21.fungi.leaves.PCA.6.dpi,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Bd21.fungi.leaves.PCA.6.dpi',
                xlim = c(-25, 25),
                ylim = c(0,50),
                FCcutoff=5,
                pCutoff = 0.00001,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/Bd21.fungi.leaves.PCA.6.dpi.strict.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)
#Bd21.3.leaves.bacteria.XT.2.dpi
Bd21.3.leaves.bacteria.XT.2.dpi <- results(ddsTC, contrast=list("groupBd21.3.bacteria.leaves.XT.2.dpi", "groupBd21.3.bacteria.leaves.Mock.2.dpi"))
EnhancedVolcano(Bd21.3.leaves.bacteria.XT.2.dpi,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Bd21.3.leaves.bacteria.XT.2.dpi',
                xlim = c(-25, 25),
                ylim = c(0,50),
                FCcutoff=1,
                pCutoff = 0.05,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/Bd21.3.leaves.bacteria.XT.2.dpi.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)
#strict:
EnhancedVolcano(Bd21.3.leaves.bacteria.XT.2.dpi,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Bd21.3.leaves.bacteria.XT.2.dpi',
                xlim = c(-25, 25),
                ylim = c(0,50),
                FCcutoff=5,
                pCutoff = 0.00001,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/Bd21.3.leaves.bacteria.XT.2.dpi_strict.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)

#Bd21.3.leaves.bacteria.XT.4.dpi
Bd21.3.leaves.bacteria.XT.4.dpi <- results(ddsTC, contrast=list("groupBd21.3.bacteria.leaves.XT.4.dpi", "groupBd21.3.bacteria.leaves.Mock.4.dpi"))
EnhancedVolcano(Bd21.3.leaves.bacteria.XT.4.dpi,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Bd21.3.leaves.bacteria.XT.4.dpi',
                xlim = c(-25, 25),
                ylim = c(0,50),
                FCcutoff=1,
                pCutoff = 0.05,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/Bd21.3.leaves.bacteria.XT.4.dpi.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)
#strict:
EnhancedVolcano(Bd21.3.leaves.bacteria.XT.4.dpi,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Bd21.3.leaves.bacteria.XT.4.dpi',
                xlim = c(-25, 25),
                ylim = c(0,50),
                FCcutoff=5,
                pCutoff = 0.00001,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/Bd21.3.leaves.bacteria.XT.4.dpi_strict.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)
#Bd21.3.leaves.virus.BSMV.3.dpi
Bd21.3.leaves.virus.BSMV.3.dpi <- results(ddsTC, contrast=list("groupBd21.3.virus.leaves.BSMV.3.dpi", "groupBd21.3.virus.leaves.Mock.3.dpi"))
EnhancedVolcano(Bd21.3.leaves.virus.BSMV.3.dpi,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Bd21.3.leaves.virus.BSMV.3.dpi',
                xlim = c(-25, 25),
                ylim = c(0,50),
                FCcutoff=1,
                pCutoff = 0.05,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/Bd21.3.leaves.virus.BSMV.3.dpi.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)
#Bd21.3.leaves.virus.BSMV.3.dpi.strict
EnhancedVolcano(Bd21.3.leaves.virus.BSMV.3.dpi,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Bd21.3.leaves.virus.BSMV.3.dpi',
                xlim = c(-25, 25),
                ylim = c(0,50),
                FCcutoff=5,
                pCutoff = 0.00001,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/Bd21.3.leaves.virus.BSMV.3.dpi.strict.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)
#Bd21.3.leaves.virus.BSMV.7.dpi
Bd21.3.leaves.virus.BSMV.7.dpi <- results(ddsTC, contrast=list("groupBd21.3.virus.leaves.BSMV.7.dpi", "groupBd21.3.virus.leaves.Mock.7.dpi"))
EnhancedVolcano(Bd21.3.leaves.virus.BSMV.7.dpi,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Bd21.3.leaves.virus.BSMV.7.dpi',
                xlim = c(-25, 25),
                ylim = c(0,50),
                FCcutoff=1,
                pCutoff = 0.05,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/Bd21.3.leaves.virus.BSMV.7.dpi.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)
#Bd21.3.leaves.virus.BSMV.7.dpi.strict
EnhancedVolcano(Bd21.3.leaves.virus.BSMV.7.dpi,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Bd21.3.leaves.virus.BSMV.7.dpi',
                xlim = c(-25, 25),
                ylim = c(0,50),
                FCcutoff=5,
                pCutoff = 0.00001,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/Bd21.3.leaves.virus.BSMV.7.dpi.strict.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)
#Bd21.3.leaves.virus.BSMV.14.dpi
Bd21.3.leaves.virus.BSMV.14.dpi <- results(ddsTC, contrast=list("groupBd21.3.virus.leaves.BSMV.14.dpi", "groupBd21.3.virus.leaves.Mock.14.dpi"))
EnhancedVolcano(Bd21.3.leaves.virus.BSMV.14.dpi,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Bd21.3.leaves.virus.BSMV.14.dpi',
                xlim = c(-25, 25),
                ylim = c(0,50),
                FCcutoff=1,
                pCutoff = 0.05,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/Bd21.3.leaves.virus.BSMV.14.dpi.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)
#Bd21.3.leaves.virus.BSMV.14.dpi.strict
EnhancedVolcano(Bd21.3.leaves.virus.BSMV.14.dpi,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Bd21.3.leaves.virus.BSMV.14.dpi',
                xlim = c(-25, 25),
                ylim = c(0,50),
                FCcutoff=5,
                pCutoff = 0.00001,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/Bd21.3.leaves.virus.BSMV.14.dpi.strict.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)
#Bd21.3.leaves.virus.BSMV.3.dpi
Bd21.3.leaves.virus.MMMV.3.dpi <- results(ddsTC, contrast=list("groupBd21.3.virus.leaves.MMMV.3.dpi", "groupBd21.3.virus.leaves.Mock.3.dpi"))
EnhancedVolcano(Bd21.3.leaves.virus.MMMV.3.dpi,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Bd21.3.leaves.virus.MMMV.3.dpi',
                xlim = c(-25, 25),
                ylim = c(0,50),
                FCcutoff=1,
                pCutoff = 0.05,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/Bd21.3.leaves.virus.MMMV.3.dpi.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)
#Bd21.3.leaves.virus.BSMV.3.dpi.strict
EnhancedVolcano(Bd21.3.leaves.virus.MMMV.3.dpi,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Bd21.3.leaves.virus.MMMV.3.dpi',
                xlim = c(-25, 25),
                ylim = c(0,50),
                FCcutoff=5,
                pCutoff = 0.00005,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/Bd21.3.leaves.virus.MMMV.3.dpi.strict.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)
#Bd21.3.leaves.virus.BSMV.7.dpi
Bd21.3.leaves.virus.MMMV.7.dpi <- results(ddsTC, contrast=list("groupBd21.3.virus.leaves.MMMV.7.dpi", "groupBd21.3.virus.leaves.Mock.7.dpi"))
EnhancedVolcano(Bd21.3.leaves.virus.MMMV.7.dpi,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Bd21.3.leaves.virus.MMMV.7.dpi',
                xlim = c(-25, 25),
                ylim = c(0,50),
                FCcutoff=1,
                pCutoff = 0.05,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/Bd21.3.leaves.virus.MMMV.7.dpi.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)
#Bd21.3.leaves.virus.BSMV.7.dpi.strict
EnhancedVolcano(Bd21.3.leaves.virus.MMMV.7.dpi,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Bd21.3.leaves.virus.MMMV.7.dpi',
                xlim = c(-25, 25),
                ylim = c(0,50),
                FCcutoff=5,
                pCutoff = 0.00001,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/Bd21.3.leaves.virus.MMMV.7.dpi.strict.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)
#Bd21.3.leaves.virus.BSMV.14.dpi
Bd21.3.leaves.virus.MMMV.14.dpi <- results(ddsTC, contrast=list("groupBd21.3.virus.leaves.MMMV.14.dpi", "groupBd21.3.virus.leaves.Mock.14.dpi"))
EnhancedVolcano(Bd21.3.leaves.virus.MMMV.14.dpi,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Bd21.3.leaves.virus.BSMV.14.dpi',
                xlim = c(-25, 25),
                ylim = c(0,50),
                FCcutoff=1,
                pCutoff = 0.05,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/Bd21.3.leaves.virus.MMMV.14.dpi.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)
EnhancedVolcano(Bd21.3.leaves.virus.MMMV.14.dpi,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Bd21.3.leaves.virus.BSMV.14.dpi',
                xlim = c(-25, 25),
                ylim = c(0,50),
                FCcutoff=5,
                pCutoff = 0.00001,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/Bd21.3.leaves.virus.MMMV.14.dpi.strict.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)
#Bd21.3.leaves.virus.WSMV.3.dpi
Bd21.3.leaves.virus.WSMV.3.dpi <- results(ddsTC, contrast=list("groupBd21.3.virus.leaves.WSMV.3.dpi", "groupBd21.3.virus.leaves.Mock.3.dpi"))
EnhancedVolcano(Bd21.3.leaves.virus.WSMV.3.dpi,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Bd21.3.leaves.virus.WSMV.3.dpi',
                xlim = c(-25, 25),
                ylim = c(0,50),
                FCcutoff=1,
                pCutoff = 0.05,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/Bd21.3.leaves.virus.WSMV.3.dpi.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)
#Bd21.3.leaves.virus.WSMV.3.dpi.strict
EnhancedVolcano(Bd21.3.leaves.virus.WSMV.3.dpi,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Bd21.3.leaves.virus.WSMV.3.dpi',
                xlim = c(-25, 25),
                ylim = c(0,50),
                FCcutoff=5,
                pCutoff = 0.00001,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/Bd21.3.leaves.virus.WSMV.3.dpi.strict.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)
#Bd21.3.leaves.virus.WSMV.7.dpi
Bd21.3.leaves.virus.WSMV.7.dpi <- results(ddsTC, contrast=list("groupBd21.3.virus.leaves.WSMV.7.dpi", "groupBd21.3.virus.leaves.Mock.7.dpi"))
EnhancedVolcano(Bd21.3.leaves.virus.WSMV.7.dpi,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Bd21.3.leaves.virus.WSMV.7.dpi',
                xlim = c(-25, 25),
                ylim = c(0,50),
                FCcutoff=1,
                pCutoff = 0.05,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/Bd21.3.leaves.virus.WSMV.7.dpi.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)
#Bd21.3.leaves.virus.WSMV.7.dpi.strict
EnhancedVolcano(Bd21.3.leaves.virus.WSMV.7.dpi,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Bd21.3.leaves.virus.WSMV.7.dpi',
                xlim = c(-25, 25),
                ylim = c(0,50),
                FCcutoff=5,
                pCutoff = 0.00001,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/Bd21.3.leaves.virus.WSMV.7.dpi.strict.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)
#Bd21.3.leaves.virus.BSMV.14.dpi
Bd21.3.leaves.virus.WSMV.14.dpi <- results(ddsTC, contrast=list("groupBd21.3.virus.leaves.WSMV.14.dpi", "groupBd21.3.virus.leaves.Mock.14.dpi"))
EnhancedVolcano(Bd21.3.leaves.virus.WSMV.14.dpi,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Bd21.3.leaves.virus.WSMV.14.dpi',
                xlim = c(-25, 25),
                ylim = c(0,50),
                FCcutoff=1,
                pCutoff = 0.05,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/Bd21.3.leaves.virus.WSMV.14.dpi.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)
#Bd21.3.leaves.virus.BSMV.14.dpi.strict
EnhancedVolcano(Bd21.3.leaves.virus.WSMV.14.dpi,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Bd21.3.leaves.virus.WSMV.14.dpi',
                xlim = c(-25, 25),
                ylim = c(0,50),
                FCcutoff=5,
                pCutoff = 0.00001,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/Bd21.3.leaves.virus.WSMV.14.dpi.strict.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)
###ROOT
#Bd21.3.root.virus.BSMV.3.dpi
Bd21.3.root.virus.BSMV.3.dpi <- results(ddsTC, contrast=list("groupBd21.3.virus.root.BSMV.3.dpi", "groupBd21.3.virus.root.Mock.3.dpi"))
EnhancedVolcano(Bd21.3.root.virus.BSMV.3.dpi,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Bd21.3.root.virus.BSMV.3.dpi',
                xlim = c(-25, 25),
                ylim = c(0,50),
                FCcutoff=1,
                pCutoff = 0.05,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/Bd21.3.root.virus.BSMV.3.dpi.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)
#Bd21.3.root.virus.BSMV.3.dpi.strict
EnhancedVolcano(Bd21.3.root.virus.BSMV.3.dpi,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Bd21.3.root.virus.BSMV.3.dpi',
                xlim = c(-25, 25),
                ylim = c(0,50),
                FCcutoff=5,
                pCutoff = 0.00001,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/Bd21.3.root.virus.BSMV.3.dpi..strict.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)
#Bd21.3.root.virus.BSMV.7.dpi
Bd21.3.root.virus.BSMV.7.dpi <- results(ddsTC, contrast=list("groupBd21.3.virus.root.BSMV.7.dpi", "groupBd21.3.virus.root.Mock.7.dpi"))
EnhancedVolcano(Bd21.3.root.virus.BSMV.7.dpi,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Bd21.3.root.virus.BSMV.7.dpi',
                xlim = c(-25, 25),
                ylim = c(0,50),
                FCcutoff=1,
                pCutoff = 0.05,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/Bd21.3.root.virus.BSMV.7.dpi.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)
#Bd21.3.root.virus.BSMV.7.dpi.strict
EnhancedVolcano(Bd21.3.root.virus.BSMV.7.dpi,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Bd21.3.root.virus.BSMV.7.dpi',
                xlim = c(-25, 25),
                ylim = c(0,50),
                FCcutoff=5,
                pCutoff = 0.00001,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/Bd21.3.root.virus.BSMV.7.dpi.strict.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)
#Bd21.3.root.virus.BSMV.14.dpi
Bd21.3.root.virus.BSMV.14.dpi <- results(ddsTC, contrast=list("groupBd21.3.virus.root.BSMV.14.dpi", "groupBd21.3.virus.root.Mock.14.dpi"))
EnhancedVolcano(Bd21.3.root.virus.BSMV.14.dpi,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Bd21.3.root.virus.BSMV.14.dpi',
                xlim = c(-25, 25),
                ylim = c(0,50),
                FCcutoff=1,
                pCutoff = 0.05,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/Bd21.3.root.virus.BSMV.14.dpi.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)
#Bd21.3.root.virus.BSMV.14.dpi.strict
EnhancedVolcano(Bd21.3.root.virus.BSMV.14.dpi,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Bd21.3.root.virus.BSMV.14.dpi',
                xlim = c(-25, 25),
                ylim = c(0,50),
                FCcutoff=5,
                pCutoff = 0.00001,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/Bd21.3.root.virus.BSMV.14.dpi.strict.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)
#Bd21.3.root.virus.BSMV.3.dpi
Bd21.3.root.virus.MMMV.3.dpi <- results(ddsTC, contrast=list("groupBd21.3.virus.root.MMMV.3.dpi", "groupBd21.3.virus.root.Mock.3.dpi"))
EnhancedVolcano(Bd21.3.root.virus.MMMV.3.dpi,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Bd21.3.root.virus.MMMV.3.dpi',
                xlim = c(-25, 25),
                ylim = c(0,50),
                FCcutoff=1,
                pCutoff = 0.05,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/Bd21.3.root.virus.MMMV.3.dpi.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)
#Bd21.3.root.virus.BSMV.3.dpi.strict.

EnhancedVolcano(Bd21.3.root.virus.MMMV.3.dpi,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Bd21.3.root.virus.MMMV.3.dpi',
                xlim = c(-25, 25),
                ylim = c(0,50),
                FCcutoff=5,
                pCutoff = 0.00001,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/Bd21.3.root.virus.MMMV.3.dpi.strict.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)
#Bd21.3.root.virus.BSMV.7.dpi
Bd21.3.root.virus.MMMV.7.dpi <- results(ddsTC, contrast=list("groupBd21.3.virus.root.MMMV.7.dpi", "groupBd21.3.virus.root.Mock.7.dpi"))
EnhancedVolcano(Bd21.3.root.virus.MMMV.7.dpi,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Bd21.3.root.virus.MMMV.7.dpi',
                xlim = c(-25, 25),
                ylim = c(0,50),
                FCcutoff=1,
                pCutoff = 0.05,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/Bd21.3.root.virus.MMMV.7.dpi.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)
#Bd21.3.root.virus.BSMV.7.dpi.strict
EnhancedVolcano(Bd21.3.root.virus.MMMV.7.dpi,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Bd21.3.root.virus.MMMV.7.dpi',
                xlim = c(-25, 25),
                ylim = c(0,50),
                FCcutoff=5,
                pCutoff = 0.00001,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/Bd21.3.root.virus.MMMV.7.dpi.strict.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)

#Bd21.3.root.virus.BSMV.14.dpi
Bd21.3.root.virus.MMMV.14.dpi <- results(ddsTC, contrast=list("groupBd21.3.virus.root.MMMV.14.dpi", "groupBd21.3.virus.root.Mock.14.dpi"))
EnhancedVolcano(Bd21.3.root.virus.MMMV.14.dpi,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Bd21.3.root.virus.BSMV.14.dpi',
                xlim = c(-25, 25),
                ylim = c(0,50),
                FCcutoff=1,
                pCutoff = 0.05,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/Bd21.3.root.virus.MMMV.14.dpi.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)

#Bd21.3.root.virus.BSMV.14.dpi.strict
EnhancedVolcano(Bd21.3.root.virus.MMMV.14.dpi,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Bd21.3.root.virus.BSMV.14.dpi',
                xlim = c(-25, 25),
                ylim = c(0,50),
                FCcutoff=5,
                pCutoff = 0.00001,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/Bd21.3.root.virus.MMMV.14.dpi.strict.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)
#Bd21.3.root.virus.WSMV.3.dpi
Bd21.3.root.virus.WSMV.3.dpi <- results(ddsTC, contrast=list("groupBd21.3.virus.root.WSMV.3.dpi", "groupBd21.3.virus.root.Mock.3.dpi"))
EnhancedVolcano(Bd21.3.root.virus.WSMV.3.dpi,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Bd21.3.root.virus.WSMV.3.dpi',
                xlim = c(-25, 25),
                ylim = c(0,50),
                FCcutoff=1,
                pCutoff = 0.05,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/Bd21.3.root.virus.WSMV.3.dpi.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)
#Bd21.3.root.virus.WSMV.3.dpi.strict
EnhancedVolcano(Bd21.3.root.virus.WSMV.3.dpi,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Bd21.3.root.virus.WSMV.3.dpi',
                xlim = c(-25, 25),
                ylim = c(0,50),
                FCcutoff=5,
                pCutoff = 0.00001,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/Bd21.3.root.virus.WSMV.3.dpi.strict.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)
#Bd21.3.root.virus.WSMV.7.dpi
Bd21.3.root.virus.WSMV.7.dpi <- results(ddsTC, contrast=list("groupBd21.3.virus.root.WSMV.7.dpi", "groupBd21.3.virus.root.Mock.7.dpi"))
EnhancedVolcano(Bd21.3.root.virus.WSMV.7.dpi,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Bd21.3.root.virus.WSMV.7.dpi',
                xlim = c(-25, 25),
                ylim = c(0,50),
                FCcutoff=1,
                pCutoff = 0.05,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/Bd21.3.root.virus.WSMV.7.dpi.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)
#Bd21.3.root.virus.WSMV.7.dpi.strict
EnhancedVolcano(Bd21.3.root.virus.WSMV.7.dpi,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Bd21.3.root.virus.WSMV.7.dpi',
                xlim = c(-25, 25),
                ylim = c(0,50),
                FCcutoff=5,
                pCutoff = 0.00001,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/Bd21.3.root.virus.WSMV.7.dpi.strict.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)

#Bd21.3.root.virus.BSMV.14.dpi
Bd21.3.root.virus.WSMV.14.dpi <- results(ddsTC, contrast=list("groupBd21.3.virus.root.WSMV.14.dpi", "groupBd21.3.virus.root.Mock.14.dpi"))
EnhancedVolcano(Bd21.3.root.virus.WSMV.14.dpi,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Bd21.3.root.virus.WSMV.14.dpi',
                xlim = c(-25, 25),
                ylim = c(0,50),
                FCcutoff=1,
                pCutoff = 0.05,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/Bd21.3.root.virus.WSMV.14.dpi.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)
#Bd21.3.root.virus.BSMV.14.dpi.strict
EnhancedVolcano(Bd21.3.root.virus.WSMV.14.dpi,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Bd21.3.root.virus.WSMV.14.dpi',
                xlim = c(-25, 25),
                ylim = c(0,50),
                FCcutoff=5,
                pCutoff = 0.00001,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/Bd21.3.root.virus.WSMV.14.dpi.strict.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)

#Define the DEG as foldchanges >1 and adjustP-value <0.05
#Bd21.fungi.leaves.PCA.2.dpi
Bd21.fungi.leaves.PCA.2.dpi <- Bd21.fungi.leaves.PCA.2.dpi[!is.na(Bd21.fungi.leaves.PCA.2.dpi$padj),]
sigres_Bd21.fungi.leaves.PCA.2.dpi_up <- Bd21.fungi.leaves.PCA.2.dpi [(Bd21.fungi.leaves.PCA.2.dpi$padj<0.05 & Bd21.fungi.leaves.PCA.2.dpi $log2FoldChange>= 1),]
nb_Bd21.fungi.leaves.PCA.2.dpi_up <- nrow(sigres_Bd21.fungi.leaves.PCA.2.dpi_up)
# 1157
write.table(x = sigres_Bd21.fungi.leaves.PCA.2.dpi_up,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.fungi.leaves.PCA.2.dpi_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

sigres_Bd21.fungi.leaves.PCA.2.dpi_down <- Bd21.fungi.leaves.PCA.2.dpi[(Bd21.fungi.leaves.PCA.2.dpi$padj<0.05 & Bd21.fungi.leaves.PCA.2.dpi$log2FoldChange<= -1),]
nb_Bd21.fungi.leaves.PCA.2.dpi_down <- nrow(sigres_Bd21.fungi.leaves.PCA.2.dpi_down)
#1401
write.table(x = sigres_Bd21.fungi.leaves.PCA.2.dpi_down,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.fungi.leaves.PCA.2.dpi_down.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
sigres_Bd21.fungi.leaves.PCA.2.dpi <- Bd21.fungi.leaves.PCA.2.dpi[(Bd21.fungi.leaves.PCA.2.dpi$padj<0.05 & abs(Bd21.fungi.leaves.PCA.2.dpi$log2FoldChange) >= 1 ),]
nb_sigres_Bd21.fungi.leaves.PCA.2.dpi_DEG <- nrow(sigres_Bd21.fungi.leaves.PCA.2.dpi)
#2558
write.table(x = sigres_Bd21.fungi.leaves.PCA.2.dpi,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_sigres_Bd21.fungi.leaves.PCA.2.dpi_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
#strict: logfoldchanges >5 and the P-value <0.00005
sigres_Bd21.fungi.leaves.PCA.2.dpi_strict <- Bd21.fungi.leaves.PCA.2.dpi[(Bd21.fungi.leaves.PCA.2.dpi$padj<0.00001 & abs(Bd21.fungi.leaves.PCA.2.dpi$log2FoldChange) >= 5 ),]
nrow(sigres_Bd21.fungi.leaves.PCA.2.dpi_strict)
#29
sigres_Bd21.fungi.leaves.PCA.2.dpi_up_strict <- Bd21.fungi.leaves.PCA.2.dpi[(Bd21.fungi.leaves.PCA.2.dpi$padj<0.00001 & Bd21.fungi.leaves.PCA.2.dpi$log2FoldChange >= 5 ),]
write.table(x = sigres_Bd21.fungi.leaves.PCA.2.dpi_strict,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_sigres_Bd21.fungi.leaves.PCA.2.dpi_DEG_strict.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
write.table(x = sigres_Bd21.fungi.leaves.PCA.2.dpi_up_strict,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.fungi.leaves.PCA.2.dpi_up_strict.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
#Bd21.fungi.leaves.PCA.4.dpi
Bd21.fungi.leaves.PCA.4.dpi <- Bd21.fungi.leaves.PCA.4.dpi[!is.na(Bd21.fungi.leaves.PCA.4.dpi$padj),]
sigres_Bd21.fungi.leaves.PCA.4.dpi_up <- Bd21.fungi.leaves.PCA.4.dpi [(Bd21.fungi.leaves.PCA.4.dpi $padj<0.05 & Bd21.fungi.leaves.PCA.4.dpi $log2FoldChange>= 1),]
nb_Bd21.fungi.leaves.PCA.4.dpi_up <- nrow(sigres_Bd21.fungi.leaves.PCA.4.dpi_up)
# 653
write.table(x = sigres_Bd21.fungi.leaves.PCA.4.dpi_up,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.fungi.leaves.PCA.4.dpi_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
sigres_Bd21.fungi.leaves.PCA.4.dpi_down <- Bd21.fungi.leaves.PCA.4.dpi[(Bd21.fungi.leaves.PCA.4.dpi$padj<0.05 & Bd21.fungi.leaves.PCA.4.dpi$log2FoldChange<= -1),]
nb_Bd21.fungi.leaves.PCA.4.dpi_down <- nrow(sigres_Bd21.fungi.leaves.PCA.4.dpi_down)
#1631
write.table(x = sigres_Bd21.fungi.leaves.PCA.4.dpi_down,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.fungi.leaves.PCA.4.dpi_down.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
sigres_Bd21.fungi.leaves.PCA.4.dpi <- Bd21.fungi.leaves.PCA.4.dpi[(Bd21.fungi.leaves.PCA.4.dpi$padj<0.05 & abs(Bd21.fungi.leaves.PCA.4.dpi$log2FoldChange) >= 1 ),]
nb_sigres_Bd21.fungi.leaves.PCA.4.dpi_DEG <- nrow(sigres_Bd21.fungi.leaves.PCA.4.dpi)
#2284
write.table(x = sigres_Bd21.fungi.leaves.PCA.4.dpi,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_sigres_Bd21.fungi.leaves.PCA.4.dpi_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
#strict:
sigres_Bd21.fungi.leaves.PCA.4.dpi_strict <- Bd21.fungi.leaves.PCA.4.dpi[(Bd21.fungi.leaves.PCA.4.dpi$padj<0.00001 & abs(Bd21.fungi.leaves.PCA.4.dpi$log2FoldChange) >= 5 ),]
nb_sigres_Bd21.fungi.leaves.PCA.4.dpi_DEG_strict <- nrow(sigres_Bd21.fungi.leaves.PCA.4.dpi_strict)
#112
write.table(x = sigres_Bd21.fungi.leaves.PCA.4.dpi_strict,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_sigres_Bd21.fungi.leaves.PCA.4.dpi_DEG_strict.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
sigres_Bd21.fungi.leaves.PCA.4.dpi_up_strict <- Bd21.fungi.leaves.PCA.4.dpi[(Bd21.fungi.leaves.PCA.4.dpi$padj<0.00001 & Bd21.fungi.leaves.PCA.4.dpi$log2FoldChange >= 5 ),]
write.table(x = sigres_Bd21.fungi.leaves.PCA.4.dpi_up_strict,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.fungi.leaves.PCA.4.dpi_up_strict.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

#Bd21.fungi.leaves.PCA.6.dpi
Bd21.fungi.leaves.PCA.6.dpi <- Bd21.fungi.leaves.PCA.6.dpi[!is.na(Bd21.fungi.leaves.PCA.6.dpi$padj),]
sigres_Bd21.fungi.leaves.PCA.6.dpi_up <- Bd21.fungi.leaves.PCA.6.dpi [(Bd21.fungi.leaves.PCA.6.dpi $padj<0.05 & Bd21.fungi.leaves.PCA.6.dpi $log2FoldChange>= 1),]
nb_Bd21.fungi.leaves.PCA.6.dpi_up <- nrow(sigres_Bd21.fungi.leaves.PCA.6.dpi_up)
# 547
write.table(x = sigres_Bd21.fungi.leaves.PCA.6.dpi_up,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.fungi.leaves.PCA.6.dpi_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
sigres_Bd21.fungi.leaves.PCA.6.dpi_down <- Bd21.fungi.leaves.PCA.6.dpi[(Bd21.fungi.leaves.PCA.6.dpi$padj<0.05 & Bd21.fungi.leaves.PCA.6.dpi$log2FoldChange<= -1),]
nb_Bd21.fungi.leaves.PCA.6.dpi_down <- nrow(sigres_Bd21.fungi.leaves.PCA.6.dpi_down)
#1154
write.table(x = sigres_Bd21.fungi.leaves.PCA.6.dpi_down,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.fungi.leaves.PCA.6.dpi_down.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
#DEG:
sigres_Bd21.fungi.leaves.PCA.6.dpi <- Bd21.fungi.leaves.PCA.6.dpi[(Bd21.fungi.leaves.PCA.6.dpi$padj<0.05 & abs(Bd21.fungi.leaves.PCA.6.dpi$log2FoldChange) >= 1 ),]
nb_sigres_Bd21.fungi.leaves.PCA.6.dpi_DEG <- nrow(sigres_Bd21.fungi.leaves.PCA.6.dpi)
#1701
write.table(x = sigres_Bd21.fungi.leaves.PCA.6.dpi,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_sigres_Bd21.fungi.leaves.PCA.6.dpi_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
sigres_Bd21.fungi.leaves.PCA.6.dpi_strict <- Bd21.fungi.leaves.PCA.6.dpi[(Bd21.fungi.leaves.PCA.6.dpi$padj<0.00001 & abs(Bd21.fungi.leaves.PCA.6.dpi$log2FoldChange) >= 5 ),]
nb_sigres_Bd21.fungi.leaves.PCA.6.dpi_DEG_strict <- nrow(sigres_Bd21.fungi.leaves.PCA.6.dpi_strict)
#20
write.table(x = sigres_Bd21.fungi.leaves.PCA.6.dpi_strict,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_sigres_Bd21.fungi.leaves.PCA.6.dpi_DEG_strict.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
sigres_Bd21.fungi.leaves.PCA.all_strict <- union(union(rownames(sigres_Bd21.fungi.leaves.PCA.2.dpi_strict),rownames(sigres_Bd21.fungi.leaves.PCA.4.dpi_strict)),rownames(sigres_Bd21.fungi.leaves.PCA.6.dpi_strict))
write.table(x = sigres_Bd21.fungi.leaves.PCA.all_strict,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_sigres_Bd21.fungi.leaves.PCA.all_strict_gene.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
length(sigres_Bd21.fungi.leaves.PCA.all_strict )
#142
intersect(intersect(rownames(sigres_Bd21.fungi.leaves.PCA.2.dpi_strict),rownames(sigres_Bd21.fungi.leaves.PCA.4.dpi_strict)),rownames(sigres_Bd21.fungi.leaves.PCA.6.dpi_strict))
#[1] "Bradi2g09680.v3.2" "Bradi2g20840.v3.2" "Bradi5g17850.v3.2" all of those three are " pectinesterase"

#Bd21.3.leaves.bacteria.XT.2.dpi
Bd21.3.leaves.bacteria.XT.2.dpi <- Bd21.3.leaves.bacteria.XT.2.dpi[!is.na(Bd21.3.leaves.bacteria.XT.2.dpi$padj),]
sigres_Bd21.3.leaves.bacteria.XT.2.dpi_up <- Bd21.3.leaves.bacteria.XT.2.dpi [(Bd21.3.leaves.bacteria.XT.2.dpi$padj<0.05 & Bd21.3.leaves.bacteria.XT.2.dpi $log2FoldChange>= 1),]
nb_Bd21.3.leaves.bacteria.XT.2.dpi_up <- nrow(sigres_Bd21.3.leaves.bacteria.XT.2.dpi_up)
# 221
write.table(x = sigres_Bd21.3.leaves.bacteria.XT.2.dpi_up,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.leaves.bacteria.XT.2.dpi_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
sigres_Bd21.3.leaves.bacteria.XT.2.dpi_down <- Bd21.3.leaves.bacteria.XT.2.dpi[(Bd21.3.leaves.bacteria.XT.2.dpi$padj<0.05 & Bd21.3.leaves.bacteria.XT.2.dpi$log2FoldChange<= -1),]
nb_Bd21.3.leaves.bacteria.XT.2.dpi_down <- nrow(sigres_Bd21.3.leaves.bacteria.XT.2.dpi_down)
#347
write.table(x = sigres_Bd21.3.leaves.bacteria.XT.2.dpi_down,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.leaves.bacteria.XT.2.dpi_down.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
#DEGs:
sigres_Bd21.3.leaves.bacteria.XT.2.dpi <- Bd21.3.leaves.bacteria.XT.2.dpi[(Bd21.3.leaves.bacteria.XT.2.dpi$padj<0.05 & abs(Bd21.3.leaves.bacteria.XT.2.dpi$log2FoldChange) >= 1 ),]
nb_sigres_Bd21.3.leaves.bacteria.XT.2.dpi_DEG <- nrow(sigres_Bd21.3.leaves.bacteria.XT.2.dpi)
#568
write.table(x = sigres_Bd21.3.leaves.bacteria.XT.2.dpi,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.leaves.bacteria.XT.2.dpi_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
##strict
sigres_Bd21.3.leaves.bacteria.XT.2.dpi_strict <- Bd21.3.leaves.bacteria.XT.2.dpi[(Bd21.3.leaves.bacteria.XT.2.dpi$padj<0.00001 & abs(Bd21.3.leaves.bacteria.XT.2.dpi$log2FoldChange) >= 5 ),]
nb_sigres_Bd21.3.leaves.bacteria.XT.2.dpi_DEG_strict <- nrow(sigres_Bd21.3.leaves.bacteria.XT.2.dpi_strict)
#27
sigres_Bd21.3.leaves.bacteria.XT.2.dpi_up_strict <- Bd21.3.leaves.bacteria.XT.2.dpi[(Bd21.3.leaves.bacteria.XT.2.dpi$padj<0.00001 & Bd21.3.leaves.bacteria.XT.2.dpi$log2FoldChange >= 5 ),]
#15
write.table(x = sigres_Bd21.3.leaves.bacteria.XT.2.dpi_strict,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.leaves.bacteria.XT.2.dpi_DEG_strict.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

write.table(x = sigres_Bd21.3.leaves.bacteria.XT.2.dpi_up_strict,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.leaves.bacteria.XT.2.dpi_DEG_strict.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

#Bd21.3.leaves.bacteria.XT.4.dpi
Bd21.3.leaves.bacteria.XT.4.dpi <- Bd21.3.leaves.bacteria.XT.4.dpi[!is.na(Bd21.3.leaves.bacteria.XT.4.dpi$padj),]
sigres_Bd21.3.leaves.bacteria.XT.4.dpi_up <- Bd21.3.leaves.bacteria.XT.4.dpi [(Bd21.3.leaves.bacteria.XT.4.dpi $padj<0.05 & Bd21.3.leaves.bacteria.XT.4.dpi $log2FoldChange>= 1),]
nb_Bd21.3.leaves.bacteria.XT.4.dpi_up <- nrow(sigres_Bd21.3.leaves.bacteria.XT.4.dpi_up)
# 2034
write.table(x = sigres_Bd21.3.leaves.bacteria.XT.4.dpi_up,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.leaves.bacteria.XT.4.dpi_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
sigres_Bd21.3.leaves.bacteria.XT.4.dpi_down <- Bd21.3.leaves.bacteria.XT.4.dpi[(Bd21.3.leaves.bacteria.XT.4.dpi$padj<0.05 & Bd21.3.leaves.bacteria.XT.4.dpi$log2FoldChange<= -1),]
nb_Bd21.3.leaves.bacteria.XT.4.dpi_down <- nrow(sigres_Bd21.3.leaves.bacteria.XT.4.dpi_down)
#1444
write.table(x = sigres_Bd21.3.leaves.bacteria.XT.4.dpi_down,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.leaves.bacteria.XT.4.dpi_down.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
sigres_Bd21.3.leaves.bacteria.XT.4.dpi <- Bd21.3.leaves.bacteria.XT.4.dpi[(Bd21.3.leaves.bacteria.XT.4.dpi$padj<0.05 & abs(Bd21.3.leaves.bacteria.XT.4.dpi$log2FoldChange) >= 1 ),]
nb_sigres_Bd21.3.leaves.bacteria.XT.4.dpi_DEG <- nrow(sigres_Bd21.3.leaves.bacteria.XT.4.dpi)
#3478
write.table(x = sigres_Bd21.3.leaves.bacteria.XT.4.dpi,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.leaves.bacteria.XT.4.dpi_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
#strict
sigres_Bd21.3.leaves.bacteria.XT.4.dpi_strict <- Bd21.3.leaves.bacteria.XT.4.dpi[(Bd21.3.leaves.bacteria.XT.4.dpi$padj<0.00001 & abs(Bd21.3.leaves.bacteria.XT.4.dpi$log2FoldChange) >= 5 ),]
nb_sigres_Bd21.3.leaves.bacteria.XT.4.dpi_DEG_strict <- nrow(sigres_Bd21.3.leaves.bacteria.XT.4.dpi_strict)
#75
union(rownames(sigres_Bd21.3.leaves.bacteria.XT.2.dpi_strict),rownames(sigres_Bd21.3.leaves.bacteria.XT.4.dpi_strict))
#[1] "Bradi1g25157.v3.2" "Bradi1g26870.v3.2" "Bradi1g38786.v3.2" "Bradi1g59800.v3.2"
#[5] "Bradi1g61397.v3.2" "Bradi1g63370.v3.2" "Bradi1g69690.v3.2" "Bradi1g74410.v3.2"
#[9] "Bradi2g06497.v3.2" "Bradi2g08730.v3.2" "Bradi2g27140.v3.2" "Bradi2g34380.v3.2"
#[13] "Bradi2g48760.v3.2" "Bradi2g51170.v3.2" "Bradi2g52090.v3.2" "Bradi2g52490.v3.2"
#[17] "Bradi2g53610.v3.2" "Bradi2g60860.v3.2" "Bradi3g44391.v3.2" "Bradi3g48770.v3.2"
#[21] "Bradi3g52220.v3.2" "Bradi3g54950.v3.2" "Bradi3g55270.v3.2" "Bradi3g59187.v3.2"
#[25] "Bradi3g60880.v3.2" "Bradi4g38006.v3.2" "Bradi4g38965.v3.2" "Bradi1g04720.v3.2"
#[29] "Bradi1g12440.v3.2" "Bradi1g15695.v3.2" "Bradi1g21870.v3.2" "Bradi1g22940.v3.2"
#[33] "Bradi1g34727.v3.2" "Bradi1g35340.v3.2" "Bradi1g39190.v3.2" "Bradi1g44370.v3.2"
#[37] "Bradi1g53630.v3.2" "Bradi1g54890.v3.2" "Bradi1g57525.v3.2" "Bradi1g57532.v3.2"
#[41] "Bradi1g61340.v3.2" "Bradi1g68530.v3.2" "Bradi1g69330.v3.2" "Bradi1g72500.v3.2"
#[45] "Bradi1g75735.v3.2" "Bradi1g76630.v3.2" "Bradi2g05950.v3.2" "Bradi2g06990.v3.2"
#[49] "Bradi2g25280.v3.2" "Bradi2g30800.v3.2" "Bradi2g38833.v3.2" "Bradi2g44420.v3.2"
#[53] "Bradi2g46120.v3.2" "Bradi2g49781.v3.2" "Bradi2g57460.v3.2" "Bradi2g60870.v3.2"
#[57] "Bradi3g10030.v3.2" "Bradi3g10500.v3.2" "Bradi3g15956.v3.2" "Bradi3g18070.v3.2"
#[61] "Bradi3g21030.v3.2" "Bradi3g22032.v3.2" "Bradi3g22215.v3.2" "Bradi3g22240.v3.2"
#[65] "Bradi3g27500.v3.2" "Bradi3g32950.v3.2" "Bradi3g35240.v3.2" "Bradi3g35660.v3.2"
#[69] "Bradi3g47390.v3.2" "Bradi3g49010.v3.2" "Bradi3g58770.v3.2" "Bradi3g58915.v3.2"
#[73] "Bradi3g59835.v3.2" "Bradi3g59842.v3.2" "Bradi3g59853.v3.2" "Bradi4g02480.v3.2"
#[77] "Bradi4g05040.v3.2" "Bradi4g20770.v3.2" "Bradi4g28612.v3.2" "Bradi4g30880.v3.2"
#[81] "Bradi4g34722.v3.2" "Bradi4g37220.v3.2" "Bradi4g39611.v3.2" "Bradi4g39950.v3.2"
#[85] "Bradi5g03477.v3.2" "Bradi5g06523.v3.2" "Bradi5g21447.v3.2" "Bradi5g23120.v3.2"

intersect(rownames(sigres_Bd21.3.leaves.bacteria.XT.2.dpi_strict),rownames(sigres_Bd21.3.leaves.bacteria.XT.4.dpi_strict))

#[1] "Bradi1g26870.v3.2" "Bradi1g59800.v3.2" "Bradi1g63370.v3.2" "Bradi1g69690.v3.2"
#[5] "Bradi1g74410.v3.2" "Bradi2g06497.v3.2" "Bradi2g34380.v3.2" "Bradi2g48760.v3.2"
#[9] "Bradi2g53610.v3.2" "Bradi3g52220.v3.2" "Bradi3g54950.v3.2" "Bradi3g55270.v3.2"
#[13] "Bradi3g59187.v3.2" "Bradi4g38006.v3.2"


write.table(x = sigres_Bd21.3.leaves.bacteria.XT.4.dpi_strict,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.leaves.bacteria.XT.4.dpi_DEG_strict.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
#Bd21.3.leaves.virus.BSMV.3.dpi
Bd21.3.leaves.virus.BSMV.3.dpi <- Bd21.3.leaves.virus.BSMV.3.dpi[!is.na(Bd21.3.leaves.virus.BSMV.3.dpi$padj),]
sigres_Bd21.3.leaves.virus.BSMV.3.dpi_up <- Bd21.3.leaves.virus.BSMV.3.dpi [(Bd21.3.leaves.virus.BSMV.3.dpi $padj<0.05 & Bd21.3.leaves.virus.BSMV.3.dpi $log2FoldChange>= 1),]
nb_Bd21.3.leaves.virus.BSMV.3.dpi_up <- nrow(sigres_Bd21.3.leaves.virus.BSMV.3.dpi_up)
# 696
write.table(x = sigres_Bd21.3.leaves.virus.BSMV.3.dpi_up,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.leaves.virus.BSMV.3.dpi_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

sigres_Bd21.3.leaves.virus.BSMV.3.dpi_down <- Bd21.3.leaves.virus.BSMV.3.dpi[(Bd21.3.leaves.virus.BSMV.3.dpi$padj<0.05 & Bd21.3.leaves.virus.BSMV.3.dpi$log2FoldChange<= -1),]
nb_Bd21.3.leaves.virus.BSMV.3.dpi_down <- nrow(sigres_Bd21.3.leaves.virus.BSMV.3.dpi_down)
#1574
write.table(x = sigres_Bd21.3.leaves.virus.BSMV.3.dpi_down,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.leaves.virus.BSMV.3.dpi_down.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

sigres_Bd21.3.leaves.virus.BSMV.3.dpi <- Bd21.3.leaves.virus.BSMV.3.dpi[(Bd21.3.leaves.virus.BSMV.3.dpi$padj<0.05 & abs(Bd21.3.leaves.virus.BSMV.3.dpi$log2FoldChange) >= 1 ),]
nb_sigres_Bd21.3.leaves.virus.BSMV.3.dpi_DEG <- nrow(sigres_Bd21.3.leaves.virus.BSMV.3.dpi)
#1907
write.table(x = sigres_Bd21.3.leaves.virus.BSMV.3.dpi,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.leaves.virus.BSMV.3.dpi_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

sigres_Bd21.3.leaves.virus.BSMV.3.dpi_strict <- Bd21.3.leaves.virus.BSMV.3.dpi[(Bd21.3.leaves.virus.BSMV.3.dpi$padj<0.00001 & abs(Bd21.3.leaves.virus.BSMV.3.dpi$log2FoldChange) >= 5 ),]
nb_sigres_Bd21.3.leaves.virus.BSMV.3.dpi_DEG_strict <- nrow(sigres_Bd21.3.leaves.virus.BSMV.3.dpi_strict)
#11
write.table(x = sigres_Bd21.3.leaves.virus.BSMV.3.dpi_strict,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.leaves.virus.BSMV.3.dpi_DEG_strict.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
#up
Bd21.3.leaves.virus.BSMV.3.dpi[(Bd21.3.leaves.virus.BSMV.3.dpi$padj<0.00001 & Bd21.3.leaves.virus.BSMV.3.dpi$log2FoldChange >= 5 ),]
#Bradi1g24340.v3.2 5.53908333946217e-18 1.13398883667139e-15
#Bradi1g37600.v3.2 5.41387617221511e-31 1.25245711682123e-27
#Bradi3g32940.v3.2 1.82565372083121e-28 2.49171305331446e-25
#Bradi5g20718.v3.2  4.4573891867709e-07 9.00531909139808e-06
#Bd21.3.leaves.virus.BSMV.7.dpi
Bd21.3.leaves.virus.BSMV.7.dpi <- Bd21.3.leaves.virus.BSMV.7.dpi[!is.na(Bd21.3.leaves.virus.BSMV.7.dpi$padj),]
sigres_Bd21.3.leaves.virus.BSMV.7.dpi_up <- Bd21.3.leaves.virus.BSMV.7.dpi [(Bd21.3.leaves.virus.BSMV.7.dpi $padj<0.05 & Bd21.3.leaves.virus.BSMV.7.dpi $log2FoldChange>= 1),]
nb_Bd21.3.leaves.virus.BSMV.7.dpi_up <- nrow(sigres_Bd21.3.leaves.virus.BSMV.7.dpi_up)
# 929
write.table(x = sigres_Bd21.3.leaves.virus.BSMV.7.dpi_up,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.leaves.virus.BSMV.7.dpi_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

sigres_Bd21.3.leaves.virus.BSMV.7.dpi_down <- Bd21.3.leaves.virus.BSMV.7.dpi[(Bd21.3.leaves.virus.BSMV.7.dpi$padj<0.05 & Bd21.3.leaves.virus.BSMV.7.dpi$log2FoldChange<= -1),]
nb_Bd21.3.leaves.virus.BSMV.7.dpi_down <- nrow(sigres_Bd21.3.leaves.virus.BSMV.7.dpi_down)
#911
write.table(x = sigres_Bd21.3.leaves.virus.BSMV.7.dpi_down,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.leaves.virus.BSMV.7.dpi_down.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
#sigres_Bd21.3.leaves.virus.BSMV.7.dpi.DEG
sigres_Bd21.3.leaves.virus.BSMV.7.dpi <- Bd21.3.leaves.virus.BSMV.7.dpi[(Bd21.3.leaves.virus.BSMV.7.dpi$padj<0.05 & abs(Bd21.3.leaves.virus.BSMV.7.dpi$log2FoldChange) >= 1 ),]
nb_sigres_Bd21.3.leaves.virus.BSMV.7.dpi_DEG <- nrow(sigres_Bd21.3.leaves.virus.BSMV.7.dpi)
#1840
write.table(x = sigres_Bd21.3.leaves.virus.BSMV.7.dpi,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_sigres_Bd21.3.leaves.virus.BSMV.7.dpi_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
#sigres_Bd21.3.leaves.virus.BSMV.7.dpi.Strict
sigres_Bd21.3.leaves.virus.BSMV.7.dpi.strict <- Bd21.3.leaves.virus.BSMV.7.dpi[(Bd21.3.leaves.virus.BSMV.7.dpi$padj<0.00001 & abs(Bd21.3.leaves.virus.BSMV.7.dpi$log2FoldChange) >= 5 ),]
nb_sigres_Bd21.3.leaves.virus.BSMV.7.dpi_DEG.strict <- nrow(sigres_Bd21.3.leaves.virus.BSMV.7.dpi.strict)
#27
write.table(x = sigres_Bd21.3.leaves.virus.BSMV.7.dpi.strict,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_sigres_Bd21.3.leaves.virus.BSMV.7.dpi_DEG_strict.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
sigres_Bd21.3.leaves.virus.BSMV.7.dpi.strict_up <- Bd21.3.leaves.virus.BSMV.7.dpi[(Bd21.3.leaves.virus.BSMV.7.dpi$padj<0.00001 & Bd21.3.leaves.virus.BSMV.7.dpi$log2FoldChange >= 5 ),]
#26
write.table(x = sigres_Bd21.3.leaves.virus.BSMV.7.dpi.strict_up,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_sigres_Bd21.3.leaves.virus.BSMV.7.dpi_DEG.strict_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

#Bd21.3.leaves.virus.BSMV.14.dpi
Bd21.3.leaves.virus.BSMV.14.dpi <- Bd21.3.leaves.virus.BSMV.14.dpi[!is.na(Bd21.3.leaves.virus.BSMV.14.dpi$padj),]
sigres_Bd21.3.leaves.virus.BSMV.14.dpi_up <- Bd21.3.leaves.virus.BSMV.14.dpi [(Bd21.3.leaves.virus.BSMV.14.dpi $padj<0.05 & Bd21.3.leaves.virus.BSMV.14.dpi $log2FoldChange>= 1),]
nb_Bd21.3.leaves.virus.BSMV.14.dpi_up <- nrow(sigres_Bd21.3.leaves.virus.BSMV.14.dpi_up)
# 1698
write.table(x = sigres_Bd21.3.leaves.virus.BSMV.14.dpi_up,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.leaves.virus.BSMV.14.dpi_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
sigres_Bd21.3.leaves.virus.BSMV.14.dpi_down <- Bd21.3.leaves.virus.BSMV.14.dpi[(Bd21.3.leaves.virus.BSMV.14.dpi$padj<0.05 & Bd21.3.leaves.virus.BSMV.14.dpi$log2FoldChange<= -1),]
nb_Bd21.3.leaves.virus.BSMV.14.dpi_down <- nrow(sigres_Bd21.3.leaves.virus.BSMV.14.dpi_down)
#1397
write.table(x = sigres_Bd21.3.leaves.virus.BSMV.14.dpi_down,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.leaves.virus.BSMV.14.dpi_down.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

#sigres_Bd21.3.leaves.virus.BSMV.14.dpi.DEG
sigres_Bd21.3.leaves.virus.BSMV.14.dpi <- Bd21.3.leaves.virus.BSMV.14.dpi[(Bd21.3.leaves.virus.BSMV.14.dpi$padj<0.05 & abs(Bd21.3.leaves.virus.BSMV.14.dpi$log2FoldChange) >= 1 ),]
nb_sigres_Bd21.3.leaves.virus.BSMV.14.dpi_DEG <- nrow(sigres_Bd21.3.leaves.virus.BSMV.14.dpi)
#2761
write.table(x = sigres_Bd21.3.leaves.virus.BSMV.14.dpi,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.leaves.virus.BSMV.14.dpi_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

#sigres_Bd21.3.leaves.virus.BSMV.14.dpi.DEG.strict
sigres_Bd21.3.leaves.virus.BSMV.14.dpi.strict <- Bd21.3.leaves.virus.BSMV.14.dpi[(Bd21.3.leaves.virus.BSMV.14.dpi$padj<0.00001 & abs(Bd21.3.leaves.virus.BSMV.14.dpi$log2FoldChange) >= 5 ),]
nb_sigres_Bd21.3.leaves.virus.BSMV.14.dpi_DEG.strict <- nrow(sigres_Bd21.3.leaves.virus.BSMV.14.dpi.strict)
#35
write.table(x = sigres_Bd21.3.leaves.virus.BSMV.14.dpi.strict,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.leaves.virus.BSMV.14.dpi_DEG.strict.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
sigres_Bd21.3.leaves.virus.BSMV.14.dpi.strict_up <- Bd21.3.leaves.virus.BSMV.14.dpi[(Bd21.3.leaves.virus.BSMV.14.dpi$padj<0.00001 & Bd21.3.leaves.virus.BSMV.14.dpi$log2FoldChange > 5 ),]
#35 so all of the 35 are the up regulating!!

intersect(intersect(rownames(sigres_Bd21.3.leaves.virus.BSMV.14.dpi.strict),rownames(sigres_Bd21.3.leaves.virus.BSMV.7.dpi.strict)),rownames(sigres_Bd21.3.leaves.virus.BSMV.3.dpi_strict))
#Nothing!!!
###MMMV
Bd21.3.leaves.virus.MMMV.3.dpi <- Bd21.3.leaves.virus.MMMV.3.dpi[!is.na(Bd21.3.leaves.virus.MMMV.3.dpi$padj),]
sigres_Bd21.3.leaves.virus.MMMV.3.dpi_up <- Bd21.3.leaves.virus.MMMV.3.dpi [(Bd21.3.leaves.virus.MMMV.3.dpi $padj<0.05 & Bd21.3.leaves.virus.MMMV.3.dpi $log2FoldChange>= 1),]
nb_Bd21.3.leaves.virus.MMMV.3.dpi_up <- nrow(sigres_Bd21.3.leaves.virus.MMMV.3.dpi_up)
# 1432
write.table(x = sigres_Bd21.3.leaves.virus.MMMV.3.dpi_up,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.leaves.virus.MMMV.3.dpi_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
sigres_Bd21.3.leaves.virus.MMMV.3.dpi_down <- Bd21.3.leaves.virus.MMMV.3.dpi[(Bd21.3.leaves.virus.MMMV.3.dpi$padj<0.05 & Bd21.3.leaves.virus.MMMV.3.dpi$log2FoldChange<= -1),]
nb_Bd21.3.leaves.virus.MMMV.3.dpi_down <- nrow(sigres_Bd21.3.leaves.virus.MMMV.3.dpi_down)
#586
write.table(x = sigres_Bd21.3.leaves.virus.MMMV.3.dpi_down,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.leaves.virus.MMMV.3.dpi_down.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
#sigres_Bd21.3.leaves.virus.MMMV.3.dpi.DEG
sigres_Bd21.3.leaves.virus.MMMV.3.dpi <- Bd21.3.leaves.virus.MMMV.3.dpi[(Bd21.3.leaves.virus.MMMV.3.dpi$padj<0.05 & abs(Bd21.3.leaves.virus.MMMV.3.dpi$log2FoldChange) >= 1 ),]
nb_sigres_Bd21.3.leaves.virus.MMMV.3.dpi_DEG <- nrow(sigres_Bd21.3.leaves.virus.MMMV.3.dpi)
#2235
write.table (x = sigres_Bd21.3.leaves.virus.MMMV.3.dpi,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.leaves.virus.MMMV.3.dpi_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
#sigres_Bd21.3.leaves.virus.MMMV.3.dpi.DEG.strict
sigres_Bd21.3.leaves.virus.MMMV.3.dpi.strict <- Bd21.3.leaves.virus.MMMV.3.dpi[(Bd21.3.leaves.virus.MMMV.3.dpi$padj<0.00001 & abs(Bd21.3.leaves.virus.MMMV.3.dpi$log2FoldChange) >= 5 ),]
nb_sigres_Bd21.3.leaves.virus.MMMV.3.dpi_DEG.strict <- nrow(sigres_Bd21.3.leaves.virus.MMMV.3.dpi.strict)
#24
write.table(x = sigres_Bd21.3.leaves.virus.MMMV.3.dpi.strict,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.leaves.virus.MMMV.3.dpi_DEG.strict.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
sigres_Bd21.3.leaves.virus.MMMV.3.dpi.strict_up <- Bd21.3.leaves.virus.MMMV.3.dpi[(Bd21.3.leaves.virus.MMMV.3.dpi$padj<0.00001 & Bd21.3.leaves.virus.MMMV.3.dpi$log2FoldChange >= 5 ),]
write.table(x = sigres_Bd21.3.leaves.virus.MMMV.3.dpi.strict_up,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.leaves.virus.MMMV.3.dpi_DEG.strict_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
#Bd21.3.leaves.virus.MMMV.7.dpi
Bd21.3.leaves.virus.MMMV.7.dpi <- Bd21.3.leaves.virus.MMMV.7.dpi[!is.na(Bd21.3.leaves.virus.MMMV.7.dpi$padj),]
sigres_Bd21.3.leaves.virus.MMMV.7.dpi_up <- Bd21.3.leaves.virus.MMMV.7.dpi [(Bd21.3.leaves.virus.MMMV.7.dpi $padj<0.05 & Bd21.3.leaves.virus.MMMV.7.dpi $log2FoldChange>= 1),]
nb_Bd21.3.leaves.virus.MMMV.7.dpi_up <- nrow(sigres_Bd21.3.leaves.virus.MMMV.7.dpi_up)
# 2225
write.table(x = sigres_Bd21.3.leaves.virus.MMMV.7.dpi_up,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.leaves.virus.MMMV.7.dpi_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
sigres_Bd21.3.leaves.virus.MMMV.7.dpi_down <- Bd21.3.leaves.virus.MMMV.7.dpi[(Bd21.3.leaves.virus.MMMV.7.dpi$padj<0.05 & Bd21.3.leaves.virus.MMMV.7.dpi$log2FoldChange<= -1),]
nb_Bd21.3.leaves.virus.MMMV.7.dpi_down <- nrow(sigres_Bd21.3.leaves.virus.MMMV.7.dpi_down)
#2262
write.table(x = sigres_Bd21.3.leaves.virus.MMMV.7.dpi_down,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.leaves.virus.MMMV.7.dpi_down.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

#sigres_Bd21.3.leaves.virus.MMMV.7.dpi.strict
sigres_Bd21.3.leaves.virus.MMMV.7.dpi.strict <- Bd21.3.leaves.virus.MMMV.7.dpi[(Bd21.3.leaves.virus.MMMV.7.dpi$padj<0.00001 & abs(Bd21.3.leaves.virus.MMMV.7.dpi$log2FoldChange) >= 5 ),]
nb_sigres_Bd21.3.leaves.virus.MMMV.7.dpi_DEG.strict <- nrow(sigres_Bd21.3.leaves.virus.MMMV.7.dpi.strict)
#226
write.table(x = sigres_Bd21.3.leaves.virus.MMMV.7.dpi.strict,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_sigres_Bd21.3.leaves.virus.MMMV.7.dpi_DEG.strict.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

#up_strct
sigres_Bd21.3.leaves.virus.MMMV.7.dpi.strict_up <- Bd21.3.leaves.virus.MMMV.7.dpi[(Bd21.3.leaves.virus.MMMV.7.dpi$padj<0.00001 & Bd21.3.leaves.virus.MMMV.7.dpi$log2FoldChange >= 5 ),]
nrow(sigres_Bd21.3.leaves.virus.MMMV.7.dpi.strict_up)
#220
write.table(x = sigres_Bd21.3.leaves.virus.MMMV.7.dpi.strict_up,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_sigres_Bd21.3.leaves.virus.MMMV.7.dpi_DEG.strict_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
#Bd21.3.leaves.virus.MMMV.14.dpi
Bd21.3.leaves.virus.MMMV.14.dpi <- Bd21.3.leaves.virus.MMMV.14.dpi[!is.na(Bd21.3.leaves.virus.MMMV.14.dpi$padj),]
sigres_Bd21.3.leaves.virus.MMMV.14.dpi_up <- Bd21.3.leaves.virus.MMMV.14.dpi [(Bd21.3.leaves.virus.MMMV.14.dpi $padj<0.05 & Bd21.3.leaves.virus.MMMV.14.dpi $log2FoldChange>= 1),]
nb_Bd21.3.leaves.virus.MMMV.14.dpi_up <- nrow(sigres_Bd21.3.leaves.virus.MMMV.14.dpi_up)
# 3919
write.table(x = sigres_Bd21.3.leaves.virus.MMMV.14.dpi_up,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.leaves.virus.MMMV.14.dpi_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
sigres_Bd21.3.leaves.virus.MMMV.14.dpi_down <- Bd21.3.leaves.virus.MMMV.14.dpi[(Bd21.3.leaves.virus.MMMV.14.dpi$padj<0.05 & Bd21.3.leaves.virus.MMMV.14.dpi$log2FoldChange<= -1),]
nb_Bd21.3.leaves.virus.MMMV.14.dpi_down <- nrow(sigres_Bd21.3.leaves.virus.MMMV.14.dpi_down)
#1961
write.table(x = sigres_Bd21.3.leaves.virus.MMMV.14.dpi_down,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.leaves.virus.MMMV.14.dpi_down.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
#DEG
sigres_Bd21.3.leaves.virus.MMMV.14.dpi <- Bd21.3.leaves.virus.MMMV.14.dpi[(Bd21.3.leaves.virus.MMMV.14.dpi$padj<0.05 & abs(Bd21.3.leaves.virus.MMMV.14.dpi$log2FoldChange) >= 1 ),]
nb_sigres_Bd21.3.leaves.virus.MMMV.14.dpi_DEG <- nrow(sigres_Bd21.3.leaves.virus.MMMV.14.dpi)
#5880
write.table(x = sigres_Bd21.3.leaves.virus.MMMV.14.dpi,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.leaves.virus.MMMV.14.dpi_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
#DEG.strict
sigres_Bd21.3.leaves.virus.MMMV.14.dpi.strict <- Bd21.3.leaves.virus.MMMV.14.dpi[(Bd21.3.leaves.virus.MMMV.14.dpi$padj<0.00001 & abs(Bd21.3.leaves.virus.MMMV.14.dpi$log2FoldChange) >= 5 ),]
nb_sigres_Bd21.3.leaves.virus.MMMV.14.dpi_DEG.strict <- nrow(sigres_Bd21.3.leaves.virus.MMMV.14.dpi.strict)
#362
write.table(x = sigres_Bd21.3.leaves.virus.MMMV.14.dpi.strict,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.leaves.virus.MMMV.14.dpi_DEG.strict.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

sigres_Bd21.3.leaves.virus.MMMV.14.dpi.strict_up <- Bd21.3.leaves.virus.MMMV.14.dpi[(Bd21.3.leaves.virus.MMMV.14.dpi$padj<0.00001 & Bd21.3.leaves.virus.MMMV.14.dpi$log2FoldChange >= 5 ),]
write.table(x = sigres_Bd21.3.leaves.virus.MMMV.14.dpi.strict_up,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.leaves.virus.MMMV.14.dpi_DEG.strict_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

intersect(intersect(rownames(sigres_Bd21.3.leaves.virus.MMMV.14.dpi.strict),rownames(sigres_Bd21.3.leaves.virus.MMMV.7.dpi.strict)),rownames(sigres_Bd21.3.leaves.virus.MMMV.3.dpi.strict))

##WSMV
Bd21.3.leaves.virus.WSMV.3.dpi <- Bd21.3.leaves.virus.WSMV.3.dpi[!is.na(Bd21.3.leaves.virus.WSMV.3.dpi$padj),]
sigres_Bd21.3.leaves.virus.WSMV.3.dpi_up <- Bd21.3.leaves.virus.WSMV.3.dpi [(Bd21.3.leaves.virus.WSMV.3.dpi $padj<0.05 & Bd21.3.leaves.virus.WSMV.3.dpi $log2FoldChange>= 1),]
nb_Bd21.3.leaves.virus.WSMV.3.dpi_up <- nrow(sigres_Bd21.3.leaves.virus.WSMV.3.dpi_up)
# 736
write.table(x = sigres_Bd21.3.leaves.virus.WSMV.3.dpi_up,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.leaves.virus.WSMV.3.dpi_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
sigres_Bd21.3.leaves.virus.WSMV.3.dpi_down <- Bd21.3.leaves.virus.WSMV.3.dpi[(Bd21.3.leaves.virus.WSMV.3.dpi$padj<0.05 & Bd21.3.leaves.virus.WSMV.3.dpi$log2FoldChange<= -1),]
nb_Bd21.3.leaves.virus.WSMV.3.dpi_down <- nrow(sigres_Bd21.3.leaves.virus.WSMV.3.dpi_down)
#694
write.table(x = sigres_Bd21.3.leaves.virus.WSMV.3.dpi_down,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.leaves.virus.WSMV.3.dpi_down.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
##WSMV DEGs
sigres_Bd21.3.leaves.virus.WSMV.3.dpi <- Bd21.3.leaves.virus.WSMV.3.dpi[(Bd21.3.leaves.virus.WSMV.3.dpi$padj<0.05 & abs(Bd21.3.leaves.virus.WSMV.3.dpi$log2FoldChange) >= 1 ),]
nb_sigres_Bd21.3.leaves.virus.WSMV.3.dpi_DEG <- nrow(sigres_Bd21.3.leaves.virus.WSMV.3.dpi)
#1430
write.table(x = sigres_Bd21.3.leaves.virus.WSMV.3.dpi,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.leaves.virus.WSMV.3.dpi_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
##WSMV DEGs strict
sigres_Bd21.3.leaves.virus.WSMV.3.dpi.strict <- Bd21.3.leaves.virus.WSMV.3.dpi[(Bd21.3.leaves.virus.WSMV.3.dpi$padj<0.00001 & abs(Bd21.3.leaves.virus.WSMV.3.dpi$log2FoldChange) >= 5 ),]
nb_sigres_Bd21.3.leaves.virus.WSMV.3.dpi_DEG.strict <- nrow(sigres_Bd21.3.leaves.virus.WSMV.3.dpi.strict)
#3
write.table(x = sigres_Bd21.3.leaves.virus.WSMV.3.dpi.strict,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.leaves.virus.WSMV.3.dpi_DEG.strict.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
Bd21.3.leaves.virus.WSMV.3.dpi[(Bd21.3.leaves.virus.WSMV.3.dpi$padj<0.00001 & Bd21.3.leaves.virus.WSMV.3.dpi$log2FoldChange >= 5 ),]
#Bradi2g61050.v3.2

#Bd21.3.leaves.virus.WSMV.7.dpi
Bd21.3.leaves.virus.WSMV.7.dpi <- Bd21.3.leaves.virus.WSMV.7.dpi[!is.na(Bd21.3.leaves.virus.WSMV.7.dpi$padj),]
sigres_Bd21.3.leaves.virus.WSMV.7.dpi_up <- Bd21.3.leaves.virus.WSMV.7.dpi [(Bd21.3.leaves.virus.WSMV.7.dpi $padj<0.05 & Bd21.3.leaves.virus.WSMV.7.dpi $log2FoldChange>= 1),]
nb_Bd21.3.leaves.virus.WSMV.7.dpi_up <- nrow(sigres_Bd21.3.leaves.virus.WSMV.7.dpi_up)
# 613
write.table(x = sigres_Bd21.3.leaves.virus.WSMV.7.dpi_up,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.leaves.virus.WSMV.7.dpi_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
sigres_Bd21.3.leaves.virus.WSMV.7.dpi_down <- Bd21.3.leaves.virus.WSMV.7.dpi[(Bd21.3.leaves.virus.WSMV.7.dpi$padj<0.05 & Bd21.3.leaves.virus.WSMV.7.dpi$log2FoldChange<= -1),]
nb_Bd21.3.leaves.virus.WSMV.7.dpi_down <- nrow(sigres_Bd21.3.leaves.virus.WSMV.7.dpi_down)
#1375
write.table(x = sigres_Bd21.3.leaves.virus.WSMV.7.dpi_down,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.leaves.virus.WSMV.7.dpi_down.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
sigres_Bd21.3.leaves.virus.WSMV.7.dpi <- Bd21.3.leaves.virus.WSMV.7.dpi[(Bd21.3.leaves.virus.WSMV.7.dpi$padj<0.05 & abs(Bd21.3.leaves.virus.WSMV.7.dpi$log2FoldChange) >= 1 ),]
nb_sigres_Bd21.3.leaves.virus.WSMV.7.dpi_DEG <- nrow(sigres_Bd21.3.leaves.virus.WSMV.7.dpi)
#1988
write.table(x = sigres_Bd21.3.leaves.virus.WSMV.7.dpi,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_sigres_Bd21.3.leaves.virus.WSMV.7.dpi_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
#DEG.strict
sigres_Bd21.3.leaves.virus.WSMV.7.dpi.strict <- Bd21.3.leaves.virus.WSMV.7.dpi[(Bd21.3.leaves.virus.WSMV.7.dpi$padj<0.00001 & abs(Bd21.3.leaves.virus.WSMV.7.dpi$log2FoldChange) >= 5 ),]
nb_sigres_Bd21.3.leaves.virus.WSMV.7.dpi_DEG.strict <- nrow(sigres_Bd21.3.leaves.virus.WSMV.7.dpi.strict)
#8
write.table(x = sigres_Bd21.3.leaves.virus.WSMV.7.dpi.strict,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_sigres_Bd21.3.leaves.virus.WSMV.7.dpi_DEG.strict.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
Bd21.3.leaves.virus.WSMV.7.dpi[(Bd21.3.leaves.virus.WSMV.7.dpi$padj<0.00001 & Bd21.3.leaves.virus.WSMV.7.dpi$log2FoldChange < -5 ),]
#all of them are down-regulating!!!
#Bd21.3.leaves.virus.WSMV.14.dpi
Bd21.3.leaves.virus.WSMV.14.dpi <- Bd21.3.leaves.virus.WSMV.14.dpi[!is.na(Bd21.3.leaves.virus.WSMV.14.dpi$padj),]
sigres_Bd21.3.leaves.virus.WSMV.14.dpi_up <- Bd21.3.leaves.virus.WSMV.14.dpi [(Bd21.3.leaves.virus.WSMV.14.dpi $padj<0.05 & Bd21.3.leaves.virus.WSMV.14.dpi $log2FoldChange>= 1),]
nb_Bd21.3.leaves.virus.WSMV.14.dpi_up <- nrow(sigres_Bd21.3.leaves.virus.WSMV.14.dpi_up)
# 454
write.table(x = sigres_Bd21.3.leaves.virus.WSMV.14.dpi_up,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.leaves.virus.WSMV.14.dpi_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
sigres_Bd21.3.leaves.virus.WSMV.14.dpi_down <- Bd21.3.leaves.virus.WSMV.14.dpi[(Bd21.3.leaves.virus.WSMV.14.dpi$padj<0.05 & Bd21.3.leaves.virus.WSMV.14.dpi$log2FoldChange<= -1),]
nb_Bd21.3.leaves.virus.WSMV.14.dpi_down <- nrow(sigres_Bd21.3.leaves.virus.WSMV.14.dpi_down)
#695
write.table(x = sigres_Bd21.3.leaves.virus.WSMV.14.dpi_down,
           file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.leaves.virus.WSMV.14.dpi_down.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
#DEG
sigres_Bd21.3.leaves.virus.WSMV.14.dpi <- Bd21.3.leaves.virus.WSMV.14.dpi[(Bd21.3.leaves.virus.WSMV.14.dpi$padj<0.05 & abs(Bd21.3.leaves.virus.WSMV.14.dpi$log2FoldChange) >= 1 ),]
nb_sigres_Bd21.3.leaves.virus.WSMV.14.dpi_DEG <- nrow(sigres_Bd21.3.leaves.virus.WSMV.14.dpi)
#1149
write.table(x = sigres_Bd21.3.leaves.virus.WSMV.14.dpi,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.leaves.virus.WSMV.14.dpi_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

#DEG.strict
sigres_Bd21.3.leaves.virus.WSMV.14.dpi.strict <- Bd21.3.leaves.virus.WSMV.14.dpi[(Bd21.3.leaves.virus.WSMV.14.dpi$padj<0.00001 & abs(Bd21.3.leaves.virus.WSMV.14.dpi$log2FoldChange) >= 5 ),]
nb_sigres_Bd21.3.leaves.virus.WSMV.14.dpi_DEG.strict <- nrow(sigres_Bd21.3.leaves.virus.WSMV.14.dpi.strict)
#2
write.table(x = sigres_Bd21.3.leaves.virus.WSMV.14.dpi.strict,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.leaves.virus.WSMV.14.dpi_DEG.strict.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

Bd21.3.leaves.virus.WSMV.14.dpi[(Bd21.3.leaves.virus.WSMV.14.dpi$padj<0.00001 & Bd21.3.leaves.virus.WSMV.14.dpi$log2FoldChange >= 5 ),]
#Bradi1g38786.v3.2

#Root virus
#Bd21.3.root.virus.BSMV.3.dpi
Bd21.3.root.virus.BSMV.3.dpi <- Bd21.3.root.virus.BSMV.3.dpi[!is.na(Bd21.3.root.virus.BSMV.3.dpi$padj),]
sigres_Bd21.3.root.virus.BSMV.3.dpi_up <- Bd21.3.root.virus.BSMV.3.dpi [(Bd21.3.root.virus.BSMV.3.dpi $padj<0.05 & Bd21.3.root.virus.BSMV.3.dpi $log2FoldChange>= 1),]
nb_Bd21.3.root.virus.BSMV.3.dpi_up <- nrow(sigres_Bd21.3.root.virus.BSMV.3.dpi_up)
# 945
write.table(x = sigres_Bd21.3.root.virus.BSMV.3.dpi_up,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.root.virus.BSMV.3.dpi_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

sigres_Bd21.3.root.virus.BSMV.3.dpi_down <- Bd21.3.root.virus.BSMV.3.dpi[(Bd21.3.root.virus.BSMV.3.dpi$padj<0.05 & Bd21.3.root.virus.BSMV.3.dpi$log2FoldChange<= -1),]
nb_Bd21.3.root.virus.BSMV.3.dpi_down <- nrow(sigres_Bd21.3.root.virus.BSMV.3.dpi_down)
#1344
write.table(x = sigres_Bd21.3.root.virus.BSMV.3.dpi_down,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.root.virus.BSMV.3.dpi_down.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
#DEG
sigres_Bd21.3.root.virus.BSMV.3.dpi <- Bd21.3.root.virus.BSMV.3.dpi[(Bd21.3.root.virus.BSMV.3.dpi$padj<0.05 & abs(Bd21.3.root.virus.BSMV.3.dpi$log2FoldChange) >= 1 ),]
nb_sigres_Bd21.3.root.virus.BSMV.3.dpi_DEG <- nrow(sigres_Bd21.3.root.virus.BSMV.3.dpi)
#2289
write.table(x = sigres_Bd21.3.root.virus.BSMV.3.dpi,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.root.virus.BSMV.3.dpi_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
#DEG.strict
sigres_Bd21.3.root.virus.BSMV.3.dpi.strict <- Bd21.3.root.virus.BSMV.3.dpi[(Bd21.3.root.virus.BSMV.3.dpi$padj<0.00001 & abs(Bd21.3.root.virus.BSMV.3.dpi$log2FoldChange) >= 5 ),]
nb_sigres_Bd21.3.root.virus.BSMV.3.dpi_DEG.strict <- nrow(sigres_Bd21.3.root.virus.BSMV.3.dpi.strict)
#16
write.table(x = sigres_Bd21.3.root.virus.BSMV.3.dpi.strict,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.root.virus.BSMV.3.dpi_DEG.strict.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

sigres_Bd21.3.root.virus.BSMV.3.dpi.strict_up <- Bd21.3.root.virus.BSMV.3.dpi[(Bd21.3.root.virus.BSMV.3.dpi$padj<0.00001 & Bd21.3.root.virus.BSMV.3.dpi$log2FoldChange >= 5 ),]
write.table(x = sigres_Bd21.3.root.virus.BSMV.3.dpi.strict_up,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.root.virus.BSMV.3.dpi_DEG.strict_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
#Bd21.3.root.virus.BSMV.7.dpi
Bd21.3.root.virus.BSMV.7.dpi <- Bd21.3.root.virus.BSMV.7.dpi[!is.na(Bd21.3.root.virus.BSMV.7.dpi$padj),]
sigres_Bd21.3.root.virus.BSMV.7.dpi_up <- Bd21.3.root.virus.BSMV.7.dpi [(Bd21.3.root.virus.BSMV.7.dpi $padj<0.05 & Bd21.3.root.virus.BSMV.7.dpi $log2FoldChange>= 1),]
nb_Bd21.3.root.virus.BSMV.7.dpi_up_down <- nrow(sigres_Bd21.3.root.virus.BSMV.7.dpi_up)
# 967
write.table(x = sigres_Bd21.3.root.virus.BSMV.7.dpi_up,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.root.virus.BSMV.7.dpi_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
#Bd21.3.root.virus.BSMV.7.dpi
Bd21.3.root.virus.BSMV.7.dpi <- Bd21.3.root.virus.BSMV.7.dpi[!is.na(Bd21.3.root.virus.BSMV.7.dpi$padj),]
sigres_Bd21.3.root.virus.BSMV.7.dpi_up <- Bd21.3.root.virus.BSMV.7.dpi [(Bd21.3.root.virus.BSMV.7.dpi $padj<0.05 & Bd21.3.root.virus.BSMV.7.dpi $log2FoldChange>= 1),]
nb_Bd21.3.root.virus.BSMV.7.dpi_up_down <- nrow(sigres_Bd21.3.root.virus.BSMV.7.dpi_up)
# 967
write.table(x = sigres_Bd21.3.root.virus.BSMV.7.dpi_up,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.root.virus.BSMV.7.dpi_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
sigres_Bd21.3.root.virus.BSMV.7.dpi_down <- Bd21.3.root.virus.BSMV.7.dpi[(Bd21.3.root.virus.BSMV.7.dpi$padj<0.05 & Bd21.3.root.virus.BSMV.7.dpi$log2FoldChange<= -1),]
nb_Bd21.3.root.virus.BSMV.7.dpi_down <- nrow(sigres_Bd21.3.root.virus.BSMV.7.dpi_down)
#847
write.table(x = sigres_Bd21.3.root.virus.BSMV.7.dpi_down,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.root.virus.BSMV.7.dpi_down.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
#DEGs
sigres_Bd21.3.root.virus.BSMV.7.dpi <- Bd21.3.root.virus.BSMV.7.dpi[(Bd21.3.root.virus.BSMV.7.dpi$padj<0.05 & abs(Bd21.3.root.virus.BSMV.7.dpi$log2FoldChange) >= 1 ),]
nb_sigres_Bd21.3.root.virus.BSMV.7.dpi_DEG <- nrow(sigres_Bd21.3.root.virus.BSMV.7.dpi)
#1814
write.table(x = sigres_Bd21.3.root.virus.BSMV.7.dpi,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_sigres_Bd21.3.root.virus.BSMV.7.dpi_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

#DEGs.strict.
sigres_Bd21.3.root.virus.BSMV.7.dpi.strict <- Bd21.3.root.virus.BSMV.7.dpi[(Bd21.3.root.virus.BSMV.7.dpi$padj<0.00001 & abs(Bd21.3.root.virus.BSMV.7.dpi$log2FoldChange) >= 5 ),]
nb_sigres_Bd21.3.root.virus.BSMV.7.dpi_DEG.strict <- nrow(sigres_Bd21.3.root.virus.BSMV.7.dpi.strict)
#21
write.table(x = sigres_Bd21.3.root.virus.BSMV.7.dpi.strict,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_sigres_Bd21.3.root.virus.BSMV.7.dpi_DEG.strict.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
sigres_Bd21.3.root.virus.BSMV.7.dpi.strict_up <- Bd21.3.root.virus.BSMV.7.dpi[(Bd21.3.root.virus.BSMV.7.dpi$padj<0.00001 & Bd21.3.root.virus.BSMV.7.dpi$log2FoldChange >= 5 ),]
write.table(x = sigres_Bd21.3.root.virus.BSMV.7.dpi.strict_up,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_sigres_Bd21.3.root.virus.BSMV.7.dpi_DEG.strict_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
#Bd21.3.root.virus.BSMV.14.dpi
Bd21.3.root.virus.BSMV.14.dpi <- Bd21.3.root.virus.BSMV.14.dpi[!is.na(Bd21.3.root.virus.BSMV.14.dpi$padj),]
sigres_Bd21.3.root.virus.BSMV.14.dpi_up <- Bd21.3.root.virus.BSMV.14.dpi [(Bd21.3.root.virus.BSMV.14.dpi $padj<0.05 & Bd21.3.root.virus.BSMV.14.dpi $log2FoldChange>= 1),]
nb_Bd21.3.root.virus.BSMV.14.dpi_up <- nrow(sigres_Bd21.3.root.virus.BSMV.14.dpi_up)
# 2637
write.table(x = sigres_Bd21.3.root.virus.BSMV.14.dpi_up,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.root.virus.BSMV.14.dpi_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
sigres_Bd21.3.root.virus.BSMV.14.dpi_down <- Bd21.3.root.virus.BSMV.14.dpi[(Bd21.3.root.virus.BSMV.14.dpi$padj<0.05 & Bd21.3.root.virus.BSMV.14.dpi$log2FoldChange<= -1),]
nb_Bd21.3.root.virus.BSMV.14.dpi_down <- nrow(sigres_Bd21.3.root.virus.BSMV.14.dpi_down)
#1299
write.table(x = sigres_Bd21.3.root.virus.BSMV.14.dpi_down,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.root.virus.BSMV.14.dpi_down.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
#DGE
sigres_Bd21.3.root.virus.BSMV.14.dpi <- Bd21.3.root.virus.BSMV.14.dpi[(Bd21.3.root.virus.BSMV.14.dpi$padj<0.05 & abs(Bd21.3.root.virus.BSMV.14.dpi$log2FoldChange) >= 1 ),]
nb_sigres_Bd21.3.root.virus.BSMV.14.dpi_DEG <- nrow(sigres_Bd21.3.root.virus.BSMV.14.dpi)
#3936
write.table(x = sigres_Bd21.3.root.virus.BSMV.14.dpi,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.root.virus.BSMV.14.dpi_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

#DGE.strict
sigres_Bd21.3.root.virus.BSMV.14.dpi.strict <- Bd21.3.root.virus.BSMV.14.dpi[(Bd21.3.root.virus.BSMV.14.dpi$padj<0.00001 & abs(Bd21.3.root.virus.BSMV.14.dpi$log2FoldChange) >= 5 ),]
nb_sigres_Bd21.3.root.virus.BSMV.14.dpi_DEG.strict <- nrow(sigres_Bd21.3.root.virus.BSMV.14.dpi.strict)
#36
write.table(x = sigres_Bd21.3.root.virus.BSMV.14.dpi.strict,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.root.virus.BSMV.14.dpi_DEG.strict.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
sigres_Bd21.3.root.virus.BSMV.14.dpi.strict_up <- Bd21.3.root.virus.BSMV.14.dpi[(Bd21.3.root.virus.BSMV.14.dpi$padj<0.00001 & Bd21.3.root.virus.BSMV.14.dpi$log2FoldChange >= 5 ),]
write.table(x = sigres_Bd21.3.root.virus.BSMV.14.dpi.strict_up,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.root.virus.BSMV.14.dpi_DEG.strict_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
intersect(intersect(rownames(sigres_Bd21.3.root.virus.BSMV.14.dpi.strict),rownames(sigres_Bd21.3.root.virus.BSMV.7.dpi.strict)),rownames(sigres_Bd21.3.root.virus.BSMV.3.dpi.strict))
#[1] "Bradi2g60870.v3.2" "Bradi3g20980.v3.2" "Bradi3g42251.v3.2" "Bradi5g04630.v3.2"
intersect(rownames(sigres_Bd21.3.leaves.virus.BSMV.3.dpi_strict),rownames(sigres_Bd21.3.root.virus.BSMV.3.dpi.strict))
#character(0)
###MMMV
Bd21.3.root.virus.MMMV.3.dpi <- Bd21.3.root.virus.MMMV.3.dpi[!is.na(Bd21.3.root.virus.MMMV.3.dpi$padj),]
sigres_Bd21.3.root.virus.MMMV.3.dpi_up <- Bd21.3.root.virus.MMMV.3.dpi [(Bd21.3.root.virus.MMMV.3.dpi $padj<0.05 & Bd21.3.root.virus.MMMV.3.dpi $log2FoldChange>= 1),]
nb_Bd21.3.root.virus.MMMV.3.dpi_up <- nrow(sigres_Bd21.3.root.virus.MMMV.3.dpi_up)
# 1369
write.table(x = sigres_Bd21.3.root.virus.MMMV.3.dpi_up,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.root.virus.MMMV.3.dpi_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
sigres_Bd21.3.root.virus.MMMV.3.dpi_down <- Bd21.3.root.virus.MMMV.3.dpi[(Bd21.3.root.virus.MMMV.3.dpi$padj<0.05 & Bd21.3.root.virus.MMMV.3.dpi$log2FoldChange<= -1),]
nb_Bd21.3.root.virus.MMMV.3.dpi_down <- nrow(sigres_Bd21.3.root.virus.MMMV.3.dpi_down)
#770
write.table(x = sigres_Bd21.3.root.virus.MMMV.3.dpi_down,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.root.virus.MMMV.3.dpi_down.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
#DEG
sigres_Bd21.3.root.virus.MMMV.3.dpi <- Bd21.3.root.virus.MMMV.3.dpi[(Bd21.3.root.virus.MMMV.3.dpi$padj<0.05 & abs(Bd21.3.root.virus.MMMV.3.dpi$log2FoldChange) >= 1 ),]
nb_sigres_Bd21.3.root.virus.MMMV.3.dpi_DEG <- nrow(sigres_Bd21.3.root.virus.MMMV.3.dpi)
#2139
write.table(x = sigres_Bd21.3.root.virus.MMMV.3.dpi,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.root.virus.MMMV.3.dpi_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

#DEG.strict
sigres_Bd21.3.root.virus.MMMV.3.dpi.strict <- Bd21.3.root.virus.MMMV.3.dpi[(Bd21.3.root.virus.MMMV.3.dpi$padj<0.00001 & abs(Bd21.3.root.virus.MMMV.3.dpi$log2FoldChange) >= 5 ),]
nb_sigres_Bd21.3.root.virus.MMMV.3.dpi_DEG.strict <- nrow(sigres_Bd21.3.root.virus.MMMV.3.dpi.strict)
#73
write.table(x = sigres_Bd21.3.root.virus.MMMV.3.dpi.strict,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.root.virus.MMMV.3.dpi_DEG.strict.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
sigres_Bd21.3.root.virus.MMMV.3.dpi.strict_up <- Bd21.3.root.virus.MMMV.3.dpi[(Bd21.3.root.virus.MMMV.3.dpi$padj<0.00001 & Bd21.3.root.virus.MMMV.3.dpi$log2FoldChange >= 5 ),]
write.table(x = sigres_Bd21.3.root.virus.MMMV.3.dpi.strict_up,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.root.virus.MMMV.3.dpi_DEG.strict_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
intersect(rownames(sigres_Bd21.3.leaves.virus.MMMV.3.dpi.strict),rownames(sigres_Bd21.3.root.virus.MMMV.3.dpi.strict))
#character(0)

#Bd21.3.root.virus.MMMV.7.dpi
Bd21.3.root.virus.MMMV.7.dpi <- Bd21.3.root.virus.MMMV.7.dpi[!is.na(Bd21.3.root.virus.MMMV.7.dpi$padj),]
sigres_Bd21.3.root.virus.MMMV.7.dpi_up <- Bd21.3.root.virus.MMMV.7.dpi [(Bd21.3.root.virus.MMMV.7.dpi $padj<0.05 & Bd21.3.root.virus.MMMV.7.dpi $log2FoldChange>= 1),]
nb_Bd21.3.root.virus.MMMV.7.dpi_up <- nrow(sigres_Bd21.3.root.virus.MMMV.7.dpi_up)
# 1776
write.table(x = sigres_Bd21.3.root.virus.MMMV.7.dpi_up,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.root.virus.MMMV.7.dpi_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
sigres_Bd21.3.root.virus.MMMV.7.dpi_down <- Bd21.3.root.virus.MMMV.7.dpi[(Bd21.3.root.virus.MMMV.7.dpi$padj<0.05 & Bd21.3.root.virus.MMMV.7.dpi$log2FoldChange<= -1),]
nb_Bd21.3.root.virus.MMMV.7.dpi_down <- nrow(sigres_Bd21.3.root.virus.MMMV.7.dpi_down)
#1031
write.table(x = sigres_Bd21.3.root.virus.MMMV.7.dpi_down,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.root.virus.MMMV.7.dpi_down.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
sigres_Bd21.3.root.virus.MMMV.7.dpi <- Bd21.3.root.virus.MMMV.7.dpi[(Bd21.3.root.virus.MMMV.7.dpi$padj<0.05 & abs(Bd21.3.root.virus.MMMV.7.dpi$log2FoldChange) >= 1 ),]
nb_sigres_Bd21.3.root.virus.MMMV.7.dpi_DEG <- nrow(sigres_Bd21.3.root.virus.MMMV.7.dpi)
#2986
write.table(x = sigres_Bd21.3.root.virus.MMMV.7.dpi,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_sigres_Bd21.3.root.virus.MMMV.7.dpi_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
#DEG.strict
sigres_Bd21.3.root.virus.MMMV.7.dpi.strict <- Bd21.3.root.virus.MMMV.7.dpi[(Bd21.3.root.virus.MMMV.7.dpi$padj<0.00001 & abs(Bd21.3.root.virus.MMMV.7.dpi$log2FoldChange) >= 5 ),]
nb_sigres_Bd21.3.root.virus.MMMV.7.dpi_DEG.strict <- nrow(sigres_Bd21.3.root.virus.MMMV.7.dpi.strict)
#21
write.table(x = sigres_Bd21.3.root.virus.MMMV.7.dpi.strict,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_sigres_Bd21.3.root.virus.MMMV.7.dpi_DEG.strict.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
sigres_Bd21.3.root.virus.MMMV.7.dpi.strict_up <- Bd21.3.root.virus.MMMV.7.dpi[(Bd21.3.root.virus.MMMV.7.dpi$padj<0.00001 & Bd21.3.root.virus.MMMV.7.dpi$log2FoldChange >= 5 ),]
write.table(x = sigres_Bd21.3.root.virus.MMMV.7.dpi.strict_up,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_sigres_Bd21.3.root.virus.MMMV.7.dpi_DEG.strict_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
#Bd21.3.root.virus.MMMV.14.dpi
Bd21.3.root.virus.MMMV.14.dpi <- Bd21.3.root.virus.MMMV.14.dpi[!is.na(Bd21.3.root.virus.MMMV.14.dpi$padj),]
sigres_Bd21.3.root.virus.MMMV.14.dpi_up <- Bd21.3.root.virus.MMMV.14.dpi [(Bd21.3.root.virus.MMMV.14.dpi $padj<0.05 & Bd21.3.root.virus.MMMV.14.dpi $log2FoldChange>= 1),]
nb_Bd21.3.root.virus.MMMV.14.dpi_up <- nrow(sigres_Bd21.3.root.virus.MMMV.14.dpi_up)
# 4216
write.table(x = sigres_Bd21.3.root.virus.MMMV.14.dpi_up,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.root.virus.MMMV.14.dpi_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

sigres_Bd21.3.root.virus.MMMV.14.dpi_down <- Bd21.3.root.virus.MMMV.14.dpi[(Bd21.3.root.virus.MMMV.14.dpi$padj<0.05 & Bd21.3.root.virus.MMMV.14.dpi$log2FoldChange<= -1),]
nb_Bd21.3.root.virus.MMMV.14.dpi_down <- nrow(sigres_Bd21.3.root.virus.MMMV.14.dpi_down)
#4788
write.table(x = sigres_Bd21.3.root.virus.MMMV.14.dpi_down,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.root.virus.MMMV.14.dpi_down.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
#DEG
sigres_Bd21.3.root.virus.MMMV.14.dpi <- Bd21.3.root.virus.MMMV.14.dpi[(Bd21.3.root.virus.MMMV.14.dpi$padj<0.05 & abs(Bd21.3.root.virus.MMMV.14.dpi$log2FoldChange) >= 1 ),]
nb_sigres_Bd21.3.root.virus.MMMV.14.dpi_DEG <- nrow(sigres_Bd21.3.root.virus.MMMV.14.dpi)
#9004
write.table(x = sigres_Bd21.3.root.virus.MMMV.14.dpi,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.root.virus.MMMV.14.dpi_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

#DEG.strict
sigres_Bd21.3.root.virus.MMMV.14.dpi.strict <- Bd21.3.root.virus.MMMV.14.dpi[(Bd21.3.root.virus.MMMV.14.dpi$padj<0.00001 & abs(Bd21.3.root.virus.MMMV.14.dpi$log2FoldChange) >= 5 ),]
nb_sigres_Bd21.3.root.virus.MMMV.14.dpi_DEG.strict <- nrow(sigres_Bd21.3.root.virus.MMMV.14.dpi.strict)
#320
write.table(x = sigres_Bd21.3.root.virus.MMMV.14.dpi.strict,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.root.virus.MMMV.14.dpi_DEG.strict.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
intersect(intersect(rownames(sigres_Bd21.3.root.virus.MMMV.3.dpi.strict),rownames(sigres_Bd21.3.root.virus.MMMV.7.dpi.strict)),rownames(sigres_Bd21.3.root.virus.MMMV.14.dpi.strict))

##WSMV
Bd21.3.root.virus.WSMV.3.dpi <- Bd21.3.root.virus.WSMV.3.dpi[!is.na(Bd21.3.root.virus.WSMV.3.dpi$padj),]
sigres_Bd21.3.root.virus.WSMV.3.dpi_up <- Bd21.3.root.virus.WSMV.3.dpi [(Bd21.3.root.virus.WSMV.3.dpi $padj<0.05 & Bd21.3.root.virus.WSMV.3.dpi $log2FoldChange>= 1),]
nb_Bd21.3.root.virus.WSMV.3.dpi_up <- nrow(sigres_Bd21.3.root.virus.WSMV.3.dpi_up)
# 1540
write.table(x = sigres_Bd21.3.root.virus.WSMV.3.dpi_up,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.root.virus.WSMV.3.dpi_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
sigres_Bd21.3.root.virus.WSMV.3.dpi_down <- Bd21.3.root.virus.WSMV.3.dpi[(Bd21.3.root.virus.WSMV.3.dpi$padj<0.05 & Bd21.3.root.virus.WSMV.3.dpi$log2FoldChange<= -1),]
nb_Bd21.3.root.virus.WSMV.3.dpi_down <- nrow(sigres_Bd21.3.root.virus.WSMV.3.dpi_down)
#1808
write.table(x = sigres_Bd21.3.root.virus.WSMV.3.dpi_down,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.root.virus.WSMV.3.dpi_down.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
##WSMV DEG
sigres_Bd21.3.root.virus.WSMV.3.dpi <- Bd21.3.root.virus.WSMV.3.dpi[(Bd21.3.root.virus.WSMV.3.dpi$padj<0.05 & abs(Bd21.3.root.virus.WSMV.3.dpi$log2FoldChange) >= 1 ),]
nb_sigres_Bd21.3.root.virus.WSMV.3.dpi_DEG <- nrow(sigres_Bd21.3.root.virus.WSMV.3.dpi)
#3348
write.table(x = sigres_Bd21.3.root.virus.WSMV.3.dpi,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.root.virus.WSMV.3.dpi_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
##WSMV DEG.strict
sigres_Bd21.3.root.virus.WSMV.3.dpi.strict <- Bd21.3.root.virus.WSMV.3.dpi[(Bd21.3.root.virus.WSMV.3.dpi$padj<0.00001 & abs(Bd21.3.root.virus.WSMV.3.dpi$log2FoldChange) >= 5 ),]
nb_sigres_Bd21.3.root.virus.WSMV.3.dpi_DEG.strict <- nrow(sigres_Bd21.3.root.virus.WSMV.3.dpi.strict)
#24
write.table(x = sigres_Bd21.3.root.virus.WSMV.3.dpi.strict,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.root.virus.WSMV.3.dpi_DEG.strict.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
sigres_Bd21.3.root.virus.WSMV.3.dpi.strict_up <- Bd21.3.root.virus.WSMV.3.dpi[(Bd21.3.root.virus.WSMV.3.dpi$padj<0.00001 & Bd21.3.root.virus.WSMV.3.dpi$log2FoldChange >= 5 ),]
write.table(x = sigres_Bd21.3.root.virus.WSMV.3.dpi.strict_up,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.root.virus.WSMV.3.dpi_DEG.strict_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
#Bd21.3.root.virus.WSMV.7.dpi
Bd21.3.root.virus.WSMV.7.dpi <- Bd21.3.root.virus.WSMV.7.dpi[!is.na(Bd21.3.root.virus.WSMV.7.dpi$padj),]
sigres_Bd21.3.root.virus.WSMV.7.dpi_up <- Bd21.3.root.virus.WSMV.7.dpi [(Bd21.3.root.virus.WSMV.7.dpi $padj<0.05 & Bd21.3.root.virus.WSMV.7.dpi $log2FoldChange>= 1),]
nb_Bd21.3.root.virus.WSMV.7.dpi_up <- nrow(sigres_Bd21.3.root.virus.WSMV.7.dpi_up)
# 1707
write.table(x = sigres_Bd21.3.root.virus.WSMV.7.dpi_up,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.root.virus.WSMV.7.dpi_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
sigres_Bd21.3.root.virus.WSMV.7.dpi_down <- Bd21.3.root.virus.WSMV.7.dpi[(Bd21.3.root.virus.WSMV.7.dpi$padj<0.05 & Bd21.3.root.virus.WSMV.7.dpi$log2FoldChange<= -1),]
nb_Bd21.3.root.virus.WSMV.7.dpi_down <- nrow(sigres_Bd21.3.root.virus.WSMV.7.dpi_down)
#2091
write.table(x = sigres_Bd21.3.root.virus.WSMV.7.dpi_down,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.root.virus.WSMV.7.dpi_down.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
#DEG
sigres_Bd21.3.root.virus.WSMV.7.dpi <- Bd21.3.root.virus.WSMV.7.dpi[(Bd21.3.root.virus.WSMV.7.dpi$padj<0.05 & abs(Bd21.3.root.virus.WSMV.7.dpi$log2FoldChange) >= 1 ),]
nb_sigres_Bd21.3.root.virus.WSMV.7.dpi_DEG <- nrow(sigres_Bd21.3.root.virus.WSMV.7.dpi)
#3798
write.table(x = sigres_Bd21.3.root.virus.WSMV.7.dpi,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_sigres_Bd21.3.root.virus.WSMV.7.dpi_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
#DEG.strict
sigres_Bd21.3.root.virus.WSMV.7.dpi.strict <- Bd21.3.root.virus.WSMV.7.dpi[(Bd21.3.root.virus.WSMV.7.dpi$padj<0.00001 & abs(Bd21.3.root.virus.WSMV.7.dpi$log2FoldChange) >= 5 ),]
nb_sigres_Bd21.3.root.virus.WSMV.7.dpi_DEG.strict <- nrow(sigres_Bd21.3.root.virus.WSMV.7.dpi.strict)
#48
write.table(x = sigres_Bd21.3.root.virus.WSMV.7.dpi.strict,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_sigres_Bd21.3.root.virus.WSMV.7.dpi_DEG.strict.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

sigres_Bd21.3.root.virus.WSMV.7.dpi.strict_up <- Bd21.3.root.virus.WSMV.7.dpi[(Bd21.3.root.virus.WSMV.7.dpi$padj<0.00001 & Bd21.3.root.virus.WSMV.7.dpi$log2FoldChange >= 5 ),]
write.table(x = sigres_Bd21.3.root.virus.WSMV.7.dpi.strict_up,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_sigres_Bd21.3.root.virus.WSMV.7.dpi_DEG.strict_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

#Bd21.3.root.virus.WSMV.14.dpi
Bd21.3.root.virus.WSMV.14.dpi <- Bd21.3.root.virus.WSMV.14.dpi[!is.na(Bd21.3.root.virus.WSMV.14.dpi$padj),]
sigres_Bd21.3.root.virus.WSMV.14.dpi_up <- Bd21.3.root.virus.WSMV.14.dpi [(Bd21.3.root.virus.WSMV.14.dpi $padj<0.05 & Bd21.3.root.virus.WSMV.14.dpi $log2FoldChange>= 1),]
nb_Bd21.3.root.virus.WSMV.14.dpi_up <- nrow(sigres_Bd21.3.root.virus.WSMV.14.dpi_up)
# 1380
write.table(x = sigres_Bd21.3.root.virus.WSMV.14.dpi_up,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.root.virus.WSMV.14.dpi_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
sigres_Bd21.3.root.virus.WSMV.14.dpi_down <- Bd21.3.root.virus.WSMV.14.dpi[(Bd21.3.root.virus.WSMV.14.dpi$padj<0.05 & Bd21.3.root.virus.WSMV.14.dpi$log2FoldChange<= -1),]
nb_Bd21.3.root.virus.WSMV.14.dpi_down <- nrow(sigres_Bd21.3.root.virus.WSMV.14.dpi_down)
#540
write.table(x = sigres_Bd21.3.root.virus.WSMV.14.dpi_down,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.root.virus.WSMV.14.dpi_down.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
#DEG
sigres_Bd21.3.root.virus.WSMV.14.dpi <- Bd21.3.root.virus.WSMV.14.dpi[(Bd21.3.root.virus.WSMV.14.dpi$padj<0.05 & abs(Bd21.3.root.virus.WSMV.14.dpi$log2FoldChange) >= 1 ),]
nb_sigres_Bd21.3.root.virus.WSMV.14.dpi_DEG <- nrow(sigres_Bd21.3.root.virus.WSMV.14.dpi)
#1920
write.table(x = sigres_Bd21.3.root.virus.WSMV.14.dpi,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.root.virus.WSMV.14.dpi_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
#DEG.strict
sigres_Bd21.3.root.virus.WSMV.14.dpi.strict <- Bd21.3.root.virus.WSMV.14.dpi[(Bd21.3.root.virus.WSMV.14.dpi$padj<0.00001 & abs(Bd21.3.root.virus.WSMV.14.dpi$log2FoldChange) >= 5 ),]
nb_sigres_Bd21.3.root.virus.WSMV.14.dpi_DEG.strict <- nrow(sigres_Bd21.3.root.virus.WSMV.14.dpi.strict)
#28
write.table(x = sigres_Bd21.3.root.virus.WSMV.14.dpi.strict,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.root.virus.WSMV.14.dpi_DEG.strict.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

sigres_Bd21.3.root.virus.WSMV.14.dpi.strict_up <- Bd21.3.root.virus.WSMV.14.dpi[(Bd21.3.root.virus.WSMV.14.dpi$padj<0.00001 & Bd21.3.root.virus.WSMV.14.dpi$log2FoldChange >= 5 ),]
write.table(x = sigres_Bd21.3.root.virus.WSMV.14.dpi.strict_up,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/sigres_Bd21.3.root.virus.WSMV.14.dpi.strict_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
intersect(intersect(rownames(sigres_Bd21.3.root.virus.WSMV.3.dpi.strict),rownames(sigres_Bd21.3.root.virus.WSMV.7.dpi.strict)),rownames(sigres_Bd21.3.root.virus.WSMV.14.dpi.strict))

#Plot the fungi treatment along the time!!!
##GGPLOt
library(ggplot2)
fungi_PCA <- data.frame(Categories=rep(c("Up", "Down"), each=3),
                        Time=rep(c("2dpi", "4dpi", "6dpi"),2),
                        Number=c(nb_Bd21.fungi.leaves.PCA.2.dpi_up, nb_Bd21.fungi.leaves.PCA.4.dpi_up, nb_Bd21.fungi.leaves.PCA.6.dpi_up,nb_Bd21.fungi.leaves.PCA.2.dpi_down, nb_Bd21.fungi.leaves.PCA.4.dpi_down, nb_Bd21.fungi.leaves.PCA.6.dpi_down))
head(fungi_PCA)
str(fungi_PCA)
fungi_PCA$Time <- factor(fungi_PCA$Time, levels = c("2dpi", "4dpi", "6dpi"))
par(mfrow = c(1,1))
par(mar=c(5,5,5,5))

pdf(file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/Bd21.fungi.leaves.PCA_DEG.pdf",width = 8,height = 8)
p <- ggplot(data=fungi_PCA, aes(x=Time, y=Number, fill=Categories)) +
  geom_bar(stat="identity", position=position_dodge())+ylim(0,3000)

heatplot <- p + labs(title="Fungi_PCA_leaves_Bd21", 
                     x="Time (days)", y = "The number of genes")+
  scale_fill_manual(values=c('black','lightgray')) + theme_bw() +theme(
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
    axis.text.x = element_text(color = "black", size = 20),
    axis.text.y = element_text(color = "black", size = 20),  
    axis.title.x = element_text(color = "black", size = 20),
    axis.title.y = element_text(color = "black", size = 20)
  )
print(heatplot)
dev.off()

##Bacterial_XT
library(ggplot2)
bac_XT <- data.frame(Categories=rep(c("Up", "Down"), each=2),
                     Time=rep(c("2dpi", "4dpi"),2),
                     Number=c(nb_Bd21.3.leaves.bacteria.XT.2.dpi_up, nb_Bd21.3.leaves.bacteria.XT.4.dpi_up, nb_Bd21.3.leaves.bacteria.XT.2.dpi_down, nb_Bd21.3.leaves.bacteria.XT.4.dpi_down))
head(bac_XT)
str(bac_XT)
bac_XT$Time <- factor(bac_XT$Time, levels = c("2dpi", "4dpi"))
par(mfrow = c(1,1))
par(mar=c(5,5,5,5))

pdf(file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/Bd21.3.leaves.bacteria.XT_DEG.pdf",width = 8,height = 8)
p <- ggplot(data=bac_XT, aes(x=Time, y=Number, fill=Categories)) +
  geom_bar(stat="identity", position=position_dodge())+ylim(0,3000)

heatplot <- p + labs(title="Bd21.3.leaves.bacteria.XT", 
                     x="Time (days)", y = "The number of genes")+
  scale_fill_manual(values=c('black','lightgray')) + theme_bw() +theme(
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
    axis.text.x = element_text(color = "black", size = 20),
    axis.text.y = element_text(color = "black", size = 20),  
    axis.title.x = element_text(color = "black", size = 20),
    axis.title.y = element_text(color = "black", size = 20)
  )

print(heatplot)
dev.off()
##leaves.virus.BSMV
##GGPLOt
library(ggplot2)
L_BSMV <- data.frame(Categories=rep(c("Up", "Down"), each=3),
                     Time=rep(c("3dpi", "7dpi", "14dpi"),2),
                     Number=c(nb_Bd21.3.leaves.virus.BSMV.3.dpi_up, nb_Bd21.3.leaves.virus.BSMV.7.dpi_up, nb_Bd21.3.leaves.virus.BSMV.14.dpi_up,nb_Bd21.3.leaves.virus.BSMV.3.dpi_down, nb_Bd21.3.leaves.virus.BSMV.7.dpi_down, nb_Bd21.3.leaves.virus.BSMV.14.dpi_down))
head(L_BSMV)
str(L_BSMV)
L_BSMV$Time <- factor(L_BSMV$Time, levels = c("3dpi", "7dpi", "14dpi"))
par(mfrow = c(1,1))
par(mar=c(5,5,5,5))

pdf(file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/Bd21.3.leaves.virus.BSMV_DEG.pdf",width = 8,height = 8)
p <- ggplot(data=L_BSMV, aes(x=Time, y=Number, fill=Categories)) +
  geom_bar(stat="identity", position=position_dodge())+ylim(0,5000)

heatplot <- p + labs(title="Bd21.3.leaves.virus.BSMV", 
                     x="Time (days)", y = "The number of genes")+
  scale_fill_manual(values=c('black','lightgray')) + theme_bw() +theme(
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
    axis.text.x = element_text(color = "black", size = 20),
    axis.text.y = element_text(color = "black", size = 20),  
    axis.title.x = element_text(color = "black", size = 20),
    axis.title.y = element_text(color = "black", size = 20)
  )

print(heatplot)
dev.off()
##leaves.virus.MMMV
##GGPLOt
library(ggplot2)
L_MMMV <- data.frame(Categories=rep(c("Up", "Down"), each=3),
                     Time=rep(c("3dpi", "7dpi", "14dpi"),2),
                     Number=c(nb_Bd21.3.leaves.virus.MMMV.3.dpi_up, nb_Bd21.3.leaves.virus.MMMV.7.dpi_up, nb_Bd21.3.leaves.virus.MMMV.14.dpi_up,nb_Bd21.3.leaves.virus.MMMV.3.dpi_down, nb_Bd21.3.leaves.virus.MMMV.7.dpi_down, nb_Bd21.3.leaves.virus.MMMV.14.dpi_down))
head(L_MMMV)
str(L_MMMV)
L_MMMV$Time <- factor(L_MMMV$Time, levels = c("3dpi", "7dpi", "14dpi"))
par(mfrow = c(1,1))
par(mar=c(5,5,5,5))

pdf(file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/Bd21.3.leaves.virus.MMMV_DEG.pdf",width = 8,height = 8)
p <- ggplot(data=L_MMMV, aes(x=Time, y=Number, fill=Categories)) +
  geom_bar(stat="identity", position=position_dodge())+ylim(0,5000)

heatplot <- p + labs(title="Bd21.3.leaves.virus.MMMV", 
                     x="Time (days)", y = "The number of genes")+
  scale_fill_manual(values=c('black','lightgray')) + theme_bw() +theme(
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
    axis.text.x = element_text(color = "black", size = 20),
    axis.text.y = element_text(color = "black", size = 20),  
    axis.title.x = element_text(color = "black", size = 20),
    axis.title.y = element_text(color = "black", size = 20)
  )

print(heatplot)
dev.off()
##leaves.virus.WSMV
##GGPLOt
library(ggplot2)
L_WSMV <- data.frame(Categories=rep(c("Up", "Down"), each=3),
                     Time=rep(c("3dpi", "7dpi", "14dpi"),2),
                     Number=c(nb_Bd21.3.leaves.virus.WSMV.3.dpi_up, nb_Bd21.3.leaves.virus.WSMV.7.dpi_up, nb_Bd21.3.leaves.virus.WSMV.14.dpi_up,nb_Bd21.3.leaves.virus.WSMV.3.dpi_down, nb_Bd21.3.leaves.virus.WSMV.7.dpi_down, nb_Bd21.3.leaves.virus.WSMV.14.dpi_down))
head(L_WSMV)
str(L_WSMV)
L_WSMV$Time <- factor(L_WSMV$Time, levels = c("3dpi", "7dpi", "14dpi"))
par(mfrow = c(1,1))
par(mar=c(5,5,5,5))

pdf(file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/Bd21.3.leaves.virus.WSMV_DEG.pdf",width = 8,height = 8)
p <- ggplot(data=L_WSMV, aes(x=Time, y=Number, fill=Categories)) +
  geom_bar(stat="identity", position=position_dodge())+ylim(0,5000)

heatplot <- p + labs(title="Bd21.3.leaves.virus.WSMV", 
                     x="Time (days)", y = "The number of genes")+
  scale_fill_manual(values=c('black','lightgray')) + theme_bw() +theme(
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
    axis.text.x = element_text(color = "black", size = 20),
    axis.text.y = element_text(color = "black", size = 20),  
    axis.title.x = element_text(color = "black", size = 20),
    axis.title.y = element_text(color = "black", size = 20)
  )

print(heatplot)
dev.off()
##root.virus.BSMV
##GGPLOt
library(ggplot2)
R_BSMV <- data.frame(Categories=rep(c("Up", "Down"), each=3),
                     Time=rep(c("3dpi", "7dpi", "14dpi"),2),
                     Number=c(nb_Bd21.3.root.virus.BSMV.3.dpi_up, nb_Bd21.3.root.virus.BSMV.7.dpi_up, nb_Bd21.3.root.virus.BSMV.14.dpi_up,nb_Bd21.3.root.virus.BSMV.3.dpi_down, nb_Bd21.3.root.virus.BSMV.7.dpi_down, nb_Bd21.3.root.virus.BSMV.14.dpi_down))
head(R_BSMV)
str(R_BSMV)
R_BSMV$Time <- factor(R_BSMV$Time, levels = c("3dpi", "7dpi", "14dpi"))
par(mfrow = c(1,1))
par(mar=c(5,5,5,5))

pdf(file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/Bd21.3.root.virus.BSMV_DEG.pdf",width = 8,height = 8)
p <- ggplot(data=R_BSMV, aes(x=Time, y=Number, fill=Categories)) +
  geom_bar(stat="identity", position=position_dodge())+ylim(0,5000)

heatplot <- p + labs(title="Bd21.3.root.virus.BSMV", 
                     x="Time (days)", y = "The number of genes")+
  scale_fill_manual(values=c('black','lightgray')) + theme_bw() +theme(
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
    axis.text.x = element_text(color = "black", size = 20),
    axis.text.y = element_text(color = "black", size = 20),  
    axis.title.x = element_text(color = "black", size = 20),
    axis.title.y = element_text(color = "black", size = 20)
  )

print(heatplot)
dev.off()


##root.virus.MMMV
##GGPLOt
library(ggplot2)
R_MMMV <- data.frame(Categories=rep(c("Up", "Down"), each=3),
                     Time=rep(c("3dpi", "7dpi", "14dpi"),2),
                     Number=c(nb_Bd21.3.root.virus.MMMV.3.dpi_up, nb_Bd21.3.root.virus.MMMV.7.dpi_up, nb_Bd21.3.root.virus.MMMV.14.dpi_up,nb_Bd21.3.root.virus.MMMV.3.dpi_down, nb_Bd21.3.root.virus.MMMV.7.dpi_down, nb_Bd21.3.root.virus.MMMV.14.dpi_down))
head(R_MMMV)
str(R_MMMV)
R_MMMV$Time <- factor(R_MMMV$Time, levels = c("3dpi", "7dpi", "14dpi"))
par(mfrow = c(1,1))
par(mar=c(5,5,5,5))

pdf(file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/Bd21.3.root.virus.MMMV_DEG.pdf",width = 8,height = 8)
p <- ggplot(data=R_MMMV, aes(x=Time, y=Number, fill=Categories)) +
  geom_bar(stat="identity", position=position_dodge())+ylim(0,5000)

heatplot <- p + labs(title="Bd21.3.root.virus.MMMV", 
                     x="Time (days)", y = "The number of genes")+
  scale_fill_manual(values=c('black','lightgray')) + theme_bw() +theme(
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
    axis.text.x = element_text(color = "black", size = 20),
    axis.text.y = element_text(color = "black", size = 20),  
    axis.title.x = element_text(color = "black", size = 20),
    axis.title.y = element_text(color = "black", size = 20)
  )

print(heatplot)
dev.off()



##root.virus.WSMV
##GGPLOt
library(ggplot2)
R_WSMV <- data.frame(Categories=rep(c("Up", "Down"), each=3),
                     Time=rep(c("3dpi", "7dpi", "14dpi"),2),
                     Number=c(nb_Bd21.3.root.virus.WSMV.3.dpi_up, nb_Bd21.3.root.virus.WSMV.7.dpi_up, nb_Bd21.3.root.virus.WSMV.14.dpi_up,nb_Bd21.3.root.virus.WSMV.3.dpi_down, nb_Bd21.3.root.virus.WSMV.7.dpi_down, nb_Bd21.3.root.virus.WSMV.14.dpi_down))
head(R_WSMV)
str(R_WSMV)
R_WSMV$Time <- factor(R_WSMV$Time, levels = c("3dpi", "7dpi", "14dpi"))
par(mfrow = c(1,1))
par(mar=c(5,5,5,5))

pdf(file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/Bd21.3.root.virus.WSMV_DEG.pdf",width = 8,height = 8)
p <- ggplot(data=R_WSMV, aes(x=Time, y=Number, fill=Categories)) +
  geom_bar(stat="identity", position=position_dodge())+ylim(0,5000)

heatplot <- p + labs(title="Bd21.3.root.virus.WSMV", 
                     x="Time (days)", y = "The number of genes")+
  scale_fill_manual(values=c('black','lightgray')) + theme_bw() +theme(
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
    axis.text.x = element_text(color = "black", size = 20),
    axis.text.y = element_text(color = "black", size = 20),  
    axis.title.x = element_text(color = "black", size = 20),
    axis.title.y = element_text(color = "black", size = 20)
  )

print(heatplot)
dev.off()
#Fungi:
fungi.leaves.PCA_gene_list <- union(union(rownames(sigres_Bd21.fungi.leaves.PCA.2.dpi),rownames(sigres_Bd21.fungi.leaves.PCA.4.dpi)),rownames(sigres_Bd21.fungi.leaves.PCA.6.dpi))
length(fungi.leaves.PCA_gene_list)
#4375

#Bacteria:
bacteria.leaves.XT_gene_list <- union(rownames(sigres_Bd21.3.leaves.bacteria.XT.2.dpi),rownames(sigres_Bd21.3.leaves.bacteria.XT.4.dpi))
length(bacteria.leaves.XT_gene_list)
#3809
#virus-leaves-BSMV
Bd21.3.leaves.virus.BSMV_gene_list <- union(union(rownames(sigres_Bd21.3.leaves.virus.BSMV.3.dpi),rownames(sigres_Bd21.3.leaves.virus.BSMV.7.dpi)),rownames(sigres_Bd21.3.leaves.virus.BSMV.14.dpi))
length(Bd21.3.leaves.virus.BSMV_gene_list)
# 5526
Bd21.3.leaves.virus.MMMV_gene_list <- union(union(rownames(sigres_Bd21.3.leaves.virus.MMMV.3.dpi),rownames(sigres_Bd21.3.leaves.virus.MMMV.7.dpi)),rownames(sigres_Bd21.3.leaves.virus.MMMV.14.dpi))
length(Bd21.3.leaves.virus.MMMV_gene_list)
#8924
Bd21.3.leaves.virus.WSMV_gene_list <- union(union(rownames(sigres_Bd21.3.leaves.virus.WSMV.3.dpi),rownames(sigres_Bd21.3.leaves.virus.WSMV.7.dpi)),rownames(sigres_Bd21.3.leaves.virus.WSMV.14.dpi))
length(Bd21.3.leaves.virus.WSMV_gene_list)
#3673
#root-BSMV
Bd21.3.root.virus.BSMV_gene_list <- union(union(rownames(sigres_Bd21.3.root.virus.BSMV.3.dpi),rownames(sigres_Bd21.3.root.virus.BSMV.7.dpi)),rownames(sigres_Bd21.3.root.virus.BSMV.14.dpi))
length(Bd21.3.root.virus.BSMV_gene_list)
#6118
#root-MMMV
Bd21.3.root.virus.MMMV_gene_list <- union(union(rownames(sigres_Bd21.3.root.virus.MMMV.3.dpi),rownames(sigres_Bd21.3.root.virus.MMMV.7.dpi)),rownames(sigres_Bd21.3.root.virus.MMMV.14.dpi))
length(Bd21.3.root.virus.MMMV_gene_list)
#10833
#root-WSMV
Bd21.3.root.virus.WSMV_gene_list <- union(union(rownames(sigres_Bd21.3.root.virus.WSMV.3.dpi),rownames(sigres_Bd21.3.root.virus.WSMV.7.dpi)),rownames(sigres_Bd21.3.root.virus.WSMV.14.dpi))
length(Bd21.3.root.virus.WSMV_gene_list)
#6322
#top gene list:
top_gene_list <- union(union(union(union(union(union(union(fungi.leaves.PCA_gene_list,bacteria.leaves.XT_gene_list),
                                                     Bd21.3.leaves.virus.BSMV_gene_list),
                                               Bd21.3.leaves.virus.MMMV_gene_list),
                                         Bd21.3.leaves.virus.WSMV_gene_list),
                                   Bd21.3.root.virus.BSMV_gene_list),
                             Bd21.3.root.virus.MMMV_gene_list),
                       Bd21.3.root.virus.WSMV_gene_list)

length(top_gene_list)
#17204
write.table(x = top_gene_list,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/top_DEG_17204_norma_pythogen_gene_list",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

write.table(x = fungi.leaves.PCA_gene_list,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/fungi.leaves.PCA_combined_gene_list",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
write.table(x = bacteria.leaves.XT_gene_list,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/bacteria.leaves.XT_combined_gene_list",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
write.table(x = Bd21.3.leaves.virus.BSMV_gene_list,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/Bd21.3.leaves.virus.BSMV_combined_gene_list",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

write.table(x = Bd21.3.leaves.virus.MMMV_gene_list,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/Bd21.3.leaves.virus.MMMV_combined_gene_list",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

write.table(x = Bd21.3.leaves.virus.WSMV_gene_list,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/Bd21.3.leaves.virus.WSMV_combined_gene_list",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

write.table(x = Bd21.3.root.virus.BSMV_gene_list,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/Bd21.3.root.virus.BSMV_combined_gene_list",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

write.table(x = Bd21.3.root.virus.MMMV_gene_list,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/Bd21.3.root.virus.MMMV_combined_gene_list",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

write.table(x = Bd21.3.root.virus.WSMV_gene_list,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/Bd21.3.root.virus.WSMV_combined_gene_list",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
#plot the genes with logfolderchanges >1
####reorder the columns
pdf(file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/heatmap_17204_DGE_pathogene_endophyte_noUB_del1rep_reordered.pdf",width = 15,height = 15)

par(mar=c(1,1,1,1))

palette <- colorRampPalette(c("red","white","blue"))(256)
reorder_rld <- assay(rld)[,c(4,5,6,1,2,3,10,11,12,7,8,9,16,17,18,13,14,15,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,67,100,101,102)]
head(reorder_rld)
heatmap.2( reorder_rld[row.names(assay(rld)) %in% top_gene_list,],  scale="row", Colv=FALSE, 
           trace="none", dendrogram='row',
           col = palette,
           lhei = c(0.8,7),
           ColSideColors = c( "Bd21.fungi.Mock"="#cc561e","Bd21.fungi.PCA" ="#aa2b1d", "Bd21-3.bacteria.Mock" = "#ef8d32","Bd21-3.bacteria.XT" = "#beca5c", "Bd21-3.virus.Mock"="#75cfb8","Bd21-3.virus.BSMV"="#bbdfc8","Bd21-3.virus.MMMV"="#f0e5d8","Bd21-3.virus.WSMV"="#ffc478")[
             colData(rld)$treatment] )
dev.off()

####reorder the columns
pdf(file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/only_pythogen/heatmap_17204_DGE_pathogene_endophyte_noUB_del1rep_reordered#.pdf",width = 15,height = 15)

par(mar=c(1,1,1,1))
reorder_rld <- assay(rld)[,c(4,5,6,10,11,12,16,17,18,1,2,3,7,8,9,13,14,15,19,20,21,25,26,27,22,23,24,28,29,30,31,32,33,43,44,45,55,56,57,34,35,36,46,47,48,58,59,60,37,38,39,49,50,51,61,62,63,40,41,42,52,53,54,64,65,66,68,69,70,80,81,82,92,93,94,71,72,73,83,84,85,95,96,97,74,75,76,86,87,88,98,99,67,77,78,79,89,90,91,100,101,102)]

head(reorder_rld)
head(reorder_rld)
heatmap.2( reorder_rld[row.names(assay(rld)) %in% top_gene_list,],  scale="row", Colv=FALSE, 
           trace="none", dendrogram='row',
           col = palette,
           lhei = c(0.8,7),
           ColSideColors = c( "Bd21.fungi.Mock"="#cc561e","Bd21.fungi.PCA" ="#aa2b1d", "Bd21-3.bacteria.Mock" = "#ef8d32","Bd21-3.bacteria.XT" = "#beca5c", "Bd21-3.virus.Mock"="#75cfb8","Bd21-3.virus.BSMV"="#bbdfc8","Bd21-3.virus.MMMV"="#f0e5d8","Bd21-3.virus.WSMV"="#ffc478")[
             colData(rld)$treatment] )
dev.off()
###
###

