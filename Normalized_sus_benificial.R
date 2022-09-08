library(ggplot2)
library(edgeR)
library(DESeq2)
library(tidyr)
library(ggrepel)
library(EnhancedVolcano)

#source("https://bioconductor.org/biocLite.R")
#biocLite(c("BatchQC"))
#library(BatchQC)

counts <- read.delim("/global/projectb/scratch/llei2019/CSP_Kranthi/sus_beni/final_susceptible_benificial_count_formal.txt", row.names = 1)
head(counts)
nrow(counts)
#[1] 32439
metaData <- read.delim("/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/meta_data_pathogene_endophyte_noUB##.txt",header = T)
head(metaData)
nrow(metaData)
metaData$header <- paste0(metaData$accessions, ".", metaData$experiment, ".", metaData$tissue, ".", metaData$treatment, ".", metaData$timepoint, ".", metaData$replicate)
head(metaData)
#140
#32439*140
#extract the column I needed!!!
col.num <- which( colnames(counts) %in% metaData$code)
row.num <- which( metaData$code %in% colnames(counts) )
counts <- counts[, col.num]
mata <- metaData[row.num,]
ncol(counts)
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
#3      3      4      5      6      7      8     10     12     14     16     19     21 
#39%    40%    41%    42%    43%    44%    45%    46%    47%    48%    49%    50%    51% 
#25     28     32     37     42     47     53     59     66     74     82     91    101 
#52%    53%    54%    55%    56%    57%    58%    59%    60%    61%    62%    63%    64% 
#111    123    135    147    161    176    192    209    227    246    266    288    312 
#65%    66%    67%    68%    69%    70%    71%    72%    73%    74%    75%    76%    77% 
#336    362    390    420    452    486    521    560    600    644    691    742    796 
#78%    79%    80%    81%    82%    83%    84%    85%    86%    87%    88%    89%    90% 
#854    918    986   1060   1142   1231   1331   1440   1562   1698   1853   2030   2234 
#91%    92%    93%    94%    95%    96%    97%    98%    99%   100% 
#  2473   2758   3108   3554   4134   4940   6139   8198  13139 808320 
count.cutoff = 3 # 60% of the data were above 3
bioreplicates.cutoff = 3# 73% of the data were above 5
## RAW COUNTS ##
keep <- rowSums(counts >= count.cutoff) >= bioreplicates.cutoff
counts <- counts[keep, ]
#?rowSums
nrow(counts)
#[1] 30836
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
#0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 
#14%          15%          16%          17%          18%          19%          20% 
#0.000000e+00 0.000000e+00 2.088135e-02 2.704418e-02 3.767557e-02 4.875274e-02 6.062461e-02 
#21%          22%          23%          24%          25%          26%          27% 
#7.446943e-02 9.259164e-02 1.139515e-01 1.377538e-01 1.677677e-01 1.972760e-01 2.324466e-01 
#28%          29%          30%          31%          32%          33%          34% 
#2.747118e-01 3.230489e-01 3.786185e-01 4.430295e-01 5.116571e-01 5.902406e-01 6.787061e-01 
#35%          36%          37%          38%          39%          40%          41% 
#7.862126e-01 8.976637e-01 1.021254e+00 1.157722e+00 1.312024e+00 1.478077e+00 1.659300e+00 
#42%          43%          44%          45%          46%          47%          48% 
#1.858440e+00 2.072302e+00 2.304924e+00 2.561939e+00 2.838302e+00 3.134015e+00 3.449701e+00 
#49%          50%          51%          52%          53%          54%          55% 
#3.787995e+00 4.148342e+00 4.540261e+00 4.960696e+00 5.408520e+00 5.882479e+00 6.390628e+00 
#56%          57%          58%          59%          60%          61%          62% 
#6.925657e+00 7.496779e+00 8.098892e+00 8.739935e+00 9.411569e+00 1.011553e+01 1.086065e+01 
#63%          64%          65%          66%          67%          68%          69% 
#1.163764e+01 1.246521e+01 1.333265e+01 1.425319e+01 1.522584e+01 1.625773e+01 1.735134e+01 
#70%          71%          72%          73%          74%          75%          76% 
#1.851351e+01 1.972625e+01 2.102325e+01 2.240940e+01 2.388827e+01 2.544864e+01 2.710944e+01 
#77%          78%          79%          80%          81%          82%          83% 
#2.888123e+01 3.077392e+01 3.282549e+01 3.500796e+01 3.737435e+01 3.992891e+01 4.271011e+01 
#84%          85%          86%          87%          88%          89%          90% 
#4.574081e+01 4.908396e+01 5.278115e+01 5.688596e+01 6.150794e+01 6.675811e+01 7.292201e+01 
#91%          92%          93%          94%          95%          96%          97% 
#8.008590e+01 8.861300e+01 9.910958e+01 1.123239e+02 1.295806e+02 1.533872e+02 1.889171e+02 
#98%          99%         100% 
#2.505080e+02 3.949451e+02 2.800450e+04 

count.cutoff = 1 #63% of the data is above 1
bioreplicates.cutoff = 3
keep <- rowSums(normalized.counts >= count.cutoff) >= bioreplicates.cutoff
#abiotic_shoot_count(keep)
normalized.counts <- normalized.counts[keep, ]
nrow(normalized.counts)
#[1] 24867
nrow(counts)
#[1] 30836
#keep the datapoint pass the threshold
counts <- counts[row.names(normalized.counts),]
nrow(counts)
#24867
grep("Bradi4g39317", rownames(counts))#check if this NBSLRR gene members got filtered from this link: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3368180/pdf/CFG2012-418208.pdf
#[1] 21898
#it suggested that the filtering seems good!
#21794/36927=0.6318141 only keep around 63% genes


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
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/normalized_VST_pathogene_endophyte_noUB.txt",
            sep = "\t",
            eol = "\n",
            col.names = TRUE,
            row.names = TRUE
)

ddsTC <- DESeq(dds)
#heatmap:
library("pheatmap")
rld <- vst( ddsTC )
pdf(file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/PCA_pathogene_endophyte_noUB#.pdf",width = 15,height = 15)
par(mar=c(1,1,1,1))
plotPCA(rld, intgroup="group")
dev.off()
###
pcaData <- plotPCA(rld, intgroup = c(  "tissue", "treatment", "timepoint"),returnData = TRUE) # vsd and plotPCA are part of DESeq2 package, nothing with my example below. 
head(pcaData)
percentVar <- round(100 * attr(pcaData, "percentVar")) 
pdf(file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/PCA_pathogene_endophyte_noUB###.pdf",width = 20,height = 15)
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
###
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
head(pcaData)
percentVar <- round(100 * attr(pcaData, "percentVar")) 
pdf(file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/PCA_pathogene_endophyte_noUB_pc1_pc3#.pdf",width = 20,height = 15)
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
pdf(file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/heatmap_1000_DGE_pathogene_endophyte_noUB#.pdf",width = 15,height = 15)

par(mar=c(1,1,1,1))
#heatmap(assay(rld)[ topVarGenes, ],
#        labRow = "",
#        scale = "row")
palette <- colorRampPalette(c("red","white","blue"))(256)

heatmap.2( assay(rld)[ topVarGenes, ], scale="row", Colv=FALSE, 
           trace="none", dendrogram='row', 
           col = palette,
           lhei = c(0.8,7),
           #labCol = colData(rld)$group,
           ColSideColors = c( "Bd21.fungi.Mock"="#cc561e","Bd21.fungi.PCA" ="#aa2b1d", "Bd21-3.bacteria.Mock" = "#ef8d32","Bd21-3.bacteria.XT" = "#beca5c", "Bd21-3.virus.Mock"="#75cfb8","Bd21-3.virus.BSMV"="#bbdfc8","Bd21-3.virus.MMMV"="#f0e5d8","Bd21-3.virus.WSMV"="#ffc478","Bd21-3.Endophytes.Mock"="#822659","Bd21-3.Endophytes.SV"="#f8a1d1")[
             colData(rld)$treatment ] )
dev.off()

####So I need to drop out the leaves.MMMV.14.dpi.3 since it is a outlier!!!Rerun eveything!
counts <- read.delim("/global/projectb/scratch/llei2019/CSP_Kranthi/sus_beni/final_susceptible_benificial_count_formal.txt", row.names = 1)
head(counts)
nrow(counts)
#[1] 32439
metaData <- read.delim("/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/meta_data_pathogene_endophyte_noUB_del1rep.txt",header = T)
head(metaData)
nrow(metaData)
#139
metaData$header <- paste0(metaData$accessions, ".", metaData$experiment, ".", metaData$tissue, ".", metaData$treatment, ".", metaData$timepoint, ".", metaData$replicate)
head(metaData)
#140
#32439*140
#extract the column I needed!!!
col.num <- which( colnames(counts) %in% metaData$code)
row.num <- which( metaData$code %in% colnames(counts) )
counts <- counts[, col.num]
mata <- metaData[row.num,]
ncol(counts)
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
#0%       1%       2%       3%       4%       5%       6%       7%       8%       9% 
#  0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0 
#10%      11%      12%      13%      14%      15%      16%      17%      18%      19% 
#  0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0 
#20%      21%      22%      23%      24%      25%      26%      27%      28%      29% 
#  1.0      1.0      1.0      1.0      2.0      2.0      3.0      3.0      4.0      5.0 
#30%      31%      32%      33%      34%      35%      36%      37%      38%      39% 
#  6.0      7.0      8.0     10.0     12.0     14.0     16.0     19.0     22.0     25.0 
#40%      41%      42%      43%      44%      45%      46%      47%      48%      49% 
#  28.0     32.0     37.0     42.0     47.0     53.0     60.0     67.0     74.0     83.0 
#50%      51%      52%      53%      54%      55%      56%      57%      58%      59% 
#  92.0    101.0    112.0    123.0    135.0    148.0    162.0    177.0    193.0    209.0 
#60%      61%      62%      63%      64%      65%      66%      67%      68%      69% 
#  228.0    247.0    267.0    289.0    313.0    337.0    364.0    392.0    422.0    454.0 
#70%      71%      72%      73%      74%      75%      76%      77%      78%      79% 
#  487.0    523.0    561.0    602.0    646.0    693.0    744.0    799.0    857.0    921.0 
#80%      81%      82%      83%      84%      85%      86%      87%      88%      89% 
#  989.0   1064.0   1146.0   1236.0   1335.0   1445.0   1567.0   1704.0   1859.0   2037.0 
#90%      91%      92%      93%      94%      95%      96%      97%      98%      99% 
#  2242.0   2481.0   2767.0   3118.0   3566.0   4148.0   4956.0   6159.0   8225.6  13181.0 
#100% 
#808320.0 
bioreplicates.cutoff = 3# 74% of the data were above 5
## RAW COUNTS ##
keep <- rowSums(counts >= count.cutoff) >= bioreplicates.cutoff
counts <- counts[keep, ]
#?rowSums
nrow(counts)
#[1] 30832
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
#0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 
#14%          15%          16%          17%          18%          19%          20% 
#0.000000e+00 0.000000e+00 2.108002e-02 2.712440e-02 3.767557e-02 4.931901e-02 6.118878e-02 
#21%          22%          23%          24%          25%          26%          27% 
#7.510165e-02 9.269302e-02 1.147764e-01 1.382971e-01 1.684764e-01 1.987786e-01 2.329776e-01 
#28%          29%          30%          31%          32%          33%          34% 
#2.757244e-01 3.240707e-01 3.797217e-01 4.447657e-01 5.136370e-01 5.922564e-01 6.827636e-01 
#35%          36%          37%          38%          39%          40%          41% 
#7.874655e-01 8.998487e-01 1.024934e+00 1.161327e+00 1.315663e+00 1.481655e+00 1.664951e+00 
#42%          43%          44%          45%          46%          47%          48% 
#1.862221e+00 2.078537e+00 2.309480e+00 2.566555e+00 2.842442e+00 3.139883e+00 3.454486e+00 
#49%          50%          51%          52%          53%          54%          55% 
#3.792254e+00 4.153472e+00 4.545903e+00 4.965553e+00 5.411758e+00 5.888226e+00 6.394538e+00 
#56%          57%          58%          59%          60%          61%          62% 
#6.932108e+00 7.499923e+00 8.104434e+00 8.743146e+00 9.414588e+00 1.011846e+01 1.086334e+01 
#63%          64%          65%          66%          67%          68%          69% 
#1.163993e+01 1.246726e+01 1.333375e+01 1.425523e+01 1.522743e+01 1.625792e+01 1.735092e+01 
#70%          71%          72%          73%          74%          75%          76% 
#1.851255e+01 1.972425e+01 2.102202e+01 2.240791e+01 2.388574e+01 2.544249e+01 2.710302e+01 
#77%          78%          79%          80%          81%          82%          83% 
#2.887498e+01 3.076448e+01 3.281108e+01 3.499392e+01 3.735620e+01 3.990809e+01 4.268920e+01 
#84%          85%          86%          87%          88%          89%          90% 
#4.571806e+01 4.906392e+01 5.275421e+01 5.685613e+01 6.147764e+01 6.672042e+01 7.287916e+01 
#91%          92%          93%          94%          95%          96%          97% 
#8.004742e+01 8.858212e+01 9.907465e+01 1.122865e+02 1.295659e+02 1.533661e+02 1.889152e+02 
#98%          99%         100% 
#2.505223e+02 3.950147e+02 2.800450e+04 

count.cutoff = 1 #63% of the data is above 1
bioreplicates.cutoff = 3
keep <- rowSums(normalized.counts >= count.cutoff) >= bioreplicates.cutoff
#abiotic_shoot_count(keep)
normalized.counts <- normalized.counts[keep, ]
nrow(normalized.counts)
#[1] 24852
nrow(counts)
#[1] 30836
#keep the datapoint pass the threshold
counts <- counts[row.names(normalized.counts),]
nrow(counts)
#24852
grep("Bradi4g39317", rownames(counts))#check if this NBSLRR gene members got filtered from this link: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3368180/pdf/CFG2012-418208.pdf
#[1] 21884
#it suggested that the filtering seems good!
#21884/32439=0.6746201 only keep around 63% genes

## VST instead of voom
ncol(counts)
#139
head(counts)
colnames(counts)
#creat a group 
design <- model.matrix(~0 + group,metaData)
dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, colData = metaData, design = design)
head(counts(dds))
expr.vst <- assay(DESeq2::vst(dds))
head(expr.vst)
#str(expr.vst)
expr.vst <- as.data.frame(expr.vst)
#This file is for WGCNA
write.table(x = expr.vst,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/normalized_VST_final_frezze_tissue_formal.txt",
            sep = "\t",
            eol = "\n",
            col.names = TRUE,
            row.names = TRUE
)

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
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/normalized_VST_pathogene_endophyte_noUB.txt",
            sep = "\t",
            eol = "\n",
            col.names = TRUE,
            row.names = TRUE
)

ddsTC <- DESeq(dds)
head(ddsTC)
#heatmap:
library("pheatmap")
rld <- vst( ddsTC )
head(rld)
pdf(file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/PCA_pathogene_endophyte_noUB_del1rep.pdf",width = 15,height = 15)
par(mar=c(1,1,1,1))
plotPCA(rld, intgroup="group")
dev.off()
###
pcaData <- plotPCA(rld, intgroup = c(  "tissue", "treatment", "timepoint"),returnData = TRUE) # vsd and plotPCA are part of DESeq2 package, nothing with my example below. 
head(pcaData)
percentVar <- round(100 * attr(pcaData, "percentVar")) 
pdf(file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/PCA_pathogene_endophyte_noUB_del1rep.pdf",width = 20,height = 15)
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
###
#Then we can plot the other PCs 
# The function is the basically the same as https://github.com/mikelove/DESeq2/blob/master/R/plots.R.
# Inspired by https://www.biostars.org/p/243695/. 
# I added two variables, pp1 and pp2 to let user chose which principle to plot. pp1 and pp2 only take integer values.
#Test PC1 vs PC3
plotPCA.DESeqTransform = function(pp1=2, pp2=3, object, intgroup="condition", ntop=500, returnData=FALSE)
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
  
  ggplot(data=d, aes_string(x="PC2", y="PC3", color="group")) + geom_point(size=3) + 
    xlab(paste0("PC2: ",round(percentVar[1] * 100),"% variance")) +
    ylab(paste0("PC3: ",round(percentVar[2] * 100),"% variance"))
}

pcaData <- plotPCA.DESeqTransform(object=rld, intgroup = c(  "tissue", "treatment", "timepoint"),returnData = TRUE) # vsd and plotPCA are part of DESeq2 package, nothing with my example below. 
head(pcaData)
head(pcaData)
percentVar <- round(100 * attr(pcaData, "percentVar")) 
pdf(file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/PCA_pathogene_endophyte_noUB_pc2_pc3_del1rep.pdf",width = 20,height = 15)
par(mar=c(1,1,1,1))
ggplot(pcaData, aes(x = PC1, y = PC2, color = factor(treatment), shape = factor(tissue))) + 
  geom_point(size =6, aes(fill=factor(treatment))) + 
  geom_point(size =6) + 
  scale_shape_manual(values=c(21,22,23,24)) + 
  #scale_alpha_manual(values=c("Bd21"=0, "Bd21-3"=1)) + 
  xlab(paste0("PC2: ", percentVar[1], "% variance")) + 
  ylab(paste0("PC3: ", percentVar[2], "% variance")) + 
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
pdf(file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/heatmap_1000_DGE_pathogene_endophyte_noUB_del1rep.pdf",width = 15,height = 15)

par(mar=c(1,1,1,1))
#heatmap(assay(rld)[ topVarGenes, ],
#        labRow = "",
#        scale = "row")
palette <- colorRampPalette(c("red","white","blue"))(256)

heatmap.2( assay(rld)[ topVarGenes, ], scale="row", Colv=FALSE, 
           trace="none", dendrogram='row', 
           col = palette,
           lhei = c(0.8,7),
           #labCol = colData(rld)$group,
           ColSideColors = c( "Bd21.fungi.Mock"="#cc561e","Bd21.fungi.PCA" ="#aa2b1d", "Bd21-3.bacteria.Mock" = "#ef8d32","Bd21-3.bacteria.XT" = "#beca5c", "Bd21-3.virus.Mock"="#75cfb8","Bd21-3.virus.BSMV"="#bbdfc8","Bd21-3.virus.MMMV"="#f0e5d8","Bd21-3.virus.WSMV"="#ffc478","Bd21-3.Endophytes.Mock"="#822659","Bd21-3.Endophytes.SV"="#f8a1d1")[
             colData(rld)$treatment ] )
dev.off()
####reorder the columns
pdf(file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/heatmap_1000_DGE_pathogene_endophyte_noUB_del1rep_reordered.pdf",width = 15,height = 15)

par(mar=c(1,1,1,1))
#heatmap(assay(rld)[ topVarGenes, ],
#        labRow = "",
#        scale = "row")
palette <- colorRampPalette(c("red","white","blue"))(256)
reorder_rld <- assay(rld)[,c(4,5,6,1,2,3,10,11,12,7,8,9,16,17,18,13,14,15,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,67,101,102,68,103,104,105,106,113,114,115,107,108,109,116,117,118,110,111,112,119,120,121,122,123,124,131,132,133,125,126,127,134,135,136,128,129,130,137,138,139)]
head(reorder_rld)
heatmap.2( reorder_rld[ topVarGenes, ], scale="row", Colv=FALSE, 
           trace="none", dendrogram='row',
           col = palette,
           lhei = c(0.8,7),
           ColSideColors = c( "Bd21.fungi.Mock"="#cc561e","Bd21.fungi.PCA" ="#aa2b1d", "Bd21-3.bacteria.Mock" = "#ef8d32","Bd21-3.bacteria.XT" = "#beca5c", "Bd21-3.virus.Mock"="#75cfb8","Bd21-3.virus.BSMV"="#bbdfc8","Bd21-3.virus.MMMV"="#f0e5d8","Bd21-3.virus.WSMV"="#ffc478","Bd21-3.Endophytes.Mock"="#822659","Bd21-3.Endophytes.SV"="#f8a1d1")[
             colData(rld)$treatment] )
dev.off()
###
##cold versus freeze leaf
colData(ddsTC)
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
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/Bd21.fungi.leaves.PCA.2.dpi.pdf",
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
                pCutoff = 0.0000000001,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/Bd21.fungi.leaves.PCA.2.dpi.strict.pdf",
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
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/Bd21.fungi.leaves.PCA.4.dpi.pdf",
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
                pCutoff = 0.0000000001,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/Bd21.fungi.leaves.PCA.4.dpi.strict.pdf",
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
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/Bd21.fungi.leaves.PCA.6.dpi.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)
resultsNames(ddsTC)
#strict:
EnhancedVolcano(Bd21.fungi.leaves.PCA.6.dpi,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Bd21.fungi.leaves.PCA.6.dpi',
                xlim = c(-25, 25),
                ylim = c(0,50),
                FCcutoff=5,
                pCutoff = 0.0000000001,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/Bd21.fungi.leaves.PCA.6.dpi.strict.pdf",
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
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/Bd21.3.leaves.bacteria.XT.2.dpi.pdf",
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
                  pCutoff = 0.0000000001,
                  pointSize = 1.5,
                  colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/Bd21.3.leaves.bacteria.XT.2.dpi_strict.pdf",
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
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/Bd21.3.leaves.bacteria.XT.4.dpi.pdf",
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
                pCutoff = 0.0000000001,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/Bd21.3.leaves.bacteria.XT.4.dpi_strict.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)
resultsNames(ddsTC)

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
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/Bd21.3.leaves.virus.BSMV.3.dpi.pdf",
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
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/Bd21.3.leaves.virus.BSMV.7.dpi.pdf",
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
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/Bd21.3.leaves.virus.BSMV.14.dpi.pdf",
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
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/Bd21.3.leaves.virus.MMMV.3.dpi.pdf",
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
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/Bd21.3.leaves.virus.MMMV.7.dpi.pdf",
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
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/Bd21.3.leaves.virus.MMMV.14.dpi.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)

##
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
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/Bd21.3.leaves.virus.WSMV.3.dpi.pdf",
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
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/Bd21.3.leaves.virus.WSMV.7.dpi.pdf",
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
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/Bd21.3.leaves.virus.WSMV.14.dpi.pdf",
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
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/Bd21.3.root.virus.BSMV.3.dpi.pdf",
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
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/Bd21.3.root.virus.BSMV.7.dpi.pdf",
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
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/Bd21.3.root.virus.BSMV.14.dpi.pdf",
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
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/Bd21.3.root.virus.MMMV.3.dpi.pdf",
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
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/Bd21.3.root.virus.MMMV.7.dpi.pdf",
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
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/Bd21.3.root.virus.MMMV.14.dpi.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)

##
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
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/Bd21.3.root.virus.WSMV.3.dpi.pdf",
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
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/Bd21.3.root.virus.WSMV.7.dpi.pdf",
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
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/Bd21.3.root.virus.WSMV.14.dpi.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)
resultsNames(ddsTC)

#Endophytes-root
#Bd21.3.Endophytes.root.SV.3dpi
Bd21.3.Endophytes.root.SV.3dpi <- results(ddsTC, contrast=list("groupBd21.3.Endophytes.root.SV.3dpi", "groupBd21.3.Endophytes.root.Mock.3dpi"))
EnhancedVolcano(Bd21.3.Endophytes.root.SV.3dpi,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Bd21.3.Endophytes.root.SV.3dpi',
                xlim = c(-25, 25),
                ylim = c(0,50),
                FCcutoff=1,
                pCutoff = 0.05,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/Bd21.3.Endophytes.root.SV.3dpi.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)

Bd21.3.Endophytes.root.SV.6dpi <- results(ddsTC, contrast=list("groupBd21.3.Endophytes.root.SV.6dpi", "groupBd21.3.Endophytes.root.Mock.6dpi"))
EnhancedVolcano(Bd21.3.Endophytes.root.SV.6dpi,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Bd21.3.Endophytes.root.SV.6dpi',
                xlim = c(-25, 25),
                ylim = c(0,50),
                FCcutoff=1,
                pCutoff = 0.05,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/Bd21.3.Endophytes.root.SV.6dpi.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)

Bd21.3.Endophytes.root.SV.13dpi <- results(ddsTC, contrast=list("groupBd21.3.Endophytes.root.SV.13dpi", "groupBd21.3.Endophytes.root.Mock.13dpi"))
EnhancedVolcano(Bd21.3.Endophytes.root.SV.13dpi,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Bd21.3.Endophytes.root.SV.13dpi',
                xlim = c(-25, 25),
                ylim = c(0,50),
                FCcutoff=1,
                pCutoff = 0.05,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/Bd21.3.Endophytes.root.SV.13dpi.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)

#Bd21.3.shoot.virus.WSMV.3.dpi
Bd21.3.shoot.virus.WSMV.3.dpi <- results(ddsTC, contrast=list("groupBd21.3.virus.shoot.WSMV.3.dpi", "groupBd21.3.virus.shoot.Mock.3.dpi"))
EnhancedVolcano(Bd21.3.shoot.virus.WSMV.3.dpi,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Bd21.3.shoot.virus.WSMV.3.dpi',
                xlim = c(-25, 25),
                ylim = c(0,50),
                FCcutoff=1,
                pCutoff = 0.05,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/Bd21.3.shoot.virus.WSMV.3.dpi.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)

#Bd21.3.shoot.virus.WSMV.7.dpi
Bd21.3.shoot.virus.WSMV.7.dpi <- results(ddsTC, contrast=list("groupBd21.3.virus.shoot.WSMV.7.dpi", "groupBd21.3.virus.shoot.Mock.7.dpi"))
EnhancedVolcano(Bd21.3.shoot.virus.WSMV.7.dpi,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Bd21.3.shoot.virus.WSMV.7.dpi',
                xlim = c(-25, 25),
                ylim = c(0,50),
                FCcutoff=1,
                pCutoff = 0.05,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/Bd21.3.shoot.virus.WSMV.7.dpi.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)

#Bd21.3.shoot.virus.BSMV.14.dpi
Bd21.3.shoot.virus.WSMV.14.dpi <- results(ddsTC, contrast=list("groupBd21.3.virus.shoot.WSMV.14.dpi", "groupBd21.3.virus.shoot.Mock.14.dpi"))
EnhancedVolcano(Bd21.3.shoot.virus.WSMV.14.dpi,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Bd21.3.shoot.virus.WSMV.14.dpi',
                xlim = c(-25, 25),
                ylim = c(0,50),
                FCcutoff=1,
                pCutoff = 0.05,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/Bd21.3.shoot.virus.WSMV.14.dpi.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)
resultsNames(ddsTC)

#Endophytes-shoot
#Bd21.3.Endophytes.shoot.SV.3dpi
Bd21.3.Endophytes.shoot.SV.3dpi <- results(ddsTC, contrast=list("groupBd21.3.Endophytes.shoot.SV.3dpi", "groupBd21.3.Endophytes.shoot.Mock.3dpi"))
EnhancedVolcano(Bd21.3.Endophytes.shoot.SV.3dpi,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Bd21.3.Endophytes.shoot.SV.3dpi',
                xlim = c(-25, 25),
                ylim = c(0,50),
                FCcutoff=1,
                pCutoff = 0.05,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/Bd21.3.Endophytes.shoot.SV.3dpi.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)

Bd21.3.Endophytes.shoot.SV.6dpi <- results(ddsTC, contrast=list("groupBd21.3.Endophytes.shoot.SV.6dpi", "groupBd21.3.Endophytes.shoot.Mock.6dpi"))
EnhancedVolcano(Bd21.3.Endophytes.shoot.SV.6dpi,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Bd21.3.Endophytes.shoot.SV.6dpi',
                xlim = c(-25, 25),
                ylim = c(0,50),
                FCcutoff=1,
                pCutoff = 0.05,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/Bd21.3.Endophytes.shoot.SV.6dpi.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)

Bd21.3.Endophytes.shoot.SV.13dpi <- results(ddsTC, contrast=list("groupBd21.3.Endophytes.shoot.SV.13dpi", "groupBd21.3.Endophytes.shoot.Mock.13dpi"))
EnhancedVolcano(Bd21.3.Endophytes.shoot.SV.13dpi,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Bd21.3.Endophytes.shoot.SV.13dpi',
                xlim = c(-25, 25),
                ylim = c(0,50),
                FCcutoff=1,
                pCutoff = 0.05,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/Bd21.3.Endophytes.shoot.SV.13dpi.pdf",
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
# 936
write.table(x = sigres_Bd21.fungi.leaves.PCA.2.dpi_up,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.fungi.leaves.PCA.2.dpi_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

sigres_Bd21.fungi.leaves.PCA.2.dpi_down <- Bd21.fungi.leaves.PCA.2.dpi[(Bd21.fungi.leaves.PCA.2.dpi$padj<0.05 & Bd21.fungi.leaves.PCA.2.dpi$log2FoldChange<= -1),]
nb_Bd21.fungi.leaves.PCA.2.dpi_down <- nrow(sigres_Bd21.fungi.leaves.PCA.2.dpi_down)
#1256
write.table(x = sigres_Bd21.fungi.leaves.PCA.2.dpi_down,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.fungi.leaves.PCA.2.dpi_down.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)


sigres_Bd21.fungi.leaves.PCA.2.dpi <- Bd21.fungi.leaves.PCA.2.dpi[(Bd21.fungi.leaves.PCA.2.dpi$padj<0.05 & abs(Bd21.fungi.leaves.PCA.2.dpi$log2FoldChange) >= 1 ),]
nb_sigres_Bd21.fungi.leaves.PCA.2.dpi_DEG <- nrow(sigres_Bd21.fungi.leaves.PCA.2.dpi)
#2192
sigres_Bd21.fungi.leaves.PCA.2.dpi_strict <- Bd21.fungi.leaves.PCA.2.dpi[(Bd21.fungi.leaves.PCA.2.dpi$padj<0.0000000001 & abs(Bd21.fungi.leaves.PCA.2.dpi$log2FoldChange) >= 5 ),]
Bd21.fungi.leaves.PCA.2.dpi[(Bd21.fungi.leaves.PCA.2.dpi$padj<0.0000000001 & Bd21.fungi.leaves.PCA.2.dpi$log2FoldChange >= 5 ),]

nb_sigres_Bd21.fungi.leaves.PCA.2.dpi_DEG_strict <- nrow(sigres_Bd21.fungi.leaves.PCA.2.dpi_strict)
write.table(x = sigres_Bd21.fungi.leaves.PCA.2.dpi_strict,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_sigres_Bd21.fungi.leaves.PCA.2.dpi_DEG_strict.txt",
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
# 455
write.table(x = sigres_Bd21.fungi.leaves.PCA.4.dpi_up,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.fungi.leaves.PCA.4.dpi_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

sigres_Bd21.fungi.leaves.PCA.4.dpi_down <- Bd21.fungi.leaves.PCA.4.dpi[(Bd21.fungi.leaves.PCA.4.dpi$padj<0.05 & Bd21.fungi.leaves.PCA.4.dpi$log2FoldChange<= -1),]
nb_Bd21.fungi.leaves.PCA.4.dpi_down <- nrow(sigres_Bd21.fungi.leaves.PCA.4.dpi_down)
#1510
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
#1965
write.table(x = sigres_Bd21.fungi.leaves.PCA.4.dpi,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_sigres_Bd21.fungi.leaves.PCA.4.dpi_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
sigres_Bd21.fungi.leaves.PCA.4.dpi_strict <- Bd21.fungi.leaves.PCA.4.dpi[(Bd21.fungi.leaves.PCA.4.dpi$padj<0.0000000001 & abs(Bd21.fungi.leaves.PCA.4.dpi$log2FoldChange) >= 5 ),]
nb_sigres_Bd21.fungi.leaves.PCA.4.dpi_DEG_strict <- nrow(sigres_Bd21.fungi.leaves.PCA.4.dpi_strict)
write.table(x = sigres_Bd21.fungi.leaves.PCA.4.dpi_strict,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_sigres_Bd21.fungi.leaves.PCA.4.dpi_DEG_strict.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
Bd21.fungi.leaves.PCA.4.dpi[(Bd21.fungi.leaves.PCA.4.dpi$padj<0.0000000001 & Bd21.fungi.leaves.PCA.4.dpi$log2FoldChange >= 5 ),]

#Bd21.fungi.leaves.PCA.6.dpi
Bd21.fungi.leaves.PCA.6.dpi <- Bd21.fungi.leaves.PCA.6.dpi[!is.na(Bd21.fungi.leaves.PCA.6.dpi$padj),]
sigres_Bd21.fungi.leaves.PCA.6.dpi_up <- Bd21.fungi.leaves.PCA.6.dpi [(Bd21.fungi.leaves.PCA.6.dpi $padj<0.05 & Bd21.fungi.leaves.PCA.6.dpi $log2FoldChange>= 1),]
nb_Bd21.fungi.leaves.PCA.6.dpi_up <- nrow(sigres_Bd21.fungi.leaves.PCA.6.dpi_up)
# 547
write.table(x = sigres_Bd21.fungi.leaves.PCA.6.dpi_up,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.fungi.leaves.PCA.6.dpi_up.txt",
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
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.fungi.leaves.PCA.6.dpi_down.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)


sigres_Bd21.fungi.leaves.PCA.6.dpi <- Bd21.fungi.leaves.PCA.6.dpi[(Bd21.fungi.leaves.PCA.6.dpi$padj<0.05 & abs(Bd21.fungi.leaves.PCA.6.dpi$log2FoldChange) >= 1 ),]
nb_sigres_Bd21.fungi.leaves.PCA.6.dpi_DEG <- nrow(sigres_Bd21.fungi.leaves.PCA.6.dpi)
#1701
write.table(x = sigres_Bd21.fungi.leaves.PCA.6.dpi,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_sigres_Bd21.fungi.leaves.PCA.6.dpi_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

sigres_Bd21.fungi.leaves.PCA.6.dpi_strict <- Bd21.fungi.leaves.PCA.6.dpi[(Bd21.fungi.leaves.PCA.6.dpi$padj<0.0000000001 & abs(Bd21.fungi.leaves.PCA.6.dpi$log2FoldChange) >= 5 ),]
nb_sigres_Bd21.fungi.leaves.PCA.6.dpi_DEG_strict <- nrow(sigres_Bd21.fungi.leaves.PCA.6.dpi_strict)
write.table(x = sigres_Bd21.fungi.leaves.PCA.6.dpi_strict,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_sigres_Bd21.fungi.leaves.PCA.6.dpi_DEG_strict.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

sigres_Bd21.fungi.leaves.PCA.all_strict <- union(union(rownames(sigres_Bd21.fungi.leaves.PCA.2.dpi_strict),rownames(sigres_Bd21.fungi.leaves.PCA.4.dpi_strict)),rownames(sigres_Bd21.fungi.leaves.PCA.6.dpi_strict))
write.table(x = sigres_Bd21.fungi.leaves.PCA.all_strict,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_sigres_Bd21.fungi.leaves.PCA.all_strict_gene.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
Bd21.fungi.leaves.PCA.6.dpi[(Bd21.fungi.leaves.PCA.6.dpi$padj<0.0000000001 & Bd21.fungi.leaves.PCA.6.dpi$log2FoldChange >= 5 ),]

intersect(intersect(rownames(sigres_Bd21.fungi.leaves.PCA.2.dpi_strict),rownames(sigres_Bd21.fungi.leaves.PCA.4.dpi_strict)),rownames(sigres_Bd21.fungi.leaves.PCA.6.dpi_strict))

#Bd21.3.leaves.bacteria.XT.2.dpi
head(Bd21.3.leaves.bacteria.XT.2.dpi)
Bd21.3.leaves.bacteria.XT.2.dpi <- Bd21.3.leaves.bacteria.XT.2.dpi[!is.na(Bd21.3.leaves.bacteria.XT.2.dpi$padj),]
sigres_Bd21.3.leaves.bacteria.XT.2.dpi_up <- Bd21.3.leaves.bacteria.XT.2.dpi [(Bd21.3.leaves.bacteria.XT.2.dpi$padj<0.05 & Bd21.3.leaves.bacteria.XT.2.dpi $log2FoldChange>= 1),]
nb_Bd21.3.leaves.bacteria.XT.2.dpi_up <- nrow(sigres_Bd21.3.leaves.bacteria.XT.2.dpi_up)
# 152
write.table(x = sigres_Bd21.3.leaves.bacteria.XT.2.dpi_up,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.leaves.bacteria.XT.2.dpi_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

sigres_Bd21.3.leaves.bacteria.XT.2.dpi_down <- Bd21.3.leaves.bacteria.XT.2.dpi[(Bd21.3.leaves.bacteria.XT.2.dpi$padj<0.05 & Bd21.3.leaves.bacteria.XT.2.dpi$log2FoldChange<= -1),]
nb_Bd21.3.leaves.bacteria.XT.2.dpi_down <- nrow(sigres_Bd21.3.leaves.bacteria.XT.2.dpi_down)
#188
write.table(x = sigres_Bd21.3.leaves.bacteria.XT.2.dpi_down,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.leaves.bacteria.XT.2.dpi_down.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)


sigres_Bd21.3.leaves.bacteria.XT.2.dpi <- Bd21.3.leaves.bacteria.XT.2.dpi[(Bd21.3.leaves.bacteria.XT.2.dpi$padj<0.05 & abs(Bd21.3.leaves.bacteria.XT.2.dpi$log2FoldChange) >= 1 ),]
nb_sigres_Bd21.3.leaves.bacteria.XT.2.dpi_DEG <- nrow(sigres_Bd21.3.leaves.bacteria.XT.2.dpi)
#340
write.table(x = sigres_Bd21.3.leaves.bacteria.XT.2.dpi,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.leaves.bacteria.XT.2.dpi_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
##strict
sigres_Bd21.3.leaves.bacteria.XT.2.dpi_strict <- Bd21.3.leaves.bacteria.XT.2.dpi[(Bd21.3.leaves.bacteria.XT.2.dpi$padj<0.0000000001 & abs(Bd21.3.leaves.bacteria.XT.2.dpi$log2FoldChange) >= 5 ),]
nb_sigres_Bd21.3.leaves.bacteria.XT.2.dpi_DEG_strict <- nrow(sigres_Bd21.3.leaves.bacteria.XT.2.dpi)
Bd21.3.leaves.bacteria.XT.2.dpi[(Bd21.3.leaves.bacteria.XT.2.dpi$padj<0.0000000001 & Bd21.3.leaves.bacteria.XT.2.dpi$log2FoldChange >= 5 ),]

write.table(x = sigres_Bd21.3.leaves.bacteria.XT.2.dpi_strict,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.leaves.bacteria.XT.2.dpi_DEG_strict.txt",
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
# 1850
write.table(x = sigres_Bd21.3.leaves.bacteria.XT.4.dpi_up,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.leaves.bacteria.XT.4.dpi_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

sigres_Bd21.3.leaves.bacteria.XT.4.dpi_down <- Bd21.3.leaves.bacteria.XT.4.dpi[(Bd21.3.leaves.bacteria.XT.4.dpi$padj<0.05 & Bd21.3.leaves.bacteria.XT.4.dpi$log2FoldChange<= -1),]
nb_Bd21.3.leaves.bacteria.XT.4.dpi_down <- nrow(sigres_Bd21.3.leaves.bacteria.XT.4.dpi_down)
#1318
write.table(x = sigres_Bd21.3.leaves.bacteria.XT.4.dpi_down,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.leaves.bacteria.XT.4.dpi_down.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)


sigres_Bd21.3.leaves.bacteria.XT.4.dpi <- Bd21.3.leaves.bacteria.XT.4.dpi[(Bd21.3.leaves.bacteria.XT.4.dpi$padj<0.05 & abs(Bd21.3.leaves.bacteria.XT.4.dpi$log2FoldChange) >= 1 ),]
nb_sigres_Bd21.3.leaves.bacteria.XT.4.dpi_DEG <- nrow(sigres_Bd21.3.leaves.bacteria.XT.4.dpi)
#3168
write.table(x = sigres_Bd21.3.leaves.bacteria.XT.4.dpi,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.leaves.bacteria.XT.4.dpi_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
#strict
sigres_Bd21.3.leaves.bacteria.XT.4.dpi_strict <- Bd21.3.leaves.bacteria.XT.4.dpi[(Bd21.3.leaves.bacteria.XT.4.dpi$padj<0.0000000001 & abs(Bd21.3.leaves.bacteria.XT.4.dpi$log2FoldChange) >= 5 ),]
nb_sigres_Bd21.3.leaves.bacteria.XT.4.dpi_DEG_strict <- nrow(sigres_Bd21.3.leaves.bacteria.XT.4.dpi)
rownames(Bd21.3.leaves.bacteria.XT.4.dpi[(Bd21.3.leaves.bacteria.XT.4.dpi$padj<0.0000000001 & Bd21.3.leaves.bacteria.XT.4.dpi$log2FoldChange >= 5 ),])
union(rownames(sigres_Bd21.3.leaves.bacteria.XT.2.dpi_strict),rownames(sigres_Bd21.3.leaves.bacteria.XT.4.dpi_strict))

intersect(rownames(sigres_Bd21.3.leaves.bacteria.XT.2.dpi_strict),rownames(sigres_Bd21.3.leaves.bacteria.XT.4.dpi_strict))


write.table(x = sigres_Bd21.3.leaves.bacteria.XT.4.dpi_strict,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.leaves.bacteria.XT.4.dpi_DEG_strict.txt",
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
# 596
write.table(x = sigres_Bd21.3.leaves.virus.BSMV.3.dpi_up,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.leaves.virus.BSMV.3.dpi_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

sigres_Bd21.3.leaves.virus.BSMV.3.dpi_down <- Bd21.3.leaves.virus.BSMV.3.dpi[(Bd21.3.leaves.virus.BSMV.3.dpi$padj<0.05 & Bd21.3.leaves.virus.BSMV.3.dpi$log2FoldChange<= -1),]
nb_Bd21.3.leaves.virus.BSMV.3.dpi_down <- nrow(sigres_Bd21.3.leaves.virus.BSMV.3.dpi_down)
#1311
write.table(x = sigres_Bd21.3.leaves.virus.BSMV.3.dpi_down,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.leaves.virus.BSMV.3.dpi_down.txt",
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
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.leaves.virus.BSMV.3.dpi_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
#Bd21.3.leaves.virus.BSMV.7.dpi
Bd21.3.leaves.virus.BSMV.7.dpi <- Bd21.3.leaves.virus.BSMV.7.dpi[!is.na(Bd21.3.leaves.virus.BSMV.7.dpi$padj),]
sigres_Bd21.3.leaves.virus.BSMV.7.dpi_up <- Bd21.3.leaves.virus.BSMV.7.dpi [(Bd21.3.leaves.virus.BSMV.7.dpi $padj<0.05 & Bd21.3.leaves.virus.BSMV.7.dpi $log2FoldChange>= 1),]
nb_Bd21.3.leaves.virus.BSMV.7.dpi_up <- nrow(sigres_Bd21.3.leaves.virus.BSMV.7.dpi_up)
# 682
write.table(x = sigres_Bd21.3.leaves.virus.BSMV.7.dpi_up,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.leaves.virus.BSMV.7.dpi_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

sigres_Bd21.3.leaves.virus.BSMV.7.dpi_down <- Bd21.3.leaves.virus.BSMV.7.dpi[(Bd21.3.leaves.virus.BSMV.7.dpi$padj<0.05 & Bd21.3.leaves.virus.BSMV.7.dpi$log2FoldChange<= -1),]
nb_Bd21.3.leaves.virus.BSMV.7.dpi_down <- nrow(sigres_Bd21.3.leaves.virus.BSMV.7.dpi_down)
#724
write.table(x = sigres_Bd21.3.leaves.virus.BSMV.7.dpi_down,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.leaves.virus.BSMV.7.dpi_down.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)


sigres_Bd21.3.leaves.virus.BSMV.7.dpi <- Bd21.3.leaves.virus.BSMV.7.dpi[(Bd21.3.leaves.virus.BSMV.7.dpi$padj<0.05 & abs(Bd21.3.leaves.virus.BSMV.7.dpi$log2FoldChange) >= 1 ),]
nb_sigres_Bd21.3.leaves.virus.BSMV.7.dpi_DEG <- nrow(sigres_Bd21.3.leaves.virus.BSMV.7.dpi)
#1406
write.table(x = sigres_Bd21.3.leaves.virus.BSMV.7.dpi,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_sigres_Bd21.3.leaves.virus.BSMV.7.dpi_DEG.txt",
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
# 1544
write.table(x = sigres_Bd21.3.leaves.virus.BSMV.14.dpi_up,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.leaves.virus.BSMV.14.dpi_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

sigres_Bd21.3.leaves.virus.BSMV.14.dpi_down <- Bd21.3.leaves.virus.BSMV.14.dpi[(Bd21.3.leaves.virus.BSMV.14.dpi$padj<0.05 & Bd21.3.leaves.virus.BSMV.14.dpi$log2FoldChange<= -1),]
nb_Bd21.3.leaves.virus.BSMV.14.dpi_down <- nrow(sigres_Bd21.3.leaves.virus.BSMV.14.dpi_down)
#1217
write.table(x = sigres_Bd21.3.leaves.virus.BSMV.14.dpi_down,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.leaves.virus.BSMV.14.dpi_down.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)


sigres_Bd21.3.leaves.virus.BSMV.14.dpi <- Bd21.3.leaves.virus.BSMV.14.dpi[(Bd21.3.leaves.virus.BSMV.14.dpi$padj<0.05 & abs(Bd21.3.leaves.virus.BSMV.14.dpi$log2FoldChange) >= 1 ),]
nb_sigres_Bd21.3.leaves.virus.BSMV.14.dpi_DEG <- nrow(sigres_Bd21.3.leaves.virus.BSMV.14.dpi)
#2761
write.table(x = sigres_Bd21.3.leaves.virus.BSMV.14.dpi,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.leaves.virus.BSMV.14.dpi_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
###MMMV
Bd21.3.leaves.virus.MMMV.3.dpi <- Bd21.3.leaves.virus.MMMV.3.dpi[!is.na(Bd21.3.leaves.virus.MMMV.3.dpi$padj),]
sigres_Bd21.3.leaves.virus.MMMV.3.dpi_up <- Bd21.3.leaves.virus.MMMV.3.dpi [(Bd21.3.leaves.virus.MMMV.3.dpi $padj<0.05 & Bd21.3.leaves.virus.MMMV.3.dpi $log2FoldChange>= 1),]
nb_Bd21.3.leaves.virus.MMMV.3.dpi_up <- nrow(sigres_Bd21.3.leaves.virus.MMMV.3.dpi_up)
# 1140
write.table(x = sigres_Bd21.3.leaves.virus.MMMV.3.dpi_up,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.leaves.virus.MMMV.3.dpi_up.txt",
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
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.leaves.virus.MMMV.3.dpi_down.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)


sigres_Bd21.3.leaves.virus.MMMV.3.dpi <- Bd21.3.leaves.virus.MMMV.3.dpi[(Bd21.3.leaves.virus.MMMV.3.dpi$padj<0.05 & abs(Bd21.3.leaves.virus.MMMV.3.dpi$log2FoldChange) >= 1 ),]
nb_sigres_Bd21.3.leaves.virus.MMMV.3.dpi_DEG <- nrow(sigres_Bd21.3.leaves.virus.MMMV.3.dpi)
#1726
write.table(x = sigres_Bd21.3.leaves.virus.MMMV.3.dpi,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.leaves.virus.MMMV.3.dpi_DEG.txt",
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
# 2111
write.table(x = sigres_Bd21.3.leaves.virus.MMMV.7.dpi_up,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.leaves.virus.MMMV.7.dpi_up.txt",
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
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.leaves.virus.MMMV.7.dpi_down.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)


sigres_Bd21.3.leaves.virus.MMMV.7.dpi <- Bd21.3.leaves.virus.MMMV.7.dpi[(Bd21.3.leaves.virus.MMMV.7.dpi$padj<0.05 & abs(Bd21.3.leaves.virus.MMMV.7.dpi$log2FoldChange) >= 1 ),]
nb_sigres_Bd21.3.leaves.virus.MMMV.7.dpi_DEG <- nrow(sigres_Bd21.3.leaves.virus.MMMV.7.dpi)
#4373
write.table(x = sigres_Bd21.3.leaves.virus.MMMV.7.dpi,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_sigres_Bd21.3.leaves.virus.MMMV.7.dpi_DEG.txt",
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
# 3913
write.table(x = sigres_Bd21.3.leaves.virus.MMMV.14.dpi_up,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.leaves.virus.MMMV.14.dpi_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

sigres_Bd21.3.leaves.virus.MMMV.14.dpi_down <- Bd21.3.leaves.virus.MMMV.14.dpi[(Bd21.3.leaves.virus.MMMV.14.dpi$padj<0.05 & Bd21.3.leaves.virus.MMMV.14.dpi$log2FoldChange<= -1),]
nb_Bd21.3.leaves.virus.MMMV.14.dpi_down <- nrow(sigres_Bd21.3.leaves.virus.MMMV.14.dpi_down)
#1859
write.table(x = sigres_Bd21.3.leaves.virus.MMMV.14.dpi_down,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.leaves.virus.MMMV.14.dpi_down.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)


sigres_Bd21.3.leaves.virus.MMMV.14.dpi <- Bd21.3.leaves.virus.MMMV.14.dpi[(Bd21.3.leaves.virus.MMMV.14.dpi$padj<0.05 & abs(Bd21.3.leaves.virus.MMMV.14.dpi$log2FoldChange) >= 1 ),]
nb_sigres_Bd21.3.leaves.virus.MMMV.14.dpi_DEG <- nrow(sigres_Bd21.3.leaves.virus.MMMV.14.dpi)
#5772
write.table(x = sigres_Bd21.3.leaves.virus.MMMV.14.dpi,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.leaves.virus.MMMV.14.dpi_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
##WSMV
Bd21.3.leaves.virus.WSMV.3.dpi <- Bd21.3.leaves.virus.WSMV.3.dpi[!is.na(Bd21.3.leaves.virus.WSMV.3.dpi$padj),]
sigres_Bd21.3.leaves.virus.WSMV.3.dpi_up <- Bd21.3.leaves.virus.WSMV.3.dpi [(Bd21.3.leaves.virus.WSMV.3.dpi $padj<0.05 & Bd21.3.leaves.virus.WSMV.3.dpi $log2FoldChange>= 1),]
nb_Bd21.3.leaves.virus.WSMV.3.dpi_up <- nrow(sigres_Bd21.3.leaves.virus.WSMV.3.dpi_up)
# 595
write.table(x = sigres_Bd21.3.leaves.virus.WSMV.3.dpi_up,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.leaves.virus.WSMV.3.dpi_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

sigres_Bd21.3.leaves.virus.WSMV.3.dpi_down <- Bd21.3.leaves.virus.WSMV.3.dpi[(Bd21.3.leaves.virus.WSMV.3.dpi$padj<0.05 & Bd21.3.leaves.virus.WSMV.3.dpi$log2FoldChange<= -1),]
nb_Bd21.3.leaves.virus.WSMV.3.dpi_down <- nrow(sigres_Bd21.3.leaves.virus.WSMV.3.dpi_down)
#492
write.table(x = sigres_Bd21.3.leaves.virus.WSMV.3.dpi_down,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.leaves.virus.WSMV.3.dpi_down.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

##WSMV
sigres_Bd21.3.leaves.virus.WSMV.3.dpi <- Bd21.3.leaves.virus.WSMV.3.dpi[(Bd21.3.leaves.virus.WSMV.3.dpi$padj<0.05 & abs(Bd21.3.leaves.virus.WSMV.3.dpi$log2FoldChange) >= 1 ),]
nb_sigres_Bd21.3.leaves.virus.WSMV.3.dpi_DEG <- nrow(sigres_Bd21.3.leaves.virus.WSMV.3.dpi)
#1087
write.table(x = sigres_Bd21.3.leaves.virus.WSMV.3.dpi,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.leaves.virus.WSMV.3.dpi_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
#Bd21.3.leaves.virus.WSMV.7.dpi
Bd21.3.leaves.virus.WSMV.7.dpi <- Bd21.3.leaves.virus.WSMV.7.dpi[!is.na(Bd21.3.leaves.virus.WSMV.7.dpi$padj),]
sigres_Bd21.3.leaves.virus.WSMV.7.dpi_up <- Bd21.3.leaves.virus.WSMV.7.dpi [(Bd21.3.leaves.virus.WSMV.7.dpi $padj<0.05 & Bd21.3.leaves.virus.WSMV.7.dpi $log2FoldChange>= 1),]
nb_Bd21.3.leaves.virus.WSMV.7.dpi_up <- nrow(sigres_Bd21.3.leaves.virus.WSMV.7.dpi_up)
# 497
write.table(x = sigres_Bd21.3.leaves.virus.WSMV.7.dpi_up,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.leaves.virus.WSMV.7.dpi_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

sigres_Bd21.3.leaves.virus.WSMV.7.dpi_down <- Bd21.3.leaves.virus.WSMV.7.dpi[(Bd21.3.leaves.virus.WSMV.7.dpi$padj<0.05 & Bd21.3.leaves.virus.WSMV.7.dpi$log2FoldChange<= -1),]
nb_Bd21.3.leaves.virus.WSMV.7.dpi_down <- nrow(sigres_Bd21.3.leaves.virus.WSMV.7.dpi_down)
#1126
write.table(x = sigres_Bd21.3.leaves.virus.WSMV.7.dpi_down,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.leaves.virus.WSMV.7.dpi_down.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)


sigres_Bd21.3.leaves.virus.WSMV.7.dpi <- Bd21.3.leaves.virus.WSMV.7.dpi[(Bd21.3.leaves.virus.WSMV.7.dpi$padj<0.05 & abs(Bd21.3.leaves.virus.WSMV.7.dpi$log2FoldChange) >= 1 ),]
nb_sigres_Bd21.3.leaves.virus.WSMV.7.dpi_DEG <- nrow(sigres_Bd21.3.leaves.virus.WSMV.7.dpi)
#1623
write.table(x = sigres_Bd21.3.leaves.virus.WSMV.7.dpi,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_sigres_Bd21.3.leaves.virus.WSMV.7.dpi_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
#Bd21.3.leaves.virus.WSMV.14.dpi
Bd21.3.leaves.virus.WSMV.14.dpi <- Bd21.3.leaves.virus.WSMV.14.dpi[!is.na(Bd21.3.leaves.virus.WSMV.14.dpi$padj),]
sigres_Bd21.3.leaves.virus.WSMV.14.dpi_up <- Bd21.3.leaves.virus.WSMV.14.dpi [(Bd21.3.leaves.virus.WSMV.14.dpi $padj<0.05 & Bd21.3.leaves.virus.WSMV.14.dpi $log2FoldChange>= 1),]
nb_Bd21.3.leaves.virus.WSMV.14.dpi_up <- nrow(sigres_Bd21.3.leaves.virus.WSMV.14.dpi_up)
# 2437
write.table(x = sigres_Bd21.3.leaves.virus.WSMV.14.dpi_up,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.leaves.virus.WSMV.14.dpi_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

sigres_Bd21.3.leaves.virus.WSMV.14.dpi_down <- Bd21.3.leaves.virus.WSMV.14.dpi[(Bd21.3.leaves.virus.WSMV.14.dpi$padj<0.05 & Bd21.3.leaves.virus.WSMV.14.dpi$log2FoldChange<= -1),]
nb_Bd21.3.leaves.virus.WSMV.14.dpi_down <- nrow(sigres_Bd21.3.leaves.virus.WSMV.14.dpi_down)
#788
write.table(x = sigres_Bd21.3.leaves.virus.WSMV.14.dpi_down,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.leaves.virus.WSMV.14.dpi_down.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)


sigres_Bd21.3.leaves.virus.WSMV.14.dpi <- Bd21.3.leaves.virus.WSMV.14.dpi[(Bd21.3.leaves.virus.WSMV.14.dpi$padj<0.05 & abs(Bd21.3.leaves.virus.WSMV.14.dpi$log2FoldChange) >= 1 ),]
nb_sigres_Bd21.3.leaves.virus.WSMV.14.dpi_DEG <- nrow(sigres_Bd21.3.leaves.virus.WSMV.14.dpi)
#3225
write.table(x = sigres_Bd21.3.leaves.virus.WSMV.14.dpi,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.leaves.virus.WSMV.14.dpi_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
#Root virus
#Bd21.3.root.virus.BSMV.3.dpi
Bd21.3.root.virus.BSMV.3.dpi <- Bd21.3.root.virus.BSMV.3.dpi[!is.na(Bd21.3.root.virus.BSMV.3.dpi$padj),]
sigres_Bd21.3.root.virus.BSMV.3.dpi_up <- Bd21.3.root.virus.BSMV.3.dpi [(Bd21.3.root.virus.BSMV.3.dpi $padj<0.05 & Bd21.3.root.virus.BSMV.3.dpi $log2FoldChange>= 1),]
nb_Bd21.3.root.virus.BSMV.3.dpi_up <- nrow(sigres_Bd21.3.root.virus.BSMV.3.dpi_up)
# 809
write.table(x = sigres_Bd21.3.root.virus.BSMV.3.dpi_up,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.root.virus.BSMV.3.dpi_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

sigres_Bd21.3.root.virus.BSMV.3.dpi_down <- Bd21.3.root.virus.BSMV.3.dpi[(Bd21.3.root.virus.BSMV.3.dpi$padj<0.05 & Bd21.3.root.virus.BSMV.3.dpi$log2FoldChange<= -1),]
nb_Bd21.3.root.virus.BSMV.3.dpi_down <- nrow(sigres_Bd21.3.root.virus.BSMV.3.dpi_down)
#1110
write.table(x = sigres_Bd21.3.root.virus.BSMV.3.dpi_down,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.root.virus.BSMV.3.dpi_down.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)


sigres_Bd21.3.root.virus.BSMV.3.dpi <- Bd21.3.root.virus.BSMV.3.dpi[(Bd21.3.root.virus.BSMV.3.dpi$padj<0.05 & abs(Bd21.3.root.virus.BSMV.3.dpi$log2FoldChange) >= 1 ),]
nb_sigres_Bd21.3.root.virus.BSMV.3.dpi_DEG <- nrow(sigres_Bd21.3.root.virus.BSMV.3.dpi)
#1919
write.table(x = sigres_Bd21.3.root.virus.BSMV.3.dpi,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.root.virus.BSMV.3.dpi_DEG.txt",
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
# 774
write.table(x = sigres_Bd21.3.root.virus.BSMV.7.dpi_up,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.root.virus.BSMV.7.dpi_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

sigres_Bd21.3.root.virus.BSMV.7.dpi_down <- Bd21.3.root.virus.BSMV.7.dpi[(Bd21.3.root.virus.BSMV.7.dpi$padj<0.05 & Bd21.3.root.virus.BSMV.7.dpi$log2FoldChange<= -1),]
nb_Bd21.3.root.virus.BSMV.7.dpi_down <- nrow(sigres_Bd21.3.root.virus.BSMV.7.dpi_down)
#514
write.table(x = sigres_Bd21.3.root.virus.BSMV.7.dpi_down,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.root.virus.BSMV.7.dpi_down.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)


sigres_Bd21.3.root.virus.BSMV.7.dpi <- Bd21.3.root.virus.BSMV.7.dpi[(Bd21.3.root.virus.BSMV.7.dpi$padj<0.05 & abs(Bd21.3.root.virus.BSMV.7.dpi$log2FoldChange) >= 1 ),]
nb_sigres_Bd21.3.root.virus.BSMV.7.dpi_DEG <- nrow(sigres_Bd21.3.root.virus.BSMV.7.dpi)
#1288
write.table(x = sigres_Bd21.3.root.virus.BSMV.7.dpi,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_sigres_Bd21.3.root.virus.BSMV.7.dpi_DEG.txt",
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
# 2405
write.table(x = sigres_Bd21.3.root.virus.BSMV.14.dpi_up,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.root.virus.BSMV.14.dpi_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

sigres_Bd21.3.root.virus.BSMV.14.dpi_down <- Bd21.3.root.virus.BSMV.14.dpi[(Bd21.3.root.virus.BSMV.14.dpi$padj<0.05 & Bd21.3.root.virus.BSMV.14.dpi$log2FoldChange<= -1),]
nb_Bd21.3.root.virus.BSMV.14.dpi_down <- nrow(sigres_Bd21.3.root.virus.BSMV.14.dpi_down)
#1191
write.table(x = sigres_Bd21.3.root.virus.BSMV.14.dpi_down,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.root.virus.BSMV.14.dpi_down.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)


sigres_Bd21.3.root.virus.BSMV.14.dpi <- Bd21.3.root.virus.BSMV.14.dpi[(Bd21.3.root.virus.BSMV.14.dpi$padj<0.05 & abs(Bd21.3.root.virus.BSMV.14.dpi$log2FoldChange) >= 1 ),]
nb_sigres_Bd21.3.root.virus.BSMV.14.dpi_DEG <- nrow(sigres_Bd21.3.root.virus.BSMV.14.dpi)
#3596
write.table(x = sigres_Bd21.3.root.virus.BSMV.14.dpi,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.root.virus.BSMV.14.dpi_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
###MMMV
Bd21.3.root.virus.MMMV.3.dpi <- Bd21.3.root.virus.MMMV.3.dpi[!is.na(Bd21.3.root.virus.MMMV.3.dpi$padj),]
sigres_Bd21.3.root.virus.MMMV.3.dpi_up <- Bd21.3.root.virus.MMMV.3.dpi [(Bd21.3.root.virus.MMMV.3.dpi $padj<0.05 & Bd21.3.root.virus.MMMV.3.dpi $log2FoldChange>= 1),]
nb_Bd21.3.root.virus.MMMV.3.dpi_up <- nrow(sigres_Bd21.3.root.virus.MMMV.3.dpi_up)
# 1232
write.table(x = sigres_Bd21.3.root.virus.MMMV.3.dpi_up,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.root.virus.MMMV.3.dpi_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

sigres_Bd21.3.root.virus.MMMV.3.dpi_down <- Bd21.3.root.virus.MMMV.3.dpi[(Bd21.3.root.virus.MMMV.3.dpi$padj<0.05 & Bd21.3.root.virus.MMMV.3.dpi$log2FoldChange<= -1),]
nb_Bd21.3.root.virus.MMMV.3.dpi_down <- nrow(sigres_Bd21.3.root.virus.MMMV.3.dpi_down)
#549
write.table(x = sigres_Bd21.3.root.virus.MMMV.3.dpi_down,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.root.virus.MMMV.3.dpi_down.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)


sigres_Bd21.3.root.virus.MMMV.3.dpi <- Bd21.3.root.virus.MMMV.3.dpi[(Bd21.3.root.virus.MMMV.3.dpi$padj<0.05 & abs(Bd21.3.root.virus.MMMV.3.dpi$log2FoldChange) >= 1 ),]
nb_sigres_Bd21.3.root.virus.MMMV.3.dpi_DEG <- nrow(sigres_Bd21.3.root.virus.MMMV.3.dpi)
#1781
write.table(x = sigres_Bd21.3.root.virus.MMMV.3.dpi,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.root.virus.MMMV.3.dpi_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
#Bd21.3.root.virus.MMMV.7.dpi
Bd21.3.root.virus.MMMV.7.dpi <- Bd21.3.root.virus.MMMV.7.dpi[!is.na(Bd21.3.root.virus.MMMV.7.dpi$padj),]
sigres_Bd21.3.root.virus.MMMV.7.dpi_up <- Bd21.3.root.virus.MMMV.7.dpi [(Bd21.3.root.virus.MMMV.7.dpi $padj<0.05 & Bd21.3.root.virus.MMMV.7.dpi $log2FoldChange>= 1),]
nb_Bd21.3.root.virus.MMMV.7.dpi_up <- nrow(sigres_Bd21.3.root.virus.MMMV.7.dpi_up)
# 1776
write.table(x = sigres_Bd21.3.root.virus.MMMV.7.dpi_up,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.root.virus.MMMV.7.dpi_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

sigres_Bd21.3.root.virus.MMMV.7.dpi_down <- Bd21.3.root.virus.MMMV.7.dpi[(Bd21.3.root.virus.MMMV.7.dpi$padj<0.05 & Bd21.3.root.virus.MMMV.7.dpi$log2FoldChange<= -1),]
nb_Bd21.3.root.virus.MMMV.7.dpi_down <- nrow(sigres_Bd21.3.root.virus.MMMV.7.dpi_down)
#822
write.table(x = sigres_Bd21.3.root.virus.MMMV.7.dpi_down,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.root.virus.MMMV.7.dpi_down.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)


sigres_Bd21.3.root.virus.MMMV.7.dpi <- Bd21.3.root.virus.MMMV.7.dpi[(Bd21.3.root.virus.MMMV.7.dpi$padj<0.05 & abs(Bd21.3.root.virus.MMMV.7.dpi$log2FoldChange) >= 1 ),]
nb_sigres_Bd21.3.root.virus.MMMV.7.dpi_DEG <- nrow(sigres_Bd21.3.root.virus.MMMV.7.dpi)
#2598
write.table(x = sigres_Bd21.3.root.virus.MMMV.7.dpi,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_sigres_Bd21.3.root.virus.MMMV.7.dpi_DEG.txt",
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
# 4199
write.table(x = sigres_Bd21.3.root.virus.MMMV.14.dpi_up,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.root.virus.MMMV.14.dpi_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

sigres_Bd21.3.root.virus.MMMV.14.dpi_down <- Bd21.3.root.virus.MMMV.14.dpi[(Bd21.3.root.virus.MMMV.14.dpi$padj<0.05 & Bd21.3.root.virus.MMMV.14.dpi$log2FoldChange<= -1),]
nb_Bd21.3.root.virus.MMMV.14.dpi_down <- nrow(sigres_Bd21.3.root.virus.MMMV.14.dpi_down)
#4715
write.table(x = sigres_Bd21.3.root.virus.MMMV.14.dpi_down,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.root.virus.MMMV.14.dpi_down.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)


sigres_Bd21.3.root.virus.MMMV.14.dpi <- Bd21.3.root.virus.MMMV.14.dpi[(Bd21.3.root.virus.MMMV.14.dpi$padj<0.05 & abs(Bd21.3.root.virus.MMMV.14.dpi$log2FoldChange) >= 1 ),]
nb_sigres_Bd21.3.root.virus.MMMV.14.dpi_DEG <- nrow(sigres_Bd21.3.root.virus.MMMV.14.dpi)
#8914
write.table(x = sigres_Bd21.3.root.virus.MMMV.14.dpi,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.root.virus.MMMV.14.dpi_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
##WSMV
Bd21.3.root.virus.WSMV.3.dpi <- Bd21.3.root.virus.WSMV.3.dpi[!is.na(Bd21.3.root.virus.WSMV.3.dpi$padj),]
sigres_Bd21.3.root.virus.WSMV.3.dpi_up <- Bd21.3.root.virus.WSMV.3.dpi [(Bd21.3.root.virus.WSMV.3.dpi $padj<0.05 & Bd21.3.root.virus.WSMV.3.dpi $log2FoldChange>= 1),]
nb_Bd21.3.root.virus.WSMV.3.dpi_up <- nrow(sigres_Bd21.3.root.virus.WSMV.3.dpi_up)
# 1397
write.table(x = sigres_Bd21.3.root.virus.WSMV.3.dpi_up,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.root.virus.WSMV.3.dpi_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

sigres_Bd21.3.root.virus.WSMV.3.dpi_down <- Bd21.3.root.virus.WSMV.3.dpi[(Bd21.3.root.virus.WSMV.3.dpi$padj<0.05 & Bd21.3.root.virus.WSMV.3.dpi$log2FoldChange<= -1),]
nb_Bd21.3.root.virus.WSMV.3.dpi_down <- nrow(sigres_Bd21.3.root.virus.WSMV.3.dpi_down)
#1616
write.table(x = sigres_Bd21.3.root.virus.WSMV.3.dpi_down,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.root.virus.WSMV.3.dpi_down.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

##WSMV
sigres_Bd21.3.root.virus.WSMV.3.dpi <- Bd21.3.root.virus.WSMV.3.dpi[(Bd21.3.root.virus.WSMV.3.dpi$padj<0.05 & abs(Bd21.3.root.virus.WSMV.3.dpi$log2FoldChange) >= 1 ),]
nb_sigres_Bd21.3.root.virus.WSMV.3.dpi_DEG <- nrow(sigres_Bd21.3.root.virus.WSMV.3.dpi)
#3013
write.table(x = sigres_Bd21.3.root.virus.WSMV.3.dpi,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.root.virus.WSMV.3.dpi_DEG.txt",
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
# 1590
write.table(x = sigres_Bd21.3.root.virus.WSMV.7.dpi_up,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.root.virus.WSMV.7.dpi_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

sigres_Bd21.3.root.virus.WSMV.7.dpi_down <- Bd21.3.root.virus.WSMV.7.dpi[(Bd21.3.root.virus.WSMV.7.dpi$padj<0.05 & Bd21.3.root.virus.WSMV.7.dpi$log2FoldChange<= -1),]
nb_Bd21.3.root.virus.WSMV.7.dpi_down <- nrow(sigres_Bd21.3.root.virus.WSMV.7.dpi_down)
#1852
write.table(x = sigres_Bd21.3.root.virus.WSMV.7.dpi_down,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.root.virus.WSMV.7.dpi_down.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)


sigres_Bd21.3.root.virus.WSMV.7.dpi <- Bd21.3.root.virus.WSMV.7.dpi[(Bd21.3.root.virus.WSMV.7.dpi$padj<0.05 & abs(Bd21.3.root.virus.WSMV.7.dpi$log2FoldChange) >= 1 ),]
nb_sigres_Bd21.3.root.virus.WSMV.7.dpi_DEG <- nrow(sigres_Bd21.3.root.virus.WSMV.7.dpi)
#3442
write.table(x = sigres_Bd21.3.root.virus.WSMV.7.dpi,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_sigres_Bd21.3.root.virus.WSMV.7.dpi_DEG.txt",
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
# 1154
write.table(x = sigres_Bd21.3.root.virus.WSMV.14.dpi_up,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.root.virus.WSMV.14.dpi_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

sigres_Bd21.3.root.virus.WSMV.14.dpi_down <- Bd21.3.root.virus.WSMV.14.dpi[(Bd21.3.root.virus.WSMV.14.dpi$padj<0.05 & Bd21.3.root.virus.WSMV.14.dpi$log2FoldChange<= -1),]
nb_Bd21.3.root.virus.WSMV.14.dpi_down <- nrow(sigres_Bd21.3.root.virus.WSMV.14.dpi_down)
#421
write.table(x = sigres_Bd21.3.root.virus.WSMV.14.dpi_down,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.root.virus.WSMV.14.dpi_down.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)


sigres_Bd21.3.root.virus.WSMV.14.dpi <- Bd21.3.root.virus.WSMV.14.dpi[(Bd21.3.root.virus.WSMV.14.dpi$padj<0.05 & abs(Bd21.3.root.virus.WSMV.14.dpi$log2FoldChange) >= 1 ),]
nb_sigres_Bd21.3.root.virus.WSMV.14.dpi_DEG <- nrow(sigres_Bd21.3.root.virus.WSMV.14.dpi)
#1575
write.table(x = sigres_Bd21.3.root.virus.WSMV.14.dpi,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.root.virus.WSMV.14.dpi_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
#Bd21.3.Endophytes.root.SV.3dpi
Bd21.3.Endophytes.root.SV.3dpi <- Bd21.3.Endophytes.root.SV.3dpi[!is.na(Bd21.3.Endophytes.root.SV.3dpi$padj),]
sigres_Bd21.3.Endophytes.root.SV.3dpi_up <- Bd21.3.Endophytes.root.SV.3dpi [(Bd21.3.Endophytes.root.SV.3dpi $padj<0.05 & Bd21.3.Endophytes.root.SV.3dpi $log2FoldChange>= 1),]
nb_Bd21.3.Endophytes.root.SV.3dpi_up <- nrow(sigres_Bd21.3.Endophytes.root.SV.3dpi_up)
# 10
write.table(x = sigres_Bd21.3.Endophytes.root.SV.3dpi_up,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.Endophytes.root.SV.3dpi_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

sigres_Bd21.3.Endophytes.root.SV.3dpi_down <- Bd21.3.Endophytes.root.SV.3dpi[(Bd21.3.Endophytes.root.SV.3dpi$padj<0.05 & Bd21.3.Endophytes.root.SV.3dpi$log2FoldChange<= -1),]
nb_Bd21.3.Endophytes.root.SV.3dpi_down <- nrow(sigres_Bd21.3.Endophytes.root.SV.3dpi_down)
#41
write.table(x = sigres_Bd21.3.Endophytes.root.SV.3dpi_down,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.Endophytes.root.SV.3dpi_down.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
sigres_Bd21.3.Endophytes.root.SV.3dpi <- Bd21.3.Endophytes.root.SV.3dpi[(Bd21.3.Endophytes.root.SV.3dpi$padj<0.05 & abs(Bd21.3.Endophytes.root.SV.3dpi$log2FoldChange) >= 1 ),]
nb_sigres_Bd21.3.Endophytes.root.SV.3dpi_DEG <- nrow(sigres_Bd21.3.Endophytes.root.SV.3dpi)
#51
write.table(x = sigres_Bd21.3.Endophytes.root.SV.3dpi,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_sigres_Bd21.3.Endophytes.root.SV.3dpi_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
#Bd21.3.Endophytes.root.SV.6dpi
Bd21.3.Endophytes.root.SV.6dpi <- Bd21.3.Endophytes.root.SV.6dpi[!is.na(Bd21.3.Endophytes.root.SV.6dpi$padj),]
sigres_Bd21.3.Endophytes.root.SV.6dpi_up <- Bd21.3.Endophytes.root.SV.6dpi [(Bd21.3.Endophytes.root.SV.6dpi $padj<0.05 & Bd21.3.Endophytes.root.SV.6dpi $log2FoldChange>= 1),]
nb_Bd21.3.Endophytes.root.SV.6dpi_up <- nrow(sigres_Bd21.3.Endophytes.root.SV.6dpi_up)
# 54
write.table(x = sigres_Bd21.3.Endophytes.root.SV.6dpi_up,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.Endophytes.root.SV.6dpi_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

sigres_Bd21.3.Endophytes.root.SV.6dpi_down <- Bd21.3.Endophytes.root.SV.6dpi[(Bd21.3.Endophytes.root.SV.6dpi$padj<0.05 & Bd21.3.Endophytes.root.SV.6dpi$log2FoldChange<= -1),]
nb_Bd21.3.Endophytes.root.SV.6dpi_down <- nrow(sigres_Bd21.3.Endophytes.root.SV.6dpi_down)
#13
write.table(x = sigres_Bd21.3.Endophytes.root.SV.6dpi_down,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.Endophytes.root.SV.6dpi_down.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)


sigres_Bd21.3.Endophytes.root.SV.6dpi <- Bd21.3.Endophytes.root.SV.6dpi[(Bd21.3.Endophytes.root.SV.6dpi$padj<0.05 & abs(Bd21.3.Endophytes.root.SV.6dpi$log2FoldChange) >= 1 ),]
nb_sigres_Bd21.3.Endophytes.root.SV.6dpi_DEG <- nrow(sigres_Bd21.3.Endophytes.root.SV.6dpi)
#67
write.table(x = sigres_Bd21.3.Endophytes.root.SV.6dpi,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_sigres_Bd21.3.Endophytes.root.SV.6dpi_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
#Bd21.3.Endophytes.root.SV.13dpi
Bd21.3.Endophytes.root.SV.13dpi <- Bd21.3.Endophytes.root.SV.13dpi[!is.na(Bd21.3.Endophytes.root.SV.13dpi$padj),]
sigres_Bd21.3.Endophytes.root.SV.13dpi_up <- Bd21.3.Endophytes.root.SV.13dpi [(Bd21.3.Endophytes.root.SV.13dpi $padj<0.05 & Bd21.3.Endophytes.root.SV.13dpi $log2FoldChange>= 1),]
nb_Bd21.3.Endophytes.root.SV.13dpi_up <- nrow(sigres_Bd21.3.Endophytes.root.SV.13dpi_up)
# 12
write.table(x = sigres_Bd21.3.Endophytes.root.SV.13dpi_up,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.Endophytes.root.SV.13dpi_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

sigres_Bd21.3.Endophytes.root.SV.13dpi_down <- Bd21.3.Endophytes.root.SV.13dpi[(Bd21.3.Endophytes.root.SV.13dpi$padj<0.05 & Bd21.3.Endophytes.root.SV.13dpi$log2FoldChange<= -1),]
nb_Bd21.3.Endophytes.root.SV.13dpi_down <- nrow(sigres_Bd21.3.Endophytes.root.SV.13dpi_down)
#169
write.table(x = sigres_Bd21.3.Endophytes.root.SV.13dpi_down,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.Endophytes.root.SV.13dpi_down.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)


sigres_Bd21.3.Endophytes.root.SV.13dpi <- Bd21.3.Endophytes.root.SV.13dpi[(Bd21.3.Endophytes.root.SV.13dpi$padj<0.05 & abs(Bd21.3.Endophytes.root.SV.13dpi$log2FoldChange) >= 1 ),]
nb_sigres_Bd21.3.Endophytes.root.SV.13dpi_DEG <- nrow(sigres_Bd21.3.Endophytes.root.SV.13dpi)
#181
write.table(x = sigres_Bd21.3.Endophytes.root.SV.13dpi,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_sigres_Bd21.3.Endophytes.root.SV.13dpi_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

#Bd21.3.Endophytes.shoot.SV.3dpi
Bd21.3.Endophytes.shoot.SV.3dpi <- Bd21.3.Endophytes.shoot.SV.3dpi[!is.na(Bd21.3.Endophytes.shoot.SV.3dpi$padj),]
sigres_Bd21.3.Endophytes.shoot.SV.3dpi_up <- Bd21.3.Endophytes.shoot.SV.3dpi [(Bd21.3.Endophytes.shoot.SV.3dpi $padj<0.05 & Bd21.3.Endophytes.shoot.SV.3dpi $log2FoldChange>= 1),]
nb_Bd21.3.Endophytes.shoot.SV.3dpi_up <- nrow(sigres_Bd21.3.Endophytes.shoot.SV.3dpi_up)
# 12
write.table(x = sigres_Bd21.3.Endophytes.shoot.SV.3dpi_up,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.Endophytes.shoot.SV.3dpi_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

sigres_Bd21.3.Endophytes.shoot.SV.3dpi_down <- Bd21.3.Endophytes.shoot.SV.3dpi[(Bd21.3.Endophytes.shoot.SV.3dpi$padj<0.05 & Bd21.3.Endophytes.shoot.SV.3dpi$log2FoldChange<= -1),]
nb_Bd21.3.Endophytes.shoot.SV.3dpi_down <- nrow(sigres_Bd21.3.Endophytes.shoot.SV.3dpi_down)
#14
write.table(x = sigres_Bd21.3.Endophytes.shoot.SV.3dpi_down,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.Endophytes.shoot.SV.3dpi_down.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
sigres_Bd21.3.Endophytes.shoot.SV.3dpi <- Bd21.3.Endophytes.shoot.SV.3dpi[(Bd21.3.Endophytes.shoot.SV.3dpi$padj<0.05 & abs(Bd21.3.Endophytes.shoot.SV.3dpi$log2FoldChange) >= 1 ),]
nb_sigres_Bd21.3.Endophytes.shoot.SV.3dpi_DEG <- nrow(sigres_Bd21.3.Endophytes.shoot.SV.3dpi)
#26
write.table(x = sigres_Bd21.3.Endophytes.shoot.SV.3dpi,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_sigres_Bd21.3.Endophytes.shoot.SV.3dpi_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
#Bd21.3.Endophytes.shoot.SV.6dpi
Bd21.3.Endophytes.shoot.SV.6dpi <- Bd21.3.Endophytes.shoot.SV.6dpi[!is.na(Bd21.3.Endophytes.shoot.SV.6dpi$padj),]
sigres_Bd21.3.Endophytes.shoot.SV.6dpi_up <- Bd21.3.Endophytes.shoot.SV.6dpi [(Bd21.3.Endophytes.shoot.SV.6dpi $padj<0.05 & Bd21.3.Endophytes.shoot.SV.6dpi $log2FoldChange>= 1),]
nb_Bd21.3.Endophytes.shoot.SV.6dpi_up <- nrow(sigres_Bd21.3.Endophytes.shoot.SV.6dpi_up)
# 4
write.table(x = sigres_Bd21.3.Endophytes.shoot.SV.6dpi_up,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.Endophytes.shoot.SV.6dpi_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

sigres_Bd21.3.Endophytes.shoot.SV.6dpi_down <- Bd21.3.Endophytes.shoot.SV.6dpi[(Bd21.3.Endophytes.shoot.SV.6dpi$padj<0.05 & Bd21.3.Endophytes.shoot.SV.6dpi$log2FoldChange<= -1),]
nb_Bd21.3.Endophytes.shoot.SV.6dpi_down <- nrow(sigres_Bd21.3.Endophytes.shoot.SV.6dpi_down)
#9
write.table(x = sigres_Bd21.3.Endophytes.shoot.SV.6dpi_down,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.Endophytes.shoot.SV.6dpi_down.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)


sigres_Bd21.3.Endophytes.shoot.SV.6dpi <- Bd21.3.Endophytes.shoot.SV.6dpi[(Bd21.3.Endophytes.shoot.SV.6dpi$padj<0.05 & abs(Bd21.3.Endophytes.shoot.SV.6dpi$log2FoldChange) >= 1 ),]
nb_sigres_Bd21.3.Endophytes.shoot.SV.6dpi_DEG <- nrow(sigres_Bd21.3.Endophytes.shoot.SV.6dpi)
#13
write.table(x = sigres_Bd21.3.Endophytes.shoot.SV.6dpi,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_sigres_Bd21.3.Endophytes.shoot.SV.6dpi_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
#Bd21.3.Endophytes.shoot.SV.13dpi
Bd21.3.Endophytes.shoot.SV.13dpi <- Bd21.3.Endophytes.shoot.SV.13dpi[!is.na(Bd21.3.Endophytes.shoot.SV.13dpi$padj),]
sigres_Bd21.3.Endophytes.shoot.SV.13dpi_up <- Bd21.3.Endophytes.shoot.SV.13dpi [(Bd21.3.Endophytes.shoot.SV.13dpi $padj<0.05 & Bd21.3.Endophytes.shoot.SV.13dpi $log2FoldChange>= 1),]
nb_Bd21.3.Endophytes.shoot.SV.13dpi_up <- nrow(sigres_Bd21.3.Endophytes.shoot.SV.13dpi_up)
# 9
write.table(x = sigres_Bd21.3.Endophytes.shoot.SV.13dpi_up,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.Endophytes.shoot.SV.13dpi_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

sigres_Bd21.3.Endophytes.shoot.SV.13dpi_down <- Bd21.3.Endophytes.shoot.SV.13dpi[(Bd21.3.Endophytes.shoot.SV.13dpi$padj<0.05 & Bd21.3.Endophytes.shoot.SV.13dpi$log2FoldChange<= -1),]
nb_Bd21.3.Endophytes.shoot.SV.13dpi_down <- nrow(sigres_Bd21.3.Endophytes.shoot.SV.13dpi_down)
#5
write.table(x = sigres_Bd21.3.Endophytes.shoot.SV.13dpi_down,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_Bd21.3.Endophytes.shoot.SV.13dpi_down.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)


sigres_Bd21.3.Endophytes.shoot.SV.13dpi <- Bd21.3.Endophytes.shoot.SV.13dpi[(Bd21.3.Endophytes.shoot.SV.13dpi$padj<0.05 & abs(Bd21.3.Endophytes.shoot.SV.13dpi$log2FoldChange) >= 1 ),]
nb_sigres_Bd21.3.Endophytes.shoot.SV.13dpi_DEG <- nrow(sigres_Bd21.3.Endophytes.shoot.SV.13dpi)
#14
write.table(x = sigres_Bd21.3.Endophytes.shoot.SV.13dpi,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/sigres_sigres_Bd21.3.Endophytes.shoot.SV.13dpi_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
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

pdf(file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/Bd21.fungi.leaves.PCA_DEG.pdf",width = 8,height = 8)
p <- ggplot(data=fungi_PCA, aes(x=Time, y=Number, fill=Categories)) +
  geom_bar(stat="identity", position=position_dodge())+ylim(0,2000)

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

pdf(file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/Bd21.3.leaves.bacteria.XT_DEG.pdf",width = 8,height = 8)
p <- ggplot(data=bac_XT, aes(x=Time, y=Number, fill=Categories)) +
  geom_bar(stat="identity", position=position_dodge())+ylim(0,2000)

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

pdf(file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/Bd21.3.leaves.virus.BSMV_DEG.pdf",width = 8,height = 8)
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

pdf(file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/Bd21.3.leaves.virus.MMMV_DEG.pdf",width = 8,height = 8)
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

pdf(file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/Bd21.3.leaves.virus.WSMV_DEG.pdf",width = 8,height = 8)
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

pdf(file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/Bd21.3.root.virus.BSMV_DEG.pdf",width = 8,height = 8)
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

pdf(file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/Bd21.3.root.virus.MMMV_DEG.pdf",width = 8,height = 8)
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

pdf(file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/Bd21.3.root.virus.WSMV_DEG.pdf",width = 8,height = 8)
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

##Endophytes.shoot.SV
##GGPLOt
library(ggplot2)
S_SV <- data.frame(Categories=rep(c("Up", "Down"), each=3),
                   Time=rep(c("3dpi", "6dpi", "13dpi"),2),
                   Number=c(nb_Bd21.3.Endophytes.shoot.SV.3dpi_up, nb_Bd21.3.Endophytes.shoot.SV.6dpi_up, 
                            nb_Bd21.3.Endophytes.shoot.SV.13dpi_up,
                            nb_Bd21.3.Endophytes.shoot.SV.3dpi_down, 
                            nb_Bd21.3.Endophytes.shoot.SV.6dpi_down, 
                            nb_Bd21.3.Endophytes.shoot.SV.13dpi_down))
head(S_SV)
str(S_SV)
S_SV$Time <- factor(S_SV$Time, levels = c("3dpi", "6dpi", "13dpi"))
par(mfrow = c(1,1))
par(mar=c(5,5,5,5))

pdf(file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/Bd21.3.Endophytes.shoot.SV_DEG.pdf",width = 8,height = 8)
p <- ggplot(data=S_SV, aes(x=Time, y=Number, fill=Categories)) +
  geom_bar(stat="identity", position=position_dodge())+ylim(0,60)

heatplot <- p + labs(title="Bd21.3.Endophytes.shoot.SV", 
                     x="Time (hours)", y = "The number of genes")+
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

##Bd21.3.Endophytes.root
##GGPLOt
library(ggplot2)
R_SV <- data.frame(Categories=rep(c("Up", "Down"), each=3),
                   Time=rep(c("3dpi", "6dpi", "13dpi"),2),
                   Number=c(nb_Bd21.3.Endophytes.root.SV.3dpi_up, nb_Bd21.3.Endophytes.root.SV.6dpi_up, 
                            nb_Bd21.3.Endophytes.root.SV.13dpi_up,
                            nb_Bd21.3.Endophytes.root.SV.3dpi_down, 
                            nb_Bd21.3.Endophytes.root.SV.6dpi_down, 
                            nb_Bd21.3.Endophytes.root.SV.13dpi_down))
head(R_SV)
str(R_SV)
R_SV$Time <- factor(R_SV$Time, levels = c("3dpi", "6dpi", "13dpi"))
par(mfrow = c(1,1))
par(mar=c(5,5,5,5))

pdf(file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/Bd21.3.Endophytes.root.SV_DEG.pdf",width = 8,height = 8)
p <- ggplot(data=R_SV, aes(x=Time, y=Number, fill=Categories)) +
  geom_bar(stat="identity", position=position_dodge())+ylim(0,200)

heatplot <- p + labs(title="Bd21.3.Endophytes.root.SV", 
                     x="Time (hours)", y = "The number of genes")+
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
#3953

#Bacteria:
bacteria.leaves.XT_gene_list <- union(rownames(sigres_Bd21.3.leaves.bacteria.XT.2.dpi),rownames(sigres_Bd21.3.leaves.bacteria.XT.4.dpi))
length(bacteria.leaves.XT_gene_list)
#3365
#virus-leaves-BSMV
Bd21.3.leaves.virus.BSMV_gene_list <- union(union(rownames(sigres_Bd21.3.leaves.virus.BSMV.3.dpi),rownames(sigres_Bd21.3.leaves.virus.BSMV.7.dpi)),rownames(sigres_Bd21.3.leaves.virus.BSMV.14.dpi))
length(Bd21.3.leaves.virus.BSMV_gene_list)
#4962
Bd21.3.leaves.virus.MMMV_gene_list <- union(union(rownames(sigres_Bd21.3.leaves.virus.MMMV.3.dpi),rownames(sigres_Bd21.3.leaves.virus.MMMV.7.dpi)),rownames(sigres_Bd21.3.leaves.virus.MMMV.14.dpi))
length(Bd21.3.leaves.virus.MMMV_gene_list)
#8734
Bd21.3.leaves.virus.WSMV_gene_list <- union(union(rownames(sigres_Bd21.3.leaves.virus.WSMV.3.dpi),rownames(sigres_Bd21.3.leaves.virus.WSMV.7.dpi)),rownames(sigres_Bd21.3.leaves.virus.WSMV.14.dpi))
length(Bd21.3.leaves.virus.WSMV_gene_list)
#5086
#root-BSMV
Bd21.3.root.virus.BSMV_gene_list <- union(union(rownames(sigres_Bd21.3.root.virus.BSMV.3.dpi),rownames(sigres_Bd21.3.root.virus.BSMV.7.dpi)),rownames(sigres_Bd21.3.root.virus.BSMV.14.dpi))
length(Bd21.3.root.virus.BSMV_gene_list)
#5470
#root-MMMV
Bd21.3.root.virus.MMMV_gene_list <- union(union(rownames(sigres_Bd21.3.root.virus.MMMV.3.dpi),rownames(sigres_Bd21.3.root.virus.MMMV.7.dpi)),rownames(sigres_Bd21.3.root.virus.MMMV.14.dpi))
length(Bd21.3.root.virus.MMMV_gene_list)
#10616
#root-WSMV
Bd21.3.root.virus.WSMV_gene_list <- union(union(rownames(sigres_Bd21.3.root.virus.WSMV.3.dpi),rownames(sigres_Bd21.3.root.virus.WSMV.7.dpi)),rownames(sigres_Bd21.3.root.virus.WSMV.14.dpi))
length(Bd21.3.root.virus.WSMV_gene_list)
#5837

#Endophytes.root.SV
Bd21.3.Endophytes.root.SV_gene_list <- union(union(rownames(sigres_Bd21.3.Endophytes.root.SV.3dpi),rownames(sigres_Bd21.3.Endophytes.root.SV.6dpi)),rownames(sigres_Bd21.3.Endophytes.root.SV.13dpi))
length(Bd21.3.Endophytes.root.SV_gene_list)
#265
#shoot
Bd21.3.Endophytes.shoot.SV_gene_list <- union(union(rownames(sigres_Bd21.3.Endophytes.shoot.SV.3dpi),rownames(sigres_Bd21.3.Endophytes.shoot.SV.6dpi)),rownames(sigres_Bd21.3.Endophytes.shoot.SV.13dpi))
length(Bd21.3.Endophytes.shoot.SV_gene_list)
#40

#top gene list:
top_gene_list <- union(union(union(union(union(union(union(union(union(fungi.leaves.PCA_gene_list,bacteria.leaves.XT_gene_list),
                                                                 Bd21.3.leaves.virus.BSMV_gene_list),
                                                           Bd21.3.leaves.virus.MMMV_gene_list),
                                                           Bd21.3.leaves.virus.WSMV_gene_list),
                                                     Bd21.3.root.virus.BSMV_gene_list),
                                               Bd21.3.root.virus.MMMV_gene_list),
                                         Bd21.3.root.virus.WSMV_gene_list),
                                   Bd21.3.Endophytes.root.SV_gene_list),
                             Bd21.3.Endophytes.shoot.SV_gene_list)
length(top_gene_list)
#17346

write.table(x = top_gene_list,
            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/top_DEG_17346_gene_list",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

#plot the genes with logfolderchanges >1
####reorder the columns
pdf(file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/heatmap_17346_DGE_pathogene_endophyte_noUB_del1rep_reordered.pdf",width = 15,height = 15)

par(mar=c(1,1,1,1))

palette <- colorRampPalette(c("red","white","blue"))(256)
reorder_rld <- assay(rld)[,c(4,5,6,1,2,3,10,11,12,7,8,9,16,17,18,13,14,15,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,67,101,102,68,103,104,105,106,113,114,115,107,108,109,116,117,118,110,111,112,119,120,121,122,123,124,131,132,133,125,126,127,134,135,136,128,129,130,137,138,139)]
head(reorder_rld)
heatmap.2( reorder_rld[row.names(assay(rld)) %in% top_gene_list,],  scale="row", Colv=FALSE, 
           trace="none", dendrogram='row',
           col = palette,
           lhei = c(0.8,7),
           ColSideColors = c( "Bd21.fungi.Mock"="#cc561e","Bd21.fungi.PCA" ="#aa2b1d", "Bd21-3.bacteria.Mock" = "#ef8d32","Bd21-3.bacteria.XT" = "#beca5c", "Bd21-3.virus.Mock"="#75cfb8","Bd21-3.virus.BSMV"="#bbdfc8","Bd21-3.virus.MMMV"="#f0e5d8","Bd21-3.virus.WSMV"="#ffc478","Bd21-3.Endophytes.Mock"="#822659","Bd21-3.Endophytes.SV"="#f8a1d1")[
             colData(rld)$treatment] )
dev.off()
###


pdf(file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/heatmap_4297_DGE_freeze_syl_novargen#.pdf",width = 15,height = 15)
par(mar=c(1,1,1,1))
#heatmap(assay(rld)[ topVarGenes, ],
#        labRow = "",
#        scale = "row")
palette <- colorRampPalette(c("red","white","blue"))(256)
head(assay(rld))
reorder_rld <- assay(rld)[,c(1,2,3,7,8,9,13,14,15,4,5,6,10,11,12,16,17,18)]
head(reorder_rld)
heatmap.2( reorder_rld[row.names(assay(rld)) %in% top_gene_list,], scale="row", Colv=FALSE, 
           trace="none", dendrogram='row',
           #Rowv=FALSE,
           #Colv=FALSE,
           col = palette,
           lhei = c(0.8,7),
           ColSideColors = c( leaf ="green", crown="gold")[
             colData(rld)$tissue ] )
heatmap.2( assay(rld)[row.names(assay(rld)) %in% top_gene_list,], scale="row", Colv=FALSE, 
           trace="none", dendrogram='row',
           #Rowv=FALSE,
           #Colv=FALSE,
           col = palette,
           lhei = c(0.8,7),
           ColSideColors = c( leaf ="green", crown="gold")[
             colData(rld)$tissue ] )
dev.off()
#plot the genes with logfolderchanges >2
pdf(file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/heatmap_9475_DGE_freeze_syl_novargen.pdf",width = 15,height = 15)
par(mar=c(1,1,1,1))
#heatmap(assay(rld)[ topVarGenes, ],
#        labRow = "",
#        scale = "row")
palette <- colorRampPalette(c("red","white","blue"))(256)

heatmap.2( reorder_rld[row.names(assay(rld)) %in% top_gene_list,], scale="row", Colv=FALSE, 
           trace="none", dendrogram='row',
           #Rowv=FALSE,
           #Colv=FALSE,
           col = palette,
           lhei = c(0.8,7),
           ColSideColors = c( cold="#1149F6", freeze="#11C5F6", recovery="gray")[
             colData(rld)$treatment ] )
dev.off()

###rlog normalization
pdf(file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/heatmap_110_strictDGE_cold_syl_rlog.pdf",width = 15,height = 15)

par(mar=c(1,1,1,1))
#heatmap(assay(rld)[ topVarGenes, ],
#        labRow = "",
#        scale = "row")
palette <- colorRampPalette(c("red","white","blue"))(256)

heatmap.2( assay(rlog_ds)[row.names(assay(rlog_ds)) %in% top_gene_list,], scale="row", Colv=FALSE, 
           trace="none", dendrogram='row', 
           col = palette,
           lhei = c(0.8,7),
           ColSideColors = c( CK="gray", cold="#11C5F6", freeze="#1149F6")[
             colData(rld)$treatment ] )
dev.off()
#rlog_ds from the plot, rlog is not great for normalization!!!!!!


pdf(file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/heatmap_110_strictDGE_cold_syl_logcount.pdf",width = 15,height = 15)

par(mar=c(1,1,1,1))
##Test another way
#selected <- rownames(sig)
output <- counts(dds,normalized=TRUE)[rownames(dds) %in% top_gene_list,]
output <- log2(output)
head(output)
#output[is.na(output)] = 0
#output_no_NA <- output[!is.na(resOrdered$padj),]
#na.dist <- function(x) {
#  t.dist <- dist(x)
#  t.dist <- as.matrix(t.dist)
#  t.limit <- 1.1*max(t.dist,na.rm=T)
#  t.dist[is.na(t.dist)] <- t.limit
##  t.dist <- as.dist(t.dist)
# return(t.dist)
#}
dist_no_na <- function(mat) {
  edist <- dist(mat)
  edist[which(is.na(edist))] <- max(edist, na.rm=TRUE) * 1.1 
  return(edist)
}
output[!is.finite(output)] <- 0

palette <- colorRampPalette(c("red","white","blue"))(256)
#docs <- dist( as.matrix(output), method = "euclidean")
#hclust_dist<- as.dist(docs)
#hclust_dist[is.na(hclust_dist)] <- 0
#hclust_dist[is.nan(hclust_dist)] <- 0
#sum(is.infinite(hclust_dist))  # THIS SHOULD BE 0
#h <- hclust(hclust_dist, "ward.D2")

heatmap.2(output,Colv=FALSE,
          distfun=dist_no_na,
          col = palette, scale="row", 
          dendrogram="row",
          trace="none",margin=c(4,6), 
          cexRow=0.5, cexCol=1, keysize=1 )
dev.off()
#check
Brasy8G092900.v1.1
counts[(row.names(counts) == "Brasy8G092900.v1.1"),]
#head(count)
assay(rlog_ds)[row.names(assay(rlog_ds)) == "Brasy8G092900.v1.1",]

###only focus on the root samples:
root <- counts[,c(31:66,122:139)]#only extract the root and see
colnames(root) <- c("Bd21-3.virus.root.Mock.3.dpi.1","Bd21-3.virus.root.Mock.3.dpi.2","Bd21-3.virus.root.Mock.3.dpi.3",
                    "Bd21-3.virus.root.BSMV.3.dpi.1","Bd21-3.virus.root.BSMV.3.dpi.2","Bd21-3.virus.root.BSMV.3.dpi.3",
                    "Bd21-3.virus.root.MMMV.3.dpi.1","Bd21-3.virus.root.MMMV.3.dpi.2","Bd21-3.virus.root.MMMV.3.dpi.3",
                    "Bd21-3.virus.root.WSMV.3.dpi.1","Bd21-3.virus.root.WSMV.3.dpi.2","Bd21-3.virus.root.WSMV.3.dpi.3",
                    "Bd21-3.virus.root.Mock.7.dpi.1","Bd21-3.virus.root.Mock.7.dpi.2","Bd21-3.virus.root.Mock.7.dpi.3",
                    "Bd21-3.virus.root.BSMV.7.dpi.1","Bd21-3.virus.root.BSMV.7.dpi.2","Bd21-3.virus.root.BSMV.7.dpi.3",
                    "Bd21-3.virus.root.MMMV.7.dpi.1","Bd21-3.virus.root.MMMV.7.dpi.2","Bd21-3.virus.root.MMMV.7.dpi.3",
                    "Bd21-3.virus.root.WSMV.7.dpi.1","Bd21-3.virus.root.WSMV.7.dpi.2","Bd21-3.virus.root.WSMV.7.dpi.3",
                    "Bd21-3.virus.root.Mock.14.dpi.1","Bd21-3.virus.root.Mock.14.dpi.2","Bd21-3.virus.root.Mock.14.dpi.3",
                    "Bd21-3.virus.root.BSMV.14.dpi.1","Bd21-3.virus.root.BSMV.14.dpi.2","Bd21-3.virus.root.BSMV.14.dpi.3",
                    "Bd21-3.virus.root.MMMV.14.dpi.1","Bd21-3.virus.root.MMMV.14.dpi.2","Bd21-3.virus.root.MMMV.14.dpi.3",
                    "Bd21-3.virus.root.WSMV.14.dpi.1","Bd21-3.virus.root.WSMV.14.dpi.2","Bd21-3.virus.root.WSMV.14.dpi.3",
                    "Bd21-3.Endophytes.root.Mock.3.dpi.1","Bd21-3.Endophytes.root.Mock.3.dpi.2","Bd21-3.Endophytes.root.Mock.3.dpi.3",
                    "Bd21-3.Endophytes.root.Mock.6.dpi.1","Bd21-3.Endophytes.root.Mock.6.dpi.2","Bd21-3.Endophytes.root.Mock.6.dpi.3",
                    "Bd21-3.Endophytes.root.Mock.13.dpi.1","Bd21-3.Endophytes.root.Mock.13.dpi.2","Bd21-3.Endophytes.root.Mock.13.dpi.3",
                    "Bd21-3.Endophytes.root.SV.3.dpi.1","Bd21-3.Endophytes.root.SV.3.dpi.2","Bd21-3.Endophytes.root.SV.3.dpi.3",
                    "Bd21-3.Endophytes.root.SV.6.dpi.1","Bd21-3.Endophytes.root.SV.6.dpi.2","Bd21-3.Endophytes.root.SV.6.dpi.3",
                    "Bd21-3.Endophytes.root.SV.13.dpi.1","Bd21-3.Endophytes.root.SV.13.dpi.2","Bd21-3.Endophytes.root.SV.13.dpi.3")
head(root)
metaData_root <- metaData[(metaData$tissue == "root"),]
rownames(metaData_root)
design_root <- model.matrix(~0 + group,metaData_root)
dds_root <- DESeq2::DESeqDataSetFromMatrix(countData = root, colData = metaData_root, design = design_root)
ddsTC_root <- DESeq(dds_root)
head(ddsTC)

rld_root <- vst( ddsTC_root )
head(rld_root)
pcaData_root <- plotPCA(rld_root, intgroup = c( "treatment", "timepoint"),returnData = TRUE) # vsd and plotPCA are part of DESeq2 package, nothing with my example below. 
head(pcaData_root)
percentVar_root <- round(100 * attr(pcaData_root, "percentVar")) 
pdf(file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/root_virus_endophyte.pdf",width = 20,height = 15)
par(mar=c(1,1,1,1))
ggplot(pcaData_root, aes(x = PC1, y = PC2, color = factor(treatment), shape = factor(timepoint))) + 
  geom_point(size =6, aes(fill=factor(treatment))) + 
  geom_point(size =6) + 
  scale_shape_manual(values=c(21,22,23,24,11,12)) + 
  xlab(paste0("PC1: ", percentVar_root[1], "% variance")) + 
  ylab(paste0("PC2: ", percentVar_root[2], "% variance")) + 
  ggtitle("PCA of all genes, no covariate adjusted")+
  theme(text = element_text(size=25))
dev.off()
#Then we can plot the other PCs 
# The function is the basically the same as https://github.com/mikelove/DESeq2/blob/master/R/plots.R.
# Inspired by https://www.biostars.org/p/243695/. 
# I added two variables, pp1 and pp2 to let user chose which principle to plot. pp1 and pp2 only take integer values.
#Test PC1 vs PC3
plotPCA.DESeqTransform = function(pp1=2, pp2=3, object, intgroup="condition", ntop=500, returnData=FALSE)
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
  
  ggplot(data=d, aes_string(x="PC2", y="PC3", color="group")) + geom_point(size=3) + 
    xlab(paste0("PC2: ",round(percentVar[1] * 100),"% variance")) +
    ylab(paste0("PC3: ",round(percentVar[2] * 100),"% variance"))
}

pcaData_root <- plotPCA.DESeqTransform(object=rld_root, intgroup = c("treatment", "timepoint"),returnData = TRUE) # vsd and plotPCA are part of DESeq2 package, nothing with my example below. 
head(pcaData_root)
head(pcaData_root)
pcaData_root$timepoint
percentVar_root <- round(100 * attr(pcaData_root, "percentVar")) 
pdf(file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/root_virus_endophyte_pc2_pc3_del1rep.pdf",width = 20,height = 15)
par(mar=c(1,1,1,1))
ggplot(pcaData_root, aes(x = PC1, y = PC2, color = factor(treatment), shape = factor(timepoint))) + 
  geom_point(size =6, aes(fill=factor(treatment))) + 
  geom_point(size =6) + 
  scale_shape_manual(values=c(21,22,23,24,11)) + 
  xlab(paste0("PC2: ", percentVar_root[2], "% variance")) + 
  ylab(paste0("PC3: ", percentVar_root[3], "% variance")) + 
  ggtitle("PCA of all genes, no covariate adjusted")+
  theme(text = element_text(size=25))
dev.off()

head(root)
colnames(counts)
colnames(root)
leaves <- counts[,c(67:100,102:122)]#only extract the leaves and shoot
head(leaves)
colnames(leaves)

metaData_leaves <- metaData[c(67:100,102:122),]
rownames (metaData_leaves)
design_leaves <- model.matrix(~0 + group,metaData_leaves)
dds_leaves <- DESeq2::DESeqDataSetFromMatrix(countData = leaves, colData = metaData_leaves, design = design_leaves)
ddsTC_leaves <- DESeq(dds_leaves)
head(ddsTC_leaves)

rld_leaves <- vst( ddsTC_leaves )
head(rld_leaves)
pcaData_leaves <- plotPCA(rld_leaves, intgroup = c( "treatment", "timepoint"),returnData = TRUE) # vsd and plotPCA are part of DESeq2 package, nothing with my example below. 
head(pcaData_leaves)
percentVar_leaves <- round(100 * attr(pcaData_leaves, "percentVar")) 
pdf(file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/leaves_virus_endophyte.pdf",width = 20,height = 15)
par(mar=c(1,1,1,1))
ggplot(pcaData_leaves, aes(x = PC1, y = PC2, color = factor(treatment), shape = factor(timepoint))) + 
  geom_point(size =6, aes(fill=factor(treatment))) + 
  geom_point(size =6) + 
  scale_shape_manual(values=c(21,22,23,24,11,12)) + 
  xlab(paste0("PC1: ", percentVar_leaves[1], "% variance")) + 
  ylab(paste0("PC2: ", percentVar_leaves[2], "% variance")) + 
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

pcaData_leaves <- plotPCA.DESeqTransform(object=rld_leaves, intgroup = c("treatment", "timepoint"),returnData = TRUE) # vsd and plotPCA are part of DESeq2 package, nothing with my example below. 
head(pcaData_leaves)
head(pcaData_leaves)
pcaData_leaves$timepoint
percentVar_leaves <- round(100 * attr(pcaData_leaves, "percentVar")) 
pdf(file = "/global/projectb/scratch/llei2019/CSP_Kranthi/clustering/leaves_virus_endophyte_pc1_pc3_del1rep.pdf",width = 20,height = 15)
par(mar=c(1,1,1,1))
ggplot(pcaData_leaves, aes(x = PC1, y = PC2, color = factor(treatment), shape = factor(timepoint))) + 
  geom_point(size =6, aes(fill=factor(treatment))) + 
  geom_point(size =6) + 
  scale_shape_manual(values=c(21,22,23,24,11)) + 
  xlab(paste0("PC1: ", percentVar_leaves[1], "% variance")) + 
  ylab(paste0("PC3: ", percentVar_leaves[3], "% variance")) + 
  ggtitle("PCA of all genes, no covariate adjusted")+
  theme(text = element_text(size=25))
dev.off()
