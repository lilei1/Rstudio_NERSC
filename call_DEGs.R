#Written by Li Lei 2022-02-28
#This script is for calling the DGE for tissue-specific expression, especially the Bsylvaticum!!!!
#John is interested in if callus have its own specific expression and then deceide that if we need to prepare the 
#callus to do the gene annotations

source("https://bioconductor.org/biocLite.R")
biocLite("edgeR", dependencies=TRUE, INSTALL_opts = c('--no-lock'))
#biocLite("httr", dependencies=TRUE, INSTALL_opts = c('--no-lock'))
biocLite("Hmisc")
biocLite("DESeq2")

library(ggplot2)
library(edgeR)
library(DESeq2)
library(tidyr)
library(ggrepel)
library(EnhancedVolcano)

#packageVersion("rlang")
counts <- read.delim("/global/cfs/cdirs/plantbox/llei/bscratch/llei2019/Brachypodium/Sylvaticum/RNAseq/counts.txt", row.names = 1)
head(counts)
nrow(counts)
#36927
metaData <- read.delim("/global/cfs/cdirs/plantbox/llei/bscratch/llei2019/Github/Check_callus_gene_exp/Bsyl/data/tissues_meta.txt",header = T)
head(metaData)
nrow(metaData)
#9
#extract the column I needed!!!
col.num <- which( colnames(counts) %in% metaData$lib_code)
row.num <- which( metaData$lib_code %in% colnames(counts) )
counts <- counts[, col.num]
mata <- metaData[row.num,]
ncol(counts)
head(counts)
colnames(counts) <- metaData$code_name
head(counts)
counts$avg_callus <- (counts$`jv-amy-callus-1`+counts$`jv-amy-callus-2`+counts$`jv-amy-callus-3`)/3
head(counts)
callus_specific <- counts[(counts$avg_callus >0 & 
                          counts$`jv-seedling-leaf-1` == 0 & 
                          counts$`jv-seedling-crown-2` == 0 &
                          counts$`jv-seedling-root-1` == 0 &
                          counts$`jv-seedling-leaf-2` == 0 &
                          counts$`jv-hyrdo-root-3` == 0 &
                          counts $`jv-palea-1` == 0 &
                          counts$`jv-stamen-1` == 0 &
                          counts$`jv-pistil-1` == 0 &
                          counts$`jv-hyrdo-shoot-1` == 0 &
                          counts$`jv-hyrdo-shoot-2` == 0 &
                          counts$`jv-hyrdo-shoot-3` == 0 &
                          counts$`jv-hyrdo-root-1` == 0 &
                          counts$`jv-hyrdo-root-2` == 0 &
                          counts$`mature-leaf-1` == 0 &
                          counts$`mature-leaf-2` == 0 &
                          counts$`mature-leaf-3` == 0 &
                          counts$`mature-root-1` == 0 &
                          counts$`mature-root-2` == 0 &
                          counts$`mature-root-3` == 0 &
                          counts$`mature-crown-2` == 0 &
                          counts$`mature-crown-3` == 0 &
                          counts$`mature-crown-4` == 0 &
                          counts$`mature-nod-2` == 0 &
                          counts$`mature-nod-4` == 0 &
                          counts$`mature-internod-tip-1` == 0 &
                          counts$`mature-internod-tip-2` == 0 &
                          counts$`mature-internod-mid-2` == 0 &
                          counts$`mature-internod-mid-3` == 0 &
                          counts$`mature-internod-base-1` == 0 &
                          counts$`mature-internod-base-3` == 0 &
                          counts$`mature-internod-base-4` == 0),]
head(callus_specific)
nrow(callus_specific)
#125
nrow(callus_specific[(callus_specific$avg_callus >1),])
#74

write.table(x = callus_specific,
            file = "/global/cfs/cdirs/plantbox/llei/bscratch/llei2019/Github/Check_callus_gene_exp/Bsyl/results/callus_specific_exp_gene_125.txt",
            sep = "\t",
            eol = "\n",
            col.names = TRUE,
            row.names = TRUE
)
