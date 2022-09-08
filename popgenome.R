install.packages("PopGenome")
library(PopGenome)
#inputFiles <- readData("/global/projectb/scratch/llei2019/uk/FileTransfer-NQvYhGqSJEY4ATP9/test_chr1H/vcf", format="VCF", gffpath="/global/projectb/scratch/llei2019/uk/FileTransfer-NQvYhGqSJEY4ATP9/test_chr1H/gff", big.data = TRUE, FAST = TRUE)
#GENOME.class <- concatenate.regions(inputFiles)
GENOME.class <- readData("/global/projectb/scratch/llei2019/uk/FileTransfer-NQvYhGqSJEY4ATP9/test_chr1H/vcf", format="VCF", big.data = TRUE, FAST = TRUE)
get.sum.data(GENOME.class)
show.slots(GENOME.class)
GENOME.class@n.sites
#146077751 This is the last SNP position
GENOME.class@n.biallelic.sites
#33636
 
#GENOME.class@region.data
#GENOME.class <- set.synnonsyn(inputFiles, ref.chr="/global/projectb/scratch/llei2019/uk/FileTransfer-NQvYhGqSJEY4ATP9/test_chr1H/fasta/chr1H_part1.fasta")
#GENOME.class@region.data@CodingSNPS  

cultivars <- scan("~/bscratch/uk/FileTransfer-NQvYhGqSJEY4ATP9/populationLists/cultivars.txt", what="character")
landraces <- scan("~/bscratch/uk/FileTransfer-NQvYhGqSJEY4ATP9/populationLists/landraces.txt", what="character")
wildBarleys <- scan("~/bscratch/uk/FileTransfer-NQvYhGqSJEY4ATP9/populationLists/wildBarleys.txt", what="character")
GENOME.class <- set.populations(GENOME.class,list(cultivars, landraces, wildBarleys, diploid=FALSE))
names(GENOME.class@populations)[1] <- "cultivars"
names(GENOME.class@populations)[2] <- "landraces"
names(GENOME.class@populations)[3] <- "wildBarleys"

GENOME.class <- F_ST.stats(GENOME.class, mode = "nucleotide", detail = TRUE, FAST = TRUE)
slide10k <- sliding.window.transform(GENOME.class,width=10000,jump=10000, type=2, whole.data=FALSE)
length(slide10k@region.names)
#[1] 14607
slide10k <- F_ST.stats(slide10k, mode = "nucleotide", detail = TRUE)
#fst_slide <- sliding.window.transform(GENOME.class,width=10,jump=10, type=1, whole.data=FALSE)
#fst_slide <- sliding.window.transform(GENOME.class,width=10,jump=10, type=1, whole.data=FALSE)
#fst_slide <- F_ST.stats(fst_slide, mode = "nucleotide", detail = TRUE)
fst_slide <- sliding.window.transform(GENOME.class,width=1,jump=1, type=1, whole.data=FALSE)
fst_slide <- F_ST.stats(fst_slide, mode = "nucleotide", detail = TRUE)
#
#fst_slide@nucleotide.F_ST
pairwise.FST <- t(fst_slide@nuc.F_ST.pairwise)
fst_df <- as.data.frame(pairwise.FST)
fst_df <- fst_df[,c(1,2,4)]
colnames(fst_df) <- c("cultivars_vs_landraces", "cultivars_vs_wildbarleys", "landraces_vs_wildbarleys")
head(fst_df)
fst_df$bin <- 1:length(fst_slide@region.names)
nrow(fst_df)
# set chromosome size
chr1part <- 146077751
# set window size and window jump
window_size <- 10
window_jump <- 10
# use seq to find the start points of each window
window_start <- seq(from = 1, to = chr1part, by = window_jump)
# add the size of the window to each start point 
window_stop <- window_start + window_size
windows <- data.frame(start = window_start, stop = window_stop, 
                      mid = window_start + (window_stop-window_start)/2)
mid = window_start + (window_stop-window_start)/2
head(windows)
ggplot(fst_df, aes(mid/10^6, fst_df[,3])) +
  geom_point(size = 0.5)  +
  ggtitle("Pairwise Fst -landraces_vs_wildbarleys") +
  xlab("bp") +
  ylab("Pairwise Fst") +
  theme_bw(base_size = 20) + ylim(0,1)

a <- ggplot(fst_df, aes(mid/10^6, value, colour = stat)) + geom_line()
a <- a + facet_grid(stat~., scales = "free_y")
a <- a + xlab("Position (Mb)")
a + theme_light() + theme(legend.position = "none")




GENOME.class <- diversity.stats(GENOME.class, pi = TRUE, keep.site.info = TRUE)
get.diversity(GENOME.class)
GENOME.class@nuc.diversity.within 
get.diversity(GENOME.class)[[1]]
get.diversity(GENOME.class)[[2]]
get.diversity(GENOME.class)[[3]]

 gffFile<- readData("~/bscratch/uk/FileTransfer-NQvYhGqSJEY4ATP9/gff_1H", format="gff", big.data = TRUE, FAST = TRUE)

