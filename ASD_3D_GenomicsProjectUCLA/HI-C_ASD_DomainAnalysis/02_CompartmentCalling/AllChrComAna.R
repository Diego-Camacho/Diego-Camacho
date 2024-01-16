# 	
# module load R/3.6.0 # default 3.4.0 does not have barplot(formula)

#R
## -------------------Start R script here -------------------
args = commandArgs(trailingOnly=TRUE)
SETWD = args[1]
OUT = args[2]
res=args[3]
str(res)
file.names=args[4]

library("HiTC")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

## To install these libraries:
# The first time running this section, use qrsh -l to request a node, and run this on that node instead of a login node. Cuz a login node would never install the package.
#source("http://bioconductor.org/biocLite.R")
#biocLite("BSgenome.Hsapiens.UCSC.hg19")
## Personal library created: ~/R/x86_64-pc-linux-gnu-library/3.2
#biocLite("HiTC")
## You will need to install into personal library.
## If prompted to install updates for other packages, type in no.

setwd(SETWD)

## Create GenomicRanges object of hg19 gene positions
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
transcript <- transcriptsBy(txdb, "gene")
gene <- reduce(transcript)
# Combine the first transcript of all genes into one GenomicRanges object
genes=unlist(gene)

## create function of compartment calling and ploting function
compartment <- function(x,OUT,Base,chr) {
  #pm <- getPearsonMap(x)
  #print("Plotting pearson correlation matrix ...")
  #pdf(paste0(OUT,Base,"_",chr,"_pm.pdf"))
  #mapC(HTClist(pm), maxrange=1, col.pos=c("black","red"), col.neg=c("black","blue"))
  #dev.off()
  print("Start PCA ...")
  #pr <- pca.hic(pm,normPerExpected = T, npc=2,asGRangesList = T）#,gene.gr = genes)
  genes_chr=genes[seqnames(genes)==chr]
  pr <- pca.hic(x,normPerExpected = T, npc=1,asGRangesList = T, gene.gr = genes_chr)
  pr_df <- as.data.frame(pr)
  pr_df$mid=(pr_df$start+pr_df$end-1)/2
  pr_df$score=ifelse(is.na(pr_df$score),0,pr_df$score) # replace NA with 0 to make space for plot
  pdf(paste0(OUT,Base,"_",chr,"_pc_genegr.pdf"))
  barplot(score~mid,data = pr_df,col="darkblue",border = NA)
  #pc1 <- pr_df[pr_df$group_name=="PC1",]
  #barplot(score~mid,data = pc1,col="darkblue",border = NA)
  #pc2 <- pr_df[pr_df$group_name=="PC2",]
  #barplot(score~mid,data = pc2,col="darkblue",border = NA)
  dev.off()
  #print("Writing PC1 and PC2 ...")
  #write.table(pc1, paste0(OUT,Base,"_",chr,"_pc1_genegr.txt"))
  #write.table(pc2, paste0(OUT,Base,"_",chr,"_pc2_genegr.txt"))
  print("Writing PC1 ...")
  write.table(pr_df, paste0(OUT,Base,"_",chr,"_pc1_genegr.txt"))
}

## If I need to read in pc1 later:
#pc1_n = read.table("B4337_NeuNn_chr1_pc1.txt")

## Import Hi-C data, call compartments for each chromosome
#file.names <- dir(path=".", pattern ="F*_rawdata_500000_iced.matrix")
#file.names <- dir(path=".", pattern ="B*_iced.matrix")
for (j in 1:length(file.names)) {
  print(file.names[j])
  #Base=substr(file.names[j],1, nchar(file.names[j])-nchar("_rawdata_100000_CAIC_sparse.matrix") )
  #Base=substr(file.names[j],1, nchar(file.names[j])-nchar("_rawdata__iced.matrix")-nchar(res))
  Base=substr(file.names[j],1, nchar(file.names[j])-nchar("_iced_.matrix")-nchar(res))
  iced=importC(file.names[j],xgi.bed=paste0("rawdata_",res,"_abs.bed"))
  ##chr1
  #detail(iced$chr1chr1)
  x = extractRegion(iced$chr1chr1, chr="chr1",from = 1, to = 249250621)
  chr="chr1"
  compartment(x,OUT,Base,chr)
  ##chr2
  #detail(iced$chr2chr2)
  x = extractRegion(iced$chr2chr2, chr="chr2",from = 1, to = 243199373)
  chr="chr2"
  compartment(x,OUT,Base,chr)
  ##chr3
  #detail(iced$chr3chr3)
  x = extractRegion(iced$chr3chr3, chr="chr3",from = 1, to = 198022430)
  chr="chr3"
  compartment(x,OUT,Base,chr)
  ##chr4
  #detail(iced$chr4chr4)
  x = extractRegion(iced$chr4chr4, chr="chr4",from = 1, to = 191154276)
  chr="chr4"
  compartment(x,OUT,Base,chr)
  ##chr5
  #detail(iced$chr5chr5)
  x = extractRegion(iced$chr5chr5, chr="chr5",from = 1, to = 180915260)
  chr="chr5"
  compartment(x,OUT,Base,chr)
  ##chr6
  #detail(iced$chr6chr6)
  x = extractRegion(iced$chr6chr6, chr="chr6",from = 1, to = 171115067)
  chr="chr6"
  compartment(x,OUT,Base,chr)
  ##chr7
  #detail(iced$chr7chr7)
  x = extractRegion(iced$chr7chr7, chr="chr7",from = 1, to = 159138663)
  chr="chr7"
  compartment(x,OUT,Base,chr)
  ##chr8
  #detail(iced$chr8chr8)
  x = extractRegion(iced$chr8chr8, chr="chr8",from = 1, to = 146364022)
  chr="chr8"
  compartment(x,OUT,Base,chr)
  ##chr9
  #detail(iced$chr9chr9)
  x = extractRegion(iced$chr9chr9, chr="chr9",from = 1, to = 141213431)
  chr="chr9"
  compartment(x,OUT,Base,chr)
  ##chr10
  #detail(iced$chr10chr10)
  x = extractRegion(iced$chr10chr10, chr="chr10",from = 1, to = 135534747)
  chr="chr10"
  compartment(x,OUT,Base,chr)
  ##chr11
  #detail(iced$chr11chr11)
  x = extractRegion(iced$chr11chr11, chr="chr11",from = 1, to = 135006516)
  chr="chr11"
  compartment(x,OUT,Base,chr)
  ##chr12
  #detail(iced$chr12chr12)
  x = extractRegion(iced$chr12chr12, chr="chr12",from = 1, to = 133851895)
  chr="chr12"
  compartment(x,OUT,Base,chr)
  ##chr13
  #detail(iced$chr13chr13)
  x = extractRegion(iced$chr13chr13, chr="chr13",from = 1, to = 115169878)
  chr="chr13"
  compartment(x,OUT,Base,chr)
  ##chr14
  #detail(iced$chr14chr14)
  x = extractRegion(iced$chr14chr14, chr="chr14",from = 1, to = 107349540)
  chr="chr14"
  compartment(x,OUT,Base,chr)
  ##chr15
  #detail(iced$chr15chr15)
  x = extractRegion(iced$chr15chr15, chr="chr15",from = 1, to = 102531392)
  chr="chr15"
  compartment(x,OUT,Base,chr)
  ##chr16
  #detail(iced$chr16chr16)
  x = extractRegion(iced$chr16chr16, chr="chr16",from = 1, to = 90354753)
  chr="chr16"
  compartment(x,OUT,Base,chr)
  ##chr17
  #detail(iced$chr17chr17)
  x = extractRegion(iced$chr17chr17, chr="chr17",from = 1, to = 81195210)
  chr="chr17"
  compartment(x,OUT,Base,chr)
  ##chr18
  #detail(iced$chr18chr18)
  x = extractRegion(iced$chr18chr18, chr="chr18",from = 1, to = 78077248)
  chr="chr18"
  try(compartment(x,OUT,Base,chr))
  ##chr19
  #detail(iced$chr19chr19)
  x = extractRegion(iced$chr19chr19, chr="chr19",from = 1, to = 59128983)
  chr="chr19"
  compartment(x,OUT,Base,chr)
  ##chr20
  #detail(iced$chr20chr20)
  x = extractRegion(iced$chr20chr20, chr="chr20",from = 1, to = 63025520)
  chr="chr20"
  compartment(x,OUT,Base,chr)
  ##chr21
  #detail(iced$chr21chr21)
  x = extractRegion(iced$chr21chr21, chr="chr21",from = 1, to = 48129895)
  chr="chr21"
  compartment(x,OUT,Base,chr)
  ##chr22
  #detail(iced$chr22chr22)
  x = extractRegion(iced$chr22chr22, chr="chr22",from = 1, to = 51304566)
  chr="chr22"
  compartment(x,OUT,Base,chr)
  ##chrX
  #detail(iced$chrXchrX)
  x = extractRegion(iced$chrXchrX, chr="chrX",from = 1, to = 155270560)
  chr="chrX"
  compartment(x,OUT,Base,chr)
  ##chrY
  #detail(iced$chrYchrY)
  x = extractRegion(iced$chrYchrY, chr="chrY",from = 1, to = 59373566)
  chr="chrY"
  try(compartment(x,OUT,Base,chr))
}

## -------------------End R script here -------------------

