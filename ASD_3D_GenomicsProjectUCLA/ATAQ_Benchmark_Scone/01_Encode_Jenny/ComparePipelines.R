library(ChIPseeker)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)

library(GenomicRanges)
#Change paths to select the peaks you desire

Jenny_peaks <- read.delim("/u/home/d/dr2camac/project-geschwind/Diego_Internship_Summary/ATAC-seq/data/2GGP-1-ATAC_S10_L001_peaks.narrowPeak", header=FALSE)

Encode_peaks <- read.delim("/u/home/d/dr2camac/project-geschwind/Diego_Internship_Summary/ATAC-seq/data/rep1-pr1_vs_rep1-pr2.idr0.05.bfilt.narrowPeak", header=FALSE)

Jenny_gr=GRanges(seqnames = Jenny_peaks$V1, ranges = IRanges(start = Jenny_peaks$V2,end = Jenny_peaks$V3))

Encode_gr=GRanges(seqnames = Encode_peaks$V1, ranges = IRanges(start = Encode_peaks$V2,end = Encode_peaks$V3))



EncAtt=annotatePeak(Encode_gr,TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene)
JenAtt=annotatePeak(Jenny_gr,TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene)
OveAtt=annotatePeak(gr3,TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene)

par(mfrow=c(1,2))
plotAnnoPie(EncAtt)
plotAnnoPie(JenAtt)

