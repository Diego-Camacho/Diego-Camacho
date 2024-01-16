#This is how to create a DIffBind object from at metadata table, with the specific colnames. This specific script is for generating the counts from a GR object contaning the consensus reads.

library(DiffBind)

ASC=read.delim("MetadataJennPrunned.txt")
ATAC=dba(sampleSheet=ASC)
ATAC$config$RunParallel=FALSE

cons=readRDS("ConsPeaksJennyPrunned.Rdata")

ATAC <- dba.count(ATAC, peaks=cons,score=DBA_SCORE_READS)

dba.save(ATAC,"atac_counts_Jenny")
write.table( ATAC$binding,"MatrixCountRawCJenny.txt")
