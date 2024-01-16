library(DiffBind)
library(pheatmap)
a=dba.load("atac_counts_Jenny")

h=dba.plotHeatmap(a,plot=FALSE)
h
Meta=read.table("MetadataJennPrunned.txt",header=TRUE)
Meta
Diagnosis=as.factor(Meta$Condition)
Batch=as.factor(Meta$Treatment - 1)
Age=cale(as.numeric(Meta$Factor))

Annot=data.frame(Batch,Diagnosis,Age)

rownames(Annot)=rownames(h)


#pdf("../results/03_ClusterSamples_Consensus/SampleClusterReadsPrunnedJenny.pdf")
#pheatmap(h,annotation=Annot)

