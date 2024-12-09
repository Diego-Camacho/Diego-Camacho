library(GO.db)
library(org.Mm.eg.db)
library(clusterProfiler)
library(AnnotationHub)
library(ensembldb)
library(tidyverse)
library(DEXSeq)
library(ggplot2)
library(gprofiler2)

#Loading of files and context per sample. 
#Using DEXSeq for exon differential usage

countFiles<-  list.files(pattern="count", full.names=TRUE)
flattenedFile<- list.files(".",pattern = "MM",full.names = T)
p<-rep(c(TRUE,FALSE),16)
Names<-c("Hpp_Hi_Sho_r1","Str_Hi_Sho_r1","Hpp_Hi_Sho_r2","Str_Hi_Sho_r2",
         "Hpp_Vh_Sho_r1","Str_Vh_Sho_r1","Hpp_Vh_Sho_r2","Str_Vh_Sho_r2",
         "Hpp_Hi_Con_r1","Str_Hi_Con_r1","Hpp_Hi_Con_r2","Str_Hi_Con_r2",
         "Hpp_Vh_Con_r1","Str_Vh_Con_r1","Hpp_Vh_Con_r2","Str_Vh_Con_r2",
         "Hpp_Vh_Con_r3","Str_Vh_Con_r3","Hpp_Vh_Con_r4","Str_Vh_Con_r4",
         "Hpp_Hi_Con_r3","Str_Hi_Con_r3","Hpp_Hi_Con_r4","Str_Hi_Con_r4",
         "Hpp_Vh_Sho_r3","Str_Vh_Sho_r3","Hpp_Vh_Sho_r4","Str_Vh_Sho_r4",
         "Hpp_Hi_Sho_r3","Str_Hi_Sho_r3","Hpp_Hi_Sho_r4","Str_Hi_Sho_r4")
Names<-Names[p]
n<-c(TRUE,TRUE,FALSE,FALSE,FALSE,FALSE,TRUE,TRUE,TRUE,TRUE,FALSE,FALSE,FALSE,FALSE,TRUE,TRUE)
Names<-Names[n]
condition<-sapply(strsplit(Names,"_"),"[[",2)

countFiles<-countFiles[p][n]
sampleTable = data.frame(
   row.names = Names,
   condition=condition)
sampleTable


dxd = DEXSeqDataSetFromHTSeq(
   countFiles,
   sampleData=sampleTable,
   design= ~ sample + exon + condition:exon,
   flattenedfile=flattenedFile )


dxd=estimateSizeFactors(dxd)
sampleAnnotation(dxd)

dxd = estimateDispersions(dxd)
plotDispEsts( dxd )

dxd = testForDEU(dxd)
dxd = estimateExonFoldChanges(dxd)
dxr1 = DEXSeqResults( dxd )
dxr1

#Visualazing results

plotMA(dxr1,cex=0.8)

#Saving significant results

DEUsage<-dxr1[which(dxr1$padj<0.05),]
DEUGenes<-dxr1$groupID[which(tapply( dxr1$padj < 0.05, dxr1$groupID, any )== TRUE)]
write.table(DEUGenes,"DifExonUsageGenesHS_VC.txt",quote = FALSE,row.names = FALSE, col.names = FALSE)


#Annotating regults with gprofiler 2

gostres <- gost(query = DEUGenesAlone, 
                organism = "mmusculus", ordered_query = TRUE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                measure_underrepresentation = FALSE, evcodes = FALSE, 
                user_threshold = 0.05, correction_method = "g_SCS", 
                domain_scope = "annotated", custom_bg = NULL, 
                numeric_ns = "", sources = NULL, as_short_link = FALSE)
gostplot(gostres, capped = TRUE, interactive = TRUE)


#BDNF exon usage plot
plotDEXSeq( dxr1, "ENSMUSG00000048482", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
