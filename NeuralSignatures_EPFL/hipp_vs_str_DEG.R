#Downloading, indexing and coutning of aligned reads

#gdown https://drive.google.com/######
#for bam in *.bam; do base=$(echo $bam | cut -f1,2 -d"."); echo -e "....$base......\n"; samtools sort $bam > ${base}.sorted.bam; samtools index ${base}.sorted.bam > ${base}.sorted.idx.bam; done 
#featureCounts -M -t gene -T 4 -s 0 --largestOverlap -O -a Mus_musculus.GRCm38.85.gtf -o Results/Raw_LO_Ocounts_85unstrgenemut.txt  Aligments/*sorted.bam 
#featureCounts -M -t gene -T 4 -s 0 --largestOverlap -O -a Mus_musculus.GRCm38.100.gtf -o Results/Raw_LO_Ocounts_100unstrgenemut.txt  Aligments/*sorted.bam  


#Data wrangling for edgeR analysis 
library(edgeR)
library(pheatmap)


Raw_100 <- read.table("C:/Users/HP/Desktop/Graff_Burns/Neuroepigenetics and Memory/DifExpAna/Raw_LO_Ocounts_100unstrgenemut.txt", comment.char="#",header = T)
Raw_85 <- read.table("C:/Users/HP/Desktop/Graff_Burns/Neuroepigenetics and Memory/DifExpAna/Raw_LO_Ocounts_85unstrgenemut.txt", comment.char="#",header = T)

Names<-c("Hpp_Hi_Sho_r1","Str_Hi_Sho_r1","Hpp_Hi_Sho_r2","Str_Hi_Sho_r2",
         "Hpp_Vh_Sho_r1","Str_Vh_Sho_r1","Hpp_Vh_Sho_r2","Str_Vh_Sho_r2",
         "Hpp_Hi_Con_r1","Str_Hi_Con_r1","Hpp_Hi_Con_r2","Str_Hi_Con_r2",
         "Hpp_Vh_Con_r1","Str_Vh_Con_r1","Hpp_Vh_Con_r2","Str_Vh_Con_r2",
         "Hpp_Vh_Con_r3","Str_Vh_Con_r3","Hpp_Vh_Con_r4","Str_Vh_Con_r4",
         "Hpp_Hi_Con_r3","Str_Hi_Con_r3","Hpp_Hi_Con_r4","Str_Hi_Con_r4",
         "Hpp_Vh_Sho_r3","Str_Vh_Sho_r3","Hpp_Vh_Sho_r4","Str_Vh_Sho_r4",
         "Hpp_Hi_Sho_r3","Str_Hi_Sho_r3","Hpp_Hi_Sho_r4","Str_Hi_Sho_r4")

colnames(Raw_100)[7:length(colnames(Raw_100))]<- Names
colnames(Raw_85)[7:length(colnames(Raw_85))]<- Names

row.names(Raw_100)<-Raw_100$Geneid
Raw_100<-Raw_100[7:length(colnames(Raw_100))]

row.names(Raw_85)<-Raw_85$Geneid
Raw_85<-Raw_85[7:length(colnames(Raw_85))]

head(Raw_85)

keep100<-rowSums(cpm(Raw_100)>=1)>=2
keep85<-rowSums(cpm(Raw_85)>=1)>=2

table(keep100)
Counts100<-Raw_100[keep100,]
Counts85<-Raw_85[keep85,]

dim(Counts100)[1]

# Expression summary in log2 CPM

cpm_log <- cpm(Counts85,log = T)
medial_log<- apply(cpm_log, 1, median)
blue_red<- colorRampPalette(c("blue","red"))
hist(medial_log,main="Median expression in CPM",xlab="Median log2",ylab = "Frequency",col = blue_red(20))

pheatmap(cor(cpm_log))

Names2<- Names[!Names=="Str_Vh_Con_r1"]
Counts85<-Counts85[,Names2]


#creating edgeR object

group <-factor(sub("_r.+","",colnames(Counts85)))
group


dge85<- DGEList(counts = Counts85, group = group)
dge85<-calcNormFactors(dge85)
dge85$samples


#Md plots per sample
par(mfrow =c(1,2))
for (i in c(1:31)){
  plotMD(cpm(dge85,log = TRUE),column = i)
  grid(col = "blue")
  abline(h=0,col="red",lty=2,lwd=2)
}


#PCA 

cols<-rainbow(length(levels(group)))
names(cols)<-levels(group)

plotMDS(dge85, col=cols[dge85$samples$group], labels=colnames(dge85$counts), main="MDS counts85")

design85<- model.matrix(~0+dge85$samples$group)
colnames(design85)<-levels(dge85$samples$group)
design85

#Dispersion estimate
dge85<- estimateDisp(dge85,design = design85,robust = T)
plotBCV(dge85)

#Differential expression analysis 3 conditions to assess
contrast<-makeContrasts(
  "Activated_effect"="Hpp_Hi_Sho - Hpp_Vh_Con",
  "Activated_primed"="Hpp_Hi_Sho - Hpp_Hi_Con",
  "Activated_nonprimed"="Hpp_Hi_Sho - Hpp_Vh_Sho",
  levels = dge85$design
)
contrast

contHS_VC<- contrast[,1]

contHS_HC<- contrast[,2]

contHS_VS<- contrast[,3]

glm<-glmQLFit(dge85,dispersion = dge85$trended.dispersion,robust = T)
glm2<-glmQLFit(dge85,dispersion = dge85$tagwise.dispersion)
glm3<-glmQLFit(dge85,dispersion = dge85$common.dispersion)

LFTHS_VC<-glmQLFTest(glm,contrast =  contHS_VC)
LFTHS_HC<-glmQLFTest(glm2,contrast =  contHS_HC)
LFTHS_VS<-glmQLFTest(glm,contrast =  contHS_VS)

dt.HS_VC<-decideTestsDGE(LFTHS_VC,adjust.method = "BH",p.value = 0.05,lfc = 0)
deGenes<- rownames(LFTHS_VC)[dt.HS_VC !=0]
plotSmear(LFTHS_VC,de.tags = deGenes

dt.HS_HC<-decideTestsDGE(LFTHS_HC,adjust.method = "BH",p.value = 0.05,lfc = 0)
deGenes<- rownames(LFTHS_HC)[dt.HS_HC !=0]
plotSmear(LFTHS_HC,de.tags = deGenes)

dt.HS_VS<-decideTestsDGE(LFTHS_VS,adjust.method = "BH",p.value = 0.05,lfc = 0)
deGenes<- rownames(LFTHS_VS)[dt.HS_VS !=0]
plotSmear(LFTHS_VS,de.tags = deGenes)

#GEtting gene names for significant results       
library(biomaRt)

DGE.HS_VC<-topTags(LFTHS_VC,n=Inf)$table
mart <- useDataset(dataset = "mmusculus_gene_ensembl", useMart("ensembl"))
Genes<- getBM(filters= "ensembl_gene_id", attributes=
c("ensembl_gene_id","external_gene_name"),values=sort(row.names(DGE.HS_VC)),mart= mart)

DGE.HS_VC$ensembl_gene_id<-row.names(DGE.HS_VC)
DGE.HS_VC<-merge(DGE.HS_VC,Genes,by="ensembl_gene_id")
write.csv(DGE.HS_VC,"DGE.HS_VC.csv")
head(DGE.HS_VC)

DGE.HS_HC<-topTags(LFTHS_HC,n=Inf)$table
DGE.HS_HC$ensembl_gene_id<-row.names(DGE.HS_HC)
DGE.HS_HC<-merge(DGE.HS_HC,Genes,by="ensembl_gene_id")
write.csv(DGE.HS_HC,"DGE.HS_HC.csv")
head(DGE.HS_HC)

DGE.HS_VS<-topTags(LFTHS_VS,n=Inf)$table
DGE.HS_VS$ensembl_gene_id<-row.names(DGE.HS_VS)
DGE.HS_VS<-merge(DGE.HS_VS,Genes,by="ensembl_gene_id")
write.csv(DGE.HS_VS,"DGE.HS_VS.csv")
head(DGE.HS_VS)

#Ontoly annotation
          
library(GO.db)
library(org.Mm.eg.db)
library(clusterProfiler)
library(AnnotationHub)
library(ensembldb)
library(tidyverse)

UpDGE.HS_VC<-DGE.HS_VC[DGE.HS_VC$logFC>=1 & DGE.HS_VC$FDR <0.05,]$ensembl_gene_id
UpDGE.HS_HC<-DGE.HS_HC[DGE.HS_HC$logFC>=1 & DGE.HS_HC$FDR <0.05,]$ensembl_gene_id
UpDGE.HS_VS<-DGE.HS_VS[DGE.HS_VS$logFC>=1 & DGE.HS_VS$FDR <0.05,]$ensembl_gene_id
DownDGE.HS_VC<-DGE.HS_VC[DGE.HS_VC$logFC<=-1 & DGE.HS_VC$FDR <0.05,]$ensembl_gene_id
DownDGE.HS_VS<-DGE.HS_VS[DGE.HS_VS$logFC<=-1 & DGE.HS_VS$FDR <0.05,]$ensembl_gene_id


no_dup<-which(duplicated(UpDGE.HS_VC)== FALSE)
UpDGE.HS_VC<-UpDGE.HS_VC[no_dup]

allOE_genes<- as.character(row.names(Counts85))
egoHS_VC<-enrichGO(gene = as.character(UpDGE.HS_VC),
                   universe = allOE_genes,
                   keyType = "ENSEMBL",
                   OrgDb =org.Mm.eg.db,
                   ont = "BP",
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.05,
                   readable = T)

Summ.ego.Up.HS_VC<-data.frame(egoHS_VC)

dotplot(object=egoHS_VC,x="count",orderBy="BgRatio",showCategory=50) + theme(axis.text.y = element_text(size = 7.7))+ ggtitle("Up regulated GO enrichment HS_VC.")

#Correlation bewtween tissues 
          
df<-data.frame(DGE.HppHi_Vh$logFC,DGE.StrHi_Vh$logFC)
colnames(df)<-c("Hippocampus","Striatum")

ggplot(df,aes(x=Hippocampus,y=Striatum)) + geom_point() + geom_smooth(method = "lm",formula =  y~x) + ggtitle("LFC correlation Str vs Hpp")
