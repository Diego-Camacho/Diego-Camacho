#Loading Hu et al 2021 compartment calls

Hu=read.csv("project-geschwind/Diego_Internship_Summary/Hi-C/data/3Dchromatin_NeuN.plus.neg_Compartments.csv")

#taking out NA's
Hu$NeuN.[is.na(Hu$NeuN.)]=0
Hu$NeuN..1[is.na(Hu$NeuN..1)]=0


library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
transcript <- transcriptsBy(txdb, "gene")
gene <- reduce(transcript)
genes=unlist(gene)


#Computing gene density for Hu et all calls as HiCT
Hu_bed=GRanges(seqnames=Hu$chr,ranges=IRanges(start=Hu$start,end=Hu$end))
gene_dens=countOverlaps(Hu,genes)


#Focusing now NeuNplus cells for example purposes

corCC=NULL
corSC=NULL

for (chr in unique(Hu$chr)){
	tmp=Hu[Hu$chr==chr,]
	bed=GRanges(seqnames=tmp$chr,ranges=IRanges(start=tmp$start,end=tmp$end))
	gene.density=countOverlaps(bed,genes)
	pscore=tmp$NeuN..1 #positive ones
	gd.pos <- sum(gene.density[which(pscore>=0)])
    gd.neg <- sum(gene.density[which(pscore<0)])

    if (gd.pos > gd.neg)
        cc <- ifelse(pscore>0, "A", "B")
    else{
        cc <- ifelse(pscore<0, "A", "B")
        pscore <- -pscore}
    corCC=c(corCC,cc)
    corSC=c(corSC,pscore)

}
#filtering out bins with no score

HuPos=Hu[corSC!0,1:3]
corCC=corCC[corSC!=0]
corSC=corSC[corSC!=0]

strand=NULL
strand[corCC=="A"]="+"
strand[corCC=="B"]="-"

HuNeuNPlus=GRanges(seqnames=HuPos$chr,strand=strand,ranges=IRanges(start=HuPos$start,end=HuPos$end))

#Asses overlaps with our data using strand as sencod level factor for overlap as A or B compartment

On=data.frame(findOverlaps(HuNeuNPlus,OurNeuNPlus))
#~15k 
Off=data.frame(findOverlaps(HuNeuNPlus,OurNeuNPlus,ignore.strand=TRUE))
#~23k

#Re assesment of UnShared bins.
#Subjecthits stands for Hu bins
#Extracting the bins that are shared at range levels but not at A/B compartment level.
#~9k ranges

UnSharedIdx=Off$queryHits[! Off$queryHits %in% On$queryHits]
UnSharedScores=corSC[UnSharedIdx]

UnSharedScores=corSC[UnSharedIdx]

UnSharedGenDens=countOverlaps(HuNeuNPlusUnShared,genes)


#Gen density is no correctly correlated with the positive bins from Hu et al, see the Gviz plot for further evidence, flipp the signs.

sum(UnSharedGenDens[UnSharedScores>0])
#[1]5344

sum(UnSharedGenDens[UnSharedScores<0])
#[1]5979

strand2=NULL
strand2[UnSharedScores>0]="-"
strand2[UnSharedScores<0]="+"

HuNeuNPlusCorrected=GRanges(seqnames=HuPos$chr[UnSharedIdx],strand=strand2,ranges=IRanges(start=HuPos$start[UnSharedIdx],end=HuPos$end[UnSharedIdx]))

findOverlaps(HuNeuNPlusCorrected,OurNeuNPlus)

#Plus 9886 overlaped bins yield= 25581 total overlap


#Once the changes are made, craft and quantify the true overlaps:
#The corrected Hu values is on the data folder

p=draw.pairwise.venn(27417, 26218, 15695+9886 , category = c( "Diego compartment calls \n both A and B","Hu et al compartment calls \n both A and B"), lty = rep("blank", 
    2), fill = c("light blue", "red"), alpha = rep(0.5, 2),scaled = TRUE)

o=draw.pairwise.venn(27417, 25958, 9297+9178 + 6819 , category = c( "Diego compartment calls \n both A and B","Hu et al compartment calls \n both A and B"), lty = rep("blank", 
    2), fill = c("light blue", "red"), alpha = rep(0.5, 2),scaled = TRUE)


