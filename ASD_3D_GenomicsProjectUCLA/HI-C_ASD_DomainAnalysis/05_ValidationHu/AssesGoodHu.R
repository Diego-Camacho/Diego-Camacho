library(GenomicRanges)
load("my_work_space.RData")
#Now that the scores have been re assessed and flipped the strand would be not taked into account.
Hu=read.table("/u/home/d/dr2camac/project-geschwind/Diego_Internship_Summary/Hi-C/data/Hu2021CompartmentCallsCorrectedScores.txt")

#taking out NA's
Hu$NeuN.[is.na(Hu$NeuN.)]=0
Hu$NeuN..1[is.na(Hu$NeuN..1)]=0

HuNeg=Hu[ Hu$NeuN.!=0,1:4]
HuPos=Hu[ Hu$NeuN..1!=0,c(1:3,5)]


HuNeuNPlus=GRanges(seqnames=HuPos$chr,strand=strandPos,ranges=IRanges(start=HuPos$start,end=HuPos$end))

HuNeuNNega=GRanges(seqnames=HuNeg$chr,strand=strandNeg,ranges=IRanges(start=HuNeg$start,end=HuNeg$end))


setwd("/u/project/geschwind/dr2camac/Hi-C/CompartmentAna/100kSepCA")



chr=c(1:22,"X")

#chr=c(3)
Metadata=read.csv("../HiC_QC.txt")

listFiles=list.files(pattern=paste0("_chr",chr[1],"_pc1_genegr.txt"))[-46]

tags=(unlist(strsplit(listFiles[1],"_"))[1:3])


sam=paste(tags[1],tags[2],tags[3],sep="_")

df=read.csv(listFiles[1], sep="")["score"][,1]

Age=Metadata$Age[Metadata$ID == tags[1]][1]
Region=tags[2]

for (i in 1:48){
	p=read.csv(listFiles[i +1], sep="")["score"][,1]
	df=cbind(df,p)
	tags=(unlist(strsplit(listFiles[i+1],"_"))[1:3])
	s=paste(tags[1],tags[2],tags[3],sep="_")
	sam=c(sam,s)
	a=Metadata$Age[Metadata$ID == tags[1]][1]
	r=tags[2]
	Age=c(Age,a)
	Region=c(Region,r)
}

colnames(df)=sam

Gen=df

for (i in 1:22){

listFiles=list.files(pattern=paste0("_chr",chr[i+1],"_pc1_genegr.txt"))[-46]
print(i+1)

#tags=(unlist(strsplit(listFiles[1],"_"))[1:3])


#sam=paste(tags[1],tags[2],tags[3],sep="_")

df1=read.csv(listFiles[1], sep="")["score"][,1]

#Age=Metadata$Age[Metadata$ID == tags[1]][1]
#Region=tags[2]

for (j in 1:48){
	p=read.csv(listFiles[j +1], sep="")["score"][,1]
	df1=cbind(df1,p)}


dim(df1)

Gen=rbind(Gen,df1)
}

rownames(Gen)=rownames(Our)


Gen=Gen[rowSums(Gen)!=0, c(grep("124",colnames(Gen)),grep("395",colnames(Gen)),grep("105",colnames(Gen)))]

OurPos=Gen[,grep("POS",colnames(Gen))]
OurPos=OurPos[]
OurNeg=Gen[,grep("NEG",colnames(Gen))]



tmp=unlist(strsplit(rownames(OurPos),"_"))

OurPosGr=GRanges(seqnames=tmp[seq(1,length(tmp),3)],ranges=IRanges(start=as.numeric(tmp[seq(2,length(tmp),3)]),end=as.numeric(tmp[seq(3,length(tmp),3)])))

tmp=unlist(strsplit(rownames(OurNeg),"_"))
OurPosNe=GRanges(seqnames=tmp[seq(1,length(tmp),3)],ranges=IRanges(start=as.numeric(tmp[seq(2,length(tmp),3)]),end=as.numeric(tmp[seq(3,length(tmp),3)])))

findOverlaps(OurPosGr,HuNeuNPlus)
findOverlaps(OurPosNg,HuNeuNNega)

##pngs VennDiagrams

p=draw.pairwise.venn(27417, 26218, 25581 , category = c( "Diego compartment calls \n both A and B","Hu et al compartment calls \n both A and B"), lty = rep("blank", 
    2), fill = c("light blue", "red"), alpha = rep(0.5, 2),scaled = TRUE)

o=draw.pairwise.venn(27417, 25958, 25294 , category = c( "Diego compartment calls \n both A and B","Hu et al compartment calls \n both A and B"), lty = rep("blank", 
    2), fill = c("light blue", "red"), alpha = rep(0.5, 2),scaled = TRUE)




