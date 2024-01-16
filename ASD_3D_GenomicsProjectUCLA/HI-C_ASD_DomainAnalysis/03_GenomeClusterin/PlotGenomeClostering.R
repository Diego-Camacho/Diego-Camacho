library(pheatmap)

setwd("/u/home/d/dr2camac/project-geschwind/Diego_Internship_Summary/Hi-C/data/100kCompartments")
#setwd("u/home/d/dr2camac/project-geschwind/Diego_Internship_Summary/Hi-C/data/250kCompartments")



chr=c(1:22,"X","Y")


Metadata=("../HiC_QC.txt")

listFiles=list.files(pattern=paste0("_chr",chr[1],"_pc1_genegr.txt"))

tags=(unlist(strsplit(listFiles[1],"_"))[1:3])


sam=paste(tags[1],tags[2],tags[3],sep="_")

df=read.csv(listFiles[1], sep="")["score"][,1]

Age=Metadata$Age[Metadata$ID == tags[1]][1]
Region=tags[2]

for (i in 1:49){
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

for (i in 1:23){

listFiles=list.files(pattern=paste0("_chr",chr[i+1],"_pc1_genegr.txt"))
print(i+1)

#tags=(unlist(strsplit(listFiles[1],"_"))[1:3])


#sam=paste(tags[1],tags[2],tags[3],sep="_")

df1=read.csv(listFiles[1], sep="")["score"][,1]

#Age=Metadata$Age[Metadata$ID == tags[1]][1]
#Region=tags[2]

for (j in 1:49){
	p=read.csv(listFiles[j +1], sep="")["score"][,1]
	df1=cbind(df1,p)}


dim(df1)

Gen=rbind(Gen,df1)
}




C=cor(Gen)

Annot=data.frame(cbind(Age,Region))
rownames(Annot)=colnames(C)


pdf("ClusterSamplesHiCGenome250kb.pdf",width=10,height=8))
pheatmap(C,annotation=Annot)
