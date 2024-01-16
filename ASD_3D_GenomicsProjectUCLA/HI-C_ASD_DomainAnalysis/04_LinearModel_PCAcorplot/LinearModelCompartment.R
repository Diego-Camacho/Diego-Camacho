##Set the desired chromosomes, resolution and load QC values. 


setwd("/u/home/d/dr2camac/project-geschwind/Diego_Internship_Summary/Hi-C/data/100kCompartments")
#setwd("u/home/d/dr2camac/project-geschwind/Diego_Internship_Summary/Hi-C/data/250kCompartments")


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


Metadata=read.table("../Metadata_QC.txt")

NeuN=as.factor(Metadata$NeuN)
Age=as.numeric(Metadata$Age)
Region=as.factor(Metadata$Region)
Depth=as.numeric(Metadata$R1)
Valid=as.numeric(Metadata$Valid)
Dup=as.numeric(Metadata$Dup)
Cis_short=as.numeric(Metadata$cis_short)
Batch=as.factor(Metadata$Batch)
Sex=as.factor(Metadata$Sex)
Trans=as.numeric(Metadata$trans)
Cis_long=as.numeric(Metadata$cis_long)

#######################
## PCA FOR VARIABLES ##
#######################

#Transpose Gen table to have the samples as variables. 
Gen_PC=prcomp(t(Gen))
summary(Gen_PC)
#Selected the first 16 PC that explaing 90% of variance

QC_PC=prcom(data.frame(Depth,Valid,Dup,Cis_short,Cis_long,Trans))
summary(QC_PC)
#Selected the first 2 PC that explaing 90% of variance

Corr_DF=data.frame(Gen_PC$x[,1:16],as.numeric(NeuN),Age,as.numeric(Region),Batch,as.numeric(Sex),Depth,Valid,Dup,Cis_short,Cis_long,Trans,QC_PC$x[,1:2])

###change names to indicate the % of variance that is explained

library(corrplot)

#pdf("CorrPlot.pdf")

corrplot(cor(Corr_DF),"ellipse",'grey50')


NeuN=NeuN[-46]
Age=Age[-46]
Region=Region[-46]
Sex=Sex[-46]

Beta = p = matrix(nrow = nrow(Gen), ncol =5)

for (i in 1:nrow(Gen)) {
Score=Gen[i,]
dft=data.frame(Score,NeuN,Age,Region,Sex)
Res = lm(Score ~ NeuN + Age + Region + Sex, dft)
q=summary(Res)
Beta[i,] =q$coefficients[,1][2:6]
p[i,] = q$coefficients[,4][2:6]
}


###drop NAs
nap=is.na(p[,1])
Beta=data.frame(Beta[!nap,])
p=data.frame(p[!nap,])


colnames(Beta) =c("NeuN","Age2","Age3","RegionSTR","Sex")
colnames(p)=c("NeuN","Age2","Age3","RegionSTR","Sex")


hist(p$Age)
#~ 15,000 <0.2
fdr_Age = p.adjust(p$Age,"fdr")
#~ 1,100 < 0.05

hist(p$NeuN)
#~ 20,000 <0.2
fdr_NeuN = p.adjust(p$NeuN,"fdr")
#~ 10,00 < 0.05

#see histograms attached.
