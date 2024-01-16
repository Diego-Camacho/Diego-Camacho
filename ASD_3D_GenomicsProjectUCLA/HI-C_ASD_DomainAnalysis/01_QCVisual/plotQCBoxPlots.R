library(ggplot2)
trialRanges <- read.csv("/u/home/d/dr2camac/project-geschwind/Diego_Internship_Summary/Hi-C/data/trialRanges.xlsx", header=FALSE)
 # View(trialRanges)
  
  HiC_QC= t(trialRanges)
  colnames(HiC_QC)=c("valid_interaction",HiC_QC[1,][-1])
  
  HiC_QC=HiC_QC[2:101,]
  HiC_QC=data.frame(HiC_QC)
  
  HiDS <- read.table("/u/home/d/dr2camac/project-geschwind/Diego_Internship_Summary/Hi-C/data/HiDS.txt", quote="\"", comment.char="")
  HI_tissue <- read.table("/u/home/d/dr2camac/project-geschwind/Diego_Internship_Summary/Hi-C/data/HI_tissue.tx", quote="\"", comment.char="")
  Hi_R1= read.table("/u/home/d/dr2camac/project-geschwind/Diego_Internship_Summary/Hi-C/data/total_R1", quote="\"", comment.char="")
  Hi_NeuN <- read.table("/u/home/d/dr2camac/project-geschwind/Diego_Internship_Summary/Hi-C/data/NeuNstat.txt", quote="\"", comment.char="")
  
  
  HiC_QC$ID=HiDS$V1
  HiC_QC$Tissue=HI_tissue$V1
  HiC_QC$R1=Hi_R1$V1
  HiC_QC$NeuN=toupper(Hi_NeuN$V1)
  HiC_QC$Age=rep("0",100)
  
   HiC_QC$Age[which(HiC_QC$ID== "HSB131")]= "0-5yo"
   HiC_QC$Age[which(HiC_QC$ID== "HSB143")]= "0-5yo"
   HiC_QC$Age[which(HiC_QC$ID== "HSB121")]= "0-5yo"
   HiC_QC$Age[which(HiC_QC$ID== "HSB120")]= "5-10yo"
   HiC_QC$Age[which(HiC_QC$ID== "HSB124")]= "10-18yo"
   HiC_QC$Age[which(HiC_QC$ID== "HSB105")]= "10-18yo"
   HiC_QC$Age[which(HiC_QC$ID== "HSB404")]= "5-10yo"
   HiC_QC$Age[which(HiC_QC$ID== "HSB141")]= "5-10yo"
   HiC_QC$Age[which(HiC_QC$ID== "HSB395")]= "10-18yo"
   
  HiC_QC$Group=paste(HiC_QC$Age,HiC_QC$Tissue,sep = "_")

setwd("diego/UCLA/Hi-C/")

pdf("/u/home/d/dr2camac/project-geschwind/Diego_Internship_Summary/Hi-C/results/01_QCVisual", width =10, height=8)


ggplot(HiC_QC,aes(x=Group,y=as.numeric(HiC_QC$valid_interaction)/as.numeric(HiC_QC$R1), fill=Group)) + geom_boxplot() + geom_jitter(color="black", size=0.4, alpha=0.9,width = 1) + xlab("") + ylab("Valid interactions %") + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())

ggplot(HiC_QC,aes(x=Group,y=as.numeric(HiC_QC$valid_interaction_rmdup)/as.numeric(HiC_QC$R1), fill=Group)) + geom_boxplot() + geom_jitter(color="black", size=0.4, alpha=0.9,width = 1) + xlab("") + ylab("Valid interactions_rmdup %") + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())

ggplot(HiC_QC,aes(x=Group,y=as.numeric(HiC_QC$trans_interaction)/as.numeric(HiC_QC$R1), fill=Group)) + geom_boxplot() + geom_jitter(color="black", size=0.4, alpha=0.9,width = 1) + xlab("") + ylab("Trans interactions") + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())

ggplot(HiC_QC,aes(x=Group,y=as.numeric(HiC_QC$cis_shortRange)/as.numeric(HiC_QC$R1), fill=Group)) + geom_boxplot() + geom_jitter(color="black", size=0.4, alpha=0.9,width = 1) + xlab("") + ylab("Cis-short interactions") + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())

ggplot(HiC_QC,aes(x=Group,y=as.numeric(HiC_QC$cis_longRange)/as.numeric(HiC_QC$R1), fill=Group)) + geom_boxplot() + geom_jitter(color="black", size=0.4, alpha=0.9,width = 1) + xlab("") + ylab("Cis-long interactions") + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())

ggplot(HiC_QC,aes(x=Group,y=(as.numeric(HiC_QC$valid_interaction)-as.numeric(HiC_QC$valid_interaction_rmdup))/as.numeric(HiC_QC$R1), fill=Group)) + geom_boxplot() + geom_jitter(color="black", size=0.4, alpha=0.9,width = 1) + xlab("") + ylab("Dup interactions") + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())

write.table(HiC_QC,"/u/home/d/dr2camac/project-geschwind/Diego_Internship_Summary/Hi-C/data/Hi_C_QC_Group.txt")

dev.off()

