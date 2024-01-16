args<- commandArgs(trailingOnly=TRUE)
file=args[1]

setwd("/u/project/geschwind/dr2camac/Hi-C/TADCalling")

s=unlist(strsplit(file,"_"))
res=s[4]
name=paste(s[1],s[2],s[3],s[4],sep="_")

print(name)
library(GENOVA)

sample=load_contacts(
signal_path = file,
indices_path = paste0("POS_",res,"_abs.bed"),
sample_name = "trial",
)

Isul=insulation_score(sample,window=11)

write.table(Isul$insula_score,paste0("../InsulationScore/",name,"_Insul_w11.txt"))
