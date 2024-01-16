args<- commandArgs(trailingOnly=TRUE)
file <- args[1]
chr<- as.character(args[2])
res= as.character(args[3])
window=args[4]

setwd("/u/project/geschwind/dr2camac/Hi-C/TADCalling")

library(HiCDCPlus)
print(paste0("hg19_",res,"kb_GATC_",chr,"_bintolen.txt.gz"))
gi_list<-generate_bintolen_gi_list(
  bintolen_path=paste0("hg19_",res,"kb_GATC_",chr,"_bintolen.txt.gz"))

gi_list=add_hicpro_matrix_counts(gi_list,absfile_path=paste0("/u/project/geschwind/dr2camac/Hi-C/TADCalling/POS_",res,"000_abs.bed"),matrixfile_path=file)


gi_list<-HiCDCPlus(gi_list) #HiCDCPlus_parallel runs in parallel across ncores


tads<-gi_list_topdom(gi_list,window.size = window)

print("saving tads")

tags=(unlist(strsplit(basename(file),"_"))[1:3])


sam=paste0(tags[1],"_",tags[2])

print(args)

#print(paste0(sam,"TADds_chr",chr,"_",res,"000_win5.rds"))
name=paste0(sam,"/TADds_chr_",chr,res,"win_",as.character(window),".RDS")
print(name)
saveRDS(tads,name)





