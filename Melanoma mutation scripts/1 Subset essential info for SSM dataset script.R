
## Description: This script exists simply to reduce the overall size of the simple somatic mutation dataset by keeping only essential info (e.g. genome coordinates, mutation type,etc)

out<-"output_directory/"

SSM<-read.table(file = "file_path/simple_somatic_mutation.open.MELA-AU.tsv.gz",header = T,sep = "\t")
SSM<-SSM[,c(1:17)]
SSM<-SSM[which(SSM$mutation_type =="single base substitution"),]

write.csv(SSM, paste0(out, "SSM_essential_info.MELA_AU.csv"))


## END ##