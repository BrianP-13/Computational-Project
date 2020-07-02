library(GenomicRanges)

## Description: This scripts consolidates all processed XRseq measurements for each leasion types into one master table.
# Tables will contain averaged replicate and stranded XRseq measurements


out<-"output_directory/"
path<-"path_to_final_XRseq_files/"


######
## Add repair data for 64PP lesion type

# Import 1Mb binned UV lesion tables
MT_path<-"file_path/"
MT<-read.csv(paste0(MT_path,"UVlesion_mastertable_1Mb.csv"),header = T)
IP64_MT<-MT[,c(2:6,8)]

# Import merged XRseq data for 64PP and add to table
PatternList<-c("5min","20min","1h","2h","4h")

for (i in 1:5) {
  load(file = paste0(path, "/NHF1_64_merged_strands_FINAL/XR64_",  PatternList[i], "_Rep&Strand_Mean_1Mb.RData"))
  IP64_MT<-cbind(IP64_MT, Pooled[,6])
  colnames(IP64_MT)[6+i]<-as.character(PatternList[i])
}
save(IP64_MT,file = paste0(out,"IP64_lesion_XRseq_mastertable.RData"))


######
## Add repair data for CPD lesion type

# Import 1Mb binned UV lesion tables
MT_path<-"file_path/"
MT<-read.csv(paste0(MT_path,"UVlesion_mastertable_1Mb.csv"),header = T)
CPD_MT<-MT[,c(2:6,7)]

# Import merged XRseq data for 64PP and add to table
PatternList<-c("1h","4h","8h","16h","1d","2d")

for (i in 1:6) {
  load(file = paste0(path, "/NHF1_CPD_merged_strands_FINAL/CPD_",  PatternList[i], "_Rep&Strand_Mean_1Mb.RData"))
  CPD_MT<-cbind(CPD_MT, Pooled[,6])
  colnames(CPD_MT)[6+i]<-as.character(PatternList[i])
}
save(CPD_MT,file = paste0(out,"CPD_lesion_XRseq_mastertable.RData"))


## END ##