library(GenomicRanges)

## Description: This script averages together repair signals from plus and minus strands.


out<-"output_directory/"


## CPD FILES

# paths to 1Mb binned CPD XRseq files
path<-file.path("path_to_merged_replicate_XRseq_files//")
PatternList<-c("1h","4h","8h","16h","1d","2d")

for (i in 1:length(PatternList)) {
  FileNames<-list.files(path, pattern = PatternList[i])
  
  load(file = paste0(path, FileNames[1]))
  Minus<-Pooled
  load(file = paste0(path, FileNames[2]))
  Plus<-Pooled
  
  Pooled<-cbind(Plus, Minus[,6])
  Pooled$Plus_Minus<-0
  Pooled$Plus_Minus<-rowMeans(Pooled[,c(6:7)], dims = 1)
  Pooled<-Pooled[,c(1:5,8)]
  colnames(Pooled)[6]<-"dStrand_Mean"
  save(Pooled ,file = paste0(out, "CPD_",  PatternList[i], "_Rep&Strand_Mean_1Mb.RData"))
}


## 64PP files ##

# paths to 1Mb binned 64PP XRseq files
path<-file.path("path_to_merged_replicate_XRseq_files//")
PatternList<-c("5min","20min","1h","2h","4h")

for (i in 1:length(PatternList)) {
  FileNames<-list.files(path, pattern = PatternList[i])
  
  load(file = paste0(path, FileNames[1]))
  Minus<-Pooled
  load(file = paste0(path, FileNames[2]))
  Plus<-Pooled
  
  Pooled<-cbind(Plus, Minus[,6])
  Pooled$Plus_Minus<-0
  Pooled$Plus_Minus<-rowMeans(Pooled[,c(6:7)], dims = 1)
  Pooled<-Pooled[,c(1:5,8)]
  colnames(Pooled)[6]<-"dStrand_Mean"
  
  save(Pooled ,file = paste0(out, "XR64_",  PatternList[i], "_Rep&Strand_Mean_1Mb.RData"))
}

## END ##