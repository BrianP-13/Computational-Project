library(rtracklayer)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)


## Description: This is the next step in measuring repair in 1Mb regions. This script averages together both replicates for each time point together.



out<-"output_directory/"


######
## CPD FILES

# paths to 1Mb binned CPD XRseq files
path<-file.path("path_to_binned_XRseq_files//")
PatternList<-c("1h","4h","8h","16h","1d","2d")

for (i in 1:length(PatternList)) {
  FileNames<-list.files(path, pattern = PatternList[i])
  # Rep1+Rep2: Minus strands
  load(file = paste0(path, FileNames[1]))
  Rep1<-ba
  load(file = paste0(path, FileNames[3]))
  Rep2<-ba
  
  Pooled<-cbind(Rep1, Rep2[,6])
  Pooled$R1_R2<-0
  Pooled$R1_R2<-rowMeans(Pooled[,c(6:7)], dims = 1)
  Pooled<-Pooled[,c(1:5,8)]
  colnames(Pooled)[6]<-paste0("mean_Rep")
  
  save(Pooled,file = paste0(out, "CPD_",  PatternList[i], "_mergedReps_Minus_1Mb.RData"))
  
  # Rep1+Rep2: Plus strands
  load(file = paste0(path, FileNames[2]))
  Rep1<-ba
  load(file = paste0(path, FileNames[4]))
  Rep2<-ba
  
  Pooled<-cbind(Rep1, Rep2[,6])
  Pooled$R1_R2<-0
  Pooled$R1_R2<-rowMeans(Pooled[,c(6:7)], dims = 1)
  Pooled<-Pooled[,c(1:5,8)]
  colnames(Pooled)[6]<-"mean_Rep"
  
  save(Pooled,file = paste0(out, "CPD_",  PatternList[i], "_mergedReps_Plus_1Mb.RData"))
}


######
## 64PP FILES

# paths to 1Mb binned CPD XRseq files
path<-file.path("path_to_binned_XRseq_files//")
PatternList<-c("5min","20min","1h","2h","4h")

for (i in 1:length(PatternList)) {
  FileNames<-list.files(path, pattern = PatternList[i])
  
  # Rep1+Rep2: Minus strands
  load(file = paste0(path, FileNames[1]))
  Rep1<-ba
  load(file = paste0(path, FileNames[3]))
  Rep2<-ba
  
  Pooled<-cbind(Rep1, Rep2[,6])
  Pooled$R1_R2<-0
  Pooled$R1_R2<-rowMeans(Pooled[,c(6:7)], dims = 1)
  Pooled<-Pooled[,c(1:5,8)]
  colnames(Pooled)[6]<-paste0("mean_Rep")
  
  save(Pooled,file = paste0(out, "XR64PP_",  PatternList[i], "_mergedReps_Minus_1Mb.RData"))
  
  # Rep1+Rep2: Plus strands
  load(file = paste0(path, FileNames[2]))
  Rep1<-ba
  load(file = paste0(path, FileNames[4]))
  Rep2<-ba
  
  Pooled<-cbind(Rep1, Rep2[,6])
  Pooled$R1_R2<-0
  Pooled$R1_R2<-rowMeans(Pooled[,c(6:7)], dims = 1)
  Pooled<-Pooled[,c(1:5,8)]
  colnames(Pooled)[6]<-paste0("mean_Rep")
  
  save(Pooled,file = paste0(out, "XR64PP_",  PatternList[i], "_mergedReps_Plus_1Mb.RData"))
}

## END ##