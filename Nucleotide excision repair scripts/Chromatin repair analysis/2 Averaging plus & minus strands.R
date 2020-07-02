library(GenomicRanges)

## Description: This script combines together repair signals from plus and minus strands that were binned for IMR90 chromatin state regions


out<-"output_directory/"



######
## paths to 64PP XRseq bw files
path<-"path_to_merged_replicate_XRseq_files/"
PatternList<-c("5min","20min","1h","2h","4h")

for (i in 1:5) {
  load(file = paste0(path, "XR64PP_",  PatternList[i], "_R1_R2_PLUS_core15.RData"))
  load(file = paste0(path, "XR64PP_",  PatternList[i], "_R1_R2_MINUS_core15.RData"))
  
  Merged<-Pooled_plus[,c(1:6)]
  Merged$Plus_ave<-Pooled_plus[,9]
  Merged$Minus_ave<-Pooled_minus[,9]
  Merged$Strand_mean<-rowMeans(Merged[,c(7:8)],dims = 1)
  
  save(Merged,file = paste0(out, "XR64PP_",  PatternList[i], "_PLUS&MINUS_core15.RData"))
}


######
## paths to cpd XRseq bw files
path<-"path_to_merged_replicate_XRseq_files/"
PatternList<-c("1h","4h","8h","16h","1d","2d")

for (i in 1:6) {
  load(file = paste0(path, "XRcpd_",  PatternList[i], "_R1_R2_PLUS_core15.RData"))
  load(file = paste0(path, "XRcpd_",  PatternList[i], "_R1_R2_MINUS_core15.RData"))
  
  Merged<-Pooled_plus[,c(1:6)]
  Merged$Plus_ave<-Pooled_plus[,9]
  Merged$Minus_ave<-Pooled_minus[,9]
  Merged$Strand_mean<-rowMeans(Merged[,c(7:8)],dims = 1)
  
  save(Merged,file = paste0(out, "XRcpd_",  PatternList[i], "_PLUS&MINUS_core15.RData"))
}


## END ##