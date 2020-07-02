library(ggplot2)
library(data.table)

## Description: This script creates a master table with all repair time points for each UV lesion type


out<-"output_directory/"


######
## 64PP
path<-"path_to_final_XRseq_files/"
PatternList<-c("5min","20min","1h","2h","4h")

# Create table with chromatin state coordinates
load(file = paste0(path, "XR64PP_",  PatternList[1], "_PLUS&MINUS_core15.RData"))
XR_table<-Merged[,c(1:6)]

# Add the averaged NER signal for each time point
for (i in 1:5) { 
  load(file = paste0(path, "XR64PP_",  PatternList[i], "_PLUS&MINUS_core15.RData"))
  XR_table<-cbind(XR_table, Merged$Strand_mean)
  colnames(XR_table)[6+i]<-PatternList[i]
}
save(XR_table, file = paste0(out,"XR64PP_PLUS&MINUS_ALL_TIME_POINTS_core15.RData"))


######
## CPD
path<-"path_to_final_XRseq_files/"
PatternList<-c("1h","4h","8h","16h","1d","2d")

# Create table with chromatin state coordinates
load(file = paste0(path, "XRcpd_",  PatternList[1], "_PLUS&MINUS_core15.RData"))
XR_table<-Merged[,c(1:6)]

# Add the averaged NER signal for each time point

for (i in 1:6) { 
  load(file = paste0(path, "XRcpd_",  PatternList[i], "_PLUS&MINUS_core15.RData"))
  XR_table<-cbind(XR_table, Merged$Strand_mean)
  colnames(XR_table)[6+i]<-PatternList[i]
}
save(XR_table, file = paste0(out,"XRcpd_PLUS&MINUS_ALL_TIME_POINTS_core15.RData"))


## END ##