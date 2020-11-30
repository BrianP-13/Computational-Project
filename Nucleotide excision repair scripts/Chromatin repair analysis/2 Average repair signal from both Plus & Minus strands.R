library(GenomicRanges)


### Description: Average together the total excision repair from BOTH strands ###

## Strategy:

# Step 1: Load data sets containing repair signal in chromatin states for IMR90 (both plus and minus strand repair datasets)
# Step 2: Prepare forloop script to average repair rates from both strands.

out<-"output_directory/"

######
## 64PP repair time points
path<-"path_to_merged_replicate_XRseq_files/"
PatternList<-c("5min","20min","1h","2h","4h")

for (i in 1:5) {
  load(file = paste0(path, "XR64PP_",  PatternList[i], "_R1_R2_PLUS_core15.RData"))
  Pooled_plus_Chrom15<-as.data.frame(Pooled_plus_Chrom15)
  
  load(file = paste0(path, "XR64PP_",  PatternList[i], "_R1_R2_MINUS_core15.RData"))
  Pooled_minus_Chrom15<-as.data.frame(Pooled_minus_Chrom15)
  
  Mean_strand_repair<-Pooled_plus_Chrom15[,1:6]
  Mean_strand_repair$XR64_plus_mean<-Pooled_plus_Chrom15$XR64_plus_mean
  Mean_strand_repair$XR64_minus_mean<-Pooled_minus_Chrom15$XR64_minus_mean
  Mean_strand_repair$XR64_strand_mean<-rowMeans(Mean_strand_repair[,7:8])
  colnames(Mean_strand_repair)[9]<-paste0('XR64_strand_mean_', PatternList[i])
  
  Mean_strand_repair<-makeGRangesFromDataFrame(Mean_strand_repair, keep.extra.columns = T)
  save(Mean_strand_repair,file = paste0(out, "XR64PP_Strand_mean_",  PatternList[i], "_core15.RData"))
}


######
## CPD repair
path<-"path_to_merged_replicate_XRseq_files/"
PatternList<-c("1h","4h","8h","16h","1d","2d")

for (i in 1:6) {
  load(file = paste0(path, "XRCPD_",  PatternList[i], "_R1_R2_PLUS_core15.RData"))
  Pooled_plus_Chrom15<-as.data.frame(Pooled_plus_Chrom15)
  
  load(file = paste0(path, "XRCPD_",  PatternList[i], "_R1_R2_MINUS_core15.RData"))
  Pooled_minus_Chrom15<-as.data.frame(Pooled_minus_Chrom15)
  
  Mean_strand_repair<-Pooled_plus_Chrom15[,1:6]
  Mean_strand_repair$XRCPD_plus_mean<-Pooled_plus_Chrom15$XRCPD_plus_mean
  Mean_strand_repair$XRCPD_minus_mean<-Pooled_minus_Chrom15$XRCPD_minus_mean
  Mean_strand_repair$XRCPD_strand_mean<-rowMeans(Mean_strand_repair[,7:8])
  colnames(Mean_strand_repair)[9]<-paste0('XRCPD_strand_mean_', PatternList[i])
  
  Mean_strand_repair<-makeGRangesFromDataFrame(Mean_strand_repair, keep.extra.columns = T)
  save(Mean_strand_repair,file = paste0(out, "XRCPD_Strand_mean_",  PatternList[i], "_core15.RData"))
}


### END ###
